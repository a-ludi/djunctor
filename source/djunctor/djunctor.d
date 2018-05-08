/**
    This is the main algorithm of this package.

    Copyright: © 2018 Arne Ludwig <arne.ludwig@posteo.de>
    License: Subject to the terms of the MIT license, as written in the
             included LICENSE file.
    Authors: Arne Ludwig <arne.ludwig@posteo.de>
*/
module djunctor.djunctor;

import djunctor.commandline : Options;
import djunctor.util.fasta : parseFasta, parseFastaRecord, parsePacBioHeader, reverseComplement;
import djunctor.util.log;
import djunctor.util.math : ceil, floor, mean, median;
import djunctor.util.range : Comparator;
import djunctor.util.string : indent;
import core.exception : AssertError;
import std.algorithm : all, any, canFind, chunkBy, each, equal, filter, find,
    fold, group, isSorted, joiner, map, max, maxElement, min, sort, sum, swap,
    SwapStrategy;
import std.array : appender, Appender, array;
import std.container : BinaryHeap, heapify, make;
import std.conv;
import std.exception : assertNotThrown, assertThrown;
import std.format : format, formattedWrite;
import std.math : abs, floor, sgn;
import std.range : assumeSorted, chain, chunks, ElementType, enumerate, iota,
    isForwardRange, only, refRange, retro, slide, SortedRange, tail, take, walkLength, zip;
import std.stdio : File, write, writeln;
import std.string : outdent;
import std.typecons : Flag, No, tuple, Tuple, Yes;
import vibe.data.json : Json, toJson = serializeToJson;

version (unittest)
{
    import djunctor.util.testing : MockCallable;
    import djunctor.dazzler : origGetLocalAlignments = getLocalAlignments,
        origGetMappings = getMappings;
    import djunctor.dazzler : buildDamFile, getConsensus, getFastaEntries,
        getNumContigs, getTracePointDistance, attachTracePoints;

    MockCallable!(AlignmentContainer!(AlignmentChain[]), const string, const Options) getLocalAlignments;
    MockCallable!(origGetMappings!Options) getMappings;
}
else
{
    import djunctor.dazzler : buildDamFile, getConsensus, getFastaEntries,
        getLocalAlignments, getMappings, getNumContigs, getTracePointDistance,
        attachTracePoints;
}

/// General container for alignment data.
template AlignmentContainer(R)
{
    struct AlignmentContainer
    {
        static immutable orders = tuple("a2b", "b2a");
        R a2b;
        R b2a;

        static auto getOrdered(T)(in string order, T thingA, T thingB) pure nothrow
        {
            if (order == "a2b")
                return tuple(thingA, thingB);
            else if (order == "b2a")
                return tuple(thingB, thingA);
            else
                assert(0, "illegal value");
        }

        unittest
        {
            assert(getOrdered("a2b", 0, 1) == tuple(0, 1));
            assert(getOrdered("b2a", 0, 1) == tuple(1, 0));
            assertThrown!AssertError(getOrdered("foo", 0, 1));
        }

        static auto getOrdered(string order, T)(T thingA, T thingB) pure nothrow
                if (order == "a2b" || order == "b2a")
        {
            static if (order == "a2b")
                return tuple(thingA, thingB);
            else
                return tuple(thingB, thingA);
        }

        unittest
        {
            assert(getOrdered!"a2b"(0, 1) == tuple(0, 1));
            assert(getOrdered!"b2a"(0, 1) == tuple(1, 0));
            assert(!__traits(compiles, getOrdered!"foo"(0, 1)));
        }
    }
}

/**
    Holds a chain of local alignments that form a compound alignment. An AlignmentChain should
    contain at least one element.
*/
struct AlignmentChain
{
    static struct LocalAlignment
    {
        static struct Locus
        {
            size_t begin;
            size_t end;
        }

        static struct TracePoint
        {
            size_t numDiffs;
            size_t numBasePairs;
        }

        Locus contigA;
        Locus contigB;
        size_t numDiffs;
        TracePoint[] tracePoints;
    }

    static struct Contig
    {
        size_t id;
        size_t length;
    }

    static enum Complement
    {
        no = false,
        yes = true,
    }

    static enum Order : int
    {
        none,
        ref2read,
        read2ref,
    }

    static immutable maxScore = 2 ^^ 16;

    size_t id;
    Contig contigA;
    Contig contigB;
    Complement complement;
    LocalAlignment[] localAlignments;
    size_t tracePointDistance;
    Order order;

    invariant
    {
        assert(localAlignments.length >= 1, "empty chain is forbidden");
        foreach (la; localAlignments)
        {
            assert(0 <= la.contigA.begin && la.contigA.begin < la.contigA.end
                    && la.contigA.end <= contigA.length, "non-sense alignment of contigA");
            assert(0 <= la.contigB.begin && la.contigB.begin < la.contigB.end
                    && la.contigB.end <= contigB.length, "non-sense alignment of contigB");

            assert(tracePointDistance == 0 || la.tracePoints.length > 0,
                    "missing trace points");
            if (tracePointDistance > 0)
            {
                size_t traceLength = la.tracePoints.map!"a.numBasePairs".sum;
                size_t traceDiffs = la.tracePoints.map!"a.numDiffs".sum;

                assert(la.numDiffs == traceDiffs, "missing trace points");
                assert(la.contigB.end - la.contigB.begin == traceLength,
                        "trace distance does not match alignment");
            }
        }
    }

    unittest
    {
        with (Complement) with (LocalAlignment)
            {
                auto acZeroLength = AlignmentChain(0, Contig(1, 10), Contig(1, 10), no, []);
                auto ac1 = AlignmentChain(1, Contig(1, 10), Contig(1, 10), no,
                        [LocalAlignment(Locus(1, 1), Locus(1, 9), 0)]);
                auto ac2 = AlignmentChain(2, Contig(1, 10), Contig(1, 10), no,
                        [LocalAlignment(Locus(1, 11), Locus(1, 9), 0)]);
                auto ac3 = AlignmentChain(3, Contig(1, 10), Contig(1, 10), no,
                        [LocalAlignment(Locus(1, 9), Locus(1, 1), 0)]);
                auto ac4 = AlignmentChain(4, Contig(1, 10), Contig(1, 10), no,
                        [LocalAlignment(Locus(1, 9), Locus(1, 11), 0)]);
                auto acFine = AlignmentChain(5, Contig(1, 10), Contig(1, 10),
                        no, [LocalAlignment(Locus(1, 9), Locus(1, 9), 0)]);

                assertThrown!AssertError(acZeroLength.totalLength);
                assertThrown!AssertError(ac1.totalLength);
                assertThrown!AssertError(ac2.totalLength);
                assertThrown!AssertError(ac3.totalLength);
                assertThrown!AssertError(ac4.totalLength);
                assertNotThrown!AssertError(acFine.totalLength);
            }
    }

    @property ref const(LocalAlignment) first() const pure nothrow
    {
        return localAlignments[0];
    }

    unittest
    {
        with (Complement) with (LocalAlignment)
            {
                auto firstLA = LocalAlignment(Locus(1, 9), Locus(1, 9), 0);
                auto otherLA = LocalAlignment(Locus(9, 10), Locus(9, 10), 0);
                auto ac = AlignmentChain(0, Contig(1, 10), Contig(1, 10), no, [firstLA, otherLA]);

                assert(ac.first == firstLA);
            }
    }

    @property ref const(LocalAlignment) last() const pure nothrow
    {
        return localAlignments[$ - 1];
    }

    unittest
    {
        with (Complement) with (LocalAlignment)
            {
                auto lastLA = LocalAlignment(Locus(1, 9), Locus(1, 9), 0);
                auto otherLA = LocalAlignment(Locus(9, 10), Locus(9, 10), 0);
                auto ac = AlignmentChain(0, Contig(1, 10), Contig(1, 10), no, [otherLA, lastLA]);

                assert(ac.last == lastLA);
            }
    }

    /**
        Returns true if the aligned read `contigB` (with extensions on either
        end) is fully contained in the reference `contigA`.

        According to the following 'definition' `contigA` is fully contained
        in `contigB` iff `x >= 0` and `y <= l_a`.
        ---
                0                  x   ua      va  y                        la
        contigA |------------------+---+-+---+-+---+-------------------------|
                                  / / /  | | |  \ \ \
                                 / / /   | | |   \ \ \
                                / / /    | | |    \ \ \
                               / / /     | | |     \ \ \
        contigB               |---+------+---+------+---|
                              0   ub                vb lb

        x = ua - (ub - 0) = ua - ub
        y = va + (lb - vb)
        ---
    */
    bool isFullyContained()
    {
        if (first.contigB.begin > first.contigA.begin)
        {
            // x < 0; return early to avoid negative numbers in unsigned integers
            return false;
        }

        auto x = first.contigA.begin - first.contigB.begin;
        auto y = last.contigA.end + contigB.length - last.contigB.end;

        return 0 <= x && y < contigA.length;
    }

    unittest
    {
        with (Complement) with (LocalAlignment)
            {
                auto la = LocalAlignment(Locus(30, 35), Locus(5, 10), 1);
                auto ac = AlignmentChain(0, Contig(1, 50), Contig(1, 15), no, [la]);

                // read with extension align an contigA from 25 to 40
                assert(ac.isFullyContained);
            }

        with (Complement) with (LocalAlignment)
            {
                auto la1 = LocalAlignment(Locus(10, 20), Locus(5, 10), 1);
                auto la2 = LocalAlignment(Locus(30, 40), Locus(5, 10), 1);
                auto ac = AlignmentChain(0, Contig(1, 50), Contig(1, 15), no, [la1, la2]);

                // read with extension align an contigA from 5 to 45
                assert(ac.isFullyContained);
            }

        with (Complement) with (LocalAlignment)
            {
                auto la = LocalAlignment(Locus(0, 10), Locus(5, 10), 1);
                auto ac = AlignmentChain(0, Contig(1, 50), Contig(1, 15), no, [la]);

                // read with extension align an contigA from -5 to 15
                assert(!ac.isFullyContained);
            }

        with (Complement) with (LocalAlignment)
            {
                auto la = LocalAlignment(Locus(40, 50), Locus(5, 10), 1);
                auto ac = AlignmentChain(0, Contig(1, 50), Contig(1, 15), no, [la]);

                // read with extension align an contigA from 35 to 55
                assert(!ac.isFullyContained);
            }

        with (Complement) with (LocalAlignment)
            {
                auto la1 = LocalAlignment(Locus(0, 20), Locus(5, 10), 1);
                auto la2 = LocalAlignment(Locus(30, 50), Locus(5, 10), 1);
                auto ac = AlignmentChain(0, Contig(1, 50), Contig(1, 15), no, [la1, la2]);

                // read with extension align an contigA from -5 to 55
                assert(!ac.isFullyContained);
            }
    }

    @property size_t totalLength() const pure
    {
        return last.contigA.end - first.contigA.begin;
    }

    unittest
    {
        with (Complement) with (LocalAlignment)
            {
                auto la1 = LocalAlignment(Locus(1, 3), Locus(1, 3), 1);
                auto la2 = LocalAlignment(Locus(5, 10), Locus(5, 10), 2);
                auto ac = AlignmentChain(0, Contig(1, 10), Contig(1, 10), no, [la1, la2]);

                assert(ac.totalLength == 9);
            }
    }

    @property size_t totalDiffs() const pure
    {
        return localAlignments.map!"a.numDiffs".sum;
    }

    unittest
    {
        with (Complement) with (LocalAlignment)
            {
                auto la1 = LocalAlignment(Locus(1, 3), Locus(1, 3), 1);
                auto la2 = LocalAlignment(Locus(5, 10), Locus(5, 10), 2);
                auto ac = AlignmentChain(0, Contig(1, 10), Contig(1, 10), no, [la1, la2]);

                assert(ac.totalDiffs == 3);
            }
    }

    @property size_t totalGapLength() const pure
    {
        // dfmt off
        return localAlignments
            .chunks(2)
            .map!(las => las.length < 2 ? 0 : las[1].contigA.begin - las[0].contigA.end)
            .sum;
        // dfmt on
    }

    unittest
    {
        with (Complement) with (LocalAlignment)
            {
                auto la1 = LocalAlignment(Locus(1, 3), Locus(1, 3), 1);
                auto la2 = LocalAlignment(Locus(5, 10), Locus(5, 10), 2);
                auto ac = AlignmentChain(0, Contig(1, 10), Contig(1, 10), no, [la1, la2]);

                assert(ac.totalGapLength == 2);
            }
    }

    @property size_t numMatchingBps() const pure
    {
        return totalLength - (totalDiffs + totalGapLength);
    }

    unittest
    {
        with (Complement) with (LocalAlignment)
            {
                auto la1 = LocalAlignment(Locus(1, 3), Locus(1, 3), 1);
                auto la2 = LocalAlignment(Locus(5, 10), Locus(5, 10), 2);
                auto ac = AlignmentChain(0, Contig(1, 10), Contig(1, 10), no, [la1, la2]);

                assert(ac.numMatchingBps == 9 - (3 + 2));
            }
    }

    @property size_t score() const pure
    {
        return numMatchingBps * maxScore / totalLength;
    }

    unittest
    {
        with (Complement) with (LocalAlignment)
            {
                auto la1 = LocalAlignment(Locus(1, 3), Locus(1, 3), 1);
                auto la2 = LocalAlignment(Locus(5, 10), Locus(5, 10), 2);
                auto ac = AlignmentChain(0, Contig(1, 10), Contig(1, 10), no, [la1, la2]);

                assert(ac.score == 4 * maxScore / 9);
            }
    }

    int compareIds(ref const AlignmentChain other) const pure nothrow
    {
        long idCompare = this.contigA.id - other.contigA.id;
        if (idCompare != 0)
            return cast(int) sgn(idCompare);

        idCompare = this.contigB.id - other.contigB.id;
        if (idCompare != 0)
            return cast(int) sgn(idCompare);

        return 0;
    }

    unittest
    {
        with (Complement) with (LocalAlignment)
            {
                auto la = LocalAlignment(Locus(0, 1), Locus(0, 1), 1);
                // dfmt off
                auto acs = [
                    AlignmentChain(0, Contig(1, 10), Contig(1, 10), no, [la]),
                    AlignmentChain(1, Contig(1, 10), Contig(2, 10), no, [la]),
                    AlignmentChain(2, Contig(2, 10), Contig(1, 10), no, [la]),
                    AlignmentChain(3, Contig(2, 10), Contig(2, 10), no, [la]),
                ];
                // dfmt on

                foreach (i; 0 .. acs.length)
                    foreach (j; 0 .. acs.length)
                    {
                        auto compareValue = acs[i].opCmp(acs[j]);
                        alias errorMessage = (cmp) => format!"expected %s and %s to compare %s 0 but got %d"(
                                acs[i], acs[j], cmp, compareValue);

                        if (i < j)
                            assert(compareValue < 0, errorMessage("<"));
                        else if (i > j)
                            assert(compareValue > 0, errorMessage(">"));
                        else
                            assert(compareValue == 0, errorMessage("=="));
                    }
            }
    }

    int opCmp(ref const AlignmentChain other) const pure nothrow
    {
        const int idCompare = this.compareIds(other);
        if (idCompare != 0)
            return idCompare;

        long locusCompare = this.first.contigA.begin - other.first.contigA.begin;
        if (locusCompare != 0)
            return cast(int) sgn(locusCompare);

        locusCompare = this.first.contigB.begin - other.first.contigB.begin;
        if (locusCompare != 0)
            return cast(int) sgn(locusCompare);

        locusCompare = this.last.contigA.end - other.last.contigA.end;
        if (locusCompare != 0)
            return cast(int) sgn(locusCompare);

        locusCompare = this.last.contigB.end - other.last.contigB.end;
        if (locusCompare != 0)
            return cast(int) sgn(locusCompare);

        return 0;
    }

    unittest
    {
        // see compareIds
        with (Complement) with (LocalAlignment)
            {
                auto la = LocalAlignment(Locus(0, 1), Locus(0, 1), 1);
                // dfmt off
                auto acs = [
                    AlignmentChain(0, Contig(1, 10), Contig(1, 10), no, [la]),
                    AlignmentChain(1, Contig(1, 10), Contig(2, 10), no, [la]),
                    AlignmentChain(2, Contig(2, 10), Contig(1, 10), no, [la]),
                    AlignmentChain(3, Contig(2, 10), Contig(2, 10), no, [la]),
                ];
                // dfmt on

                foreach (i; 0 .. acs.length)
                    foreach (j; 0 .. acs.length)
                    {
                        auto compareValue = acs[i].opCmp(acs[j]);
                        alias errorMessage = (cmp) => format!"expected %s and %s to compare %s 0 but got %d"(
                                acs[i], acs[j], cmp, compareValue);

                        if (i < j)
                            assert(compareValue < 0, errorMessage("<"));
                        else if (i > j)
                            assert(compareValue > 0, errorMessage(">"));
                        else
                            assert(compareValue == 0, errorMessage("=="));
                    }
            }

        // test non-id-related comparison
        with (Complement) with (LocalAlignment)
            {
                // dfmt off
                auto acs = [
                    AlignmentChain(0, Contig(1, 10), Contig(1, 10), no, [
                        LocalAlignment(Locus(0, 2), Locus(0, 2), 1),
                        LocalAlignment(Locus(3, 5), Locus(3, 5), 1)
                    ]),
                    AlignmentChain(1, Contig(1, 10), Contig(1, 10), no, [
                        LocalAlignment(Locus(0, 2), Locus(1, 2), 1),
                        LocalAlignment(Locus(3, 5), Locus(3, 5), 1)
                    ]),
                    AlignmentChain(2, Contig(1, 10), Contig(1, 10), no, [
                        LocalAlignment(Locus(1, 2), Locus(0, 2), 1),
                        LocalAlignment(Locus(3, 5), Locus(3, 5), 1)
                    ]),
                    AlignmentChain(3, Contig(1, 10), Contig(1, 10), no, [
                        LocalAlignment(Locus(1, 2), Locus(1, 2), 1),
                        LocalAlignment(Locus(3, 5), Locus(3, 5), 1)
                    ]),
                    AlignmentChain(4, Contig(1, 10), Contig(1, 10), no, [
                        LocalAlignment(Locus(1, 2), Locus(1, 2), 1),
                        LocalAlignment(Locus(3, 5), Locus(3, 6), 1)
                    ]),
                    AlignmentChain(5, Contig(1, 10), Contig(1, 10), no, [
                        LocalAlignment(Locus(1, 2), Locus(1, 2), 1),
                        LocalAlignment(Locus(3, 6), Locus(3, 5), 1)
                    ]),
                    AlignmentChain(6, Contig(1, 10), Contig(1, 10), no, [
                        LocalAlignment(Locus(1, 2), Locus(1, 2), 1),
                        LocalAlignment(Locus(3, 6), Locus(3, 6), 1)
                    ]),
                ];
                // dfmt on

                foreach (i; 0 .. acs.length)
                    foreach (j; 0 .. acs.length)
                    {
                        auto compareValue = acs[i].opCmp(acs[j]);
                        alias errorMessage = (cmp) => format!"expected %s and %s to compare %s 0 but got %d"(
                                acs[i], acs[j], cmp, compareValue);

                        if (i < j)
                            assert(compareValue < 0, errorMessage("<"));
                        else if (i > j)
                            assert(compareValue > 0, errorMessage(">"));
                        else
                            assert(compareValue == 0, errorMessage("=="));
                    }
            }
    }
}

AlignmentChain.Order reverse(in AlignmentChain.Order order) pure nothrow
{
    final switch (order)
    {
        case AlignmentChain.Order.none:
            return AlignmentChain.Order.none;
        case AlignmentChain.Order.ref2read:
            return AlignmentChain.Order.read2ref;
        case AlignmentChain.Order.read2ref:
            return AlignmentChain.Order.ref2read;
    }
}

AlignmentContainer!(AlignmentChain[]) annotateOrder(ref AlignmentContainer!(AlignmentChain[]) alignmentContainer, in AlignmentChain.Order a2bOrder) nothrow
{
    foreach (ref alignmentChain; alignmentContainer.a2b)
    {
        alignmentChain.order = a2bOrder;
    }

    auto b2aOrder = a2bOrder.reverse;
    foreach (ref alignmentChain; alignmentContainer.b2a)
    {
        alignmentChain.order = b2aOrder;
    }

    return alignmentContainer;
}

bool idsPred(in AlignmentChain ac1, in AlignmentChain ac2) pure
{
    auto cmpValue = ac1.compareIds(ac2);

    return 0 != cmpValue && cmpValue < 0;
}

unittest
{
    with (AlignmentChain) with (LocalAlignment) with (Complement)
            {
                auto la = LocalAlignment(Locus(0, 1), Locus(0, 1), 1);
                // dfmt off
                auto acs = [
                    AlignmentChain(0, Contig(1, 10), Contig(1, 10), no, [la]),
                    AlignmentChain(1, Contig(1, 10), Contig(2, 10), no, [la]),
                    AlignmentChain(2, Contig(2, 10), Contig(1, 10), no, [la]),
                    AlignmentChain(3, Contig(2, 10), Contig(2, 10), no, [la]),
                ];
                // dfmt on

                foreach (i; 0 .. acs.length)
                    foreach (j; 0 .. acs.length)
                    {
                        auto compareValue = idsPred(acs[i], acs[j]);
                        alias errorMessage = (expValue) => format!"expected idsPred(%s, %s) to be %s but got %s"(
                                acs[i], acs[j], expValue, compareValue);

                        if (i < j)
                            assert(compareValue, errorMessage(true));
                        else
                            assert(!compareValue, errorMessage(false));
                    }
            }
}

auto haveEqualIds(in AlignmentChain ac1, in AlignmentChain ac2) pure
{
    auto cmpValue = ac1.compareIds(ac2);

    return 0 == cmpValue;
}

unittest
{
    with (AlignmentChain) with (Complement) with (LocalAlignment)
            {
                // dfmt off
                auto sortedTestChains = [
                    AlignmentChain(0, Contig(1, 10), Contig(1, 10), no, [LocalAlignment(Locus(0, 10), Locus(0, 1), 0)]),
                    AlignmentChain(1, Contig(1, 10), Contig(2, 20), no, [LocalAlignment(Locus(0, 10), Locus(0, 2), 0)]),
                    AlignmentChain(2, Contig(1, 10), Contig(3, 30), no, [LocalAlignment(Locus(0, 10), Locus(0, 3), 0)]),
                    AlignmentChain(3, Contig(2, 20), Contig(1, 10), no, [LocalAlignment(Locus(0, 20), Locus(0, 4), 0)]),
                    AlignmentChain(4, Contig(2, 20), Contig(1, 10), no, [LocalAlignment(Locus(0, 20), Locus(0, 5), 0)]),
                    AlignmentChain(5, Contig(2, 20), Contig(3, 30), no, [LocalAlignment(Locus(0, 20), Locus(0, 6), 0)]),
                    AlignmentChain(6, Contig(3, 30), Contig(1, 10), no, [LocalAlignment(Locus(0, 30), Locus(0, 7), 0)]),
                    AlignmentChain(7, Contig(3, 30), Contig(2, 20), no, [LocalAlignment(Locus(0, 30), Locus(0, 8), 0)]),
                    AlignmentChain(8, Contig(3, 30), Contig(3, 30), no, [LocalAlignment(Locus(0, 30), Locus(0, 9), 0)]),
                ];

                assert(sortedTestChains.chunkBy!haveEqualIds.equal!equal([
                    sortedTestChains[0..1],
                    sortedTestChains[1..2],
                    sortedTestChains[2..3],
                    sortedTestChains[3..5],
                    sortedTestChains[5..6],
                    sortedTestChains[6..7],
                    sortedTestChains[7..8],
                    sortedTestChains[8..9],
                ]));
                // dfmt on
            }
}

auto equalIdsRange(in AlignmentChain[] acList, in size_t contigAID, in size_t contigBID) pure
{
    assert(isSorted!idsPred(acList));
    // dfmt off
    AlignmentChain needle = {
        contigA: AlignmentChain.Contig(contigAID, 1),
        contigB: AlignmentChain.Contig(contigBID, 1),
        localAlignments: [AlignmentChain.LocalAlignment(AlignmentChain.LocalAlignment.Locus(0, 1), AlignmentChain.LocalAlignment.Locus(0, 1))],
    };
    // dfmt on

    return acList.assumeSorted!idsPred.equalRange(needle);
}

unittest
{
    with (AlignmentChain) with (Complement) with (LocalAlignment)
            {
                // dfmt off
                auto sortedTestChains = [
                    AlignmentChain(0, Contig(1, 10), Contig(1, 10), no, [LocalAlignment(Locus(0, 1), Locus(0, 1), 0)]),
                    AlignmentChain(1, Contig(1, 10), Contig(2, 20), no, [LocalAlignment(Locus(0, 2), Locus(0, 2), 0)]),
                    AlignmentChain(2, Contig(1, 10), Contig(3, 30), no, [LocalAlignment(Locus(0, 3), Locus(0, 3), 0)]),
                    AlignmentChain(3, Contig(2, 20), Contig(1, 10), no, [LocalAlignment(Locus(0, 4), Locus(0, 4), 0)]),
                    AlignmentChain(4, Contig(2, 20), Contig(1, 10), no, [LocalAlignment(Locus(0, 5), Locus(0, 5), 0)]),
                    AlignmentChain(5, Contig(2, 20), Contig(3, 30), no, [LocalAlignment(Locus(0, 6), Locus(0, 6), 0)]),
                    AlignmentChain(6, Contig(3, 30), Contig(1, 10), no, [LocalAlignment(Locus(0, 7), Locus(0, 7), 0)]),
                    AlignmentChain(7, Contig(3, 30), Contig(2, 20), no, [LocalAlignment(Locus(0, 8), Locus(0, 8), 0)]),
                    AlignmentChain(8, Contig(3, 30), Contig(3, 30), no, [LocalAlignment(Locus(0, 9), Locus(0, 9), 0)]),
                ];
                // dfmt on

                assert(sortedTestChains.equalIdsRange(1, 1).equal(sortedTestChains[0 .. 1]));
                assert(sortedTestChains.equalIdsRange(2, 1).equal(sortedTestChains[3 .. 5]));
                assert(sortedTestChains.equalIdsRange(3, 1).equal(sortedTestChains[6 .. 7]));
                assert(sortedTestChains.equalIdsRange(42, 1337).equal(sortedTestChains[0 .. 0]));
            }
}

/**
    Returns true iff ac1 begins before ac2.

    **Note:** the order is relative to the orientation of the opposite contig.
*/
bool isBefore(string contig)(in AlignmentChain ac1, in AlignmentChain ac2) pure
        if (contig == "contigA" || contig == "contigB")
{
    assert(__traits(getMember, ac1, contig) == __traits(getMember, ac2, contig),
            "alignment chains do not belong to the same contig");

    return __traits(getMember, ac1.first, contig).begin < __traits(getMember,
            ac2.first, contig).begin;
}

unittest
{
    with (AlignmentChain) with (Complement) with (LocalAlignment)
            {
                // dfmt off
                auto acs = [
                    AlignmentChain(0, Contig(1, 10), Contig(1, 10), no, [LocalAlignment(Locus(1, 6), Locus(0, 1), 0)]),
                    AlignmentChain(1, Contig(1, 10), Contig(2, 10), yes, [LocalAlignment(Locus(2, 6), Locus(0, 1), 0)]),
                    AlignmentChain(2, Contig(1, 10), Contig(3, 10), no, [LocalAlignment(Locus(3, 6), Locus(0, 1), 0)]),
                    AlignmentChain(3, Contig(1, 10), Contig(4, 10), yes, [LocalAlignment(Locus(4, 6), Locus(0, 1), 0)]),
                    AlignmentChain(4, Contig(1, 10), Contig(5, 10), no, [LocalAlignment(Locus(5, 6), Locus(0, 1), 0)]),
                ];
                // dfmt on

                foreach (i; 0 .. acs.length)
                    foreach (j; 0 .. acs.length)
                    {
                        auto compareValue = isBefore!"contigA"(acs[i], acs[j]);
                        alias errorMessage = (expValue) => format!"expected idsPred(%s, %s) to be %s but got %s"(
                                acs[i], acs[j], expValue, compareValue);

                        if (i < j)
                            assert(compareValue, errorMessage(true));
                        else
                            assert(!compareValue, errorMessage(false));
                    }
            }
}

/**
    If readAlignment is a gap return true iff the first alignment begins
    before the second alignment on the reference (contigB); otherwise
    returns true.
*/
bool isInReferenceOrder(in ReadAlignment readAlignment) pure nothrow
{
    return !readAlignment.isGap || isBefore!"contigA"(readAlignment[0],
            readAlignment[1]) != readAlignment[0].complement;
}

/// Start the `djunctor` algorithm with preprocessed options.
void runWithOptions(in ref Options options)
{
    new DJunctor(options).run();
}

/**
    Alignment of a read against the reference. This is either one or two
    alignment chains which belong to the same read and one or two reference
    contig(s).
*/
alias ReadAlignment = AlignmentChain[];

/// Type of the read alignment.
static enum ReadAlignmentType
{
    front = 0,
    gap = 1,
    back = 2,
}

/**
    Returns true iff the read alignment is valid, ie. it is either an
    extension or gap.
*/
bool isValid(in ReadAlignment readAlignment) pure nothrow
{
    assert(readAlignment[0].order == AlignmentChain.Order.read2ref, "only applicable for read2ref alignments");
    return readAlignment.isExtension ^ readAlignment.isGap;
}

/**
    Get the type of the read alignment.

    See_Also: `isFrontExtension`, `isBackExtension`, `isGap`
*/
ReadAlignmentType getType(in ReadAlignment readAlignment) pure nothrow
{
    assert(readAlignment[0].order == AlignmentChain.Order.read2ref, "only applicable for read2ref alignments");
    assert(readAlignment.isValid, "invalid read alignment");

    if (readAlignment.isGap)
    {
        return ReadAlignmentType.gap;
    }
    else if (readAlignment.isFrontExtension)
    {
        return ReadAlignmentType.front;
    }
    else
    {
        return ReadAlignmentType.back;
    }
}

/**
    Returns true iff the read alignment is an extension, ie. it is a front or
    back extension.

    See_Also: `isFrontExtension`, `isBackExtension`
*/
bool isExtension(in ReadAlignment readAlignment) pure nothrow
{
    assert(readAlignment[0].order == AlignmentChain.Order.read2ref, "only applicable for read2ref alignments");
    return readAlignment.isFrontExtension ^ readAlignment.isBackExtension;
}

/**
    Returns true iff the read alignment is an front extension, ie. it is an
    extension and reaches over the front of the reference contig.

    ---
    Case 1 (complement alignment):

              0             ry lr
        ref   |--<--<--<-+-<-+--|
                         | | |
        read          |--+->-+->-->-->--|
                      0     ay         la

    Case 2 (non-complement alignment):

                      0  rx
        ref           |--+->-+->-->-->--|
                         | | |
        read  |-->-->-->-+->-+--|
              0          ax
    ---
*/
bool isFrontExtension(in ReadAlignment readAlignment) pure nothrow
{
    assert(readAlignment[0].order == AlignmentChain.Order.read2ref, "only applicable for read2ref alignments");
    if (readAlignment.length != 1)
    {
        return false;
    }

    auto alignment = readAlignment[0];

    if (alignment.complement)
    {
        auto readExtensionLength = alignment.contigA.length - alignment.last.contigA.end;
        auto referenceExtensionLength = alignment.contigB.length - alignment.last.contigB.end;

        return readExtensionLength > referenceExtensionLength;
    }
    else
    {
        auto readExtensionLength = alignment.first.contigA.begin;
        auto referenceExtensionLength = alignment.first.contigB.begin;

        return readExtensionLength > referenceExtensionLength;
    }
}

/**
    Returns true iff the read alignment is an back extension, ie. it is an
    extension and reaches over the back of the reference contig.

    ---
    Case 1 (complement alignment):

                      0  rx
        ref           |--+-<-+-<--<--<--|
                         | | |
        read  |-->-->-->-+->-+--|
              0          ax

    Case 2 (non-complement alignment):

              0             ry lr
        ref   |-->-->-->-+->-+--|
                         | | |
        read          |--+->-+->-->-->--|
                      0     ay         la
    ---
*/
bool isBackExtension(in ReadAlignment readAlignment) pure nothrow
{
    assert(readAlignment[0].order == AlignmentChain.Order.read2ref, "only applicable for read2ref alignments");
    if (readAlignment.length != 1)
    {
        return false;
    }

    auto alignment = readAlignment[0];

    if (alignment.complement)
    {
        auto readExtensionLength = alignment.first.contigA.begin;
        auto referenceExtensionLength = alignment.first.contigB.begin;

        return readExtensionLength > referenceExtensionLength;
    }
    else
    {
        auto readExtensionLength = alignment.contigA.length - alignment.last.contigA.end;
        auto referenceExtensionLength = alignment.contigB.length - alignment.last.contigB.end;

        return readExtensionLength > referenceExtensionLength;
    }
}

/**
    Returns true iff the read alignment spans a gap, ie. two alignments on
    different reference contigs with different individual extension type are
    involved.

    ---
    Case 1 (complement alignment):

              0             ry lr   0             ry lr
        ref   |--<--<--<-+-<-+--|   |--+-<-+-<-+-<-+--|
                         | | |         | | |
        read          |--+->-+->-->-->-+->-+--|
                      0     ay         la

    Case 1 (complement alignment):

              0             ry lr1  0 rx             lr2
        ref   |-->-->-->-+->-+--|   |--+->-+->-+->-+--|
                         | | |         | | |
        read          |--+->-+->-->-->-+->-+--|
                      0     ay         bx     la
    ---
*/
bool isGap(in ReadAlignment readAlignment) pure nothrow
{
    // dfmt off
    return readAlignment.length == 2 &&
        readAlignment[0].contigB.id != readAlignment[1].contigB.id &&
        (
            (readAlignment[0 .. 1].isBackExtension && readAlignment[1 .. 2].isFrontExtension) ||
            (readAlignment[0 .. 1].isFrontExtension && readAlignment[1 .. 2].isBackExtension)
        );
    // dfmt on
}

unittest
{
    with (AlignmentChain) with (LocalAlignment) with (Complement)
            {
                // dfmt off
                auto testCases = [
                    "innerAlignmentComplement": tuple(
                        [
                            AlignmentChain(
                                1,
                                Contig(1, 10),
                                Contig(1, 100),
                                no,
                                [
                                    LocalAlignment(
                                        Locus(0, 1),
                                        Locus(10, 11),
                                        0,
                                    ),
                                    LocalAlignment(
                                        Locus(9, 10),
                                        Locus(19, 20),
                                        0,
                                    ),
                                ],
                                0,
                                Order.read2ref,
                            ),
                        ],
                        No.isValid,
                        ReadAlignmentType.front,
                        No.isExtension,
                        No.isFrontExtension,
                        No.isBackExtension,
                        No.isGap,
                    ),
                    "innerAlignmentComplement": tuple(
                        [
                            AlignmentChain(
                                2,
                                Contig(1, 10),
                                Contig(1, 100),
                                yes,
                                [
                                    LocalAlignment(
                                        Locus(0, 1),
                                        Locus(10, 11),
                                        0,
                                    ),
                                    LocalAlignment(
                                        Locus(9, 10),
                                        Locus(19, 20),
                                        0,
                                    ),
                                ],
                                0,
                                Order.read2ref,
                            ),
                        ],
                        No.isValid,
                        ReadAlignmentType.front,
                        No.isExtension,
                        No.isFrontExtension,
                        No.isBackExtension,
                        No.isGap,
                    ),
                    "frontExtension": tuple(
                        [
                            AlignmentChain(
                                3,
                                Contig(1, 10),
                                Contig(1, 100),
                                no,
                                [
                                    LocalAlignment(
                                        Locus(5, 6),
                                        Locus(2, 3),
                                        0,
                                    ),
                                    LocalAlignment(
                                        Locus(9, 10),
                                        Locus(5, 6),
                                        0,
                                    ),
                                ],
                                0,
                                Order.read2ref,
                            ),
                        ],
                        Yes.isValid,
                        ReadAlignmentType.front,
                        Yes.isExtension,
                        Yes.isFrontExtension,
                        No.isBackExtension,
                        No.isGap,
                    ),
                    "frontExtensionComplement": tuple(
                        [
                            AlignmentChain(
                                4,
                                Contig(1, 10),
                                Contig(1, 100),
                                yes,
                                [
                                    LocalAlignment(
                                        Locus(0, 1),
                                        Locus(94, 95),
                                        0,
                                    ),
                                    LocalAlignment(
                                        Locus(4, 5),
                                        Locus(97, 98),
                                        0,
                                    ),
                                ],
                                0,
                                Order.read2ref,
                            ),
                        ],
                        Yes.isValid,
                        ReadAlignmentType.front,
                        Yes.isExtension,
                        Yes.isFrontExtension,
                        No.isBackExtension,
                        No.isGap,
                    ),
                    "backExtension": tuple(
                        [
                            AlignmentChain(
                                5,
                                Contig(1, 10),
                                Contig(1, 100),
                                no,
                                [
                                    LocalAlignment(
                                        Locus(0, 1),
                                        Locus(94, 95),
                                        0,
                                    ),
                                    LocalAlignment(
                                        Locus(4, 5),
                                        Locus(97, 98),
                                        0,
                                    ),
                                ],
                                0,
                                Order.read2ref,
                            ),
                        ],
                        Yes.isValid,
                        ReadAlignmentType.back,
                        Yes.isExtension,
                        No.isFrontExtension,
                        Yes.isBackExtension,
                        No.isGap,
                    ),
                    "backExtensionComplement": tuple(
                        [
                            AlignmentChain(
                                6,
                                Contig(1, 10),
                                Contig(1, 100),
                                yes,
                                [
                                    LocalAlignment(
                                        Locus(5, 6),
                                        Locus(2, 3),
                                        0,
                                    ),
                                    LocalAlignment(
                                        Locus(9, 10),
                                        Locus(5, 6),
                                        0,
                                    ),
                                ],
                                0,
                                Order.read2ref,
                            ),
                        ],
                        Yes.isValid,
                        ReadAlignmentType.back,
                        Yes.isExtension,
                        No.isFrontExtension,
                        Yes.isBackExtension,
                        No.isGap,
                    ),
                    "gap": tuple(
                        [
                            AlignmentChain(
                                7,
                                Contig(1, 10),
                                Contig(1, 100),
                                no,
                                [
                                    LocalAlignment(
                                        Locus(0, 1),
                                        Locus(94, 95),
                                        0,
                                    ),
                                    LocalAlignment(
                                        Locus(4, 5),
                                        Locus(97, 98),
                                        0,
                                    ),
                                ],
                                0,
                                Order.read2ref,
                            ),
                            AlignmentChain(
                                8,
                                Contig(1, 10),
                                Contig(2, 100),
                                no,
                                [
                                    LocalAlignment(
                                        Locus(5, 6),
                                        Locus(2, 3),
                                        0,
                                    ),
                                    LocalAlignment(
                                        Locus(9, 10),
                                        Locus(5, 6),
                                        0,
                                    ),
                                ],
                                0,
                                Order.read2ref,
                            ),
                        ],
                        Yes.isValid,
                        ReadAlignmentType.gap,
                        No.isExtension,
                        No.isFrontExtension,
                        No.isBackExtension,
                        Yes.isGap,
                    ),
                    "gapComplement": tuple(
                        [
                            AlignmentChain(
                                9,
                                Contig(1, 10),
                                Contig(1, 100),
                                yes,
                                [
                                    LocalAlignment(
                                        Locus(5, 6),
                                        Locus(2, 3),
                                        0,
                                    ),
                                    LocalAlignment(
                                        Locus(9, 10),
                                        Locus(5, 6),
                                        0,
                                    ),
                                ],
                                0,
                                Order.read2ref,
                            ),
                            AlignmentChain(
                                10,
                                Contig(1, 10),
                                Contig(2, 100),
                                yes,
                                [
                                    LocalAlignment(
                                        Locus(0, 1),
                                        Locus(94, 95),
                                        0,
                                    ),
                                    LocalAlignment(
                                        Locus(4, 5),
                                        Locus(97, 98),
                                        0,
                                    ),
                                ],
                                0,
                                Order.read2ref,
                            ),
                        ],
                        Yes.isValid,
                        ReadAlignmentType.gap,
                        No.isExtension,
                        No.isFrontExtension,
                        No.isBackExtension,
                        Yes.isGap,
                    ),
                ];
                // dfmt on

                alias getFailureMessage = (testCase, testFunction, expectedValue) => format!"expected %s(%s) to be %s"(
                        testFunction, testCase, expectedValue);

                foreach (testCase, testData; testCases)
                {
                    auto readAlignment = testData[0];

                    assert(testData[1] == isValid(readAlignment),
                            getFailureMessage(testCase, "isValid", testData[1]));
                    if (isValid(readAlignment))
                        assert(testData[2] == getType(readAlignment),
                                getFailureMessage(testCase, "getType", testData[2]));
                    else
                        assertThrown!AssertError(getType(readAlignment),
                                format!"expected getType(%s) to throw"(testCase));
                    assert(testData[3] == isExtension(readAlignment),
                            getFailureMessage(testCase, "isExtension", testData[3]));
                    assert(testData[4] == isFrontExtension(readAlignment),
                            getFailureMessage(testCase, "isFrontExtension", testData[4]));
                    assert(testData[5] == isBackExtension(readAlignment),
                            getFailureMessage(testCase, "isBackExtension", testData[5]));
                    assert(testData[6] == isGap(readAlignment),
                            getFailureMessage(testCase, "isGap", testData[6]));
                }
            }
}

/**
    Returns the (approximate) size of the insertion produced by this
    read alignment.

    See_Also: `getGapSize`, `getExtensionSize`
*/
long getInsertionSize(in ReadAlignment readAlignment) pure
{
    final switch (readAlignment.getType)
    {
    case ReadAlignmentType.front:
    case ReadAlignmentType.back:
        return readAlignment.getExtensionSize();
    case ReadAlignmentType.gap:
        return readAlignment.getGapSize();
    }
}

unittest
{
    with (AlignmentChain) with (LocalAlignment) with (Complement)
            {
                // dfmt off
                auto extensionAlignment = [
                    AlignmentChain(0, Contig(1, 50), Contig(1, 50), no, [
                        LocalAlignment(Locus(0, 15), Locus(10, 20), 0),
                        LocalAlignment(Locus(20, 40), Locus(30, 50), 1),
                    ], 0, Order.read2ref)
                ];
                auto gapAlignment = [
                    AlignmentChain(0, Contig(1, 80), Contig(1, 50), no, [
                        LocalAlignment(Locus(10, 15), Locus(30, 35), 0),
                        LocalAlignment(Locus(20, 30), Locus(40, 50), 1),
                    ], 0, Order.read2ref),
                    AlignmentChain(1, Contig(1, 80), Contig(2, 50), no, [
                        LocalAlignment(Locus(50, 55), Locus(0, 10), 2),
                        LocalAlignment(Locus(60, 70), Locus(15, 20), 3),
                    ], 0, Order.read2ref),
                ];
                // dfmt on

                assert(extensionAlignment.getInsertionSize == extensionAlignment.getExtensionSize);
                assert(gapAlignment.getInsertionSize == gapAlignment.getGapSize);
            }
}

/**
    Returns the (approximate) size of the extension constituted by alignmentsRange.
    ---
    extensionSize ~= readsSequenceLength - referenceExcess

    CASE 1 (back extension, no complement):

                   0       rx  ry  lr
        reference  |-------+---+---|
                            \ \ \
                             \ \ \
                   alignment  \ \ \
                               \ \ \
        read            |-->-->-+->-+->-->--|
                        0       ax  ay      la

        readsSequenceLength ~= la - ay
        referenceExcess ~= lr - ry

    CASE 2 (front extension, no complement):

                            0    rx  ry lr
        reference           |---+---+-------|
                               / / /
                              / / /
                  alignment  / / /
                            / / /
        read       |-->-->-+->-+->-->--|
                   0       ax  ay      la

        readsSequenceLength ~= ax - 0 = ax
        referenceExcess ~= rx - 0 = rx

    CASE 3 (back extension, complement):

                   0       rx  ry  lr
        reference  |-------+---+---|
                            \ \ \
                             \ \ \
                   alignment  \ \ \
                               \ \ \
        read            |--<--<-+-<-+-<--<--|
                        0       ax  ay      la

        readsSequenceLength ~= la - ay
        referenceExcess ~= lr - ry

    CASE 4 (front extension, complement):

                            0    rx  ry lr
        reference           |---+---+-------|
                               / / /
                              / / /
                  alignment  / / /
                            / / /
        read       |--<--<-+-<-+-<--<--|
                   0       ax  ay      la

        readsSequenceLength ~= ax - 0 = ax
        referenceExcess ~= rx - 0 = rx
    ---
*/
private long getExtensionSize(in ReadAlignment readAlignment) pure
{
    assert(readAlignment.isExtension);

    auto alignment = readAlignment[0];
    // CASE 1/3 (back extension)
    long readsSequenceLengthRightExtension = alignment.contigA.length - alignment.last.contigA.end;
    long referenceExcessRightExtension = alignment.contigB.length - alignment.last.contigB.end;
    // CASE 2/4 (front extension)
    long readsSequenceLengthLeftExtension = alignment.first.contigA.begin;
    long referenceExcessLeftExtension = alignment.first.contigB.begin;

    // dfmt off
    return max(
        // Case 1/3 (back extension)
        readsSequenceLengthRightExtension - referenceExcessRightExtension,
        // Case 2/4 (front extension)
        readsSequenceLengthLeftExtension - referenceExcessLeftExtension,
    );
    // dfmt on
}

unittest
{
    with (AlignmentChain) with (LocalAlignment) with (Complement)
            {
                // dfmt off
                auto testCases = [
                    // CASE 1 (right extension, no complement)
                    // ------------------------------------
                    //            0      10  50  50
                    // reference  |-------+---+---|
                    //                     \ \ \
                    // read         |-->-->-+->-+->-->--|
                    //              0       0  40      50
                    tuple(
                        [AlignmentChain(0, Contig(1, 50), Contig(1, 50), no, [
                            LocalAlignment(Locus(0, 15), Locus(10, 20), 0),
                            LocalAlignment(Locus(20, 40), Locus(30, 50), 1),
                        ], 0, Order.read2ref)],
                        10
                    ),
                    //            0      10  45  50
                    // reference  |-------+---+---|
                    //                     \ \ \
                    // read         |-->-->-+->-+->-->--|
                    //              0       0  40      50
                    tuple(
                        [AlignmentChain(1, Contig(1, 50), Contig(1, 50), no, [
                            LocalAlignment(Locus(0, 15), Locus(10, 20), 2),
                            LocalAlignment(Locus(20, 40), Locus(30, 45), 3),
                        ], 0, Order.read2ref)],
                        5
                    ),
                    // CASE 2 (left extension, no complement):
                    // ------------------------------------
                    //                  0   0  40      50
                    // reference        |---+---+-------|
                    //                     / / /
                    // read       |-->-->-+->-+->-->--|
                    //            0      10  50      50
                    tuple(
                        [AlignmentChain(2, Contig(1, 50), Contig(1, 50), no, [
                            LocalAlignment(Locus(10, 15), Locus(0, 20), 4),
                            LocalAlignment(Locus(20, 50), Locus(30, 40), 5),
                        ], 0, Order.read2ref)],
                        10
                    ),
                    //                  0   5  40      50
                    // reference        |---+---+-------|
                    //                     / / /
                    // read       |-->-->-+->-+->-->--|
                    //            0      10  50      50
                    tuple(
                        [AlignmentChain(3, Contig(1, 50), Contig(1, 50), no, [
                            LocalAlignment(Locus(10, 15), Locus(5, 20), 6),
                            LocalAlignment(Locus(20, 50), Locus(30, 40), 7),
                        ], 0, Order.read2ref)],
                        5
                    ),
                    // CASE 3 (right extension, complement):
                    // ------------------------------------
                    //            0      10  50  50
                    // reference  |-------+---+---|
                    //                     \ \ \
                    // read         |--<--<-+-<-+-<--<--|
                    //              0       0   40      50
                    tuple(
                        [AlignmentChain(4, Contig(1, 50), Contig(1, 50), yes, [
                            LocalAlignment(Locus(0, 15), Locus(10, 20), 8),
                            LocalAlignment(Locus(20, 40), Locus(30, 50), 9),
                        ], 0, Order.read2ref)],
                        10
                    ),
                    //            0      10  45  50
                    // reference  |-------+---+---|
                    //                     \ \ \
                    // read         |--<--<-+-<-+-<--<--|
                    //              0       0   50      50
                    tuple(
                        [AlignmentChain(5, Contig(1, 50), Contig(1, 50), yes, [
                            LocalAlignment(Locus(0, 15), Locus(10, 20), 10),
                            LocalAlignment(Locus(20, 40), Locus(30, 45), 11),
                        ], 0, Order.read2ref)],
                        5
                    ),
                    // CASE 4 (left extension, complement):
                    // ------------------------------------
                    //                  0   0  40      50
                    // reference        |---+---+-------|
                    //                     / / /
                    // read       |--<--<-+-<-+-<--<--|
                    //            0       10  50      50
                    tuple(
                        [AlignmentChain(6, Contig(1, 50), Contig(1, 50), yes, [
                            LocalAlignment(Locus(10, 15), Locus(0, 20), 12),
                            LocalAlignment(Locus(20, 50), Locus(30, 40), 13),
                        ], 0, Order.read2ref)],
                        10
                    ),
                    //                  0   5  40      50
                    // reference        |---+---+-------|
                    //                     / / /
                    // read       |--<--<-+-<-+-<--<--|
                    //            0       10  50      50
                    tuple(
                        [AlignmentChain(7, Contig(1, 50), Contig(1, 50), yes, [
                            LocalAlignment(Locus(10, 15), Locus(5, 20), 14),
                            LocalAlignment(Locus(20, 50), Locus(30, 40), 15),
                        ], 0, Order.read2ref)],
                        5
                    ),
                ];
                // dfmt on

                foreach (testCaseIdx, testCase; testCases.enumerate)
                {
                    auto gotValue = getExtensionSize(testCase[0]);
                    auto expValue = testCase[1];
                    auto errorMessage = format!"expected extension size %d but got %d for test case %d"(
                            expValue, gotValue, testCaseIdx);

                    assert(gotValue == expValue, errorMessage);
                }
            }
}

/**
    Returns the (approximate) size of the gap spanned by this read alignment.
    ---
    gapSize = readsSequenceLength - referenceExcess

    readsSequenceLength ~= a2x - a1y

    referenceExcess ~= referenceExcess1 + referenceExcess2

    referenceExcess1 ~= la - cay
    referenceExcess2 ~= cbx - 0 = cbx

    CASE 1 (no complement):

                        contig a                  contig b

                  0         cax cay la      0   cbx cby       lb
        reference |---------+---+---|       |---+---+---------|
                             \ \ \             / / /
                              \ \ \           / / /
                  alignment 1  \ \ \         / / /  alignment 2
                                \ \ \       / / /
        read            |-->-->--+->-+-->--+->-+-->-->--|
                        0        a1x a1y   a2x a2y      lr

    CASE 2 (complement):

                        contig b                  contig a

                  0         cax cay la      0   cbx cby       lb
        reference |---------+---+---|       |---+---+---------|
                             \ \ \             / / /
                              \ \ \           / / /
                  alignment 2  \ \ \         / / /  alignment 1
                                \ \ \       / / /
        read            |--<--<--+-<-+--<--+-<-+--<--<--|
                        0        a1x a1y   a2x a2y      lr
    ---
*/
private long getGapSize(in ReadAlignment readAlignment) pure
{
    assert(readAlignment.isGap);
    // dfmt off
    auto alignments = isBefore!"contigA"(readAlignment[0], readAlignment[1])
        ? tuple(readAlignment[0], readAlignment[1])
        : tuple(readAlignment[1], readAlignment[0]);
    // dfmt on
    auto firstAlignment = alignments[0];
    auto secondAlignment = alignments[1];

    assert(secondAlignment.first.contigA.begin > firstAlignment.last.contigA.end,
            format!"intersecting local alignments in (%s, %s)"(firstAlignment, secondAlignment));

    long readsSequenceLength = secondAlignment.first.contigA.begin - firstAlignment
        .last.contigA.end;
    // dfmt off
    long referenceExcess1 = firstAlignment.contigB.length - firstAlignment.last.contigB.end;
    long referenceExcess2 = secondAlignment.first.contigB.begin;
    // dfmt on
    long referenceExcess = referenceExcess1 + referenceExcess2;

    return readsSequenceLength - referenceExcess;
}

unittest
{
    with (AlignmentChain) with (LocalAlignment) with (Complement)
            {
                // dfmt off
                auto testCases = [
                    // CASE 1 (no complement)
                    // ------------------------------------
                    //           0     30  50  50     0   0  20     50
                    // reference |------+---+---|     |---+---+------|
                    //                   \ \ \           / / /
                    //       alignment 1  \ \ \         / / /  alignment 2
                    //                     \ \ \       / / /
                    // read        |-->-->--+->-+-->--+->-+-->-->--|
                    //             0       10   30   50  70       80
                    tuple(
                        [
                            AlignmentChain(0, Contig(1, 80), Contig(1, 50), no, [
                                LocalAlignment(Locus(10, 15), Locus(30, 35), 0),
                                LocalAlignment(Locus(20, 30), Locus(40, 50), 1),
                            ], 0, Order.read2ref),
                            AlignmentChain(1, Contig(1, 80), Contig(2, 50), no, [
                                LocalAlignment(Locus(50, 55), Locus(0, 10), 2),
                                LocalAlignment(Locus(60, 70), Locus(15, 20), 3),
                            ], 0, Order.read2ref),
                        ],
                        20
                    ),
                    //           0     30  45  50     0   5  20     50
                    // reference |------+---+---|     |---+---+------|
                    //                   \ \ \           / / /
                    //       alignment 1  \ \ \         / / /  alignment 2
                    //                     \ \ \       / / /
                    // read        |-->-->--+->-+-->--+->-+-->-->--|
                    //             0       10   30   50  70       80
                    tuple(
                        [
                            AlignmentChain(2, Contig(1, 80), Contig(1, 50), no, [
                                LocalAlignment(Locus(10, 15), Locus(30, 35), 4),
                                LocalAlignment(Locus(20, 30), Locus(40, 45), 5),
                            ], 0, Order.read2ref),
                            AlignmentChain(3, Contig(1, 80), Contig(2, 50), no, [
                                LocalAlignment(Locus(50, 55), Locus(5, 10), 6),
                                LocalAlignment(Locus(60, 70), Locus(15, 20), 7),
                            ], 0, Order.read2ref),
                        ],
                        10
                    ),
                    //           0     30  45  50     0   5  20     50
                    // reference |------+---+---|     |---+---+------|
                    //                   \ \ \           / / /
                    //       alignment 2  \ \ \         / / /  alignment 1
                    //                     \ \ \       / / /
                    // read        |-->-->--+->-+-->--+->-+-->-->--|
                    //             0       10   30   50  70       80
                    tuple(
                        [
                            AlignmentChain(4, Contig(1, 80), Contig(2, 50), no, [
                                LocalAlignment(Locus(50, 55), Locus(5, 10), 6),
                                LocalAlignment(Locus(60, 70), Locus(15, 20), 7),
                            ], 0, Order.read2ref),
                            AlignmentChain(5, Contig(1, 80), Contig(1, 50), no, [
                                LocalAlignment(Locus(10, 15), Locus(30, 35), 4),
                                LocalAlignment(Locus(20, 30), Locus(40, 45), 5),
                            ], 0, Order.read2ref),
                        ],
                        10
                    ),
                    // CASE 2 (complement)
                    // ------------------------------------
                    //           0     30  50  50     0   0  20     50
                    // reference |------+---+---|     |---+---+------|
                    //                   \ \ \           / / /
                    //       alignment 1  \ \ \         / / /  alignment 2
                    //                     \ \ \       / / /
                    // read        |--<--<--+-<-+--<--+-<-+--<--<--|
                    //             0       10   30   50  70       80
                    tuple(
                        [
                            AlignmentChain(6, Contig(1, 80), Contig(1, 50), yes, [
                                LocalAlignment(Locus(10, 15), Locus(30, 35), 8),
                                LocalAlignment(Locus(20, 30), Locus(40, 50), 9),
                            ], 0, Order.read2ref),
                            AlignmentChain(7, Contig(1, 80), Contig(2, 50), yes, [
                                LocalAlignment(Locus(50, 55), Locus(0, 10), 10),
                                LocalAlignment(Locus(60, 70), Locus(15, 20), 11),
                            ], 0, Order.read2ref),
                        ],
                        20
                    ),
                    //           0     30  45  50     0   5  20     50
                    // reference |------+---+---|     |---+---+------|
                    //                   \ \ \           / / /
                    //       alignment 1  \ \ \         / / /  alignment 2
                    //                     \ \ \       / / /
                    // read        |--<--<--+-<-+--<--+-<-+--<--<--|
                    //             0       10   30   50  70       80
                    tuple(
                        [
                            AlignmentChain(8, Contig(1, 80), Contig(1, 50), yes, [
                                LocalAlignment(Locus(10, 15), Locus(30, 35), 12),
                                LocalAlignment(Locus(20, 30), Locus(40, 45), 13),
                            ], 0, Order.read2ref),
                            AlignmentChain(9, Contig(1, 80), Contig(2, 50), yes, [
                                LocalAlignment(Locus(50, 55), Locus(5, 10), 14),
                                LocalAlignment(Locus(60, 70), Locus(15, 20), 15),
                            ], 0, Order.read2ref),
                        ],
                        10
                    ),
                    //           0     30  45  50     0   5  20     50
                    // reference |------+---+---|     |---+---+------|
                    //                   \ \ \           / / /
                    //       alignment 2  \ \ \         / / /  alignment 1
                    //                     \ \ \       / / /
                    // read        |--<--<--+-<-+--<--+-<-+--<--<--|
                    //             0       10   30   50  70       80
                    tuple(
                        [
                            AlignmentChain(10, Contig(1, 80), Contig(2, 50), yes, [
                                LocalAlignment(Locus(50, 55), Locus(5, 10), 14),
                                LocalAlignment(Locus(60, 70), Locus(15, 20), 15),
                            ], 0, Order.read2ref),
                            AlignmentChain(11, Contig(1, 80), Contig(1, 50), yes, [
                                LocalAlignment(Locus(10, 15), Locus(30, 35), 12),
                                LocalAlignment(Locus(20, 30), Locus(40, 45), 13),
                            ], 0, Order.read2ref),
                        ],
                        10
                    ),
                ];
                // dfmt on

                foreach (testCaseIdx, testCase; testCases.enumerate)
                {
                    auto gotValue = getGapSize(testCase[0]);
                    auto expValue = testCase[1];
                    auto errorMessage = format!"expected spanning gap size %d but got %d for test case %d"(
                            expValue, gotValue, testCaseIdx);

                    assert(gotValue == expValue, errorMessage);
                }
            }
}

size_t meanScore(in ReadAlignment readAlignment) pure
{
    return readAlignment.map!"a.score".mean;
}

/**
    A pile of read alignments belonging to the same gap/contig end.

    See_Also: Hit
*/
alias PileUp = ReadAlignment[];

PileUp[] buildPileUps(AlignmentContainer!(AlignmentChain[]) candidates) pure
{
    static bool orderByReferenceIds(T)(T a, T b) pure
    {
        return a[0].contigB.id < b[0].contigB.id && (!(a.length > 1
                && a[0].contigB.id == b[0].contigB.id) || a[1].contigB.id < b[1].contigB.id);
    }

    static size_t numContigsInvolved(T)(T readsByNumFlanksByGap)
    {
        auto readsByGap = readsByNumFlanksByGap.front;
        auto readsByFirstGap = readsByGap.front;
        auto reads = readsByFirstGap.front;

        return reads.length;
    }

    static bool isSameRead(T)(T a, T b)
    {
        return a.contigA.id == b.contigA.id;
    }

    static ReadAlignment sortGapAlignment(ReadAlignment readAlignment) pure nothrow
    {
        if (readAlignment.isGap)
        {
            auto isBefore = readAlignment[0].isBefore!"contigA"(readAlignment[1]);
            auto isComplement = readAlignment[0].complement;

            if (isComplement == isBefore)
                swap(readAlignment[0], readAlignment[1]);
        }

        return readAlignment;
    }

    // dfmt off
    auto sortedReadAlignments = candidates
        .b2a
        .chunkBy!isSameRead
        .map!array
        .filter!isValid // TODO: this might discard "good" reads which span
                        //       a whole contig extending it on both ends.
        .map!sortGapAlignment
        .array;
    // dfmt on
    sortedReadAlignments.sort!byGapSort;
    PileUp[] pileUpsByGap = [];
    PileUp lastPileUp;
    ReadAlignment lastReadAlignment;

    foreach (readAlignment; sortedReadAlignments)
    {
        if (pileUpsByGap.length == 0)
        {
            pileUpsByGap ~= [readAlignment];
            continue;
        }

        lastPileUp = pileUpsByGap[$ - 1];
        lastReadAlignment = lastPileUp[$ - 1];

        if (belongToSamePile(lastReadAlignment, readAlignment))
        {
            lastPileUp ~= readAlignment;
            pileUpsByGap[$ - 1] = lastPileUp;
        }
        else
        {
            pileUpsByGap ~= [readAlignment];
        }
    }

    debug logJsonDebug("pileStructure", pileUpsByGap.map!"a.length".map!Json.array);

    return pileUpsByGap;
}

bool belongToSamePile(in ReadAlignment a, in ReadAlignment b) pure nothrow
{
    assert(isInReferenceOrder(a) && isInReferenceOrder(b),
            "gap alignment chains must be in reference order");

    // dfmt off
    debug logJsonDebug(
        "aIsBackExtension", a.isBackExtension,
        "aIsFrontExtension", a.isFrontExtension,
        "aIsGap", a.isGap,
        "aMapAContigBIDMapJsonArray", a.map!"a.contigB.id".map!Json.array,
        "bIsBackExtension", b.isBackExtension,
        "bIsFrontExtension", b.isFrontExtension,
        "bIsGap", b.isGap,
        "bMapAContigBIDMapJsonArray", b.map!"a.contigB.id".map!Json.array,
    );
    // dfmt on

    // dfmt off
    return any(only(
        a.isBackExtension && b.isBackExtension && a[0].contigB.id == b[0].contigB.id,
        a.isBackExtension && b.isGap && a[0].contigB.id == b[0].contigB.id,
        a.isGap && b.isGap && a[0].contigB.id == b[0].contigB.id && a[1].contigB.id == b[1].contigB.id,
        a.isGap && b.isFrontExtension && a[1].contigB.id == b[0].contigB.id,
        a.isFrontExtension && b.isFrontExtension && a[0].contigB.id == b[0].contigB.id,
    ));
    // dfmt on
}

bool byGapSort(in ReadAlignment a, in ReadAlignment b) pure nothrow
{
    assert(isInReferenceOrder(a) && isInReferenceOrder(b),
            "gap alignment chains must be in reference order");

    if (a.isExtension && b.isExtension)
    {
        return byGapSortForExtensions(a, b);
    }
    else if (a.isExtension && b.isGap)
    {
        return byGapSortForExtensionAndGap(a, b, Yes.inOrder);
    }
    else if (a.isGap && b.isExtension)
    {
        return byGapSortForExtensionAndGap(b, a, No.inOrder);
    }
    else
    {
        try
        {
            debug assert(a.isGap && b.isGap, a.to!string ~ " " ~ b.to!string);
        }
        catch (Exception)
        {
        }

        return byGapSortForGaps(a, b);
    }
}

private bool byGapSortForExtensions(in ReadAlignment a, in ReadAlignment b) pure nothrow
{
    if (a[0].contigB.id != b[0].contigB.id)
    {
        return a[0].contigB.id < b[0].contigB.id;
    }
    else if (a.getType != b.getType)
    {
        return a.getType < b.getType;
    }
    else
    {
        // same contig and same extension type => same pile
        return false;
    }
}

private bool byGapSortForGaps(in ReadAlignment a, in ReadAlignment b) pure nothrow
{
    if (a[0].contigB.id != b[0].contigB.id)
    {
        return a[0].contigB.id < b[0].contigB.id;
    }
    else if (a[1].contigB.id != b[1].contigB.id)
    {
        return a[1].contigB.id < b[1].contigB.id;
    }
    else
    {
        // same begin contig, same end contig => same pile
        return false;
    }
}

private bool byGapSortForExtensionAndGap(in ReadAlignment ext,
        in ReadAlignment gap, Flag!"inOrder" inOrder) pure nothrow
{
    alias cmp = (extValue, gapValue) => inOrder ? extValue < gapValue : gapValue < extValue;
    auto gapPart = ext.isFrontExtension ? gap[0] : gap[1];

    if (ext[0].contigB.id != gap[0].contigB.id && ext[0].contigB.id != gap[1].contigB.id)
    {
        return cmp(ext[0].contigB.id, min(gap[0].contigB.id, gap[1].contigB.id));
    }
    else if (ext[0].contigB.id != gap[0].contigB.id && ext[0].contigB.id == gap[1].contigB.id)
    {
        // extension contig is end contig of gap => same pile if front extension else in reverse order
        return !inOrder;
    }
    else if (ext[0].contigB.id == gap[0].contigB.id && ext[0].contigB.id != gap[1].contigB.id)
    {
        // extension contig is begin contig of gap => same pile if back extension else in order
        return inOrder;
    }
    else
    {
        assert(0, "illegal gap: reference contigs are the same");
    }
}

unittest
{
    //               c1                       c2                       c3
    //        0    5   10   15         0    5   10   15         0    5   10   15
    // ref:   |---->---->----|.........|---->---->----|.........|---->---->----|
    // reads: .    .    .    :    .    :    .    .    :    .    :    .    .    .
    //        . #1 |---->---->--| .    :    .    .    :    .    :    .    .    .
    //        . #2 |----<----<--| .    :    .    .    :    .    :    .    .    .
    //        . #3 |---->---->----|    :    .    .    :    .    :    .    .    .
    //        .    . #4 |---->----|    :    .    .    :    .    :    .    .    .
    //        .    . #5 |---->---->---->----|    .    :    .    :    .    .    .
    //        .    . #6 |----<----<----<----|    .    :    .    :    .    .    .
    //        .    .    #7 |->---->---->---->--| .    :    .    :    .    .    .
    //        .    .    .    : #8 |---->---->--| .    :    .    :    .    .    .
    //        .    .    .    : #9 |----<----<--| .    :    .    :    .    .    .
    //        .    .    .    :#10 |---->---->----|    :    .    :    .    .    .
    //        .    .    .    :    #11 |----->----|    :    .    :    .    .    .
    //        .    .    .    :    .    :    .    .    :    .    :    .    .    .
    //        .    .    .    :    .    :#12 |---->---->--| .    :    .    .    .
    //        .    .    .    :    .    :#13 |----<----<--| .    :    .    .    .
    //        .    .    .    :    .    :#14 |---->---->----|    :    .    .    .
    //        .    .    .    :    .    :    .#15 |---->----|    :    .    .    .
    //        .    .    .    :    .    :    .#16 |---->---->---->----|    .    .
    //        .    .    .    :    .    :    .#17 |----<----<----<----|    .    .
    //        .    .    .    :    .    :    .   #18 |->---->---->---->--| .    .
    //        .    .    .    :    .    :    .    .    :#19 |---->---->--| .    .
    //        .    .    .    :    .    :    .    .    :#20 |----<----<--| .    .
    //        .    .    .    :    .    :    .    .    :#21 |---->---->----|    .
    //        .    .    .    :    .    :    .    .    :    #22 |----->----|    .
    import std.algorithm : clamp;
    import std.random : randomShuffle;

    with (AlignmentChain) with (LocalAlignment)
        {
            size_t alignmentChainId = 0;
            size_t contReadId = 0;
            ReadAlignment getDummyRead(size_t beginContigId, long beginIdx,
                    size_t endContigId, long endIdx, Complement complement)
            {
                static immutable contigLength = 16;
                static immutable gapLength = 9;
                static immutable numDiffs = 0;

                alignmentChainId += 2;
                auto readId = ++contReadId;
                // dfmt off
                auto readLength = beginContigId == endContigId
                    ? endIdx - beginIdx
                    : contigLength - beginIdx + gapLength + endIdx;
                // dfmt on
                size_t firstReadBeginIdx;
                size_t firstReadEndIdx;

                if (beginIdx < 0)
                {
                    firstReadBeginIdx = readLength - endIdx;
                    firstReadEndIdx = readLength;
                }
                else
                {
                    firstReadBeginIdx = 0;
                    firstReadEndIdx = contigLength - beginIdx;
                }

                beginIdx = clamp(beginIdx, 0, contigLength);
                endIdx = clamp(endIdx, 0, contigLength);

                if (complement)
                {
                    beginIdx = contigLength - beginIdx;
                    endIdx = contigLength - endIdx;
                    swap(beginIdx, endIdx);

                    firstReadBeginIdx = readLength - firstReadBeginIdx;
                    firstReadEndIdx = readLength - firstReadEndIdx;
                    swap(firstReadBeginIdx, firstReadEndIdx);
                }

                if (beginContigId == endContigId)
                {
                    // dfmt off
                    return [AlignmentChain(
                        alignmentChainId - 2,
                        Contig(readId, readLength),
                        Contig(beginContigId, contigLength),
                        complement,
                        [
                            LocalAlignment(
                                Locus(firstReadBeginIdx, firstReadBeginIdx + 1),
                                Locus(beginIdx, beginIdx + 1),
                                numDiffs,
                            ),
                            LocalAlignment(
                                Locus(firstReadEndIdx - 1, firstReadEndIdx),
                                Locus(endIdx - 1, endIdx),
                                numDiffs,
                            ),
                        ],
                        0,
                        Order.read2ref,
                    )];
                    // dfmt on
                }
                else
                {
                    auto secondReadBeginIdx = readLength - endIdx;
                    auto secondReadEndIdx = readLength;

                    if (complement)
                    {
                        secondReadBeginIdx = readLength - secondReadBeginIdx;
                        secondReadEndIdx = readLength - secondReadEndIdx;
                        swap(secondReadBeginIdx, secondReadEndIdx);
                    }

                    // dfmt off
                    return [
                        AlignmentChain(
                            alignmentChainId - 2,
                            Contig(readId, readLength),
                            Contig(beginContigId, contigLength),
                            complement,
                            [
                                LocalAlignment(
                                    Locus(firstReadBeginIdx, firstReadBeginIdx + 1),
                                    Locus(beginIdx, beginIdx + 1),
                                    numDiffs,
                                ),
                                LocalAlignment(
                                    Locus(firstReadEndIdx - 1, firstReadEndIdx),
                                    Locus(contigLength - 1, contigLength),
                                    numDiffs,
                                ),
                            ],
                            0,
                            Order.read2ref,
                        ),
                        AlignmentChain(
                            alignmentChainId - 1,
                            Contig(readId, readLength),
                            Contig(endContigId, contigLength),
                            complement,
                            [
                                LocalAlignment(
                                    Locus(secondReadBeginIdx, secondReadBeginIdx + 1),
                                    Locus(0, 1),
                                    numDiffs,
                                ),
                                LocalAlignment(
                                    Locus(secondReadEndIdx - 1, secondReadEndIdx),
                                    Locus(endIdx - 1, endIdx),
                                    numDiffs,
                                ),
                            ],
                            0,
                            Order.read2ref,
                        ),
                    ];
                    // dfmt on
                }
            }

            size_t c1 = 1;
            size_t c2 = 2;
            size_t c3 = 3;
            // dfmt off
            auto totalOrderingGroups = [
                [
                    getDummyRead(c1,  5, c1, 18, Complement.no),  //  #1
                    getDummyRead(c1,  5, c1, 18, Complement.yes), //  #2
                    getDummyRead(c1,  5, c1, 20, Complement.no),  //  #3
                    getDummyRead(c1, 10, c1, 20, Complement.no),  //  #4
                ],
                [
                    getDummyRead(c1, 10, c2,  5, Complement.no),  //  #5
                    getDummyRead(c1, 10, c2,  5, Complement.yes), //  #6
                    getDummyRead(c1, 13, c2,  8, Complement.no),  //  #7
                ],
                [
                    getDummyRead(c2, -5, c2,  8, Complement.no),  //  #8
                    getDummyRead(c2, -5, c2,  8, Complement.yes), //  #9
                    getDummyRead(c2, -5, c2, 10, Complement.no),  // #10
                    getDummyRead(c2, -1, c2, 10, Complement.no),  // #11
                ],
                [
                    getDummyRead(c2,  5, c2, 18, Complement.no),  // #12
                    getDummyRead(c2,  5, c2, 18, Complement.yes), // #13
                    getDummyRead(c2,  5, c2, 20, Complement.no),  // #14
                    getDummyRead(c2, 10, c2, 20, Complement.no),  // #15
                ],
                [
                    getDummyRead(c2, 10, c3,  5, Complement.no),  // #16
                    getDummyRead(c2, 10, c3,  5, Complement.yes), // #17
                    getDummyRead(c2, 13, c3,  8, Complement.no),  // #18
                ],
                [
                    getDummyRead(c3, -5, c3,  8, Complement.no),  // #19
                    getDummyRead(c3, -5, c3,  8, Complement.yes), // #20
                    getDummyRead(c3, -5, c3, 10, Complement.no),  // #21
                    getDummyRead(c3, -1, c3, 10, Complement.no),  // #22
                ],
            ];
            // dfmt on

            foreach (i, group1; totalOrderingGroups.enumerate)
            {
                foreach (j, group2; totalOrderingGroups.enumerate)
                {
                    foreach (readAlignment1; group1)
                    {
                        foreach (readAlignment2; group2)
                        {
                            size_t read1Id = readAlignment1[0].contigA.id;
                            size_t read2Id = readAlignment2[0].contigA.id;

                            assert((i != j) || !byGapSort(readAlignment1, readAlignment2),
                                    format!"read %d should not be before %d (same pile)"(read1Id,
                                        read2Id));
                            assert((i < j) == byGapSort(readAlignment1, readAlignment2),
                                    format!"read %d should be before %d"(read1Id, read2Id));
                            assert((i > j) == byGapSort(readAlignment2, readAlignment1),
                                    format!"read %d should be after %d"(read2Id, read1Id));

                            bool shouldGroup = i == j || (i + 1 == j && j != 3);
                            assert(shouldGroup == belongToSamePile(readAlignment1, readAlignment2),
                                    format!"reads %d and %d should%s be grouped"(read1Id,
                                        read2Id, (shouldGroup ? "" : " not")));
                        }
                    }
                }
            }
        }
}

/**
    Get the type of the read alignment.

    See_Also: `isFrontExtension`, `isBackExtension`, `isGap`
*/
ReadAlignmentType getType(in PileUp pileUp) pure nothrow
{
    assert(pileUp.isValid, "invalid read alignment");

    if (pileUp.isGap)
    {
        return ReadAlignmentType.gap;
    }
    else
    {
        assert(pileUp.isExtension);
        return pileUp[0].isFrontExtension
            ? ReadAlignmentType.front
            : ReadAlignmentType.back;
    }
}

bool isValid(in PileUp pileUp) pure nothrow
{
    return pileUp.isExtension ^ pileUp.isGap;
}

bool isExtension(in PileUp pileUp) pure nothrow
{
    if (pileUp[0].isFrontExtension)
    {
        return pileUp.all!(readAlignment => readAlignment.isFrontExtension);
    }
    else if (pileUp[0].isBackExtension)
    {
        return pileUp.all!(readAlignment => readAlignment.isBackExtension);
    }
    else
    {
        return false;
    }
}

bool isGap(in PileUp pileUp) pure nothrow
{
    return pileUp.any!(readAlignment => readAlignment.isGap);
}

/// Returns a pile up with the complementary order alignment chains.
PileUp getComplementaryOrder(in PileUp pileUp, AlignmentChain[] complementaryAlignmentChains) pure
{
    assert(complementaryAlignmentChains.isSorted);
    auto sortedComplementaryAlignmentChains = assumeSorted(complementaryAlignmentChains);
    auto complementaryPileUp = appender!PileUp;
    complementaryPileUp.reserve(pileUp.length);
    Appender!ReadAlignment complementaryReadAlignment;

    foreach (readAlignment; pileUp)
    {
        complementaryReadAlignment = appender!ReadAlignment;
        complementaryReadAlignment.reserve(readAlignment.length);

        foreach (alignmentChain; readAlignment)
        {
            auto needle = alignmentChain.getComplementaryOrder();
            auto complementaryCandidates = sortedComplementaryAlignmentChains.equalRange(needle);
            if (complementaryCandidates.empty)
            {
                throw new Exception("missing complementary alignment chain");
            }

            complementaryReadAlignment ~= complementaryCandidates.front;
        }

        complementaryPileUp ~= complementaryReadAlignment.data;
    }

    return complementaryPileUp.data;
}

unittest
{
    alias Complement = AlignmentChain.Complement;
    alias Order = AlignmentChain.Order;
    AlignmentChain getDummyAC(size_t contigA, size_t contigB, AlignmentChain.Complement complement, size_t contigABegin, size_t contigAEnd, size_t contigBBegin, size_t contigBEnd, AlignmentChain.Order order)
    {
        immutable contigLength = 100;
        // dfmt off
        AlignmentChain dummy = {
            contigA: AlignmentChain.Contig(contigA, contigLength),
            contigB: AlignmentChain.Contig(contigB, contigLength),
            complement: complement,
            localAlignments: [AlignmentChain.LocalAlignment(
                AlignmentChain.LocalAlignment.Locus(contigABegin, contigAEnd),
                AlignmentChain.LocalAlignment.Locus(contigBBegin, contigBEnd),
            )],
            0,
            order,
        };
        // dfmt on

        return dummy;
    }

    // dfmt off
    auto inPileUp = [
        [  // read 1 spans gap from 1 to 2
            getDummyAC(1, 1, Complement.no,   0,  10,  90, 100, Order.read2ref),
            getDummyAC(1, 2, Complement.no,  90, 100,   0,  10, Order.read2ref),
        ],
        [  // read 2 extends begin contig 3
            getDummyAC(2, 3, Complement.no,  90, 100,   0,  10, Order.read2ref),
        ],
        [  // read 3 extends end contig 3
            getDummyAC(3, 3, Complement.no,   0,  10,  90, 100, Order.read2ref),
        ],
        [  // read 4 spans gap from 4 to 5 (complement)
            getDummyAC(4, 4, Complement.yes,   0,  10,  90, 100, Order.read2ref),
            getDummyAC(4, 5, Complement.yes,  90, 100,   0,  10, Order.read2ref),
        ],
        [  // read 5 extends begin contig 6 (complement)
            getDummyAC(5, 6, Complement.yes,  90, 100,   0,  10, Order.read2ref),
        ],
        [  // read 6 extends end contig 6 (complement)
            getDummyAC(6, 6, Complement.yes,   0,  10,  90, 100, Order.read2ref),
        ],
    ];
    auto complementaryAlignmentChains = [
        getDummyAC(1, 1, Complement.no,  90, 100,   0,  10, Order.ref2read),
        getDummyAC(2, 1, Complement.no,   0,  10,  90, 100, Order.ref2read),
        getDummyAC(3, 2, Complement.no,   0,  10,  90, 100, Order.ref2read),
        getDummyAC(3, 3, Complement.no,  90, 100,   0,  10, Order.ref2read),
        getDummyAC(4, 4, Complement.yes,   0,  10,  90, 100, Order.ref2read),
        getDummyAC(5, 4, Complement.yes,  90, 100,   0,  10, Order.ref2read),
        getDummyAC(6, 5, Complement.yes,  90, 100,   0,  10, Order.ref2read),
        getDummyAC(6, 6, Complement.yes,   0,  10,  90, 100, Order.ref2read),
    ];
    // dfmt on

    auto complementaryPileUp = inPileUp.getComplementaryOrder(complementaryAlignmentChains);

    // dfmt off
    assert(complementaryPileUp == [
        [  // read 1 spans gap from 1 to 2
            getDummyAC(1, 1, Complement.no,  90, 100,   0,  10, Order.ref2read),
            getDummyAC(2, 1, Complement.no,   0,  10,  90, 100, Order.ref2read),
        ],
        [  // read 2 extends begin contig 3
            getDummyAC(3, 2, Complement.no,   0,  10,  90, 100, Order.ref2read),
        ],
        [  // read 3 extends end contig 3
            getDummyAC(3, 3, Complement.no,  90, 100,   0,  10, Order.ref2read),
        ],
        [  // read 4 spans gap from 4 to 5 (complement)
            getDummyAC(4, 4, Complement.yes,   0,  10,  90, 100, Order.ref2read),
            getDummyAC(5, 4, Complement.yes,  90, 100,   0,  10, Order.ref2read),
        ],
        [  // read 5 extends begin contig 6 (complement)
            getDummyAC(6, 5, Complement.yes,  90, 100,   0,  10, Order.ref2read),
        ],
        [  // read 6 extends end contig 6 (complement)
            getDummyAC(6, 6, Complement.yes,   0,  10,  90, 100, Order.ref2read),
        ],
    ]);
    // dfmt on
}

AlignmentChain getComplementaryOrder(in AlignmentChain alignmentChain) pure
{
    AlignmentChain complementary = {
        contigA: alignmentChain.contigB,
        contigB: alignmentChain.contigA,
        complement: alignmentChain.complement,
        localAlignments: [AlignmentChain.LocalAlignment(
            AlignmentChain.LocalAlignment.Locus(
                alignmentChain.first.contigB.begin,
                alignmentChain.last.contigB.end,
            ),
            AlignmentChain.LocalAlignment.Locus(
                alignmentChain.first.contigA.begin,
                alignmentChain.last.contigA.end,
            ),
        )],
        order: alignmentChain.order.reverse
    };
    // dfmt on

    if (complementary.complement)
    {
        static foreach (contig; ["contigA", "contigB"])
        {
            static foreach (coord; ["begin", "end"])
                mixin("complementary.localAlignments[0]."~contig~"."~coord~" = complementary."~contig~".length - complementary.localAlignments[0]."~contig~"."~coord~";");
            mixin("swap(complementary.localAlignments[0]."~contig~".begin, complementary.localAlignments[0]."~contig~".end);");
        }

    }

    return complementary;
}


/// Returns a list of pointers to all involved alignment chains.
AlignmentChain*[] getAlignmentChainRefs(PileUp pileUp) pure nothrow
{
    auto alignmentChainsAcc = appender!(AlignmentChain*[]);
    alignmentChainsAcc.reserve(pileUp.map!"a.length".length);

    foreach (ref readAlignment; pileUp)
    {
        foreach (ref alignmentChain; readAlignment)
        {
            alignmentChainsAcc ~= &alignmentChain;
        }
    }

    return alignmentChainsAcc.data;
}

///
unittest
{
    auto pileUp = [[AlignmentChain(), AlignmentChain()], [AlignmentChain()]];
    auto allAlignmentChains = pileUp.getAlignmentChainRefs();

    assert(allAlignmentChains.length == 3);
    assert(pileUp[0][0].id == 0);
    assert(pileUp[0][1].id == 0);
    assert(pileUp[1][0].id == 0);

    allAlignmentChains[0].id = 1;
    allAlignmentChains[1].id = 2;
    allAlignmentChains[2].id = 3;

    assert(pileUp[0][0].id == 1);
    assert(pileUp[0][1].id == 2);
    assert(pileUp[1][0].id == 3);
}

// dfmt off
alias CroppingSlice = Tuple!(
    size_t, "begin",
    size_t, "end",
);
alias CroppingRefPosition = Tuple!(
    size_t, "frontIdx",
    size_t, "frontContigId",
    size_t, "backIdx",
    size_t, "backContigId",
);
// dfmt on

CroppingSlice intersection(in CroppingSlice aSlice, in CroppingSlice bSlice) pure nothrow
{
    return CroppingSlice(max(aSlice[0], bSlice[0]), min(aSlice[1], bSlice[1]));
}

CroppingRefPosition getCroppingRefPositions(const PileUp pileUp, in size_t tracePointDistance) pure nothrow
{
    // dfmt off
    auto frontAlignments = pileUp
        .filter!(readAlignment => readAlignment.isFrontExtension || readAlignment.isGap)
        .map!"a[$ - 1]"
        .array;
    auto backAlignments = pileUp
        .filter!(readAlignment => readAlignment.isBackExtension || readAlignment.isGap)
        .map!"a[0]"
        .array;
    // dfmt on
    auto commonFrontTracePoint = frontAlignments.getCommonFrontTracePoint(tracePointDistance);
    auto commonBackTracePoint = backAlignments.getCommonBackTracePoint(tracePointDistance);
    size_t frontContig;
    size_t backContig;

    if (frontAlignments.length > 0 && backAlignments.length > 0)
    {
        frontContig = frontAlignments[0].contigB.id;
        backContig = backAlignments[0].contigB.id;
    }
    else if (frontAlignments.length > 0)
    {
        frontContig = backContig = frontAlignments[0].contigB.id;
    }
    else if (backAlignments.length > 0)
    {
        backContig = frontContig = backAlignments[0].contigB.id;
    }
    else
    {
        assert(0, "empty pile up");
    }

    // dfmt off
    return CroppingRefPosition(
        commonFrontTracePoint,
        frontContig,
        commonBackTracePoint,
        backContig,
    );
    // dfmt on
}

CroppingSlice getCroppingSlice(const AlignmentChain alignmentChain, in CroppingRefPosition croppingRefPosition) pure
{
    auto tracePointDistance = alignmentChain.tracePointDistance;
    auto alignmentType = [alignmentChain.getComplementaryOrder].getType;
    // dfmt off
    auto refPos = alignmentType == ReadAlignmentType.front
        ? croppingRefPosition.frontIdx
        : croppingRefPosition.backIdx;
    // dfmt on
    auto openInterval = alignmentType == ReadAlignmentType.front;
    size_t readCroppingPos;

    foreach (localAlignment; alignmentChain.localAlignments)
    {
        auto firstTracePointRefPos = ceil(localAlignment.contigA.begin, tracePointDistance);
        auto numTracePoints = localAlignment.tracePoints.length;
        assert(refPos >= firstTracePointRefPos);
        auto tracePointIdx = (refPos - firstTracePointRefPos) / tracePointDistance;

        if (tracePointIdx < numTracePoints)
        {
            auto endIdx = tracePointIdx + (openInterval ? 0 : 1);
            // dfmt off
            readCroppingPos =
                localAlignment.contigB.begin +
                localAlignment
                    .tracePoints[0 .. endIdx]
                    .map!"a.numBasePairs"
                    .sum;
            // dfmt on
            break;
        }
    }

    // dfmt off
    return alignmentType == ReadAlignmentType.front
        ? CroppingSlice(0, readCroppingPos)
        : CroppingSlice(readCroppingPos, alignmentChain.contigB.length);
    // dfmt on
}

/**
    Returns a back "common" trace points wrt. contigB. The returned
    trace point is not necessarily common to *all* alignment chains, ie. the
    alignment intervals are incrementally intersected reject those that would
    cause an empty intersection. The intervals are ordered descending by
    end point and starting point. Example:

    ---
    :....:....:....:....:....:....:....:....:
    :    :    :    :    |----+--------------|
    :    :    :    |---------+--------------|
    :    :    :    :    :  |-+---------|    :
    :    :    :    :    :   |+--------|:    :
    :    :    :    :    :|---+-------| :    :
    :    : |-----------|:    :    :    :    :
    |---------|    :    :    :    :    :    :

     :  = trace point
     +  = "common" trace point
    |-| = alignment interval
    ---
*/
size_t getCommonBackTracePoint(const AlignmentChain[] alignmentChains, in size_t tracePointDistance) pure nothrow
{
    // dfmt off
    auto alignmentSlices = alignmentChains
        .map!(ac => ac.complement
            ? tuple(
                ac.contigB.length - ac.first.contigB.begin + 0,
                ac.contigB.length - ac.last.contigB.end + 0,
            )
            : tuple(
                ac.last.contigB.end + 0,
                ac.first.contigB.begin + 0,
            ))
        .array;
    // dfmt on

    foreach (alignmentSlice; alignmentSlices)
    {
        assert(alignmentSlice[0] > alignmentSlice[1]);
    }
    // dfmt off
    debug logJsonDebug(
        "alignmentSlices", alignmentSlices.map!"[a[0], a[1]]".array.toJson,
        "contigBIds", alignmentChains.map!"a.contigB.id".array.toJson,
        "type", "back",
    );
    // dfmt on

    return getCommonBackTracePoint(alignmentSlices, tracePointDistance);
}

private size_t getCommonBackTracePoint(Tuple!(size_t, size_t)[] alignmentSlices, in size_t tracePointDistance) pure nothrow
{
    if (alignmentSlices.length == 0)
        return 0;

    alignmentSlices.sort!"a > b";
    auto commonSlice = tuple(0UL, size_t.max);
    auto lastCommonSlice = commonSlice;

    foreach (alignmentSlice; alignmentSlices)
    {
        commonSlice[0] = ceil(max(commonSlice[0], alignmentSlice[1]), tracePointDistance);
        commonSlice[1] = floor(min(commonSlice[1], alignmentSlice[0]), tracePointDistance);

        if (commonSlice[1] < commonSlice[0])
        {
            break;
        }

        lastCommonSlice = commonSlice;
    }

    return lastCommonSlice[0];
}

unittest
{
    immutable tracePointDistance = 5UL;
    // dfmt off
    // 0    5   10   15   20   25
    // :....:....:....:....:....:
    // |---------+--------------|
    //      |----+--------------|
    //         |-+---------|
    //       |---+-------|
    //          |+--------|
    auto alignmentSlices1 = [
        tuple(25UL,  0UL),
        tuple(25UL,  5UL),
        tuple(20UL,  8UL),
        tuple(18UL,  6UL),
        tuple(19UL,  9UL),
    ];
    // dfmt on
    auto commonBackTracePoint1 = getCommonBackTracePoint(alignmentSlices1, tracePointDistance);

    assert(commonBackTracePoint1 == 10);

    // dfmt off
    // 0    5   10   15   20   25
    // :....:....:....:....:....:
    // |-------------|
    //               |+----|
    //            |---+-------|
    //              |-+---------|
    auto alignmentSlices2 = [
        tuple(14UL,  0UL),
        tuple(20UL, 14UL),
        tuple(23UL, 11UL),
        tuple(25UL, 13UL),
    ];
    // dfmt on
    auto commonBackTracePoint2 = getCommonBackTracePoint(alignmentSlices2, tracePointDistance);

    assert(commonBackTracePoint2 == 15);

    Tuple!(size_t, size_t)[] alignmentSlices3;
    auto commonBackTracePoint3 = getCommonBackTracePoint(alignmentSlices3, tracePointDistance);

    assert(commonBackTracePoint3 == 0);
}

/**
    Returns a back "common" trace points wrt. contigB. The returned
    trace point is not necessarily common to *all* alignment chains, ie. the
    alignment intervals are incrementally intersected reject those that would
    cause an empty intersection. The intervals are ordered ascending by
    starting point and end point. Example:

    ---
    :....:....:....:....:....:....:....:....:
    |--------------+----|    :    :    :    :
    |--------------+---------|    :    :    :
    :    |---------+-|  :    :    :    :    :
    :    :|--------+|   :    :    :    :    :
    :    : |-------+---|:    :    :    :    :
    :    :    :    :    :|-----------| :    :
    :    :    :    :    :    :    |---------|

     :  = trace point
     +  = "common" trace point
    |-| = alignment interval
    ---
*/
size_t getCommonFrontTracePoint(const AlignmentChain[] alignmentChains, in size_t tracePointDistance) pure nothrow
{
    // dfmt off
    auto alignmentSlices = alignmentChains
        .map!(ac => ac.complement
            ? tuple(
                ac.contigB.length - ac.last.contigB.end + 0,
                ac.contigB.length - ac.first.contigB.begin + 0,
            )
            : tuple(
                ac.first.contigB.begin + 0,
                ac.last.contigB.end + 0,
            ))
        .array;
    // dfmt on

    foreach (alignmentSlice; alignmentSlices)
    {
        assert(alignmentSlice[0] < alignmentSlice[1]);
    }
    // dfmt off
    debug logJsonDebug(
        "alignmentSlices", alignmentSlices.map!"[a[0], a[1]]".array.toJson,
        "contigBIds", alignmentChains.map!"a.contigB.id".array.toJson,
        "type", "front",
    );
    // dfmt on

    return getCommonFrontTracePoint(alignmentSlices, tracePointDistance);
}

private size_t getCommonFrontTracePoint(Tuple!(size_t, size_t)[] alignmentSlices, in size_t tracePointDistance) pure nothrow
{
    if (alignmentSlices.length == 0)
        return size_t.max;

    alignmentSlices.sort!"a < b";
    auto commonSlice = tuple(0UL, size_t.max);
    auto lastCommonSlice = commonSlice;

    foreach (alignmentSlice; alignmentSlices)
    {
        commonSlice[0] = max(commonSlice[0], alignmentSlice[0]);
        commonSlice[0] = ceil(commonSlice[0], tracePointDistance);
        commonSlice[1] = min(commonSlice[1], alignmentSlice[1]);
        commonSlice[1] = floor(commonSlice[1], tracePointDistance);

        if (commonSlice[1] < commonSlice[0])
        {
            break;
        }

        lastCommonSlice = commonSlice;
    }

    return lastCommonSlice[1];
}

///
unittest
{
    immutable tracePointDistance = 5UL;
    // dfmt off
    // 0    5   10   15   20   25
    // :....:....:....:....:....:
    // |--------------+---------|
    // |--------------+----|
    //      |---------+-|
    //        |-------+---|
    //       |--------+|
    auto alignmentSlices1 = [
        tuple( 0UL, 25UL),
        tuple( 0UL, 20UL),
        tuple( 5UL, 17UL),
        tuple( 7UL, 19UL),
        tuple( 6UL, 16UL),
    ];
    // dfmt on
    auto commonFrontTracePoint1 = getCommonFrontTracePoint(alignmentSlices1, tracePointDistance);

    assert(commonFrontTracePoint1 == 15);

    // dfmt off
    // 0    5   10   15   20   25
    // :....:....:....:....:....:
    //            |-------------|
    //      |----+|
    //   |-------+---|
    // |---------+-|
    auto alignmentSlices2 = [
        tuple(11UL, 25UL),
        tuple( 5UL, 11UL),
        tuple( 2UL, 14UL),
        tuple( 0UL, 12UL),
    ];
    // dfmt on
    auto commonFrontTracePoint2 = getCommonFrontTracePoint(alignmentSlices2, tracePointDistance);

    assert(commonFrontTracePoint2 == 10);

    Tuple!(size_t, size_t)[] alignmentSlices3;
    auto commonFrontTracePoint3 = getCommonFrontTracePoint(alignmentSlices3, tracePointDistance);

    assert(commonFrontTracePoint3 == size_t.max);
}

/// This characterizes an insertion.
// dfmt off
alias Hit = Tuple!(
    CoordinateTransform.Insertion, "insertion",
    size_t, "readId",
    AlignmentChain.Complement, "complement",
    string, "dbFile",
);
// dfmt on

class DJunctor
{
    /// Stop after maxLoops `mainLoop`s at the latest.
    static immutable maxLoops = 1;
    /**
        Two scores are considered similar if the relative "error" of them is
        smaller than defaulMaxRelativeDiff.

        The relative error is defined as:

            relError(a, b) = abs(a - b)/max(a - b)

        **Implementation note:** all computations are done in integer
        arithmetic hence AlignmentChain.maxScore corresponds to 1 in the above
        equation.
    */
    static immutable maxRelativeDiff = AlignmentChain.maxScore / 20; // 5% of larger value
    /**
        Two scores are considered similar if the absolute "error" of them is
        smaller than defaulMaxAbsoluteDiff.

        The relative error is defined as:

            absError(a, b) = abs(a - b)

        **Implementation note:** all computations are done in integer
        arithmetic hence all scores are lower than or equal to
        AlignmentChain.maxScore .
    */
    static immutable maxAbsoluteDiff = AlignmentChain.maxScore / 100; // 1% wrt. score
    alias numEstimateUseless = (numCatUnsure, iteration) => (numCatUnsure - numCatUnsure / 20) / (
            iteration + 1);

    AlignmentContainer!(AlignmentChain[]) selfAlignment;
    AlignmentContainer!(AlignmentChain[]) readsAlignment;
    const Options options;
    /// Set of read ids not to be considered in further processing.
    size_t[] catUseless;
    /// Set of alignments to be considered in further processing.
    AlignmentContainer!(AlignmentChain[]) catCandidates;
    /// Set of alignments to be used for gap filling.
    Hit[] catHits;
    /// Access contigs distributed over several DBs by ID.
    DBUnion referenceDb;
    /// Transform coordinates before insertions to coordinates after them.
    CoordinateTransform coordTransform;
    alias T = coordTransform;
    size_t iteration;

    this(in ref Options options)
    {
        this.options = options;
        this.catUseless = [];
        this.catHits = [];
        this.iteration = 0;
        this.referenceDb.baseDb = options.refDb;
    }

    DJunctor run()
    {
        logJsonDiagnostic("state", "enter", "function", "run");
        // dfmt off
        this
            .init
            .mainLoop
            .finish;
        // dfmt on
        logJsonDiagnostic("state", "exit", "function", "run");

        return this;
    }

    protected DJunctor init()
    {
        logJsonDiagnostic("state", "enter", "function", "djunctor.init");
        selfAlignment = getLocalAlignments(options.refDb, options);
        readsAlignment = getMappings(options.refDb, options.readsDb, options);
        annotateOrder(readsAlignment, AlignmentChain.Order.ref2read);
        catCandidates = AlignmentContainer!(AlignmentChain[])(readsAlignment.a2b.dup,
                readsAlignment.b2a.dup);
        logJsonDiagnostic("state", "exit", "function", "djunctor.init");

        return this;
    }

    unittest
    {
        auto options = Options();
        auto djunctor = new DJunctor(options);

        djunctor.init();

        assert(getLocalAlignments.wasCalled);
        assert(djunctor.selfAlignment == getLocalAlignments.returnValue);
        getLocalAlignments.reset();
        assert(getMappings.wasCalled);
        assert(djunctor.readsAlignment == getMappings.returnValue);
        getMappings.reset();
    }

    protected DJunctor mainLoop()
    {
        logJsonDiagnostic("state", "enter", "function", "djunctor.mainLoop");
        do
        {
            filterUseless();
            findHits();

            // dfmt off
            logJsonDiagnostic(
                "iteration", iteration,
                "numUseless", catUseless.length,
                "numCandiatesA2b", catCandidates.a2b.length,
                "numCandiatesB2a", catCandidates.b2a.length,
                "numHits", catHits.length
            );
            // dfmt on

            if (catHits.length > 0)
            {
                insertHits();
            }

            ++iteration;
        }
        while (catHits.length > 0 && iteration < maxLoops);
        logJsonDiagnostic("state", "exit", "function", "djunctor.mainLoop");

        return this;
    }

    protected DJunctor filterUseless()
    {
        logJsonDiagnostic("state", "enter", "function", "djunctor.filterUseless");

        auto sizeReserve = numEstimateUseless(catCandidates.a2b.length, iteration);
        if (sizeReserve == 0)
        {
            // dfmt off
            logJsonDiagnostic(
                "state", "exit",
                "function", "djunctor.filterUseless",
                "reason", "skipping: expecting zero useless reads",
            );
            // dfmt on

            // Skip this step if we do not expect to find "useless" reads;
            // it is not harmful to consider some "useless" reads as
            // candidates (catCandidates).
            return this;
        }

        auto uselessAcc = appender!(size_t[]);
        auto candidatesAcc = appender!(AlignmentChain[]);
        uselessAcc.reserve(sizeReserve);
        alias isNotUseless = ac => !uselessAcc.data.canFind(ac.contigA.id);

        foreach (alignmentsChunk; catCandidates.a2b.chunkBy!haveEqualIds)
        {
            auto alignments = alignmentsChunk.array;
            alignments.sort!"a.score > b.score";

            // Mark reads as useless if either: (s. issue #3)
            // (a) The aligned read, ie. the alignment extended with the exceeding
            //     read sequence on either end, is fully contained in a single
            //     contig.
            // (b) It aligns in multiple locations of one contig with
            //     approximately equal quality.
            if (alignments[0].isFullyContained || (alignments.length > 1
                    && similarScore(alignments[0].score, alignments[1].score)))
            {
                uselessAcc ~= alignments[0].contigB.id;
            }
            else
            {
                candidatesAcc ~= alignments;
            }
        }
        // dfmt off
        logJsonDiagnostic(
            "numUselessExpected", sizeReserve,
            "numUselessFound", uselessAcc.data.length,
        );
        // dfmt on

        this.catUseless ~= uselessAcc.data;
        this.catCandidates.a2b = candidatesAcc.data;
        this.catCandidates.b2a = this.catCandidates.b2a.filter!isNotUseless.array;
        logJsonDiagnostic("state", "exit", "function", "djunctor.filterUseless");

        return this;
    }

    protected static bool similarScore(size_t a, size_t b) pure
    {
        auto diff = a > b ? a - b : b - a;
        auto magnitude = a > b ? a : b;

        return diff < maxAbsoluteDiff || (diff * AlignmentChain.maxScore / magnitude) < maxRelativeDiff;
    }

    unittest
    {
        immutable eps = 1;
        immutable refValueSmall = 2 * maxAbsoluteDiff;
        immutable refValueLarge = AlignmentChain.maxScore / 2;

        // test absolute part
        assert(similarScore(refValueSmall, refValueSmall));
        assert(similarScore(refValueSmall, refValueSmall + (maxAbsoluteDiff - 1)));
        assert(similarScore(refValueSmall, refValueSmall - (maxAbsoluteDiff - 1)));
        assert(!similarScore(refValueSmall, refValueSmall + (maxAbsoluteDiff + 1)));
        assert(!similarScore(refValueSmall, refValueSmall - (maxAbsoluteDiff + 1)));
        // test relative part
        assert(similarScore(refValueLarge, refValueLarge));
        assert(similarScore(refValueLarge,
                refValueLarge + refValueLarge * maxRelativeDiff / AlignmentChain.maxScore));
        assert(similarScore(refValueLarge,
                refValueLarge - refValueLarge * maxRelativeDiff / AlignmentChain.maxScore) + eps);
        assert(!similarScore(refValueLarge,
                refValueLarge + refValueLarge * 2 * maxRelativeDiff / AlignmentChain.maxScore));
        assert(!similarScore(refValueLarge,
                refValueLarge - refValueLarge * 2 * maxRelativeDiff / AlignmentChain.maxScore));
    }

    protected DJunctor findHits()
    {
        logJsonDiagnostic("state", "enter", "function", "djunctor.findHits");

        auto pileUps = buildPileUps(catCandidates);

        logFillingInfo!("findHits", "raw")(pileUps);
        logJsonDebug("pileUps", pileUps.toJson);

        foreach (pileUp; pileUps)
        {
            assert(pileUp.isValid, "ambiguous pile up");

            auto croppedDbResult = getCroppedPileUpDb(pileUp);
            auto hit = buildConsensus(croppedDbResult);

            // Mark best scoring read for later sequence insertion.
            catHits ~= hit;
            logJsonDebug("hit", hit.toJson);
        }

        logJsonDiagnostic("state", "exit", "function", "djunctor.findHits");

        return this;
    }

    protected Hit getCroppedPileUpDb(PileUp pileUp)
    {
        auto pileUpType = pileUp.getType;
        auto complementaryPileUp = pileUp.getComplementaryOrder(readsAlignment.a2b);
        fetchTracePoints(complementaryPileUp);
        auto tracePointDistance = complementaryPileUp[0][0].tracePointDistance;
        auto croppingRefPositions = pileUp.getCroppingRefPositions(tracePointDistance);
        logJsonDebug("croppingRefPositions", croppingRefPositions.toJson);
        size_t[] readIds = pileUp.map!"a[0].contigA.id".array;
        // dfmt off
        auto rawFastaEntries = zip(
            readIds,
            getFastaEntries(options.readsDb, readIds, options).map!parseFastaRecord,
        ).array.assumeSorted!"a[0] < b[0]";
        // dfmt on
        auto croppedFastaEntries = appender!(string[]);
        croppedFastaEntries.reserve(rawFastaEntries.length);
        size_t bestInsertionScore;
        size_t bestInsertionLength;
        size_t bestInsertionIdx;
        AlignmentChain.Complement bestInsertionComplement;

        foreach (croppedDbIdx, complementaryReadAlignment; complementaryPileUp)
        {
            immutable dummyFastaEntry = typeof(rawFastaEntries.front[1])();
            auto readId = complementaryReadAlignment[0].contigB.id;
            auto rawFastaEntry = rawFastaEntries.equalRange(tuple(readId, dummyFastaEntry)).front[1];
            immutable lineSep = typeof(rawFastaEntry).lineSep;
            auto readCroppingSlice = complementaryReadAlignment
                .map!(alignmentChain => alignmentChain.getCroppingSlice(croppingRefPositions))
                .fold!((aSlice, bSlice) => intersection(aSlice, bSlice));
            auto insertionLength = readCroppingSlice[1] - readCroppingSlice[0];

            // dfmt off
            debug logJsonDebug(
                "readId", readId,
                "complement", complementaryReadAlignment[0].complement,
                "rawReadCroppingSlices", complementaryReadAlignment.map!(alignmentChain => alignmentChain.getCroppingSlice(croppingRefPositions)).array.toJson,
                "readCroppingSlice", readCroppingSlice.toJson,
                "croppedReadLength", readCroppingSlice[1] - readCroppingSlice[0],
            );
            // dfmt on

            if (complementaryReadAlignment[0].complement)
            {
                auto readLength = complementaryReadAlignment[0].contigB.length;

                foreach (ref endpoint; readCroppingSlice)
                    endpoint = readLength - endpoint;
                swap(readCroppingSlice.expand);
            }
            assert(readCroppingSlice[0] < readCroppingSlice[1], "invalid/empty read cropping slice");

            size_t insertionScore = getInsertionScore(complementaryReadAlignment, pileUpType);

            if (insertionScore > bestInsertionScore)
            {
                bestInsertionScore = insertionScore;
                bestInsertionLength = insertionLength;
                bestInsertionIdx = croppedDbIdx;
                bestInsertionComplement = complementaryReadAlignment[0].complement;
            }

            auto croppedLength = readCroppingSlice[1] - readCroppingSlice[0];
            auto headerBuilder = parsePacBioHeader(rawFastaEntry.header);
            headerBuilder.qualityRegionBegin = 0;
            headerBuilder.qualityRegionBegin = croppedLength;
            auto newHeader = headerBuilder.to!string;
            auto expectedFastaLength = newHeader.length + croppedLength + croppedLength / options.fastaLineWidth + 2;
            auto croppedFastaEntry = appender!string;
            croppedFastaEntry.reserve(expectedFastaLength);

            croppedFastaEntry ~= newHeader;
            croppedFastaEntry ~= lineSep;
            // dfmt off
            auto croppedSequence = complementaryReadAlignment[0].complement
                ? rawFastaEntry[readCroppingSlice[0] .. readCroppingSlice[1]].array.reverseComplement
                : rawFastaEntry[readCroppingSlice[0] .. readCroppingSlice[1]].array;
            // dfmt on
            croppedFastaEntry ~= croppedSequence.chunks(options.fastaLineWidth).joiner(lineSep);
            croppedFastaEntry ~= lineSep;

            croppedFastaEntries ~= croppedFastaEntry.data;
        }

        string croppedPileUpDb = buildDamFile(croppedFastaEntries.data, options);

        with (CoordinateTransform) with (ReadAlignmentType)
            {
                return Hit(
                    Insertion.make(
                        Coordinate(
                            croppingRefPositions.backContigId,
                            croppingRefPositions.backIdx,
                        ),
                        Coordinate(
                            croppingRefPositions.frontContigId,
                            croppingRefPositions.frontIdx,
                        ),
                        bestInsertionLength,
                        pileUpType,
                    ),
                    bestInsertionIdx,
                    bestInsertionComplement,
                    croppedPileUpDb,
                );
            }
    }

    protected size_t getInsertionScore(in ReadAlignment readAlignment, in ReadAlignmentType pileUpType) pure
    {
        immutable shortAlignmentPenaltyMagnitude = AlignmentChain.maxScore / 512;
        immutable notSpanningPenaltyMagnitude = AlignmentChain.maxScore / 2;

        long expectedAlignmentCount = pileUpType == ReadAlignmentType.gap ? 2 : 1;
        long avgAlignmentLength = readAlignment.map!"a.totalLength".sum / expectedAlignmentCount;
        long avgAlignmentScore = readAlignment.map!"a.score".sum / expectedAlignmentCount;
        long shortAlignmentPenalty = floor(shortAlignmentPenaltyMagnitude * ((options.goodAnchorLength + 1) / avgAlignmentLength.to!float)^^2).to!size_t;
        size_t score = max(0, avgAlignmentScore - shortAlignmentPenalty);

        // dfmt off
        debug logJsonDebug(
            "expectedAlignmentCount", expectedAlignmentCount,
            "avgAlignmentLength", avgAlignmentLength,
            "avgAlignmentScore", avgAlignmentScore,
            "shortAlignmentPenalty", shortAlignmentPenalty,
            "score", score,
        );
        // dfmt on

        return score;
    }

    protected Hit buildConsensus(Hit croppedDbResult)
    {
        string[] consensusDbs = getConsensus(croppedDbResult.dbFile, options);

        debug logJsonDebug("consensusDbs", consensusDbs.map!Json.array);
        assert(consensusDbs.length > 0, "empty consensus");

        if (consensusDbs.length > 1)
        {
            // dfmt off
            logJsonWarn(
                "warn", "more consensus DBs than expected; this may lead to unexpected results",
                "consensusDbs", consensusDbs.map!Json.array,
            );
            // dfmt on
        }

        croppedDbResult.dbFile = consensusDbs[0];

        return croppedDbResult;
    }

    protected DJunctor fetchTracePoints(PileUp pileUp)
    {
        auto allAlignmentChains = pileUp.getAlignmentChainRefs();
        allAlignmentChains.sort!("*a < *b", SwapStrategy.stable);
        allAlignmentChains.attachTracePoints(options.refDb, options.readsDb,
                options.damapperOptions, options);

        return this;
    }

    protected DJunctor insertHits()
    {
        logJsonDiagnostic("state", "enter", "function", "djunctor.insertHits");
        foreach (hit; catHits)
        {
            insertHit(hit);
        }
        // dfmt off
        logJsonDiagnostic(
            "step", "djunctor.insertHits",
            "coordTransform", coordTransform.toJson,
            "coordTransformPython", coordTransform.toString(),
        );
        // dfmt on
        // Clear `catHits` for next iteration.
        catHits.length = 0;
        logJsonDiagnostic("state", "exit", "function", "djunctor.insertHits");

        return this;
    }

    protected DJunctor insertHit(in Hit hit)
    {
        logJsonDiagnostic("state", "enter", "function", "djunctor.insertHit");

        // dfmt off
        const(size_t)[] refContigIds = [
            hit.insertion.begin.contigId,
            hit.insertion.end.contigId,
        ];
        // dfmt on
        auto complement = hit.complement;
        size_t readId = hit.readId;
        auto read = (() => {
            auto fastaEntries = getFastaEntries(hit.dbFile, [readId], options);
            auto front = fastaEntries.front;
            auto parsedFastaRecord = front.parseFastaRecord;

            return parsedFastaRecord;
        })()();
        auto insertSequence = read[];
        CoordinateTransform.Insertion insertionInfo = hit.insertion;
        // Update insertion info with the actual consensus length.
        insertionInfo.length = read.length;
        CoordinateTransform.Insertion trInsertionInfo = insertionInfo;
        trInsertionInfo.begin = T(insertionInfo.begin);
        trInsertionInfo.end = T(insertionInfo.end);
        auto refContigFastaEntries = refContigIds.map!(
                contigId => getReferenceFastaEntry(T(contigId))).array;
        auto endContigLength = refContigFastaEntries[$ - 1].length;
        // TODO update headers
        string newHeader = refContigFastaEntries[0].header;
        size_t newContigLength = trInsertionInfo.totalInsertionLength(endContigLength);
        immutable lineSep = typeof(refContigFastaEntries[0]).lineSep;
        alias ReferenceSequence = typeof(refContigFastaEntries[0][0 .. 0]);
        ReferenceSequence beforeInsertSequence;
        ReferenceSequence afterInsertSequence;
        auto fastaBuilder = appender!string;
        auto expectedFastaLength = newHeader.length + newContigLength
            + newContigLength / options.fastaLineWidth + 2;
        fastaBuilder.reserve(expectedFastaLength);

        fastaBuilder ~= newHeader;
        fastaBuilder ~= lineSep;

        with (CoordinateTransform.Insertion.Type)
        {
            final switch (insertionInfo.type)
            {
            case gap:
                beforeInsertSequence = refContigFastaEntries[0][0 .. trInsertionInfo.begin.idx];
                afterInsertSequence = refContigFastaEntries[1][trInsertionInfo.end.idx .. $];
                break;
            case front:
                beforeInsertSequence = refContigFastaEntries[0][0 .. 0];
                afterInsertSequence = refContigFastaEntries[0][trInsertionInfo.end.idx .. $];
                break;
            case back:
                beforeInsertSequence = refContigFastaEntries[0][0 .. trInsertionInfo.begin.idx],
                afterInsertSequence = refContigFastaEntries[0][0 .. 0];
                break;
            case relabel:
                assert(0, "invalid insertion");
            }
        }

        auto joinedSequence = chain(beforeInsertSequence, insertSequence, afterInsertSequence);
        fastaBuilder ~= joinedSequence.chunks(options.fastaLineWidth).joiner(lineSep);
        fastaBuilder ~= lineSep;

        debug assert(expectedFastaLength >= fastaBuilder.data.length,
                format!"avoid reallocation: expected %d but used %d"(expectedFastaLength,
                    fastaBuilder.data.length));

        auto overlayDb = buildDamFile([fastaBuilder.data], options);

        T ~= insertionInfo;
        referenceDb.addAlias(T(refContigIds[0]), overlayDb, 1);

        immutable sequencePreviewLength = 500;
        // dfmt off
        logJsonDebug(
            "step", "insertHits",
            "type", insertionInfo.type.to!string,
            "contigIds", refContigIds.map!Json.array,
            "readId", readId,
            "complement", complement.to!bool,
            "insertSequence", insertSequence.array.to!string,
            "beforeInsertSequence", beforeInsertSequence.tail(sequencePreviewLength).array.to!string,
            "afterInsertSequence", afterInsertSequence.take(sequencePreviewLength).array.to!string,
            "insertionInfo", insertionInfo.toJson(),
            "transformedInsertionInfo", trInsertionInfo.toJson(),
            "overlayDb", overlayDb,
            "expectedFastaLength", expectedFastaLength,
            "newHeaderLength", newHeader.length,
            "newContigLength", newContigLength,
            "optionsFastaLineWidth", options.fastaLineWidth,
            "fastaBuilderDataLength", fastaBuilder.data.length,
        );
        // dfmt on
        debug assert(expectedFastaLength >= fastaBuilder.data.length,
                format!"avoid reallocation: expected %d but used %d"(expectedFastaLength,
                    fastaBuilder.data.length));

        logJsonDiagnostic("state", "exit", "function", "djunctor.insertHit");

        return this;
    }

    protected auto getReferenceFastaEntry(in size_t contigId) const
    {
        auto dbRef = referenceDb[contigId];

        return getReferenceFastaEntry(dbRef);
    }

    protected auto getReferenceFastaEntry(in DBUnion.DBReference dbRef) const
    {
        return parseFastaRecord(getFastaEntries(dbRef.dbFile, [dbRef.contigId + 0], options).front);
    }

    protected DJunctor finish()
    {
        logJsonDiagnostic("state", "enter", "function", "djunctor.finish");

        // dfmt off
        this
            .mergeOverlayDb
            .writeCoordTransform;
        // dfmt on

        logJsonDiagnostic("state", "exit", "function", "djunctor.finish");

        return this;
    }

    protected DJunctor mergeOverlayDb()
    {
        size_t numReferenceContigs = getNumContigs(referenceDb.baseDb, options);
        // dfmt off
        auto contigSources = iota(1, numReferenceContigs + 1)
            .map!(contigId => T(contigId))
            .map!(trContigId => tuple(trContigId, referenceDb[trContigId]))
            .group!"a[0] == b[0]"
            .map!"a[0]";
        auto fastaEntries = contigSources
            .save
            .map!"a[1]"
            .map!(dbRef => getReferenceFastaEntry(dbRef));
        immutable lineSep = typeof(fastaEntries.front).lineSep;

        fastaEntries
            .map!(fastaEntry => fastaEntry.toFasta(options.fastaLineWidth))
            .joiner(lineSep)
            .write;
        // dfmt on

        // dfmt off
        logJsonDebug(
            "numReferenceContigs", numReferenceContigs,
            "contigSources", contigSources.array.toJson,
        );
        // dfmt on

        foreach (newContigId, contigSource; contigSources.enumerate(1))
        {
            size_t oldContigId = contigSource[0];

            if (oldContigId != newContigId)
            {
                T ~= CoordinateTransform.Insertion.makeRelabel(oldContigId, newContigId);
                logJsonDebug("relabel", ["oldContigId" : oldContigId,
                        "newContigId" : newContigId].toJson);
            }
        }

        return this;
    }

    protected DJunctor writeCoordTransform()
    {
        if (!(options.coordTransform is null))
        {
            auto coordTransformFile = File(options.coordTransform, "w");

            coordTransformFile.write(
                    coordTransform.toString(CoordinateTransform.ExportLanguage.python));
        }

        return this;
    }

    private void logFillingInfo(string step, string readState)(in PileUp[] pileUpsByGap)
    {
        static auto getGapInfo(in PileUp pileUp)
        {
            bool isGap = pileUp.isGap;
            auto type = isGap ? ReadAlignmentType.gap : pileUp[0].getType;

            alias lengthFilter = p => isGap ? p.length == 2 : p.length == 1;

            // dfmt off
            return Json([
                "type": Json(type.to!string),
                "contigIds": pileUp
                    .filter!lengthFilter
                    .map!(readAlignment => readAlignment
                        .map!"a.contigB.id"
                        .map!Json
                        .array)
                    .front
                    .Json,
                "estimateLengthMean": pileUp
                    .filter!lengthFilter
                    .map!getInsertionSize
                    .array
                    .mean
                    .Json,
                "estimateLengthMedian": pileUp
                    .filter!lengthFilter
                    .map!getInsertionSize
                    .array
                    .median
                    .Json,
                "numReads": Json(pileUp.length),
            ]);
            // dfmt on
        }

        // dfmt off
        logJsonDiagnostic(
            "step", step,
            "readState", readState,
            "numGaps", pileUpsByGap.length,
            "gapInfo", pileUpsByGap.map!getGapInfo.array,
        );
        // dfmt on
    }
}

/**
                This realizes a transform which translates coordinates before insertions
                to coordinates after them. The insertions can be added sequentially. The
                transform is fully functional at any point in time and can be converted
                to a Python 2.7 script for external usage.
            */
struct CoordinateTransform
{
    static struct Coordinate
    {
        size_t contigId;
        size_t idx;
    }

    static struct Insertion
    {
        enum Type
        {
            front,
            gap,
            back,
            relabel,
        }

        Coordinate begin;
        Coordinate end;
        private size_t _length;

        @disable this(Coordinate);
        @disable this(Coordinate, Coordinate);
        private this(Coordinate begin, Coordinate end, size_t length) pure nothrow
        {
            this.begin = begin;
            this.end = end;
            this._length = length;
        }

        static Insertion makeRelabel(size_t oldContigId, size_t newContigId) pure nothrow
        {
            return Insertion(Coordinate(newContigId), Coordinate(oldContigId), 0);
        }

        static Insertion make(Coordinate begin, Coordinate end, size_t length,
                in ReadAlignmentType insertionType) pure nothrow
        {
            final switch (insertionType)
            {
            case ReadAlignmentType.front:
                begin.idx = 0;
                break;
            case ReadAlignmentType.gap:
                break;
            case ReadAlignmentType.back:
                end.idx = begin.idx + length;
                break;
            }

            return Insertion(begin, end, length);
        }

        invariant
        {
            assert(begin != end || (begin.idx == end.idx && _length > 0), "empty insertion");
            if (begin.contigId == end.contigId)
            {
                assert(end.idx >= begin.idx,
                        "insertion begin index should be before/equal to end index on same contig");
                assert(_length >= end.idx - begin.idx, "inserted sequence too small");
            }
        }

        @property length() pure const nothrow
        {
            return _length;
        }

        @property length(size_t newLength) pure nothrow
        {
            final switch (type)
            {
                case Type.back:
                    end.idx = begin.idx + newLength;
                    goto case;
                case Type.front:
                case Type.gap:
                    _length = newLength;
                    return;
                case Type.relabel:
                    return;
            }
        }

        @property isExtension() pure const nothrow
        {
            return begin.contigId == end.contigId;
        }

        @property isFrontExtension() pure const nothrow
        {
            return isExtension && !isBackExtension;
        }

        @property isBackExtension() pure const nothrow
        {
            return isExtension && (begin.idx > 0 && end.idx == begin.idx + length);
        }

        @property isGap() pure const nothrow
        {
            return !isExtension && !isRelabel;
        }

        @property isRelabel() pure const nothrow
        {
            return !isExtension && begin.idx == 0 && end.idx == 0 && length == 0;
        }

        @property type() pure const nothrow
        {
            if (isExtension)
            {
                return isFrontExtension ? Type.front : Type.back;
            }
            else if (isGap)
            {
                return Type.gap;
            }
            else
            {
                return Type.relabel;
            }
        }

        size_t totalInsertionLength(in size_t endContigLength) pure const nothrow
        {
            final switch (type)
            {
            case Type.gap:
            case Type.front:
            case Type.relabel:
                return begin.idx + length + endContigLength - end.idx;
            case Type.back:
                return begin.idx + length;
            }
        }
    }

    static enum ExportLanguage
    {
        python,
    }

    const(Insertion)[] insertions;

    /**
        Transform a coordinate according to insertions. If the designated
        contig is not present in insertions then the original coordinate is
        returned. Otherwise the new coordinate is computed as follows:

        1. Find insertion where end contig matches input. Note: if begin
           contig matches input not action is needed because neither the
           `idx` nor the `contigId` is affected.
        2. Travel from that point backwards to the beginning of the chain –
           ie. until begin contig and end contig of two insertions are
           unequal – summing all the coordinate shifts. The new reference
           contig is the last one in the (reverse) chain.

        For a single insertion and input coordinate `x` the shifted
        coordinate `T(x)` is calculated as follows:

        ---
        Case 1:  begin contig == end contig
            Case 1.1:  extend begin of contig

                   begin/end contig

                  0  ie  x
                  |---+--+-------|
                 /    |
                |-----|
                0    li

                T(x) = x
                     - ie  // (1) make relative to insertion point of end contig
                     + li  // (2) make relative to begin of insertion
                // choose ib = 0
                     = x - ie + li + ib

            Case 1.2:  extend end of contig

                        begin/end contig

                0       x ib        ie
                |-------+--+---| - - +
                           |         |
                           |---------|
                           0        li

                T(x) = x
                // choose ie = li + ib
                     = x - ie + li + ib

        Case 2:  begin contig != end contig

              begin contig       end contig

            0         ib     0  ie  x
            |----------+---| |---+--+-------|
                       |         |
                       |---------|
                       0        li

                        insertion

            T(x) = x
                 - ie  // (1) make relative to insertion point of end contig
                 + li  // (2) make relative to begin of insertion
                 + ib  // (3) make relative to begin of begin contig
        ---

        Returns: transformed coordinate.
    */
    Coordinate transform(in size_t contigId, in size_t idx) pure const
    {
        return transform(Coordinate(contigId, idx));
    }

    size_t transform(in size_t contigId) pure const
    {
        return transform(Coordinate(contigId, 0)).contigId;
    }

    /// ditto
    Coordinate transform(in Coordinate input) pure const
    {
        Coordinate result = input;

        // dfmt off
        auto insertionChain = insertions[]
            .retro
            .find!"a.end.contigId == b"(0 + input.contigId);
        // dfmt on
        size_t lastContigId = input.contigId;

        foreach (insertion; insertionChain)
        {
            if (insertion.end.contigId != lastContigId)
            {
                break; // chain end reached
            }

            // Note: reverse order of calculation to prevent negative
            // intermediate results
            result.idx += insertion.begin.idx; // (3)
            result.idx += insertion.length; // (2)
            result.idx -= insertion.end.idx; // (1)
            // Begin contig is the new reference
            result.contigId = insertion.begin.contigId;
            lastContigId = insertion.begin.contigId;
        }

        return result;
    }

    /// ditto
    alias opCall = transform;

    ///
    unittest
    {
        CoordinateTransform coordTransform;
        //           0        495 500 0   5  100   600
        // reference |----------+---| |---+--+-------|
        //                      |         |
        // insert               |---------|
        //                      0       100
        // dfmt off
        coordTransform.add(Insertion.make(
            Coordinate(1, 495),
            Coordinate(2, 5),
            100,
            ReadAlignmentType.gap,
        ));
        // dfmt on

        auto inputCoord = Coordinate(2, 100);
        // dfmt off
        auto transformedCoord = Coordinate(
            1,
            (
                100
                - 5    // coord on contig 2 relative to insertion point
                + 100  // make relative to begin of insertion
                + 495  // make relative to begin of contig 1
            )
        );
        // dfmt on
        assert(coordTransform.transform(inputCoord) == transformedCoord);
        assert(coordTransform(inputCoord) == transformedCoord);
    }

    unittest
    {
        with (ReadAlignmentType)
        {
            CoordinateTransform coordTransform;
            // dfmt off
            coordTransform.add(Insertion.make(
                Coordinate(100, 0),
                Coordinate(101, 0),
                0,
                gap,
            ));
            coordTransform.add(Insertion.make(
                Coordinate(1, 495),
                Coordinate(2, 5),
                100,
                gap,
            ));
            coordTransform.add(Insertion.make(
                Coordinate(2, 490),
                Coordinate(5, 10),
                150,
                gap,
            ));
            coordTransform.add(Insertion.make(
                Coordinate(5, 480),
                Coordinate(3, 20),
                200,
                gap,
            ));
            coordTransform.add(Insertion.make(
                Coordinate(102, 0),
                Coordinate(103, 0),
                0,
                gap,
            ));
            coordTransform.add(Insertion.make(
                Coordinate(200),
                Coordinate(200, 5),
                10,
                front,
            ));
            coordTransform.add(Insertion.make(
                Coordinate(201, 90),
                Coordinate(201),
                20,
                back,
            ));
            coordTransform.add(Insertion.make(
                Coordinate(202),
                Coordinate(202, 15),
                30,
                front,
            ));
            coordTransform.add(Insertion.make(
                Coordinate(202),
                Coordinate(202, 20),
                40,
                front,
            ));
            coordTransform.add(Insertion.make(
                Coordinate(203, 65),
                Coordinate(203),
                50,
                back,
            ));
            coordTransform.add(Insertion.make(
                Coordinate(203, 60),
                Coordinate(203),
                60,
                back,
            ));
            coordTransform.add(Insertion.makeRelabel(
                1024,
                1000,
            ));
            // dfmt on

            {
                // Case: contig not present in insertions
                auto inputCoord = Coordinate(1337, 42);
                assert(coordTransform.transform(inputCoord) == inputCoord);
            }
            {
                // Case: contig is never end of an insertions
                auto inputCoord = Coordinate(1, 64);
                assert(coordTransform.transform(inputCoord) == inputCoord);
            }
            {
                // Case: contig is end of the last insertions of a reverse chain.
                auto inputCoord = Coordinate(2, 64);
                // dfmt off
                auto transformedCoord = Coordinate(
                    1,
                    (
                          64   // last and only insertion in chain
                        - 5    // (1)
                        + 100  // (2)
                        + 495  // (3)
                    )
                );
                // dfmt on
                assert(coordTransform.transform(inputCoord) == transformedCoord);
            }
            {
                // Case: contig is end of an inner insertions, ie. there are
                //       insertions before and after in the chain.
                auto inputCoord = Coordinate(5, 64);
                // dfmt off
                auto transformedCoord = Coordinate(
                    1,
                    (
                        (
                              64   // last insertion in chain
                            - 5    // (1)
                            + 100  // (2)
                            + 495  // (3)
                        )      // first insertion in chain
                        - 10    // (1)
                        + 150  // (2)
                        + 490  // (3)
                    )

                );
                // dfmt on
                assert(coordTransform.transform(inputCoord) == transformedCoord);
            }
            {
                // Case: contig is end of the first insertions of a reverse chain.
                auto inputCoord = Coordinate(3, 64);
                // dfmt off
                auto transformedCoord = Coordinate(
                    1,
                    (
                        (
                            (
                                  64   // last insertion in chain
                                - 5    // (1)
                                + 100  // (2)
                                + 495  // (3)
                            )      // second insertion in chain
                            - 10    // (1)
                            + 150  // (2)
                            + 490  // (3)
                        )      // first insertion in chain
                        - 20    // (1)
                        + 200  // (2)
                        + 480  // (3)
                    )

                );
                // dfmt on
                assert(coordTransform.transform(inputCoord) == transformedCoord);
            }
            {
                // Case: coordinate on rejected part of end contig, ie. `x < ie`.
                auto inputCoord = Coordinate(2, 3);
                // dfmt off
                auto transformedCoord = Coordinate(
                    1,
                    (
                          3    // last and only insertion in chain
                        - 5    // (1)
                        + 100  // (2)
                        + 495  // (3)
                    )
                );
                // dfmt on
                assert(coordTransform.transform(inputCoord) == transformedCoord);
            }

            {
                // Case: simple front extension
                auto inputCoord = Coordinate(200, 3);
                // dfmt off
                auto transformedCoord = Coordinate(
                    200,
                    (
                          3    // only extension (in chain)
                        - 5    // (1)
                        + 10   // (2)
                        + 0    // (3)
                    )
                );
                // dfmt on
                assert(coordTransform.transform(inputCoord) == transformedCoord);
            }
            {
                // Case: simple back extension
                auto inputCoord = Coordinate(201, 3);
                // dfmt off
                auto transformedCoord = Coordinate(
                    201,
                    (
                          3    // only extension (in chain)
                        - 110  // (1)
                        + 20   // (2)
                        + 90   // (3)
                    )
                );
                // dfmt on
                assert(coordTransform.transform(inputCoord) == transformedCoord);
            }
            {
                // Case: double front extension
                auto inputCoord = Coordinate(202, 3);
                // dfmt off
                auto transformedCoord = Coordinate(
                    202,
                    (
                        (
                              3    // last extension (in chain)
                            - 15   // (1)
                            + 30   // (2)
                            +  0   // (3)
                        )     // first extension (in chain)
                        - 20  // (1)
                        + 40  // (2)
                        +  0  // (3)
                    )
                );
                // dfmt on
                assert(coordTransform.transform(inputCoord) == transformedCoord);
            }
            {
                // Case: double back extension
                auto inputCoord = Coordinate(203, 3);
                // dfmt off
                auto transformedCoord = Coordinate(
                    203,
                    (
                        (
                              3    // last extension (in chain)
                            - 115  // (1)
                            + 50   // (2)
                            + 65   // (3)
                        )      // first extension (in chain)
                        - 120  // (1)
                        + 60   // (2)
                        + 60   // (3)
                    )
                );
                // dfmt on
                assert(coordTransform.transform(inputCoord) == transformedCoord);
            }
            {
                // Case: relabel
                auto inputCoord = Coordinate(1024, 3);
                auto transformedCoord = Coordinate(1000, 3);

                assert(coordTransform.transform(inputCoord) == transformedCoord);
            }

            // dfmt off
            coordTransform.add(Insertion.make(
                Coordinate(203, 115),
                Coordinate(202, 5),
                50,
                gap,
            ));
            // dfmt on

            {
                // Case: double front extension + double back extension + gap spanned
                auto inputCoord = Coordinate(202, 3);
                // dfmt off
                auto transformedCoord = Coordinate(
                    203,
                    (
                        (
                            (
                                  3    // last front extension (in chain)
                                - 15   // (1)
                                + 30   // (2)
                                +  0   // (3)
                            )     // first front extension (in chain)
                            - 20  // (1)
                            + 40  // (2)
                            +  0  // (3)
                        )      // only gap (in chain)
                        -   5  // (1)
                        +  50  // (2)
                        + 115  // (3)
                    )
                );
                // dfmt on
                assert(coordTransform.transform(inputCoord) == transformedCoord);
            }
        }
    }

    /**
        Add an insertion to this transform.

        Returns: this transform.
        See_Also: `CoordinateTransform.transform`.
    */
    CoordinateTransform add(in Insertion newInsertion) pure
    {
        size_t idx = getInsertionIndex(newInsertion);

        if (idx == 0)
        {
            insertions = [newInsertion] ~ insertions;
        }
        else if (idx == insertions.length)
        {
            insertions ~= newInsertion;
        }
        else
        {
            insertions = insertions[0 .. idx] ~ [newInsertion] ~ insertions[idx .. $];
        }

        return this;
    }

    /// ditto
    void opOpAssign(string op)(in Insertion newInsertion) pure if (op == "~")
    {
        this.add(newInsertion);
    }

    unittest
    {
        with (ReadAlignmentType)
        {
            Insertion getDummyInsertion(in size_t beginContigId,
                    in size_t endContigId, in ReadAlignmentType insertionType = gap)
            {
                static immutable dummyLength = 42;

                return Insertion.make(Coordinate(beginContigId, dummyLength / 2),
                        Coordinate(endContigId, dummyLength / 2), dummyLength, insertionType);
            }

            CoordinateTransform getDummyTransform()
            {
                CoordinateTransform coordTransform;

                // dfmt off
                coordTransform.insertions = [
                    getDummyInsertion(1, 2),
                    getDummyInsertion(3, 4),
                    getDummyInsertion(9, 8),
                    getDummyInsertion(6, 5),
                ];
                // dfmt on

                return coordTransform;
            }

            {
                auto coordTransform = getDummyTransform();
                // Case 1 (Gap): newInsertion fits between two existing insertions
                // dfmt off
                coordTransform.add(getDummyInsertion(2, 3));
                assert(coordTransform.insertions == [
                    getDummyInsertion(1, 2),
                    getDummyInsertion(2, 3),
                    getDummyInsertion(3, 4),
                    getDummyInsertion(9, 8),
                    getDummyInsertion(6, 5),
                ]);
                coordTransform.add(getDummyInsertion(8, 6));
                assert(coordTransform.insertions == [
                    getDummyInsertion(1, 2),
                    getDummyInsertion(2, 3),
                    getDummyInsertion(3, 4),
                    getDummyInsertion(9, 8),
                    getDummyInsertion(8, 6),
                    getDummyInsertion(6, 5),
                ]);
                // dfmt on
            }
            {
                auto coordTransform = getDummyTransform();
                // Case 2* (Gap): newInsertion extends an existing insertion chain in the front
                // dfmt off
                coordTransform.add(getDummyInsertion(42, 1));
                assert(coordTransform.insertions == [
                    getDummyInsertion(42, 1),
                    getDummyInsertion(1, 2),
                    getDummyInsertion(3, 4),
                    getDummyInsertion(9, 8),
                    getDummyInsertion(6, 5),
                ]);
                // Case 2 (Gap): newInsertion extends an existing insertion chain in the front
                coordTransform.add(getDummyInsertion(42, 9));
                assert(coordTransform.insertions == [
                    getDummyInsertion(42, 1),
                    getDummyInsertion(1, 2),
                    getDummyInsertion(3, 4),
                    getDummyInsertion(42, 9),
                    getDummyInsertion(9, 8),
                    getDummyInsertion(6, 5),
                ]);
                // dfmt on
            }
            {
                auto coordTransform = getDummyTransform();
                // Case 2* (Extension): newInsertion extends an existing insertion chain in the front
                // dfmt off
                coordTransform.add(getDummyInsertion(1, 1, front));
                assert(coordTransform.insertions == [
                    getDummyInsertion(1, 1, front),
                    getDummyInsertion(1, 2),
                    getDummyInsertion(3, 4),
                    getDummyInsertion(9, 8),
                    getDummyInsertion(6, 5),
                ]);
                coordTransform.add(getDummyInsertion(1, 1, front));
                assert(coordTransform.insertions == [
                    getDummyInsertion(1, 1, front),
                    getDummyInsertion(1, 1, front),
                    getDummyInsertion(1, 2),
                    getDummyInsertion(3, 4),
                    getDummyInsertion(9, 8),
                    getDummyInsertion(6, 5),
                ]);
                // Case 2 (Extension): newInsertion extends an existing insertion chain in the front
                coordTransform.add(getDummyInsertion(9, 9, front));
                assert(coordTransform.insertions == [
                    getDummyInsertion(1, 1, front),
                    getDummyInsertion(1, 1, front),
                    getDummyInsertion(1, 2),
                    getDummyInsertion(3, 4),
                    getDummyInsertion(9, 9, front),
                    getDummyInsertion(9, 8),
                    getDummyInsertion(6, 5),
                ]);
                coordTransform.add(getDummyInsertion(9, 9, front));
                assert(coordTransform.insertions == [
                    getDummyInsertion(1, 1, front),
                    getDummyInsertion(1, 1, front),
                    getDummyInsertion(1, 2),
                    getDummyInsertion(3, 4),
                    getDummyInsertion(9, 9, front),
                    getDummyInsertion(9, 9, front),
                    getDummyInsertion(9, 8),
                    getDummyInsertion(6, 5),
                ]);
                // dfmt on
            }
            {
                auto coordTransform = getDummyTransform();
                // Case 3 (Gap): newInsertion extends an existing insertion chain in the back
                // dfmt off
                coordTransform.add(getDummyInsertion(4, 42));
                assert(coordTransform.insertions == [
                    getDummyInsertion(1, 2),
                    getDummyInsertion(3, 4),
                    getDummyInsertion(4, 42),
                    getDummyInsertion(9, 8),
                    getDummyInsertion(6, 5),
                ]);
                coordTransform.add(getDummyInsertion(5, 42));
                assert(coordTransform.insertions == [
                    getDummyInsertion(1, 2),
                    getDummyInsertion(3, 4),
                    getDummyInsertion(4, 42),
                    getDummyInsertion(9, 8),
                    getDummyInsertion(6, 5),
                    getDummyInsertion(5, 42),
                ]);
                // dfmt on
            }
            {
                auto coordTransform = getDummyTransform();
                // Case 3 (Extension): newInsertion extends an existing insertion chain in the back
                // dfmt off
                coordTransform.add(getDummyInsertion(4, 4, back));
                assert(coordTransform.insertions == [
                    getDummyInsertion(1, 2),
                    getDummyInsertion(3, 4),
                    getDummyInsertion(4, 4, back),
                    getDummyInsertion(9, 8),
                    getDummyInsertion(6, 5),
                ]);
                coordTransform.add(getDummyInsertion(4, 4, back));
                assert(coordTransform.insertions == [
                    getDummyInsertion(1, 2),
                    getDummyInsertion(3, 4),
                    getDummyInsertion(4, 4, back),
                    getDummyInsertion(4, 4, back),
                    getDummyInsertion(9, 8),
                    getDummyInsertion(6, 5),
                ]);
                coordTransform.add(getDummyInsertion(5, 5, back));
                assert(coordTransform.insertions == [
                    getDummyInsertion(1, 2),
                    getDummyInsertion(3, 4),
                    getDummyInsertion(4, 4, back),
                    getDummyInsertion(4, 4, back),
                    getDummyInsertion(9, 8),
                    getDummyInsertion(6, 5),
                    getDummyInsertion(5, 5, back),
                ]);
                coordTransform.add(getDummyInsertion(5, 5, back));
                assert(coordTransform.insertions == [
                    getDummyInsertion(1, 2),
                    getDummyInsertion(3, 4),
                    getDummyInsertion(4, 4, back),
                    getDummyInsertion(4, 4, back),
                    getDummyInsertion(9, 8),
                    getDummyInsertion(6, 5),
                    getDummyInsertion(5, 5, back),
                    getDummyInsertion(5, 5, back),
                ]);
                // dfmt on
            }
            {
                auto coordTransform = getDummyTransform();
                // Case 4 (Gap): open new insertion chain
                // dfmt off
                coordTransform.add(getDummyInsertion(1337, 42));
                assert(coordTransform.insertions == [
                    getDummyInsertion(1, 2),
                    getDummyInsertion(3, 4),
                    getDummyInsertion(9, 8),
                    getDummyInsertion(6, 5),
                    getDummyInsertion(1337, 42),
                ]);
                coordTransform.add(getDummyInsertion(1337, 42));
                assert(coordTransform.insertions == [
                    getDummyInsertion(1, 2),
                    getDummyInsertion(3, 4),
                    getDummyInsertion(9, 8),
                    getDummyInsertion(6, 5),
                    getDummyInsertion(1337, 42),
                    getDummyInsertion(1337, 42),
                ]);
                // dfmt on
            }
            {
                auto coordTransform = getDummyTransform();
                // Case 4 (Extension): open new insertion chain
                // dfmt off
                coordTransform.add(getDummyInsertion(1337, 1337, front));
                assert(coordTransform.insertions == [
                    getDummyInsertion(1, 2),
                    getDummyInsertion(3, 4),
                    getDummyInsertion(9, 8),
                    getDummyInsertion(6, 5),
                    getDummyInsertion(1337, 1337, front),
                ]);
                coordTransform.add(getDummyInsertion(42, 42, back));
                assert(coordTransform.insertions == [
                    getDummyInsertion(1, 2),
                    getDummyInsertion(3, 4),
                    getDummyInsertion(9, 8),
                    getDummyInsertion(6, 5),
                    getDummyInsertion(1337, 1337, front),
                    getDummyInsertion(42, 42, back),
                ]);
                // dfmt on
            }
            {
                CoordinateTransform emptyCoordTransform;

                // Case 4* (Gap): open new insertion chain
                assert(emptyCoordTransform.add(getDummyInsertion(1337, 42))
                        .insertions == [getDummyInsertion(1337, 42)]);
            }
            {
                CoordinateTransform emptyCoordTransform;

                // Case 4* (Extension): open new insertion chain
                assert(emptyCoordTransform.add(getDummyInsertion(1337, 1337,
                        front)).insertions == [getDummyInsertion(1337, 1337, front)]);
            }
            {
                CoordinateTransform emptyCoordTransform;

                // Operator Style
                emptyCoordTransform ~= getDummyInsertion(1337, 42);
                assert(emptyCoordTransform.insertions == [getDummyInsertion(1337, 42)]);
            }
        }
    }

    private size_t getInsertionIndex(in Insertion newInsertion) pure const
    {
        // dfmt off
        if (
            // Case 4*: open new insertion chain
            insertions.length == 0 ||
            // Case 2*: newInsertion extends an existing insertion chain in the front
            newInsertion.end.contigId == insertions[0].begin.contigId
        )
        {
            return 0;
        }
        // dfmt on
        // Case 3*: newInsertion extends an existing insertion chain in the back
        // Case 4**: open new insertion chain
        else if (insertions.length == 1)
        {
            // Note: this case is needed in order for `slice(2)` to always yield
            //       pairs (of size two).
            return 1;
        }

        size_t case1Result = size_t.max;
        size_t case2Result = size_t.max;
        size_t case3Result = size_t.max;

        foreach (i, insAB; insertions.slide(2).enumerate)
        {
            auto insA = insAB[0];
            auto insB = insAB[1];

            // Case 1: newInsertion fits between two existing insertions; take first result
            if (case1Result == size_t.max && newInsertion.begin.contigId == insA.end.contigId
                    && newInsertion.end.contigId == insB.begin.contigId)
            {
                case1Result = i + 1;
            }

            // Case 2: newInsertion extends an existing insertion chain in the front; take first result
            if (case2Result == size_t.max && newInsertion.begin.contigId != insA.end.contigId
                    && newInsertion.end.contigId == insB.begin.contigId)
            {
                case2Result = i + 1;
            }

            // Case 3: newInsertion extends an existing insertion chain in the back; take last result
            if (newInsertion.begin.contigId == insA.end.contigId
                    && newInsertion.end.contigId != insB.begin.contigId)
            {
                case3Result = i + 1;
            }
        }

        if (case1Result != size_t.max)
            return case1Result;
        else if (case2Result != size_t.max)
            return case2Result;
        else if (case3Result != size_t.max)
            return case3Result;

        // Case 4: open new insertion chain
        return insertions.length;
    }

    unittest
    {
        Insertion getDummyInsertion(in size_t beginContigId, in size_t endContigId)
        {
            return Insertion(Coordinate(beginContigId, 0), Coordinate(endContigId, 0), 0);
        }

        CoordinateTransform coordTransform;

        // dfmt off
        coordTransform.insertions = [
            getDummyInsertion(1, 2),
            getDummyInsertion(3, 4),
            getDummyInsertion(9, 8),
            getDummyInsertion(6, 5),
        ];
        // dfmt on
        {
            // Case 1: newInsertion fits between two existing insertions
            assert(coordTransform.getInsertionIndex(getDummyInsertion(2, 3)) == 1);
            assert(coordTransform.getInsertionIndex(getDummyInsertion(8, 6)) == 3);
        }
        {
            // Case 2*: newInsertion extends an existing insertion chain in the front
            assert(coordTransform.getInsertionIndex(getDummyInsertion(42, 1)) == 0);
            // Case 2: newInsertion extends an existing insertion chain in the front
            assert(coordTransform.getInsertionIndex(getDummyInsertion(42, 9)) == 2);
        }
        {
            // Case 3: newInsertion extends an existing insertion chain in the front
            assert(coordTransform.getInsertionIndex(getDummyInsertion(4, 42)) == 2);
            assert(coordTransform.getInsertionIndex(getDummyInsertion(5, 42)) == 4);
        }
        {
            // Case 4: open new insertion chain
            assert(coordTransform.getInsertionIndex(getDummyInsertion(1337, 42)) == 4);
            assert(coordTransform.getInsertionIndex(getDummyInsertion(42, 1337)) == 4);
        }
        {
            CoordinateTransform emptyCoordTransform;

            // Case 4*: open new insertion chain
            assert(emptyCoordTransform.getInsertionIndex(getDummyInsertion(1337, 42)) == 0);
            assert(emptyCoordTransform.getInsertionIndex(getDummyInsertion(42, 1337)) == 0);
        }
        coordTransform.insertions = [getDummyInsertion(1, 2)];
        {
            // Case 3*: newInsertion extends an existing insertion chain in the front
            assert(coordTransform.getInsertionIndex(getDummyInsertion(2, 42)) == 1);
        }
        {
            // Case 4**: open new insertion chain
            assert(coordTransform.getInsertionIndex(getDummyInsertion(1337, 42)) == 1);
        }
    }

    string toString(in ExportLanguage lang = ExportLanguage.python) pure const
    {
        immutable estimatePrefaceLength = 250;
        immutable estimateLengthPerInsertion = 100;
        immutable estimateEpilogueLength = 2000;

        auto result = appender!string;
        result.reserve(estimatePrefaceLength + estimateLengthPerInsertion
                * insertions.length + estimateEpilogueLength);

        final switch (lang)
        {
        case ExportLanguage.python:
            buildPythonString(result);
            break;
        }

        return result.data;
    }

    private void buildPythonString(out Appender!string builder) pure const
    {
        immutable preface = q"EOF
            from __future__ import print_function
            from collections import namedtuple
            from sys import argv, exit

            Coordinate = namedtuple("Coordinate", ["contig_id", "idx"])
            Insertion = namedtuple("Insertion", ["begin", "end", "length"])

            INSERTIONS = [
EOF".outdent;
        immutable epilogue = q"EOF
            ]


            def transform(contig_id, idx, insertions=INSERTIONS):
                result_contig_id = contig_id
                result_idx = idx
                last_contig_id = contig_id

                for insertion in reversed(insertions):
                    if insertion.end.contig_id != last_contig_id and \
                            last_contig_id != contig_id:
                        break  # chain end reached
                    elif insertion.end.contig_id == last_contig_id or \
                            last_contig_id != contig_id:
                        # chain begin reached or inside chain

                        result_idx += insertion.begin.idx
                        result_idx += insertion.length
                        result_idx -= insertion.end.idx
                        result_contig_id = insertion.begin.contig_id
                        last_contig_id = insertion.begin.contig_id

                return Coordinate(result_contig_id, result_idx)


            def print_help():
                print("Usage: {} [-h] CONTIG_ID IDX".format(argv[0]))
                print()
                print(
                    "Transform coordinates from input coordinates to output " +
                    "coordinates."
                )
                print("Prints `NEW_CONTIG_ID NEW_IDX` to STDOUT.")
                print()
                print("Positional arguments:")
                print(" CONTIG_ID  Contig ID in .dam before transformation")
                print(" IDX        Base index (zero-based) in sequence")
                print()
                print("Optional arguments:")
                print(" -h         Prints this help.")

            if __name__ == "__main__":
                if len(argv) != 3 or "-h" in argv:
                    print_help()
                    exit(1)

                try:
                    contig_id = int(argv[1])
                    idx = int(argv[2])

                    if contig_id < 0 or idx < 0:
                        raise Exception()
                except Exception:
                    print_help()
                    exit(2)

                transformed = transform(contig_id, idx)
                print("{} {}".format(transformed.contig_id, transformed.idx))
EOF".outdent;
        immutable insertionTemplate = q"EOF
            Insertion(
                Coordinate(%d, %d),
                Coordinate(%d, %d),
                %d
            ),
EOF".outdent.indent(4);

        builder ~= preface;
        foreach (insertion; insertions)
        {
            // dfmt off
            formattedWrite!insertionTemplate(
                builder,
                insertion.begin.contigId,
                insertion.begin.idx,
                insertion.end.contigId,
                insertion.end.idx,
                insertion.length,
            );
            // dfmt on
        }
        builder ~= epilogue;
    }
}

/// Access contigs distributed over various DBs by ID.
struct DBUnion
{
    alias DBReference = Tuple!(string, "dbFile", size_t, "contigId");

    string baseDb;
    DBReference[size_t] overlayDbs;

    /// Return the `dbFile` for `contigId`.
    DBReference opIndex(in size_t contigId) pure const
    {
        return overlayDbs.get(contigId, DBReference(baseDb, contigId));
    }

    unittest
    {
        auto dbQuery = DBUnion("baseDb");
        dbQuery.overlayDbs[42] = DBReference("overlayDb42", 1);

        assert(dbQuery[0] == DBReference("baseDb", 0));
        assert(dbQuery[7] == DBReference("baseDb", 7));
        assert(dbQuery[42] == DBReference("overlayDb42", 1));
    }

    /// Set the source for `contigId`.
    void addAlias(in size_t contigId, in string dbFile, in size_t newContigId) pure
    {
        addAlias(contigId, DBReference(dbFile, newContigId));
    }

    /// ditto
    void addAlias(in size_t contigId, in DBReference dbRef) pure
    {
        overlayDbs[contigId] = dbRef;
    }

    unittest
    {
        auto dbQuery = DBUnion("baseDb");

        assert(dbQuery[42] == DBReference("baseDb", 42));

        dbQuery.addAlias(42, DBReference("overlayDb42", 1));

        assert(dbQuery[42] == DBReference("overlayDb42", 1));
    }
}

///
unittest
{
    auto dbQuery = DBUnion("baseDb");

    dbQuery.addAlias(42, "overlayDb42", 1);

    assert(dbQuery[0].dbFile == "baseDb");
    assert(dbQuery[0].contigId == 0);
    assert(dbQuery[7].dbFile == "baseDb");
    assert(dbQuery[7].contigId == 7);
    assert(dbQuery[42].dbFile == "overlayDb42");
    assert(dbQuery[42].contigId == 1);

    dbQuery.addAlias(7, "overlayDb7", 1);

    assert(dbQuery[7].dbFile == "overlayDb7");
    assert(dbQuery[7].contigId == 1);
}
