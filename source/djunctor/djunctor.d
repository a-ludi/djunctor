/**
    This is the main algorithm of this package.

    Copyright: © 2018 Arne Ludwig <arne.ludwig@posteo.de>
    License: Subject to the terms of the MIT license, as written in the
             included LICENSE file.
    Authors: Arne Ludwig <arne.ludwig@posteo.de>
*/
module djunctor.djunctor;

import djunctor.commandline : Options;
import djunctor.util.fasta : buildFastaRecord, parseFasta, parseFastaRecord,
    reverseComplement;
import djunctor.util.log;
import djunctor.util.math : mean, median;
import djunctor.util.range : Comparator;
import djunctor.util.string : indent;
import core.exception : AssertError;
import std.algorithm : all, any, canFind, chunkBy, each, equal, filter, find,
    group, isSorted, map, max, sort, sum, swap;
import std.array : appender, Appender, array;
import std.container : BinaryHeap, heapify, make;
import std.conv;
import std.exception : assertNotThrown, assertThrown;
import std.format : format, formattedWrite;
import std.math : abs, sgn;
import std.range : assumeSorted, chunks, ElementType, enumerate, isForwardRange,
    only, retro, slide, SortedRange, walkLength;
import std.stdio : writeln;
import std.string : outdent;
import std.typecons : tuple, Tuple;
import vibe.data.json : Json, serializeToJson;

version (unittest)
{
    import djunctor.util.testing : MockCallable;
    import djunctor.dazzler : origGetLocalAlignments = getLocalAlignments,
        origGetMappings = getMappings;
    import djunctor.dazzler : buildDamFile, getConsensus, getFastaEntries;

    MockCallable!(AlignmentContainer!(AlignmentChain[]), const string, const Options) getLocalAlignments;
    MockCallable!(origGetMappings!Options) getMappings;
}
else
{
    import djunctor.dazzler : buildDamFile, getConsensus, getFastaEntries,
        getLocalAlignments, getMappings;
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
        struct Locus
        {
            size_t begin;
            size_t end;
        }

        Locus contigA;
        Locus contigB;
        size_t numDiffs;
    }

    static struct Contig
    {
        size_t id;
        size_t length;
    }

    static enum Complement
    {
        no,
        yes,
    }

    static immutable maxScore = 2 ^^ 16;

    size_t id;
    Contig contigA;
    Contig contigB;
    Complement complement;
    LocalAlignment[] localAlignments;

    invariant
    {
        assert(localAlignments.length >= 1, "empty chain is forbidden");
        assert(localAlignments.all!(la => 0 <= la.contigA.begin
                && la.contigA.begin < la.contigA.end && la.contigA.end <= contigA.length),
                "non-sense alignment of contigA");
        assert(localAlignments.all!(la => 0 <= la.contigB.begin
                && la.contigB.begin < la.contigB.end && la.contigB.end <= contigB.length),
                "non-sense alignment of contigB");
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

    @property auto ref LocalAlignment first() const pure nothrow
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

    @property auto ref LocalAlignment last() const pure nothrow
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
    Returns true iff ac1 begins befor ac2 taking complementary alignment into account.
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

/// Start the `djunctor` alorithm with preprocessed options.
void runWithOptions(in ref Options options)
{
    new DJunctor(options).run();
}

class DJunctor
{
    /// Stop after maxLoops `mainLoop`s at the latest.
    static immutable maxLoops = 1;
    /**
        Two scores are condsidered similar if the relative "error" of them is
        smaller than defaulMaxRelativeDiff.

        The relative error is defined as:

            relError(a, b) = abs(a - b)/max(a - b)

        **Implementation note:** all computations are done in integer
        arithmetic hence AlignmentChain.maxScore corresponds to 1 in the above
        equation.
    */
    static immutable maxRelativeDiff = AlignmentChain.maxScore / 20; // 5% of larger value
    /**
        Two scores are condsidered similar if the absolute "error" of them is
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
    alias MatchedAlignmentChains = AlignmentChain[];
    alias PileUp = MatchedAlignmentChains[];
    alias Hit = Tuple!(MatchedAlignmentChains, "alignments", string, "dbFile");

    AlignmentContainer!(AlignmentChain[]) selfAlignment;
    AlignmentContainer!(AlignmentChain[]) readsAlignment;
    const Options options;
    /// Set of read ids not to be considered in further processing.
    size_t[] catUseless;
    /// Set of alignments to be considered in further processing.
    AlignmentContainer!(AlignmentChain[]) catCandidates;
    /// Set of alignments to be used for gap filling.
    Hit[] catHits;
    CoordinateTransform coordTransform;
    size_t iteration;

    this(in ref Options options)
    {
        this.options = options;
        this.catUseless = [];
        this.catHits = [];
        this.iteration = 0;
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
            // (a) The aligned read, ie. the alignment extened with the exceeding
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

        auto groupedByGap = groupByGap(catCandidates);
        PileUp[] extendingReadsByContig = groupedByGap.extendingReadsByContig;
        PileUp[] spanningReadsByGap = groupedByGap.spanningReadsByGap;

        logFillingInfo!("findHits", "extension", "raw")(extendingReadsByContig);
        logFillingInfo!("findHits", "span", "raw")(spanningReadsByGap);

        alias getReadId = (reads) => reads[0].contigA.id;

        foreach (gap; spanningReadsByGap)
        {
            size_t refContig1Id = gap[0][0].contigB.id;
            size_t refContig2Id = gap[0][1].contigB.id;
            size_t[] readIds = gap.map!getReadId.array;
            string[] fastaEntries = getFastaEntries(options.readsDb, readIds, options).array;
            string pileupDb = buildDamFile(fastaEntries, options);
            string[] consensusDbs = getConsensus(pileupDb, options);

            // dfmt off
            logJsonDebug(
                "consensusDbs", consensusDbs.map!Json.array,
                "contigIds", [Json(refContig1Id), Json(refContig2Id)],
            );
            // dfmt on

            foreach (consensusDb; consensusDbs)
            {
                auto consensusAlignments = getMappings(options.refDb, consensusDb, options);
                auto consensusGroupedByGap = groupByGap(consensusAlignments);
                PileUp consensusPileUp = consensusGroupedByGap.spanningReadsByGap[0];
                consensusPileUp.sort!(Comparator!gapScore.gt);

                logFillingInfo!("findHits", "span", "consensus")([consensusPileUp]);
                // dfmt off
                logJsonDebug(
                    "contigIds", [Json(refContig1Id), Json(refContig2Id)],
                    "pileUp", consensusPileUp
                        .map!(spanningReads => Json([
                            "estimateSize": Json(estimateSize(spanningReads)),
                            "gapScore": Json(gapScore(spanningReads)),
                            "alignmentScores": Json([
                                Json(spanningReads[0].score),
                                Json(spanningReads[1].score),
                            ]),
                            "readId": Json(spanningReads[0].contigA.id),
                        ]))
                        .array,
                );
                // dfmt on

                // Mark best scoring read for later sequence insertion.
                catHits ~= Hit(consensusPileUp[0], consensusDb);
            }
        }

        logJsonDiagnostic("state", "exit", "function", "djunctor.findHits");

        return this;
    }

    protected Tuple!(PileUp[], "extendingReadsByContig", PileUp[], "spanningReadsByGap") groupByGap(
            AlignmentContainer!(AlignmentChain[]) candidates)
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

        // dfmt off
        auto readsByNumFlanks = candidates
            .b2a
            .chunkBy!isSameRead
            .map!array
            .array
            .sort!"a.length < b.length"
            .groupBy;
        auto readsByNumFlanksByGap = readsByNumFlanks
            .map!(group => group
                .array
                .sort!orderByReferenceIds
                .groupBy);
        // dfmt on
        PileUp[] extendingReadsByContig;
        PileUp[] spanningReadsByGap;

        if (!readsByNumFlanksByGap.empty && numContigsInvolved(readsByNumFlanksByGap) == 1)
        {
            extendingReadsByContig = readsByNumFlanksByGap.front.map!array.array;
            readsByNumFlanksByGap.popFront();
        }

        if (!readsByNumFlanksByGap.empty && numContigsInvolved(readsByNumFlanksByGap) == 2)
        {
            spanningReadsByGap = readsByNumFlanksByGap.front.map!array.array;
            readsByNumFlanksByGap.popFront();
        }

        return typeof(return)(extendingReadsByContig, spanningReadsByGap);
    }

    static protected long estimateSize(in MatchedAlignmentChains alignments) pure
    {
        switch (alignments.length)
        {
        case 1:
            return extensionSize(alignments[0]);
        case 2:
            return spanningGapSize(alignments[0], alignments[1]);
        default:
            throw new Exception(
                    format!"cannot estimate size for %d matched alignments"(alignments.length));
        }
    }

    /**
        Returns the (approximate) size of the gap spanned by alignmentsRange.
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
    static protected long spanningGapSize(in AlignmentChain alignment1, in AlignmentChain alignment2) pure
    {
        // dfmt off
        auto alignments = isBefore!"contigA"(alignment1, alignment2)
            ? tuple(alignment1, alignment2)
            : tuple(alignment2, alignment1);
        // dfmt on
        auto firstAlignment = alignments[0];
        auto secondAlignment = alignments[1];

        assert(secondAlignment.first.contigA.begin > firstAlignment.last.contigA.end,
                format!"intersecting local alignments in %s"(tuple(firstAlignment,
                    secondAlignment).to!string));
        long readsSequenceLength = secondAlignment.first.contigA.begin
            - firstAlignment.last.contigA.end;
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
                            AlignmentChain(0, Contig(1, 80), Contig(1, 50), no, [
                                LocalAlignment(Locus(10, 15), Locus(30, 35), 0),
                                LocalAlignment(Locus(20, 30), Locus(40, 50), 1),
                            ]),
                            AlignmentChain(1, Contig(1, 80), Contig(2, 50), no, [
                                LocalAlignment(Locus(50, 55), Locus(0, 10), 2),
                                LocalAlignment(Locus(60, 70), Locus(15, 20), 3),
                            ]),
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
                            AlignmentChain(2, Contig(1, 80), Contig(1, 50), no, [
                                LocalAlignment(Locus(10, 15), Locus(30, 35), 4),
                                LocalAlignment(Locus(20, 30), Locus(40, 45), 5),
                            ]),
                            AlignmentChain(3, Contig(1, 80), Contig(2, 50), no, [
                                LocalAlignment(Locus(50, 55), Locus(5, 10), 6),
                                LocalAlignment(Locus(60, 70), Locus(15, 20), 7),
                            ]),
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
                            AlignmentChain(4, Contig(1, 80), Contig(2, 50), no, [
                                LocalAlignment(Locus(50, 55), Locus(5, 10), 6),
                                LocalAlignment(Locus(60, 70), Locus(15, 20), 7),
                            ]),
                            AlignmentChain(5, Contig(1, 80), Contig(1, 50), no, [
                                LocalAlignment(Locus(10, 15), Locus(30, 35), 4),
                                LocalAlignment(Locus(20, 30), Locus(40, 45), 5),
                            ]),
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
                            AlignmentChain(6, Contig(1, 80), Contig(1, 50), yes, [
                                LocalAlignment(Locus(10, 15), Locus(30, 35), 8),
                                LocalAlignment(Locus(20, 30), Locus(40, 50), 9),
                            ]),
                            AlignmentChain(7, Contig(1, 80), Contig(2, 50), yes, [
                                LocalAlignment(Locus(50, 55), Locus(0, 10), 10),
                                LocalAlignment(Locus(60, 70), Locus(15, 20), 11),
                            ]),
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
                            AlignmentChain(8, Contig(1, 80), Contig(1, 50), yes, [
                                LocalAlignment(Locus(10, 15), Locus(30, 35), 12),
                                LocalAlignment(Locus(20, 30), Locus(40, 45), 13),
                            ]),
                            AlignmentChain(9, Contig(1, 80), Contig(2, 50), yes, [
                                LocalAlignment(Locus(50, 55), Locus(5, 10), 14),
                                LocalAlignment(Locus(60, 70), Locus(15, 20), 15),
                            ]),
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
                            AlignmentChain(10, Contig(1, 80), Contig(2, 50), yes, [
                                LocalAlignment(Locus(50, 55), Locus(5, 10), 14),
                                LocalAlignment(Locus(60, 70), Locus(15, 20), 15),
                            ]),
                            AlignmentChain(11, Contig(1, 80), Contig(1, 50), yes, [
                                LocalAlignment(Locus(10, 15), Locus(30, 35), 12),
                                LocalAlignment(Locus(20, 30), Locus(40, 45), 13),
                            ]),
                            10
                        ),
                    ];
                    // dfmt on

                    foreach (testCaseIdx, testCase; testCases.enumerate)
                    {
                        auto gotValue = spanningGapSize(testCase[0], testCase[1]);
                        auto expValue = testCase[2];
                        auto errorMessage = format!"expected spanning gap size %d but got %d for test case %d"(
                                expValue, gotValue, testCaseIdx);

                        assert(gotValue == expValue, errorMessage);
                    }
                }
    }

    /**
        Returns the (approximate) size of the extension constituted by alignmentsRange.
        ---
        extensionSize ~= readsSequenceLength - referenceExcess

        CASE 1 (right extension, no complement):

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

        CASE 2 (left extension, no complement):

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

        CASE 3 (right extension, complement):

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

        CASE 4 (left extension, complement):

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
    static protected long extensionSize(in AlignmentChain alignment) pure
    {
        // CASE 1/3 (right extension)
        long readsSequenceLengthRightExtension = alignment.contigA.length
            - alignment.last.contigA.end;
        long referenceExcessRightExtension = alignment.contigB.length - alignment.last.contigB.end;
        // CASE 2/4 (left extension)
        long readsSequenceLengthLeftExtension = alignment.first.contigA.begin;
        long referenceExcessLeftExtension = alignment.first.contigB.begin;

        // dfmt off
        return max(
            // Case 1/3 (right extension)
            readsSequenceLengthRightExtension - referenceExcessRightExtension,
            // Case 2/4 (left extension)
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
                            AlignmentChain(0, Contig(1, 50), Contig(1, 50), no, [
                                LocalAlignment(Locus(0, 15), Locus(10, 20), 0),
                                LocalAlignment(Locus(20, 40), Locus(30, 50), 1),
                            ]),
                            10
                        ),
                        //            0      10  45  50
                        // reference  |-------+---+---|
                        //                     \ \ \
                        // read         |-->-->-+->-+->-->--|
                        //              0       0  40      50
                        tuple(
                            AlignmentChain(1, Contig(1, 50), Contig(1, 50), no, [
                                LocalAlignment(Locus(0, 15), Locus(10, 20), 2),
                                LocalAlignment(Locus(20, 40), Locus(30, 45), 3),
                            ]),
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
                            AlignmentChain(2, Contig(1, 50), Contig(1, 50), no, [
                                LocalAlignment(Locus(10, 15), Locus(0, 20), 4),
                                LocalAlignment(Locus(20, 50), Locus(30, 40), 5),
                            ]),
                            10
                        ),
                        //                  0   5  40      50
                        // reference        |---+---+-------|
                        //                     / / /
                        // read       |-->-->-+->-+->-->--|
                        //            0      10  50      50
                        tuple(
                            AlignmentChain(3, Contig(1, 50), Contig(1, 50), no, [
                                LocalAlignment(Locus(10, 15), Locus(5, 20), 6),
                                LocalAlignment(Locus(20, 50), Locus(30, 40), 7),
                            ]),
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
                            AlignmentChain(4, Contig(1, 50), Contig(1, 50), yes, [
                                LocalAlignment(Locus(0, 15), Locus(10, 20), 8),
                                LocalAlignment(Locus(20, 40), Locus(30, 50), 9),
                            ]),
                            10
                        ),
                        //            0      10  45  50
                        // reference  |-------+---+---|
                        //                     \ \ \
                        // read         |--<--<-+-<-+-<--<--|
                        //              0       0   50      50
                        tuple(
                            AlignmentChain(5, Contig(1, 50), Contig(1, 50), yes, [
                                LocalAlignment(Locus(0, 15), Locus(10, 20), 10),
                                LocalAlignment(Locus(20, 40), Locus(30, 45), 11),
                            ]),
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
                            AlignmentChain(6, Contig(1, 50), Contig(1, 50), yes, [
                                LocalAlignment(Locus(10, 15), Locus(0, 20), 12),
                                LocalAlignment(Locus(20, 50), Locus(30, 40), 13),
                            ]),
                            10
                        ),
                        //                  0   5  40      50
                        // reference        |---+---+-------|
                        //                     / / /
                        // read       |--<--<-+-<-+-<--<--|
                        //            0       10  50      50
                        tuple(
                            AlignmentChain(7, Contig(1, 50), Contig(1, 50), yes, [
                                LocalAlignment(Locus(10, 15), Locus(5, 20), 14),
                                LocalAlignment(Locus(20, 50), Locus(30, 40), 15),
                            ]),
                            5
                        ),
                    ];
                    // dfmt on

                    foreach (testCaseIdx, testCase; testCases.enumerate)
                    {
                        auto gotValue = extensionSize(testCase[0]);
                        auto expValue = testCase[1];
                        auto errorMessage = format!"expected extension size %d but got %d for test case %d"(
                                expValue, gotValue, testCaseIdx);

                        assert(gotValue == expValue, errorMessage);
                    }
                }
    }

    static protected size_t gapScore(in MatchedAlignmentChains alignments) pure
    {
        assert(alignments.length == 2);
        auto score1 = alignments[0].score.to!long;
        auto score2 = alignments[1].score.to!long;

        return (score1 + score2) / 2 - abs(score1 - score2);
    }

    protected DJunctor insertHits()
    {
        logJsonDiagnostic("state", "enter", "function", "djunctor.insertHits");
        foreach (hit; catHits)
        {
            switch (hit.length)
            {
            case 1:
                extendContig(hit);
                break;
            case 2:
                fillGap(hit);
                break;
            default:
                throw new Exception(format!"cannot fill gap for %d alignment chains"(hit.length));
            }
        }
        // Clear `catHits` for next iteration.
        // dfmt off
        logJsonDiagnostic(
            "step", "djunctor.insertHits",
            "coordTransformInsertions", coordTransform.insertions[].map!(serializeToJson!(CoordinateTransform.Insertion)).array,
            "coordTransformPython", coordTransform.toString(),
        );
        // dfmt off
        catHits.length = 0;
        logJsonDiagnostic("state", "exit", "function", "djunctor.insertHits");

        return this;
    }

    protected DJunctor extendContig(in Hit hit)
    {
        assert(0, "unimplemented");
        //logJsonDiagnostic("state", "enter", "function", "djunctor.extendContig");
        //logJsonDiagnostic("state", "exit", "function", "djunctor.extendContig");

        //return this;
    }

    protected DJunctor fillGap(in Hit hit)
    {
        logJsonDiagnostic("state", "enter", "function", "djunctor.fillGap");
        size_t refContig1Id = hit.alignments[0].contigB.id;
        size_t refContig2Id = hit.alignments[1].contigB.id;
        auto complement = hit.alignments[0].complement;
        size_t readId = hit.alignments[0].contigA.id;
        auto read = getFastaEntries(hit.dbFile, [readId], options).front.parseFastaRecord;
        auto gapSequenceSlice = getGapSequenceSlice(hit);
        // dfmt off
        auto fastaSequence = complement
            ? read.reverseComplement[gapSequenceSlice]
            : read[gapSequenceSlice];
        // dfmt on
        // dfmt off
        logJsonDebug(
            "step", "insertHits",
            "type", "span",
            "contigIds", [Json(refContig1Id), Json(refContig2Id)],
            "readId", readId,
            "complement", complement.to!bool,
            "fastaSequence", fastaSequence.array.to!string,
        );
        // dfmt on
        with (CoordinateTransform)
        {
            size_t gapRefBegin;
            size_t gapRefEnd;

            if (hit.alignments[0].isBefore!"contigA"(hit.alignments[1]))
            {
                gapRefBegin = hit.alignments[0].last.contigB.end;
                gapRefEnd = hit.alignments[1].first.contigB.begin;
            }
            else
            {
                gapRefBegin = hit.alignments[1].last.contigB.end;
                gapRefEnd = hit.alignments[0].first.contigB.begin;
            }

            // dfmt off
            coordTransform.add(
                Coordinate(refContig1Id, gapRefBegin),
                Coordinate(refContig2Id, gapRefEnd),
                gapSequenceSlice[1] - gapSequenceSlice[0],
            );
            // dfmt on
        }

        logJsonDiagnostic("state", "exit", "function", "djunctor.fillGap");

        return this;
    }

    protected static Tuple!(int, int) getGapSequenceSlice(in Hit hit) pure
    {
        auto alignment1 = hit.alignments[0];
        auto alignment2 = hit.alignments[1];
        size_t begin;
        size_t end;

        if (alignment1.isBefore!"contigA"(alignment2))
        {
            begin = alignment1.last.contigA.end;
            end = alignment2.first.contigA.begin;
        }
        else
        {
            begin = alignment2.last.contigA.end;
            end = alignment1.first.contigA.begin;
        }

        return typeof(return)(begin.to!int, end.to!int);
    }

    protected DJunctor finish()
    {
        logJsonDiagnostic("state", "enter", "function", "djunctor.finish");
        logJsonDiagnostic("state", "exit", "function", "djunctor.finish");

        return this;
    }

    private void logFillingInfo(string step, string type, string readState)(in PileUp[] readsByContig)
            if (type == "extension" || type == "span")
    {
        immutable lengthFilter = type == "extension" ? "a.length == 1" : "a.length == 2";

        // dfmt off
        logJsonDiagnostic(
            "step", step,
            "type", type,
            "readState", readState,
            "numGaps", readsByContig.length,
            "estimateLengths", readsByContig
                .map!(readsByGap => readsByGap
                    .filter!lengthFilter
                    .map!estimateSize
                    .array
                    .mean)
                .map!Json
                .array,
            "numReads", readsByContig
                .map!"a.length"
                .map!Json
                .array,
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

    alias Insertion = Tuple!(Coordinate, "begin", Coordinate, "end", size_t, "length");

    static enum ExportLanguage
    {
        python,
    }

    const(Insertion)[] insertions;

    /**
        Transform a coordinate accordin to insertions. If the designated
        contig is not present in insertions then the original coordinate is
        returned. Otherwise the new coordinate is computed as follows:

        1. Find insertion where end contig matches input. Note: if begin
           contig matches input not action is needed becuase neither the
           `idx` nor the `contigId` is affected.
        2. Travel from that point backwards to the beginning of the chain –
           ie. until begin contig and end contig of two insertions are
           unequal – summing all the coordinate shifts. The new reference
           contig is the last one in the (reverse) chain.

        For a single insertion and input coordinate `x` the shifted
        coordinate `T(x)` is calculated as follows:

        ---
          begin contig       end contig

        0         ib  lb 0  ie  x      le
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
        coordTransform.add(
            Coordinate(1, 495),
            Coordinate(2, 5),
            100,
        );
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
    }

    unittest
    {
        CoordinateTransform coordTransform;
        // dfmt off
        coordTransform.add(
            Coordinate(100, 0),
            Coordinate(101, 0),
            0,
        );
        coordTransform.add(
            Coordinate(1, 495),
            Coordinate(2, 5),
            100,
        );
        coordTransform.add(
            Coordinate(2, 490),
            Coordinate(5, 10),
            150,
        );
        coordTransform.add(
            Coordinate(5, 480),
            Coordinate(3, 20),
            200,
        );
        coordTransform.add(
            Coordinate(102, 0),
            Coordinate(103, 0),
            0,
        );
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
    }

    /**
        Add an insertion to this transform.

        Returns: this tranform.
        See_Also: `CoordinateTransform.transform`.
    */
    CoordinateTransform add(in size_t beginContigId, in size_t beginIdx,
            in size_t endContigId, in size_t endIdx, size_t insertionLength)
    {
        return add(Coordinate(beginContigId, beginIdx), Coordinate(endContigId,
                endIdx), insertionLength);
    }

    /// ditto
    CoordinateTransform add(in Coordinate begin, in Coordinate end, size_t insertionLength)
    {
        return add(Insertion(begin, end, insertionLength));
    }

    /// ditto
    CoordinateTransform add(in Insertion newInsertion)
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

    unittest
    {
        Insertion getDummyInsertion(in size_t beginContigId, in size_t endContigId)
        {
            return Insertion(Coordinate(beginContigId, 0), Coordinate(endContigId, 0), 0);
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
            // Case 1: newInsertion fits between too existing insertions
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
            // Case 2*: newInsertion extends an existing insertion chain in the front
            // dfmt off
            coordTransform.add(getDummyInsertion(42, 1));
            assert(coordTransform.insertions == [
                getDummyInsertion(42, 1),
                getDummyInsertion(1, 2),
                getDummyInsertion(3, 4),
                getDummyInsertion(9, 8),
                getDummyInsertion(6, 5),
            ]);
            // Case 2: newInsertion extends an existing insertion chain in the front
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
            // Case 3: newInsertion extends an existing insertion chain in the front
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
            // Case 4: open new insertion chain
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
            CoordinateTransform emptyCoordTransform;

            // Case 4*: open new insertion chain
            assert(emptyCoordTransform.add(getDummyInsertion(1337, 42))
                    .insertions == [getDummyInsertion(1337, 42)]);
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
        // dfmt on
        {
            return 0;
        }
        // Case 3*: newInsertion extends an existing insertion chain in the back
        // Case 4**: open new insertion chain
        else if (insertions.length == 1)
        {
            // Note: this case is needed in order for `slice(2)` to always yield
            //       pairs (of size two).
            return 1;
        }

        foreach (i, insAB; insertions.slide(2).enumerate)
        {
            auto insA = insAB[0];
            auto insB = insAB[1];

            // dfmt off
            if (any(only(
                // Case 1: newInsertion fits between too existing insertions
                newInsertion.begin.contigId == insA.end.contigId && newInsertion.end.contigId == insB.begin.contigId,
                // Case 2: newInsertion extends an existing insertion chain in the front
                newInsertion.begin.contigId != insA.end.contigId && newInsertion.end.contigId == insB.begin.contigId,
                // Case 3: newInsertion extends an existing insertion chain in the back
                newInsertion.begin.contigId == insA.end.contigId && newInsertion.end.contigId != insB.begin.contigId,
            )))
            {
                return i + 1;
            }
            // dfmt on
        }

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
            // Case 1: newInsertion fits between too existing insertions
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
