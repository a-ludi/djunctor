/**
    This is the main algorithm of this package.

    Copyright: © 2018 Arne Ludwig <arne.ludwig@posteo.de>
    License: Subject to the terms of the MIT license, as written in the
             included LICENSE file.
    Authors: Arne Ludwig <arne.ludwig@posteo.de>
*/
module djunctor.djunctor;

import djunctor.commandline : Options;
import djunctor.util.log;
import djunctor.util.math : mean, median;
import core.exception : AssertError;
import std.algorithm : all, any, canFind, chunkBy, each, equal, filter, group,
    isSorted, map, max, sort, sum, swap;
import std.array : appender, array;
import std.container : BinaryHeap, heapify, make;
import std.conv;
import std.exception : assertNotThrown, assertThrown;
import std.format : format;
import std.math : sgn;
import std.range : assumeSorted, chunks, ElementType, enumerate, isForwardRange,
    only, SortedRange, walkLength;
import std.stdio : writeln;
import std.typecons : tuple, Tuple;

version (unittest)
{
    import djunctor.util.testing : MockCallable;
    import djunctor.dazzler : origGetLocalAlignments = getLocalAlignments,
        origGetMappings = getMappings;

    MockCallable!(AlignmentContainer!(AlignmentChain[]), const string, const Options) getLocalAlignments;
    MockCallable!(origGetMappings!Options) getMappings;
}
else
{
    import djunctor.dazzler : getLocalAlignments, getMappings;
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

    @property size_t totalLength()
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

    @property size_t totalDiffs()
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

    @property size_t totalGapLength()
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

    @property size_t numMatchingBps()
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

    @property size_t score()
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
    AlignmentContainer!(AlignmentChain[]) selfAlignment;
    AlignmentContainer!(AlignmentChain[]) readsAlignment;
    const Options options;
    /// Set of read ids not to be considered in further processing.
    size_t[] catUseless;
    /// Set of alignments to be considered in further processing.
    AlignmentContainer!(AlignmentChain[]) catCandidates;
    /// Set of alignments to be used for gap filling.
    BinaryHeap!(AlignmentChain[]) catHits;
    size_t iteration;

    this(in ref Options options)
    {
        this.options = options;
        this.catUseless = [];
        this.catHits = make!(typeof(this.catHits))();
        this.iteration = 0;
    }

    DJunctor run()
    {
        logDiagnostic("BEGIN djunctor");
        // dfmt off
        this
            .init
            .mainLoop
            .finish;
        // dfmt on
        logDiagnostic("END djunctor");

        return this;
    }

    protected DJunctor init()
    {
        logDiagnostic("BEGIN djunctor.init");
        selfAlignment = getLocalAlignments(options.refDb, options);
        readsAlignment = getMappings(options.refDb, options.readsDb, options);
        catCandidates = AlignmentContainer!(AlignmentChain[])(readsAlignment.a2b.dup,
                readsAlignment.b2a.dup);
        logDiagnostic("END djunctor.init");

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
        logDiagnostic("BEGIN djunctor.mainLoop");
        do
        {
            filterUseless();
            findHits();

            logDiagnostic("i: %d, useless: %d, candiates.a2b: %d candiates.b2a: %d, hits: %d",
                    iteration, catUseless.length,
                    catCandidates.a2b.length, catCandidates.b2a.length, catHits.length);

            if (catHits.length > 0)
            {
                fillGaps();
            }

            ++iteration;
        }
        while (catHits.length > 0 && iteration < maxLoops);
        logDiagnostic("END djunctor.mainLoop");

        return this;
    }

    protected DJunctor filterUseless()
    {
        auto sizeReserve = numEstimateUseless(catCandidates.a2b.length, iteration);
        if (sizeReserve == 0)
        {
            logDiagnostic("skipping filterUseless: expected count is zero");
            // Skip this step if we do not expect to find "useless" reads;
            // it is not harmful to consider some "useless" reads as
            // candidates (catCandidates).
            return this;
        }

        auto uselessAcc = appender!(size_t[]);
        auto candidatesAcc = appender!(AlignmentChain[]);
        uselessAcc.reserve(sizeReserve);
        alias isNotUseless = ac => !uselessAcc.data.canFind(ac.contigA.id);

        logDiagnostic("BEGIN djunctor.filterUseless");
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
        logDebug("useless expected: %d found: %d", sizeReserve, uselessAcc.data.length);
        this.catUseless ~= uselessAcc.data;
        this.catCandidates.a2b = candidatesAcc.data;
        this.catCandidates.b2a = this.catCandidates.b2a.filter!isNotUseless.array;
        logDiagnostic("END djunctor.filterUseless");

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

        logDiagnostic("BEGIN djunctor.findHits");
        // dfmt off
        auto readsByNumFlanks = catCandidates
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

        if (!readsByNumFlanksByGap.empty && numContigsInvolved(readsByNumFlanksByGap) == 1)
        {
            auto extendingReadsByGap = readsByNumFlanksByGap.front;
            readsByNumFlanksByGap.popFront();

            logDiagnostic("extending %d contigs", extendingReadsByGap.save.walkLength);

            // dfmt off
            logDebug(
                "extending size (raw averages): %s",
                extendingReadsByGap
                    .save
                    .map!(readsByGap => readsByGap
                        .map!(extendingAlignments => extensionSize(extendingAlignments[0]))
                        .array
                        .mean)
            );
            // dfmt on
        }

        if (!readsByNumFlanksByGap.empty && numContigsInvolved(readsByNumFlanksByGap) == 2)
        {
            auto spanningReadsByGap = readsByNumFlanksByGap.front;
            readsByNumFlanksByGap.popFront();

            logDiagnostic("spanning %d gaps", spanningReadsByGap.save.walkLength);

            // dfmt off
            logDebug(
                "spanning size (raw averages): %s",
                spanningReadsByGap
                    .save
                    .map!(readsByGap => readsByGap
                        .map!(spanningAlignments => spanningGapSize(spanningAlignments[0], spanningAlignments[1]))
                        .array
                        .mean)
            );
            // dfmt on

            /*
                for each gap/pile up:
                    1. determine the range of the reference that is
                       covered by the pile up plus some margin
                    2. get the corresponding subsequence of the reference
                    3. get the list of read sequences
                    4. dalign read sequences to the ref subsequence
                    5. build consensus using daccord
            */
        }

        logDiagnostic("END djunctor.findHits");

        return this;
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

    protected DJunctor fillGaps()
    {
        logDiagnostic("BEGIN djunctor.fillGaps");
        logDiagnostic("END djunctor.fillGaps");

        return this;
    }

    protected DJunctor finish()
    {
        logDiagnostic("BEGIN djunctor.finish");
        logDiagnostic("END djunctor.finish");

        return this;
    }
}
