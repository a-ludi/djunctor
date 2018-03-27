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
import std.algorithm : any, canFind, chunkBy, each, filter, group, isSorted,
    map, sort, swap;
import std.array : appender, array;
import std.container : BinaryHeap, heapify, make;
import std.conv;
import std.format : format;
import std.range : ElementType, isForwardRange, only, walkLength;
import std.stdio : writeln;
import std.typecons : tuple, Tuple;

/// General container for alignment data.
template AlignmentContainer(R)
{
    import std.typecons : tuple;

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
            import core.exception : AssertError;
            import std.exception : assertThrown;

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
            import core.exception : AssertError;
            import std.exception : assertThrown;

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

    Contig contigA;
    Contig contigB;
    Complement complement;
    LocalAlignment[] localAlignments;

    invariant
    {
        import std.algorithm : all;

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
        import core.exception : AssertError;
        import std.exception : assertNotThrown, assertThrown;

        with (Complement) with (LocalAlignment)
            {
                auto acZeroLength = AlignmentChain(Contig(1, 10), Contig(1, 10), no, []);
                auto ac1 = AlignmentChain(Contig(1, 10), Contig(1, 10), no,
                        [LocalAlignment(Locus(1, 1), Locus(1, 9), 0)]);
                auto ac2 = AlignmentChain(Contig(1, 10), Contig(1, 10), no,
                        [LocalAlignment(Locus(1, 11), Locus(1, 9), 0)]);
                auto ac3 = AlignmentChain(Contig(1, 10), Contig(1, 10), no,
                        [LocalAlignment(Locus(1, 9), Locus(1, 1), 0)]);
                auto ac4 = AlignmentChain(Contig(1, 10), Contig(1, 10), no,
                        [LocalAlignment(Locus(1, 9), Locus(1, 11), 0)]);
                auto acFine = AlignmentChain(Contig(1, 10), Contig(1, 10), no,
                        [LocalAlignment(Locus(1, 9), Locus(1, 9), 0)]);

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
                auto ac = AlignmentChain(Contig(1, 10), Contig(1, 10), no, [firstLA, otherLA]);

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
                auto ac = AlignmentChain(Contig(1, 10), Contig(1, 10), no, [otherLA, lastLA]);

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
                auto ac = AlignmentChain(Contig(1, 50), Contig(1, 15), no, [la]);

                // read with extension align an contigA from 25 to 40
                assert(ac.isFullyContained);
            }

        with (Complement) with (LocalAlignment)
            {
                auto la1 = LocalAlignment(Locus(10, 20), Locus(5, 10), 1);
                auto la2 = LocalAlignment(Locus(30, 40), Locus(5, 10), 1);
                auto ac = AlignmentChain(Contig(1, 50), Contig(1, 15), no, [la1, la2]);

                // read with extension align an contigA from 5 to 45
                assert(ac.isFullyContained);
            }

        with (Complement) with (LocalAlignment)
            {
                auto la = LocalAlignment(Locus(0, 10), Locus(5, 10), 1);
                auto ac = AlignmentChain(Contig(1, 50), Contig(1, 15), no, [la]);

                // read with extension align an contigA from -5 to 15
                assert(!ac.isFullyContained);
            }

        with (Complement) with (LocalAlignment)
            {
                auto la = LocalAlignment(Locus(40, 50), Locus(5, 10), 1);
                auto ac = AlignmentChain(Contig(1, 50), Contig(1, 15), no, [la]);

                // read with extension align an contigA from 35 to 55
                assert(!ac.isFullyContained);
            }

        with (Complement) with (LocalAlignment)
            {
                auto la1 = LocalAlignment(Locus(0, 20), Locus(5, 10), 1);
                auto la2 = LocalAlignment(Locus(30, 50), Locus(5, 10), 1);
                auto ac = AlignmentChain(Contig(1, 50), Contig(1, 15), no, [la1, la2]);

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
                auto ac = AlignmentChain(Contig(1, 10), Contig(1, 10), no, [la1, la2]);

                assert(ac.totalLength == 9);
            }
    }

    @property size_t totalDiffs()
    {
        import std.algorithm : map, sum;

        return localAlignments.map!"a.numDiffs".sum;
    }

    unittest
    {
        with (Complement) with (LocalAlignment)
            {
                auto la1 = LocalAlignment(Locus(1, 3), Locus(1, 3), 1);
                auto la2 = LocalAlignment(Locus(5, 10), Locus(5, 10), 2);
                auto ac = AlignmentChain(Contig(1, 10), Contig(1, 10), no, [la1, la2]);

                assert(ac.totalDiffs == 3);
            }
    }

    @property size_t totalGapLength()
    {
        import std.algorithm : map, sum;
        import std.range : chunks;

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
                auto ac = AlignmentChain(Contig(1, 10), Contig(1, 10), no, [la1, la2]);

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
                auto ac = AlignmentChain(Contig(1, 10), Contig(1, 10), no, [la1, la2]);

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
                auto ac = AlignmentChain(Contig(1, 10), Contig(1, 10), no, [la1, la2]);

                assert(ac.score == 4 * maxScore / 9);
            }
    }

    int compareIds(ref const AlignmentChain other) const pure nothrow
    {
        import std.math : sgn;

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
                    AlignmentChain(Contig(1, 10), Contig(1, 10), no, [la]),
                    AlignmentChain(Contig(1, 10), Contig(2, 10), no, [la]),
                    AlignmentChain(Contig(2, 10), Contig(1, 10), no, [la]),
                    AlignmentChain(Contig(2, 10), Contig(2, 10), no, [la]),
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
        import std.math : sgn;

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
                    AlignmentChain(Contig(1, 10), Contig(1, 10), no, [la]),
                    AlignmentChain(Contig(1, 10), Contig(2, 10), no, [la]),
                    AlignmentChain(Contig(2, 10), Contig(1, 10), no, [la]),
                    AlignmentChain(Contig(2, 10), Contig(2, 10), no, [la]),
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
                    AlignmentChain(Contig(1, 10), Contig(1, 10), no, [
                        LocalAlignment(Locus(0, 2), Locus(0, 2), 1),
                        LocalAlignment(Locus(3, 5), Locus(3, 5), 1)
                    ]),
                    AlignmentChain(Contig(1, 10), Contig(1, 10), no, [
                        LocalAlignment(Locus(0, 2), Locus(1, 2), 1),
                        LocalAlignment(Locus(3, 5), Locus(3, 5), 1)
                    ]),
                    AlignmentChain(Contig(1, 10), Contig(1, 10), no, [
                        LocalAlignment(Locus(1, 2), Locus(0, 2), 1),
                        LocalAlignment(Locus(3, 5), Locus(3, 5), 1)
                    ]),
                    AlignmentChain(Contig(1, 10), Contig(1, 10), no, [
                        LocalAlignment(Locus(1, 2), Locus(1, 2), 1),
                        LocalAlignment(Locus(3, 5), Locus(3, 5), 1)
                    ]),
                    AlignmentChain(Contig(1, 10), Contig(1, 10), no, [
                        LocalAlignment(Locus(1, 2), Locus(1, 2), 1),
                        LocalAlignment(Locus(3, 5), Locus(3, 6), 1)
                    ]),
                    AlignmentChain(Contig(1, 10), Contig(1, 10), no, [
                        LocalAlignment(Locus(1, 2), Locus(1, 2), 1),
                        LocalAlignment(Locus(3, 6), Locus(3, 5), 1)
                    ]),
                    AlignmentChain(Contig(1, 10), Contig(1, 10), no, [
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
                    AlignmentChain(Contig(1, 10), Contig(1, 10), no, [la]),
                    AlignmentChain(Contig(1, 10), Contig(2, 10), no, [la]),
                    AlignmentChain(Contig(2, 10), Contig(1, 10), no, [la]),
                    AlignmentChain(Contig(2, 10), Contig(2, 10), no, [la]),
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
    import std.algorithm : chunkBy, equal;
    import std.array : array;

    with (AlignmentChain) with (Complement) with (LocalAlignment)
            {
                // dfmt off
                auto sortedTestChains = [
                    AlignmentChain(Contig(1, 10), Contig(1, 10), no, [LocalAlignment(Locus(0, 10), Locus(0, 1), 0)]),
                    AlignmentChain(Contig(1, 10), Contig(2, 20), no, [LocalAlignment(Locus(0, 10), Locus(0, 2), 0)]),
                    AlignmentChain(Contig(1, 10), Contig(3, 30), no, [LocalAlignment(Locus(0, 10), Locus(0, 3), 0)]),
                    AlignmentChain(Contig(2, 20), Contig(1, 10), no, [LocalAlignment(Locus(0, 20), Locus(0, 4), 0)]),
                    AlignmentChain(Contig(2, 20), Contig(1, 10), no, [LocalAlignment(Locus(0, 20), Locus(0, 5), 0)]),
                    AlignmentChain(Contig(2, 20), Contig(3, 30), no, [LocalAlignment(Locus(0, 20), Locus(0, 6), 0)]),
                    AlignmentChain(Contig(3, 30), Contig(1, 10), no, [LocalAlignment(Locus(0, 30), Locus(0, 7), 0)]),
                    AlignmentChain(Contig(3, 30), Contig(2, 20), no, [LocalAlignment(Locus(0, 30), Locus(0, 8), 0)]),
                    AlignmentChain(Contig(3, 30), Contig(3, 30), no, [LocalAlignment(Locus(0, 30), Locus(0, 9), 0)]),
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
    import std.range : assumeSorted, SortedRange;

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
    import std.algorithm : equal;
    import std.array : array;

    with (AlignmentChain) with (Complement) with (LocalAlignment)
            {
                // dfmt off
                auto sortedTestChains = [
                    AlignmentChain(Contig(1, 10), Contig(1, 10), no, [LocalAlignment(Locus(0, 1), Locus(0, 1), 0)]),
                    AlignmentChain(Contig(1, 10), Contig(2, 20), no, [LocalAlignment(Locus(0, 2), Locus(0, 2), 0)]),
                    AlignmentChain(Contig(1, 10), Contig(3, 30), no, [LocalAlignment(Locus(0, 3), Locus(0, 3), 0)]),
                    AlignmentChain(Contig(2, 20), Contig(1, 10), no, [LocalAlignment(Locus(0, 4), Locus(0, 4), 0)]),
                    AlignmentChain(Contig(2, 20), Contig(1, 10), no, [LocalAlignment(Locus(0, 5), Locus(0, 5), 0)]),
                    AlignmentChain(Contig(2, 20), Contig(3, 30), no, [LocalAlignment(Locus(0, 6), Locus(0, 6), 0)]),
                    AlignmentChain(Contig(3, 30), Contig(1, 10), no, [LocalAlignment(Locus(0, 7), Locus(0, 7), 0)]),
                    AlignmentChain(Contig(3, 30), Contig(2, 20), no, [LocalAlignment(Locus(0, 8), Locus(0, 8), 0)]),
                    AlignmentChain(Contig(3, 30), Contig(3, 30), no, [LocalAlignment(Locus(0, 9), Locus(0, 9), 0)]),
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
                    AlignmentChain(Contig(1, 10), Contig(1, 10), no, [LocalAlignment(Locus(1, 6), Locus(0, 1), 0)]),
                    AlignmentChain(Contig(1, 10), Contig(2, 10), yes, [LocalAlignment(Locus(2, 6), Locus(0, 1), 0)]),
                    AlignmentChain(Contig(1, 10), Contig(3, 10), no, [LocalAlignment(Locus(3, 6), Locus(0, 1), 0)]),
                    AlignmentChain(Contig(1, 10), Contig(4, 10), yes, [LocalAlignment(Locus(4, 6), Locus(0, 1), 0)]),
                    AlignmentChain(Contig(1, 10), Contig(5, 10), no, [LocalAlignment(Locus(5, 6), Locus(0, 1), 0)]),
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
    alias numExpectedUseless = (numCatUnsure, iteration) => (numCatUnsure - numCatUnsure / 20) / (
            iteration + 1);
    AlignmentContainer!(AlignmentChain[]) selfAlignment;
    AlignmentContainer!(AlignmentChain[]) readsAlignment;
    const Options options;
    /// Set of read ids not to be considered in further processing.
    size_t[] catUseless;
    /// Set of alignments to be considered in further processing.
    AlignmentChain[] catUnsure;
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
        logInfo("BEGIN djunctor");
        // dfmt off
        this
            .init
            .mainLoop
            .finish;
        // dfmt on
        logInfo("END djunctor");

        return this;
    }

    protected DJunctor init()
    {
        import djunctor.dazzler : getLocalAlignments, getMappings;

        logInfo("BEGIN djunctor.init");
        selfAlignment = getLocalAlignments(options.refDb, options);
        readsAlignment = getMappings(options.refDb, options.readsDb, options);
        catUnsure = readsAlignment.a2b.dup;
        logInfo("END djunctor.init");

        return this;
    }

    protected DJunctor mainLoop()
    {
        logInfo("BEGIN djunctor.mainLoop");
        do
        {
            filterUseless();
            findHits();

            logDiagnostic("i: %d, useless: %d, unsure: %d, hits: %d",
                    iteration, catUseless.length, catUnsure.length, catHits.length);

            if (catHits.length > 0)
            {
                fillGaps();
            }

            ++iteration;
        }
        while (catHits.length > 0 && iteration < maxLoops);
        logInfo("END djunctor.mainLoop");

        return this;
    }

    protected DJunctor filterUseless()
    {
        auto sizeReserve = numExpectedUseless(catUnsure.length, iteration);
        if (sizeReserve == 0)
        {
            logDiagnostic("skipping filterUseless: expected count is zero");
            // Skip this step if we do not expect to find "useless" reads;
            // it is not harmful to consider some "useless" reads as
            // candidates (catUnsure).
            return this;
        }

        auto uselessAcc = appender!(size_t[]);
        uselessAcc.reserve(numExpectedUseless(catUnsure.length, iteration));

        logInfo("BEGIN djunctor.filterUseless");
        foreach (alignmentsChunk; catUnsure.chunkBy!haveEqualIds)
        {
            auto alignments = alignmentsChunk.array;
            // Mark reads as useless if either: (s. issue #3)
            // (a) The aligned read, ie. the alignment extened with the exceeding
            //     read sequence on either end, is fully contained in a single
            //     contig.
            // (b) It aligns in multiple locations of one contig with
            //     approximately equal quality.
            // dfmt off
            if (alignments.any!"a.isFullyContained" ||
                (alignments.length > 1 && similarScore(alignments[0].score, alignments[1].score)))
                // dfmt on
            {
                uselessAcc ~= alignments[0].contigB.id;
            }
        }
        logDebug("useless expected: %d found: %d",
                numExpectedUseless(catUnsure.length, iteration), uselessAcc.data.length);
        this.catUseless ~= uselessAcc.data;
        this.catUnsure = this.catUnsure.filter!(
                ac => !uselessAcc.data.canFind(ac.contigB.id)).array;
        logInfo("END djunctor.filterUseless");

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
        logInfo("BEGIN djunctor.findHits");
        logInfo("END djunctor.findHits");

        return this;
    }

    /**
        Returns the (approximate) size of the gap spanned by alignmentsRange.
        ---
                        contig a                  contig b

                  0         cax cay la      0   cbx cby       lb
        reference |---------+---+---|       |---+---+---------|
                             \ \ \             / / /
                              \ \ \           / / /
                  alignment 1  \ \ \         / / /  alignment 2
                                \ \ \       / / /
        read            |--------+---+-----+---+-------|
                        0        a1x a1y   a2x a2y     lr

        gap_size ~= a2x - a1y - (la - cay) - (cbx - 0)
                  = a2x - a1y - (la - cay + cbx)
        ---
    */
    protected long spanningGapSize(in AlignmentChain alignment1, in AlignmentChain alignment2) pure
    {
        auto alignments = isBefore!"contigA"(alignment1, alignment2) ? tuple(alignment1, alignment2) : tuple(
                alignment2, alignment1);
        auto firstAlignment = alignments[0];
        auto secondAlignment = alignments[1];

        assert(secondAlignment.first.contigA.begin > firstAlignment.last.contigA.end,
                tuple(firstAlignment, secondAlignment).to!string);
        auto readsSequenceLength = secondAlignment.first.contigA.begin
            - firstAlignment.last.contigA.end;
        auto referenceExcess = firstAlignment.contigB.length
            - firstAlignment.last.contigB.end + secondAlignment.first.contigB.begin;

        return readsSequenceLength.to!long - referenceExcess.to!long;
    }

    /**
        Returns the (approximate) size of the extension constituted by alignmentsRange.
        ---
        CASE 1 (right extension):

                       0       rx  ry  lr
            reference  |-------+---+---|
                                \ \ \
                                 \ \ \
                       alignment  \ \ \
                                   \ \ \
            read            |-------+---+-------|
                            0       ax  ay      la

            extension_size ~= la - ay - (lr - ry)


        CASE 2 (left extension):

                                0    rx  ry lr
            reference           |---+---+-------|
                                   / / /
                                  / / /
                      alignment  / / /
                                / / /
            read       |-------+---+-------|
                       0       ax  ay      la

            extension_size ~= ax - 0 - (rx - 0)
                            = ax - rx
        ---
    */
    protected long extensionSize(in AlignmentChain alignment) pure
    {
        if (alignment.last.contigB.end > alignment.contigB.length / 2)
        {
            // CASE 1 (right extension)
            auto readsSequenceLength = alignment.contigA.length - alignment.last.contigA.end;
            auto referenceExcess = alignment.contigB.length - alignment.last.contigB.end;

            return readsSequenceLength.to!long - referenceExcess.to!long;
        }
        else if (alignment.first.contigB.begin < alignment.contigB.length / 2)
        {
            // CASE 2 (left extension)
            return alignment.first.contigA.begin.to!long - alignment.first.contigB.begin.to!long;
        }
        else
        {
            assert(0, "alignment does not seem to extend the reference");
        }
    }

    protected DJunctor fillGaps()
    {
        logInfo("BEGIN djunctor.fillGaps");
        logInfo("END djunctor.fillGaps");

        return this;
    }

    protected DJunctor finish()
    {
        logInfo("BEGIN djunctor.finish");
        logInfo("END djunctor.finish");

        return this;
    }
}
