/**
    This is the main algorithm of this package.

    Copyright: © 2018 Arne Ludwig <arne.ludwig@posteo.de>
    License: Subject to the terms of the MIT license, as written in the
             included LICENSE file.
    Authors: Arne Ludwig <arne.ludwig@posteo.de>
*/
module djunctor.djunctor;

import djunctor.commandline : Options;
import djunctor.log;
import std.conv;
import std.format : format;
import std.stdio : writeln;

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

        long idCompare = other.contigA.id - this.contigA.id;
        if (idCompare != 0)
            return cast(int) sgn(idCompare);

        idCompare = other.contigB.id - this.contigB.id;
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
                            assert(compareValue > 0, errorMessage(">"));
                        else if (i > j)
                            assert(compareValue < 0, errorMessage("<"));
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

        long locusCompare = other.first.contigA.begin - this.first.contigA.begin;
        if (locusCompare != 0)
            return cast(int) sgn(locusCompare);

        locusCompare = other.first.contigB.begin - this.first.contigB.begin;
        if (locusCompare != 0)
            return cast(int) sgn(locusCompare);

        locusCompare = other.last.contigA.end - this.last.contigA.end;
        if (locusCompare != 0)
            return cast(int) sgn(locusCompare);

        locusCompare = other.last.contigB.end - this.last.contigB.end;
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
                            assert(compareValue > 0, errorMessage(">"));
                        else if (i > j)
                            assert(compareValue < 0, errorMessage("<"));
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
                            assert(compareValue > 0, errorMessage(">"));
                        else if (i > j)
                            assert(compareValue < 0, errorMessage("<"));
                        else
                            assert(compareValue == 0, errorMessage("=="));
                    }
            }
    }
}

bool idsPred(in AlignmentChain ac1, in AlignmentChain ac2) pure
{
    auto cmpValue = ac1.compareIds(ac2);

    return 0 != cmpValue && 0 < cmpValue;
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
    import std.algorithm : isSorted;
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

/// Start the `djunctor` alorithm with preprocessed options.
void runWithOptions(in ref Options options)
{
    import djunctor.dazzler : getLocalAlignments, getMappings;

    logInfo("starting");
    writeln(getMappings(options.refDb, options.readsDb, options));
    writeln(getLocalAlignments(options.refDb, options));
}
