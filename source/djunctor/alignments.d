/**
    Everything to handle local alignments and friends.

    Copyright: © 2018 Arne Ludwig <arne.ludwig@posteo.de>
    License: Subject to the terms of the MIT license, as written in the
             included LICENSE file.
    Authors: Arne Ludwig <arne.ludwig@posteo.de>
*/
module djunctor.alignments;

import djunctor.util.log;
import djunctor.util.math : absdiff;
import core.exception : AssertError;
import std.algorithm : all, any, chunkBy, equal, filter, isSorted, map, max,
    mean, min, sort, sum, swap, SwapStrategy;
import std.array : appender, Appender, array;
import std.conv : to;
import std.exception : assertNotThrown, assertThrown;
import std.format : format;
import std.math : sgn;
import std.range : assumeSorted, chunks, enumerate, only, retro;
import std.typecons : Flag, No, tuple, Tuple, Yes;
import vibe.data.json : toJson = serializeToJson;

debug import std.stdio : writeln;

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

    static immutable maxScore = 2 ^^ 16;

    size_t id;
    Contig contigA;
    Contig contigB;
    Complement complement;
    LocalAlignment[] localAlignments;
    size_t tracePointDistance;

    invariant
    {
        assert(localAlignments.length >= 1, "empty chain is forbidden");
        foreach (la; localAlignments)
        {
            assert(0 <= la.contigA.begin && la.contigA.begin < la.contigA.end
                    && la.contigA.end <= contigA.length, "non-sense alignment of contigA");
            assert(0 <= la.contigB.begin && la.contigB.begin < la.contigB.end
                    && la.contigB.end <= contigB.length, "non-sense alignment of contigB");

            assert(tracePointDistance == 0 || la.tracePoints.length > 0, "missing trace points");
            if (tracePointDistance > 0)
            {
                size_t traceLength = la.tracePoints.map!"a.numBasePairs".sum;
                size_t traceDiffs = la.tracePoints.map!"a.numDiffs".sum;

                // TODO remove logging and "soft" assertion if fixed in LAdump
                if (la.numDiffs != traceDiffs)// dfmt off
                    debug logJsonDebug(
                        "contigA", contigA.id,
                        "contigB", contigB.id,
                        "la.numDiffs", la.numDiffs,
                        "traceDiffs", traceDiffs,
                    );
                    // dfmt on
                assert(absdiff(la.numDiffs, traceDiffs) <= 1, "missing trace points");
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

    /// This alignment is called proper iff it starts and ends at a read boundary.
    @property bool isProper() const pure nothrow
    {
        // dfmt off
        return (
            first.contigA.begin == 0 ||
            first.contigB.begin == 0
        )
        &&
        (
            last.contigA.end == contigA.length ||
            last.contigB.end == contigB.length
        );
        // dfmt on
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
    bool isFullyContained() const
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

    @property size_t coveredBases(string contig)() const pure
    {
        return localAlignments.map!("a." ~ contig ~ ".end - a." ~ contig ~ ".begin").sum;
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
    Returns true iff ac1 begins before ac2 in "`contigA`-order".
*/
bool isBefore(string contig)(in AlignmentChain ac1, in AlignmentChain ac2) pure
        if (contig == "contigA" || contig == "contigB")
{
    assert(__traits(getMember, ac1, contig) == __traits(getMember, ac2, contig),
            "alignment chains do not belong to the same contig");
    static if (contig == "contigB")
        assert(ac1.complement == ac2.complement, "alignment chains must both align complemented or not");

    return __traits(getMember, ac1.first, contig).begin < __traits(getMember, ac2.first, contig).begin;
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
                        alias errorMessage = (expValue) => format!"expected isBefore!\"contigA\"(%s, %s) to be %s but got %s"(
                                acs[i], acs[j], expValue, compareValue);

                        if (i < j)
                            assert(compareValue, errorMessage(true));
                        else
                            assert(!compareValue, errorMessage(false));
                    }
            }
}

/// Type of the read alignment.
static enum ReadAlignmentType
{
    front = 0,
    gap = 1,
    back = 2,
}

/**
    Alignment of a read against the reference. This is either one or two
    alignment chains which belong to the same read and one or two reference
    contig(s).
*/
struct ReadAlignment
{
    AlignmentChain[2] _alignments;
    size_t _length;

    this(in AlignmentChain[] alignments...)
    {
        if (1 <= alignments.length && alignments.length <= 2)
        {
            this._length = alignments.length;

            foreach (i, alignment; alignments)
            {
                this._alignments[i] = cast(AlignmentChain) alignment;
            }
        }
        else
        {
            // dfmt off
            logJsonDiagnostic(
                "info", format!"creating invalid read alignment with %d local alignments"(alignments.length),
                "alignments", alignments.toJson,
            );
            // dfmt on

            this._length = 0;
        }
    }

    invariant
    {
        // If one of the alignments is in reverse-complement and the other not
        // then the same part of the read aligns twice which is an illegal
        // spanning alignment.
        assert(_length != 2 || _alignments[0].complement == _alignments[1].complement);
    }

    @property size_t length() pure const nothrow
    {
        return _length;
    }

    @property size_t opDollar() pure const nothrow
    {
        return _length;
    }

    inout(AlignmentChain[]) opIndex() inout pure nothrow
    {
        return _alignments[0 .. _length];
    }

    inout(AlignmentChain) opIndex(T)(T idx) inout pure nothrow
    {
        return this[][idx];
    }

    unittest
    {
        auto ac1 = AlignmentChain(1);
        auto ac2 = AlignmentChain(2);

        auto ra1 = ReadAlignment(ac1);

        assert(ra1.length == 1);
        assert(ra1[] == [ac1]);
        assert(ra1[0] == ac1);
        assert(ra1[$ - 1] == ac1);

        auto ra2 = ReadAlignment(ac1, ac2);

        assert(ra2.length == 2);
        assert(ra2[] == [ac1, ac2]);
        assert(ra2[0] == ac1);
        assert(ra2[1] == ac2);
        assert(ra2[$ - 1] == ac2);
    }

    /**
        If readAlignment is a gap return true iff the first alignment begins
        before the second alignment on the reference (`contigA`); otherwise
        returns true.
    */
    @property bool isInOrder() const pure nothrow
    {
        return !isGap || isBefore!"contigB"(_alignments[0], _alignments[1]);
    }

    void forceInOrder() pure nothrow
    {
        if (!isInOrder && isGap)
        {
            swap(_alignments[0], _alignments[1]);
        }
    }

    /**
        Returns true iff the read alignment is valid, ie. it is either an
        extension or gap.
    */
    @property bool isValid() const pure nothrow
    {
        return isExtension ^ isGap;
    }

    /**
        Get the type of the read alignment.

        See_Also: `isFrontExtension`, `isBackExtension`, `isGap`
    */
    @property ReadAlignmentType type() const pure nothrow
    {
        assert(isValid, "invalid read alignment");

        if (isGap)
        {
            return ReadAlignmentType.gap;
        }
        else if (isFrontExtension)
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
    @property bool isExtension() const pure nothrow
    {
        if (_length != 1)
        {
            return false;
        }

        return _isExtension(_alignments[0]);
    }

    private static bool _isExtension(in AlignmentChain alignment) pure nothrow
    {
        return _isFrontExtension(alignment) ^ _isBackExtension(alignment);
    }

    /**
        Returns true iff the read alignment is an front extension, ie. it is an
        extension and reaches over the front of the reference contig.

        ---
        Case 1 (complement alignment):

                          0  rx
            ref           |--+->-+->-->-->--|
                             | | |
            read  |--<--<--<-+-<-+--|
                  0          ax

        Case 2 (non-complement alignment):

                          0  rx
            ref           |--+->-+->-->-->--|
                             | | |
            read  |-->-->-->-+->-+--|
                  0          ax
        ---
    */
    @property bool isFrontExtension() const pure nothrow
    {
        if (_length != 1)
        {
            return false;
        }

        return _isFrontExtension(_alignments[0]);
    }

    private static bool _isFrontExtension(in AlignmentChain alignment) pure nothrow
    {
        auto readExtensionLength = alignment.first.contigB.begin;
        auto referenceExtensionLength = alignment.first.contigA.begin;

        return readExtensionLength > referenceExtensionLength;
    }

    /**
        Returns true iff the read alignment is an back extension, ie. it is an
        extension and reaches over the back of the reference contig.

        ---
        Case 1 (complement alignment):

                  0             ry lr
            ref   |-->-->-->-+->-+--|
                             | | |
            read          |--+-<-+-<--<--<--|
                          0     ay         la

        Case 2 (non-complement alignment):

                  0             ry lr
            ref   |-->-->-->-+->-+--|
                             | | |
            read          |--+->-+->-->-->--|
                          0     ay         la
        ---
    */
    @property bool isBackExtension() const pure nothrow
    {
        if (_length != 1)
        {
            return false;
        }

        return _isBackExtension(_alignments[0]);
    }

    private static bool _isBackExtension(in AlignmentChain alignment) pure nothrow
    {
        auto readExtensionLength = alignment.contigB.length - alignment.last.contigB.end;
        auto referenceExtensionLength = alignment.contigA.length - alignment.last.contigA.end;

        return readExtensionLength > referenceExtensionLength;
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
    @property bool isGap() const pure nothrow
    {
        // dfmt off
        return _length == 2 &&
            _alignments[0].contigA.id != _alignments[1].contigA.id &&
            _alignments[0].contigB.id == _alignments[1].contigB.id &&
            _isExtension(_alignments[0]) && _isExtension(_alignments[1]) &&
            (_isBackExtension(_alignments[0]) == _isFrontExtension(_alignments[1]));
        // dfmt on
    }

    unittest
    {
        with (AlignmentChain) with (LocalAlignment) with (Complement)
                {
                    // dfmt off
                    auto testCases = [
                        "innerAlignmentComplement": tuple(
                            ReadAlignment(
                                AlignmentChain(
                                    1,
                                    Contig(1, 100),
                                    Contig(1, 10),
                                    no,
                                    [
                                        LocalAlignment(
                                            Locus(10, 11),
                                            Locus(0, 1),
                                            0,
                                        ),
                                        LocalAlignment(
                                            Locus(19, 20),
                                            Locus(9, 10),
                                            0,
                                        ),
                                    ],
                                ),
                            ),
                            No.isValid,
                            ReadAlignmentType.front,
                            No.isExtension,
                            No.isFrontExtension,
                            No.isBackExtension,
                            No.isGap,
                        ),
                        "innerAlignmentComplement": tuple(
                            ReadAlignment(
                                AlignmentChain(
                                    2,
                                    Contig(1, 100),
                                    Contig(1, 10),
                                    yes,
                                    [
                                        LocalAlignment(
                                            Locus(10, 11),
                                            Locus(0, 1),
                                            0,
                                        ),
                                        LocalAlignment(
                                            Locus(19, 20),
                                            Locus(9, 10),
                                            0,
                                        ),
                                    ],
                                ),
                            ),
                            No.isValid,
                            ReadAlignmentType.front,
                            No.isExtension,
                            No.isFrontExtension,
                            No.isBackExtension,
                            No.isGap,
                        ),
                        "frontExtension": tuple(
                            ReadAlignment(
                                AlignmentChain(
                                    3,
                                    Contig(1, 100),
                                    Contig(1, 10),
                                    no,
                                    [
                                        LocalAlignment(
                                            Locus(2, 3),
                                            Locus(5, 6),
                                            0,
                                        ),
                                        LocalAlignment(
                                            Locus(5, 6),
                                            Locus(9, 10),
                                            0,
                                        ),
                                    ],
                                ),
                            ),
                            Yes.isValid,
                            ReadAlignmentType.front,
                            Yes.isExtension,
                            Yes.isFrontExtension,
                            No.isBackExtension,
                            No.isGap,
                        ),
                        "frontExtensionComplement": tuple(
                            ReadAlignment(
                                AlignmentChain(
                                    4,
                                    Contig(1, 100),
                                    Contig(1, 10),
                                    yes,
                                    [
                                        LocalAlignment(
                                            Locus(2, 3),
                                            Locus(5, 6),
                                            0,
                                        ),
                                        LocalAlignment(
                                            Locus(5, 6),
                                            Locus(9, 10),
                                            0,
                                        ),
                                    ],
                                ),
                            ),
                            Yes.isValid,
                            ReadAlignmentType.front,
                            Yes.isExtension,
                            Yes.isFrontExtension,
                            No.isBackExtension,
                            No.isGap,
                        ),
                        "backExtension": tuple(
                            ReadAlignment(
                                AlignmentChain(
                                    5,
                                    Contig(1, 100),
                                    Contig(1, 10),
                                    no,
                                    [
                                        LocalAlignment(
                                            Locus(94, 95),
                                            Locus(0, 1),
                                            0,
                                        ),
                                        LocalAlignment(
                                            Locus(97, 98),
                                            Locus(4, 5),
                                            0,
                                        ),
                                    ],
                                ),
                            ),
                            Yes.isValid,
                            ReadAlignmentType.back,
                            Yes.isExtension,
                            No.isFrontExtension,
                            Yes.isBackExtension,
                            No.isGap,
                        ),
                        "backExtensionComplement": tuple(
                            ReadAlignment(
                                AlignmentChain(
                                    6,
                                    Contig(1, 100),
                                    Contig(1, 10),
                                    yes,
                                    [
                                        LocalAlignment(
                                            Locus(94, 95),
                                            Locus(0, 1),
                                            0,
                                        ),
                                        LocalAlignment(
                                            Locus(97, 98),
                                            Locus(4, 5),
                                            0,
                                        ),
                                    ],
                                ),
                            ),
                            Yes.isValid,
                            ReadAlignmentType.back,
                            Yes.isExtension,
                            No.isFrontExtension,
                            Yes.isBackExtension,
                            No.isGap,
                        ),
                        "gap": tuple(
                            ReadAlignment(
                                AlignmentChain(
                                    7,
                                    Contig(1, 100),
                                    Contig(1, 10),
                                    no,
                                    [
                                        LocalAlignment(
                                            Locus(94, 95),
                                            Locus(0, 1),
                                            0,
                                        ),
                                        LocalAlignment(
                                            Locus(97, 98),
                                            Locus(4, 5),
                                            0,
                                        ),
                                    ],
                                ),
                                AlignmentChain(
                                    8,
                                    Contig(2, 100),
                                    Contig(1, 10),
                                    no,
                                    [
                                        LocalAlignment(
                                            Locus(2, 3),
                                            Locus(5, 6),
                                            0,
                                        ),
                                        LocalAlignment(
                                            Locus(5, 6),
                                            Locus(9, 10),
                                            0,
                                        ),
                                    ],
                                ),
                            ),
                            Yes.isValid,
                            ReadAlignmentType.gap,
                            No.isExtension,
                            No.isFrontExtension,
                            No.isBackExtension,
                            Yes.isGap,
                        ),
                        "gapComplement": tuple(
                            ReadAlignment(
                                AlignmentChain(
                                    9,
                                    Contig(1, 100),
                                    Contig(1, 10),
                                    yes,
                                    [
                                        LocalAlignment(
                                            Locus(2, 3),
                                            Locus(5, 6),
                                            0,
                                        ),
                                        LocalAlignment(
                                            Locus(5, 6),
                                            Locus(9, 10),
                                            0,
                                        ),
                                    ],
                                ),
                                AlignmentChain(
                                    10,
                                    Contig(2, 100),
                                    Contig(1, 10),
                                    yes,
                                    [
                                        LocalAlignment(
                                            Locus(94, 95),
                                            Locus(0, 1),
                                            0,
                                        ),
                                        LocalAlignment(
                                            Locus(97, 98),
                                            Locus(4, 5),
                                            0,
                                        ),
                                    ],
                                ),
                            ),
                            Yes.isValid,
                            ReadAlignmentType.gap,
                            No.isExtension,
                            No.isFrontExtension,
                            No.isBackExtension,
                            Yes.isGap,
                        ),
                    ];
                    // dfmt on

                    alias getFailureMessage = (testCase, testFunction, expectedValue) => format!"expected %s.%s to be %s"(
                            testCase, testFunction, expectedValue);

                    foreach (testCase, testData; testCases)
                    {
                        auto readAlignment = testData[0];

                        assert(testData[1] == readAlignment.isValid,
                                getFailureMessage(testCase, "isValid", testData[1]));
                        if (readAlignment.isValid)
                            assert(testData[2] == readAlignment.type,
                                    getFailureMessage(testCase, "type", testData[2]));
                        else
                            assertThrown!AssertError(readAlignment.type,
                                    format!"expected type(%s) to throw"(testCase));
                        assert(testData[3] == readAlignment.isExtension,
                                getFailureMessage(testCase, "isExtension", testData[3]));
                        assert(testData[4] == readAlignment.isFrontExtension,
                                getFailureMessage(testCase, "isFrontExtension", testData[4]));
                        assert(testData[5] == readAlignment.isBackExtension,
                                getFailureMessage(testCase, "isBackExtension", testData[5]));
                        assert(testData[6] == readAlignment.isGap,
                                getFailureMessage(testCase, "isGap", testData[6]));
                    }
                }
    }

    /**
        Returns the (approximate) size of the insertion produced by this
        read alignment.

        See_Also: `getGapSize`, `getExtensionSize`
    */
    long getInsertionSize() const pure
    {
        final switch (type)
        {
        case ReadAlignmentType.front:
        case ReadAlignmentType.back:
            return getExtensionSize();
        case ReadAlignmentType.gap:
            return getGapSize();
        }
    }

    unittest
    {
        with (AlignmentChain) with (LocalAlignment) with (Complement)
                {
                    // dfmt off
                    auto extensionAlignment = ReadAlignment(
                        AlignmentChain(0, Contig(1, 50), Contig(1, 50), no, [
                            LocalAlignment(Locus(0, 15), Locus(10, 20), 0),
                            LocalAlignment(Locus(20, 40), Locus(30, 50), 1),
                        ])
                    );
                    auto gapAlignment = ReadAlignment(
                        AlignmentChain(0, Contig(1, 50), Contig(1, 80), no, [
                            LocalAlignment(Locus(30, 35), Locus(10, 15), 0),
                            LocalAlignment(Locus(40, 50), Locus(20, 30), 1),
                        ]),
                        AlignmentChain(1, Contig(2, 50), Contig(1, 80), no, [
                            LocalAlignment(Locus(0, 10), Locus(50, 55), 2),
                            LocalAlignment(Locus(15, 20), Locus(60, 70), 3),
                        ]),
                    );
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
    private long getExtensionSize() const pure nothrow
    {
        assert(isExtension);

        auto alignment = _alignments[0];
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
                            ReadAlignment(AlignmentChain(0, Contig(1, 50), Contig(1, 50), no, [
                                LocalAlignment(Locus(0, 15), Locus(10, 20), 0),
                                LocalAlignment(Locus(20, 40), Locus(30, 50), 1),
                            ])),
                            10
                        ),
                        //            0      10  45  50
                        // reference  |-------+---+---|
                        //                     \ \ \
                        // read         |-->-->-+->-+->-->--|
                        //              0       0  40      50
                        tuple(
                            ReadAlignment(AlignmentChain(1, Contig(1, 50), Contig(1, 50), no, [
                                LocalAlignment(Locus(0, 15), Locus(10, 20), 2),
                                LocalAlignment(Locus(20, 40), Locus(30, 45), 3),
                            ])),
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
                            ReadAlignment(AlignmentChain(2, Contig(1, 50), Contig(1, 50), no, [
                                LocalAlignment(Locus(10, 15), Locus(0, 20), 4),
                                LocalAlignment(Locus(20, 50), Locus(30, 40), 5),
                            ])),
                            10
                        ),
                        //                  0   5  40      50
                        // reference        |---+---+-------|
                        //                     / / /
                        // read       |-->-->-+->-+->-->--|
                        //            0      10  50      50
                        tuple(
                            ReadAlignment(AlignmentChain(3, Contig(1, 50), Contig(1, 50), no, [
                                LocalAlignment(Locus(10, 15), Locus(5, 20), 6),
                                LocalAlignment(Locus(20, 50), Locus(30, 40), 7),
                            ])),
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
                            ReadAlignment(AlignmentChain(4, Contig(1, 50), Contig(1, 50), yes, [
                                LocalAlignment(Locus(0, 15), Locus(10, 20), 8),
                                LocalAlignment(Locus(20, 40), Locus(30, 50), 9),
                            ])),
                            10
                        ),
                        //            0      10  45  50
                        // reference  |-------+---+---|
                        //                     \ \ \
                        // read         |--<--<-+-<-+-<--<--|
                        //              0       0   50      50
                        tuple(
                            ReadAlignment(AlignmentChain(5, Contig(1, 50), Contig(1, 50), yes, [
                                LocalAlignment(Locus(0, 15), Locus(10, 20), 10),
                                LocalAlignment(Locus(20, 40), Locus(30, 45), 11),
                            ])),
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
                            ReadAlignment(AlignmentChain(6, Contig(1, 50), Contig(1, 50), yes, [
                                LocalAlignment(Locus(10, 15), Locus(0, 20), 12),
                                LocalAlignment(Locus(20, 50), Locus(30, 40), 13),
                            ])),
                            10
                        ),
                        //                  0   5  40      50
                        // reference        |---+---+-------|
                        //                     / / /
                        // read       |--<--<-+-<-+-<--<--|
                        //            0       10  50      50
                        tuple(
                            ReadAlignment(AlignmentChain(7, Contig(1, 50), Contig(1, 50), yes, [
                                LocalAlignment(Locus(10, 15), Locus(5, 20), 14),
                                LocalAlignment(Locus(20, 50), Locus(30, 40), 15),
                            ])),
                            5
                        ),
                    ];
                    // dfmt on

                    foreach (testCaseIdx, testCase; testCases.enumerate)
                    {
                        auto gotValue = testCase[0].getExtensionSize();
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
    private long getGapSize() const pure
    {
        assert(isGap);
        // dfmt off
        const(AlignmentChain[2]) orderedAlignments = isInOrder
            ? [_alignments[0], _alignments[1]]
            : [_alignments[1], _alignments[0]];
        // dfmt on
        auto firstAlignment = orderedAlignments[0];
        auto secondAlignment = orderedAlignments[1];

        assert(firstAlignment.last.contigB.end < secondAlignment.first.contigB.begin,
                format!"intersecting local alignments in [%s, %s]"(firstAlignment, secondAlignment));

        long readsSequenceLength = secondAlignment.first.contigB.begin - firstAlignment
            .last.contigB.end;
        // dfmt off
        long referenceExcess1 = firstAlignment.contigA.length - firstAlignment.last.contigA.end;
        long referenceExcess2 = secondAlignment.first.contigA.begin;
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
                            ReadAlignment(
                                AlignmentChain(0, Contig(1, 50), Contig(1, 80), no, [
                                    LocalAlignment(Locus(30, 35), Locus(10, 15), 0),
                                    LocalAlignment(Locus(40, 50), Locus(20, 30), 1),
                                ]),
                                AlignmentChain(1, Contig(2, 50), Contig(1, 80), no, [
                                    LocalAlignment(Locus(0, 10), Locus(50, 55), 2),
                                    LocalAlignment(Locus(15, 20), Locus(60, 70), 3),
                                ]),
                            ),
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
                            ReadAlignment(
                                AlignmentChain(2, Contig(1, 50), Contig(1, 80), no, [
                                    LocalAlignment(Locus(30, 35), Locus(10, 15), 4),
                                    LocalAlignment(Locus(40, 45), Locus(20, 30), 5),
                                ]),
                                AlignmentChain(3, Contig(2, 50), Contig(1, 80), no, [
                                    LocalAlignment(Locus(5, 10), Locus(50, 55), 6),
                                    LocalAlignment(Locus(15, 20), Locus(60, 70), 7),
                                ]),
                            ),
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
                            ReadAlignment(
                                AlignmentChain(4, Contig(2, 50), Contig(1, 80), no, [
                                    LocalAlignment(Locus(5, 10), Locus(50, 55), 6),
                                    LocalAlignment(Locus(15, 20), Locus(60, 70), 7),
                                ]),
                                AlignmentChain(5, Contig(1, 50), Contig(1, 80), no, [
                                    LocalAlignment(Locus(30, 35), Locus(10, 15), 4),
                                    LocalAlignment(Locus(40, 45), Locus(20, 30), 5),
                                ]),
                            ),
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
                            ReadAlignment(
                                AlignmentChain(6, Contig(1, 50), Contig(1, 80), yes, [
                                    LocalAlignment(Locus(30, 35), Locus(10, 15), 8),
                                    LocalAlignment(Locus(40, 50), Locus(20, 30), 9),
                                ]),
                                AlignmentChain(7, Contig(2, 50), Contig(1, 80), yes, [
                                    LocalAlignment(Locus(0, 10), Locus(50, 55), 10),
                                    LocalAlignment(Locus(15, 20), Locus(60, 70), 11),
                                ]),
                            ),
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
                            ReadAlignment(
                                AlignmentChain(8, Contig(1, 50), Contig(1, 80), yes, [
                                    LocalAlignment(Locus(30, 35), Locus(10, 15), 12),
                                    LocalAlignment(Locus(40, 45), Locus(20, 30), 13),
                                ]),
                                AlignmentChain(9, Contig(2, 50), Contig(1, 80), yes, [
                                    LocalAlignment(Locus(5, 10), Locus(50, 55), 14),
                                    LocalAlignment(Locus(15, 20), Locus(60, 70), 15),
                                ]),
                            ),
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
                            ReadAlignment(
                                AlignmentChain(10, Contig(2, 50), Contig(1, 80), yes, [
                                    LocalAlignment(Locus(5, 10), Locus(50, 55), 14),
                                    LocalAlignment(Locus(15, 20), Locus(60, 70), 15),
                                ]),
                                AlignmentChain(11, Contig(1, 50), Contig(1, 80), yes, [
                                    LocalAlignment(Locus(30, 35), Locus(10, 15), 12),
                                    LocalAlignment(Locus(40, 45), Locus(20, 30), 13),
                                ]),
                            ),
                            10
                        ),
                    ];
                    // dfmt on

                    foreach (testCaseIdx, testCase; testCases.enumerate)
                    {
                        auto gotValue = testCase[0].getGapSize();
                        auto expValue = testCase[1];
                        auto errorMessage = format!"expected spanning gap size %d but got %d for test case %d"(
                                expValue, gotValue, testCaseIdx);

                        assert(gotValue == expValue, errorMessage);
                    }
                }
    }

    double meanScore() const pure
    {
        return this[].map!"a.score".mean;
    }
}

/**
    A pile of read alignments belonging to the same gap/contig end.

    See_Also: Hit
*/
alias PileUp = ReadAlignment[];

PileUp[] buildPileUps(AlignmentChain[] candidates)
{
    static bool orderByReferenceIds(ReadAlignment a, ReadAlignment b) pure
    {
        return a[0].contigA.id < b[0].contigA.id && (!(a.length > 1
                && a[0].contigA.id == b[0].contigA.id) || a[1].contigA.id < b[1].contigA.id);
    }

    static size_t numContigsInvolved(T)(T readsByNumFlanksByGap)
    {
        auto readsByGap = readsByNumFlanksByGap.front;
        auto readsByFirstGap = readsByGap.front;
        auto reads = readsByFirstGap.front;

        return reads.length;
    }

    static bool isSameRead(AlignmentChain a, AlignmentChain b)
    {
        return a.contigB.id == b.contigB.id;
    }

    static auto forceInOrder(ReadAlignment readAlignment)
    {
        readAlignment.forceInOrder;

        return readAlignment;
    }

    candidates.sort!("a.contigB.id < b.contigB.id", SwapStrategy.stable);

    // dfmt off
    auto sortedReadAlignments = candidates
        .chunkBy!isSameRead
        .map!array
        .map!ReadAlignment
        .filter!"a.isValid"
        .map!forceInOrder
        .array;
    // dfmt on
    sortedReadAlignments.sort!(byGapSort, SwapStrategy.stable);
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

    debug logJsonDebug("pileStructure", pileUpsByGap.map!"a.length"
            .map!toJson
            .array);

    return pileUpsByGap;
}

bool belongToSamePile(in ReadAlignment a, in ReadAlignment b) pure nothrow
{
    assert(a.isInOrder && b.isInOrder,
            "gap alignment chains must be in reference order");

    // dfmt off
    debug logJsonDebug(
        "aIsBackExtension", a.isBackExtension,
        "aIsFrontExtension", a.isFrontExtension,
        "aIsGap", a.isGap,
        "aReferenceIds", a[].map!"a.contigA.id".array.toJson,
        "bIsBackExtension", b.isBackExtension,
        "bIsFrontExtension", b.isFrontExtension,
        "bIsGap", b.isGap,
        "bReferenceIds", b[].map!"a.contigA.id".array.toJson,
    );
    // dfmt on

    // dfmt off
    return any(only(
        a.isBackExtension && b.isBackExtension && a[0].contigA.id == b[0].contigA.id,
        a.isBackExtension && b.isGap && a[0].contigA.id == b[0].contigA.id,
        a.isGap && b.isFrontExtension && a[1].contigA.id == b[0].contigA.id,
        a.isGap && b.isGap && a[0].contigA.id == b[0].contigA.id && a[1].contigA.id == b[1].contigA.id,
        a.isFrontExtension && b.isFrontExtension && a[0].contigA.id == b[0].contigA.id,
    ));
    // dfmt on
}

bool byGapSort(in ReadAlignment a, in ReadAlignment b) pure nothrow
{
    assert(a.isInOrder && b.isInOrder,
            "gap alignment chains must be in reference order");

    // dfmt off
    debug logJsonDebug(
        "a", a[].toJson,
        "a.type", a.type.to!string,
        "b", b[].toJson,
        "b.type", b.type.to!string,
    );
    // dfmt on

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
    if (a[0].contigA.id != b[0].contigA.id)
    {
        return a[0].contigA.id < b[0].contigA.id;
    }
    else if (a.type != b.type)
    {
        return a.type < b.type;
    }
    else
    {
        // same contig and same extension type => same pile
        return false;
    }
}

private bool byGapSortForGaps(in ReadAlignment a, in ReadAlignment b) pure nothrow
{
    if (a[0].contigA.id != b[0].contigA.id)
    {
        return a[0].contigA.id < b[0].contigA.id;
    }
    else if (a[1].contigA.id != b[1].contigA.id)
    {
        return a[1].contigA.id < b[1].contigA.id;
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

    if (ext[0].contigA.id != gap[0].contigA.id && ext[0].contigA.id != gap[1].contigA.id)
    {
        return cmp(ext[0].contigA.id, min(gap[0].contigA.id, gap[1].contigA.id));
    }
    else if (ext[0].contigA.id != gap[0].contigA.id && ext[0].contigA.id == gap[1].contigA.id)
    {
        // extension contig is end contig of gap => same pile if front extension else in reverse order
        return !inOrder;
    }
    else if (ext[0].contigA.id == gap[0].contigA.id && ext[0].contigA.id != gap[1].contigA.id)
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

                if (beginContigId == endContigId)
                {
                    // dfmt off
                    return ReadAlignment(AlignmentChain(
                        alignmentChainId - 2,
                        Contig(beginContigId, contigLength),
                        Contig(readId, readLength),
                        complement,
                        [
                            LocalAlignment(
                                Locus(beginIdx, beginIdx + 1),
                                Locus(firstReadBeginIdx, firstReadBeginIdx + 1),
                                numDiffs,
                            ),
                            LocalAlignment(
                                Locus(endIdx - 1, endIdx),
                                Locus(firstReadEndIdx - 1, firstReadEndIdx),
                                numDiffs,
                            ),
                        ],
                    ));
                    // dfmt on
                }
                else
                {
                    auto secondReadBeginIdx = readLength - endIdx;
                    auto secondReadEndIdx = readLength;

                    // dfmt off
                    return ReadAlignment(
                        AlignmentChain(
                            alignmentChainId - 2,
                            Contig(beginContigId, contigLength),
                            Contig(readId, readLength),
                            complement,
                            [
                                LocalAlignment(
                                    Locus(beginIdx, beginIdx + 1),
                                    Locus(firstReadBeginIdx, firstReadBeginIdx + 1),
                                    numDiffs,
                                ),
                                LocalAlignment(
                                    Locus(contigLength - 1, contigLength),
                                    Locus(firstReadEndIdx - 1, firstReadEndIdx),
                                    numDiffs,
                                ),
                            ],
                        ),
                        AlignmentChain(
                            alignmentChainId - 1,
                            Contig(endContigId, contigLength),
                            Contig(readId, readLength),
                            complement,
                            [
                                LocalAlignment(
                                    Locus(0, 1),
                                    Locus(secondReadBeginIdx, secondReadBeginIdx + 1),
                                    numDiffs,
                                ),
                                LocalAlignment(
                                    Locus(endIdx - 1, endIdx),
                                    Locus(secondReadEndIdx - 1, secondReadEndIdx),
                                    numDiffs,
                                ),
                            ],
                        ),
                    );
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
                            size_t read1Id = readAlignment1[0].contigB.id;
                            size_t read2Id = readAlignment2[0].contigB.id;

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
        // dfmt off
        return pileUp[0].isFrontExtension
            ? ReadAlignmentType.front
            : ReadAlignmentType.back;
        // dfmt on
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
    auto pileUp = [ReadAlignment(AlignmentChain(), AlignmentChain()), ReadAlignment(AlignmentChain())];
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
