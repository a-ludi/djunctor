/**
    Defines a Region and common operation with these. A Region is a set of
    tagged intervals where differently tagged intervals are distinct.

    Copyright: © 2018 Arne Ludwig <arne.ludwig@posteo.de>
    License: Subject to the terms of the MIT license, as written in the
             included LICENSE file.
    Authors: Arne Ludwig <arne.ludwig@posteo.de>
*/
module djunctor.util.region;

import std.algorithm : all, cache, chunkBy, cmp, filter, fold, map, max, min,
    sort, sum;
import std.array : appender, array, join;
import std.exception : assertThrown;
import std.format : format;
import std.range : assumeSorted, chunks, ElementType, InputRange,
    inputRangeObject, isInputRange, only;
import std.traits : isNumeric, Unqual;
import std.stdio : writeln;

/// Returns the type of the property `tag` of `T`.
template TagType(T)
{
    private T instance;

    alias TagType = typeof(instance.tag);
}

/// Checks if T has a property `tag` implicitly convertible to `Tag` – if given.
template isTaggable(T)
{
    enum isTaggable = is(TagType!T);
}

unittest
{
    struct Taggable
    {
        int tag;
    }

    static assert(isTaggable!Taggable);
    static assert(!isTaggable!int);
}

unittest
{
    struct Taggable
    {
        int tag;
    }

    static assert(isTaggable!Taggable);
    static assert(!isTaggable!int);
}

/// Thrown if two operands require the same tag but different were provided.
static class MismatchingTagsException(Tag) : Exception
{
    const(Tag[2]) tags;

    this(in Tag tagA, in Tag tagB)
    {
        this.tags = [tagA, tagB];

        super(format!"mismatching interval tags: %s"(this.tags));
    }
}

/*
    Throws an exception if tags do not match.

    Throws: MismatchingTagsException if tags do not match.
*/
void enforceMatchingTags(Taggable)(in Taggable taggableA, in Taggable taggableB) pure
        if (isTaggable!Taggable)
{
    alias Tag = TagType!Taggable;

    if (taggableA.tag != taggableB.tag)
    {
        throw new MismatchingTagsException!Tag(taggableA.tag, taggableB.tag);
    }
}

/**
    A Region is a set of tagged intervals where differently tagged intervals are distinct.
*/
struct Region(Number, Tag, string tagAlias = null, Tag emptyTag = Tag.init)
{
    static assert(isNumeric!Number, "interval limits must be numeric");

    /**
        This is a right-open interval `[begin, end)` tagged with `tag`.
        If `tagAlias` is given then the tag may be access as a property of
        that name.
    */
    static struct TaggedInterval
    {
        Tag tag = emptyTag;
        Number begin;
        Number end;
        static if (!(tagAlias is null) && tagAlias != "tag")
        {
            mixin("alias " ~ tagAlias ~ " = tag;");
        }

        invariant
        {
            assert(begin <= end, "begin must be less than or equal to end");
        }

        /// Returns the size of this interval.
        @property Number size() pure const nothrow
        {
            return this.end - this.begin;
        }

        ///
        unittest
        {
            assert(Region!(int, int).TaggedInterval().size == 0);
            assert(TaggedInterval(1, 10, 20).size == 10);
            assert(TaggedInterval(2, 20, 40).size == 20);
            assert(TaggedInterval(3, 30, 60).size == 30);
        }

        /// Returns true iff the interval is empty. An interval is empty iff
        /// `begin == end`.
        @property bool empty() pure const nothrow
        {
            return begin >= end;
        }

        ///
        unittest
        {
            assert(TaggedInterval().empty);
            assert(!TaggedInterval(1, 10, 20).empty);
            assert(!TaggedInterval(2, 20, 40).empty);
            assert(TaggedInterval(3, 60, 60).empty);
        }

        /**
            Returns the convex hull of the intervals.

            Throws: MismatchingTagsException if `tag`s differ.
        */
        static TaggedInterval convexHull(in TaggedInterval[] intervals...) pure
        {
            if (intervals.length == 0)
            {
                return TaggedInterval();
            }

            TaggedInterval convexHullInterval = intervals[0];

            foreach (interval; intervals[1 .. $])
            {
                enforceMatchingTags(convexHullInterval, interval);

                if (convexHullInterval.empty)
                {
                    convexHullInterval = interval;
                }
                else if (interval.empty)
                {
                    continue;
                }
                else
                {
                    convexHullInterval.begin = min(convexHullInterval.begin, interval.begin);
                    convexHullInterval.end = max(convexHullInterval.end, interval.end);
                }
            }

            // dfmt off
            return convexHullInterval.empty
                ? TaggedInterval()
                : convexHullInterval;
            // dfmt on
        }

        ///
        unittest
        {
            alias R = Region!(int, int);
            alias TI = R.TaggedInterval;

            assert(TI.convexHull(TI(0, 10, 20), TI(0, 0, 5)) == TI(0, 0, 20));
            assert(TI.convexHull(TI(0, 10, 20), TI(0, 5, 15)) == TI(0, 5, 20));
            assert(TI.convexHull(TI(0, 10, 20), TI(0, 12, 18)) == TI(0, 10, 20));
            assert(TI.convexHull(TI(0, 10, 20), TI(0, 10, 20)) == TI(0, 10, 20));
            assert(TI.convexHull(TI(0, 10, 20), TI(0, 15, 25)) == TI(0, 10, 25));
            assert(TI.convexHull(TI(0, 10, 20), TI(0, 25, 30)) == TI(0, 10, 30));
            assertThrown!(MismatchingTagsException!int)(TI.convexHull(TI(0, 10, 20), TI(1, 25, 30)));
        }

        /// Returns the intersection of both intervals; empty if `tag`s differ.
        TaggedInterval opBinary(string op)(in TaggedInterval other) const pure nothrow
                if (op == "&")
        {
            if (this.tag != other.tag)
            {
                return TaggedInterval();
            }

            auto newBegin = max(this.begin, other.begin);
            auto newEnd = min(this.end, other.end);

            if (newBegin > newEnd)
            {
                return TaggedInterval();
            }

            // dfmt off
            return TaggedInterval(
                tag,
                newBegin,
                newEnd,
            );
            // dfmt on
        }

        /// ditto
        TaggedInterval opOpAssign(string op)(in TaggedInterval other) if (op == "&")
        {
            {
                auto tmp = this & other;

                this.tag = tmp.tag;
                this.begin = tmp.begin;
                this.end = tmp.end;

                return this;
            }
        }
        ///
        unittest
        {
            alias R = Region!(int, int);
            alias TI = R.TaggedInterval;

            assert((TI(0, 10, 20) & TI(0, 0, 5)).empty);
            assert((TI(0, 10, 20) & TI(0, 5, 15)) == TI(0, 10, 15));
            assert((TI(0, 10, 20) & TI(0, 12, 18)) == TI(0, 12, 18));
            assert((TI(0, 10, 20) & TI(0, 10, 20)) == TI(0, 10, 20));
            assert((TI(0, 10, 20) & TI(0, 15, 25)) == TI(0, 15, 20));
            assert((TI(0, 10, 20) & TI(0, 25, 30)).empty);
            assert((TI(0, 10, 20) & TI(1, 25, 30)).empty);
        }

        /// Returns the difference of both intervals.
        Region opBinary(string op)(in TaggedInterval other) const if (op == "-")
        {
            auto intersection = this & other;

            if (intersection.empty)
            {
                TaggedInterval thisCopy = this;

                return Region(this.empty ? [] : [thisCopy]);
            }

            // dfmt off
            return Region(only(
                TaggedInterval(
                    tag,
                    this.begin,
                    intersection.begin,
                ),
                TaggedInterval(
                    tag,
                    intersection.end,
                    this.end,
                ),
            ).filter!"!a.empty".array);
            // dfmt on
        }

        ///
        unittest
        {
            alias R = Region!(int, int);
            alias TI = R.TaggedInterval;

            assert(TI(0, 10, 20) - TI(0, 0, 5) == R([TI(0, 10, 20)]));
            assert(TI(0, 10, 20) - TI(0, 5, 15) == R([TI(0, 15, 20)]));
            assert(TI(0, 10, 20) - TI(0, 12, 18) == R([TI(0, 10, 12), TI(0, 18, 20)]));
            assert(TI(0, 10, 20) - TI(0, 10, 20) == R([]));
            assert(TI(0, 10, 20) - TI(0, 15, 25) == R([TI(0, 10, 15)]));
            assert(TI(0, 10, 20) - TI(0, 25, 30) == R([TI(0, 10, 20)]));
            assert(TI(0, 10, 20) - TI(1, 25, 30) == R([TI(0, 10, 20)]));
        }

        int opCmp(in TaggedInterval other) const pure nothrow
        {
            // dfmt off
            return cmp(
                only(this.tag, this.begin, this.end),
                only(other.tag, other.begin, other.end),
            );
            // dfmt on
        }

        ///
        unittest
        {
            alias R = Region!(int, int);
            alias TI = R.TaggedInterval;

            assert(TI(0, 10, 20) > TI(0, 0, 5));
            assert(TI(0, 10, 20) > TI(0, 5, 15));
            assert(TI(0, 10, 20) < TI(0, 12, 18));
            assert(TI(0, 10, 20) < TI(0, 15, 25));
            assert(TI(0, 10, 20) == TI(0, 10, 20));
            assert(TI(0, 10, 20) < TI(0, 25, 30));
            assert(TI(0, 10, 20) < TI(1, 25, 30));
        }

        /// Returns true iff the tagged intervals intersect.
        bool intersects(in TaggedInterval other) const pure nothrow
        {
            return !(this & other).empty;
        }

        ///
        unittest
        {
            alias R = Region!(int, int);
            alias TI = R.TaggedInterval;

            assert(!TI(0, 10, 20).intersects(TI(0, 0, 5)));
            assert(TI(0, 10, 20).intersects(TI(0, 5, 15)));
            assert(TI(0, 10, 20).intersects(TI(0, 12, 18)));
            assert(TI(0, 10, 20).intersects(TI(0, 15, 25)));
            assert(TI(0, 10, 20).intersects(TI(0, 10, 20)));
            assert(!TI(0, 10, 20).intersects(TI(0, 25, 30)));
            assert(!TI(0, 10, 20).intersects(TI(1, 25, 30)));
        }

        /// Returns true iff the tagged intervals do not intersect and `this < other`.
        bool isStrictlyBefore(in TaggedInterval other) const pure nothrow
        {
            return !this.intersects(other) && this < other;
        }

        ///
        unittest
        {
            alias R = Region!(int, int);
            alias TI = R.TaggedInterval;

            assert(!TI(0, 10, 20).isStrictlyBefore(TI(0, 0, 5)));
            assert(!TI(0, 10, 20).isStrictlyBefore(TI(0, 5, 15)));
            assert(!TI(0, 10, 20).isStrictlyBefore(TI(0, 12, 18)));
            assert(!TI(0, 10, 20).isStrictlyBefore(TI(0, 15, 25)));
            assert(!TI(0, 10, 20).isStrictlyBefore(TI(0, 10, 20)));
            assert(TI(0, 10, 20).isStrictlyBefore(TI(0, 25, 30)));
            assert(TI(0, 10, 20).isStrictlyBefore(TI(1, 25, 30)));
        }

        /// Returns true iff the tagged intervals do not intersect and `this > other`.
        bool isStrictlyAfter(in TaggedInterval other) const pure nothrow
        {
            return !this.intersects(other) && this > other;
        }

        ///
        unittest
        {
            alias R = Region!(int, int);
            alias TI = R.TaggedInterval;

            assert(TI(0, 10, 20).isStrictlyAfter(TI(0, 0, 5)));
            assert(!TI(0, 10, 20).isStrictlyAfter(TI(0, 5, 15)));
            assert(!TI(0, 10, 20).isStrictlyAfter(TI(0, 12, 18)));
            assert(!TI(0, 10, 20).isStrictlyAfter(TI(0, 15, 25)));
            assert(!TI(0, 10, 20).isStrictlyAfter(TI(0, 10, 20)));
            assert(!TI(0, 10, 20).isStrictlyAfter(TI(0, 25, 30)));
            assert(!TI(0, 10, 20).isStrictlyAfter(TI(1, 25, 30)));
        }
    }

    ///
    unittest
    {
        static immutable emptyTag = 42;
        alias R = Region!(int, int, "bucketId", emptyTag);
        alias TI = R.TaggedInterval;

        TI emptyInterval;

        // Default constructor produces empty interval.
        assert((emptyInterval).empty);
        assert(emptyInterval.tag == emptyTag);

        auto ti1 = TI(1, 0, 10);

        // The tag can be aliased:
        assert(ti1.tag == ti1.bucketId);

        auto ti2 = TI(1, 5, 15);

        // Tagged intervals with the same tag behave like regular intervals:
        assert((ti1 & ti2) == TI(1, 5, 10));

        auto ti3 = TI(2, 0, 10);

        // Tagged intervals with different tags are distinct:
        assert((ti1 & ti3).empty);
    }

    TaggedInterval[] _intervals;

    this(TaggedInterval[] intervals)
    {
        this._intervals = intervals;
        this.normalize();
    }

    ///
    unittest
    {
        alias R = Region!(int, int);
        alias TI = R.TaggedInterval;

        auto region = R([TI(0, 15, 20), TI(0, 5, 10), TI(0, 0, 10)]);

        // Intervals get implicitly normalized.
        assert(region.intervals == [TI(0, 0, 10), TI(0, 15, 20)]);
    }

    this(TaggedInterval interval)
    {
        this._intervals = [interval];
    }

    ///
    unittest
    {
        alias R = Region!(int, int);
        alias TI = R.TaggedInterval;

        R region = TI(0, 15, 20);

        assert(region.intervals == [TI(0, 15, 20)]);
    }

    this(Tag tag, Number begin, Number end)
    {
        this._intervals = [TaggedInterval(tag, begin, end)];
    }

    ///
    unittest
    {
        alias R = Region!(int, int);
        alias TI = R.TaggedInterval;

        auto region = R(0, 15, 20);

        assert(region.intervals == [TI(0, 15, 20)]);
    }

    this(this)
    {
        this._intervals = this._intervals.dup;
    }

    unittest
    {
        alias R = Region!(int, int);
        alias TI = R.TaggedInterval;

        auto tis = [TI(0, 15, 20)];
        auto region = R(tis);
        auto regionDup = region;

        assert(region == regionDup);

        // Changing the original does not affect duplicate
        region |= R(TI(0, 25, 30));

        assert(region != regionDup);
    }

    /// Return a list of the tagged intervals in this region.
    @property const(TaggedInterval)[] intervals() const pure nothrow
    {
        return _intervals;
    }

    ///
    unittest
    {
        alias R = Region!(int, int);
        alias TI = R.TaggedInterval;

        R emptyRegion;
        auto region1 = R(0, 5, 10);
        auto region2 = R([TI(0, 5, 10), TI(0, 15, 20)]);

        assert(emptyRegion.intervals == []);
        assert(region1.intervals == [TI(0, 5, 10)]);
        assert(region2.intervals == [TI(0, 5, 10), TI(0, 15, 20)]);
    }

    /// Returns the size of this region.
    Number size() pure const nothrow
    {
        return _intervals.map!"a.size".sum;
    }

    ///
    unittest
    {
        alias R = Region!(int, int);
        alias TI = R.TaggedInterval;

        R emptyRegion;
        auto region1 = R(0, 0, 10);
        auto region2 = R([TI(0, 0, 10), TI(0, 20, 30)]);
        auto region3 = R([TI(0, 0, 20), TI(0, 10, 30)]);

        assert(emptyRegion.size == 0);
        assert(region1.size == 10);
        assert(region2.size == 20);
        assert(region3.size == 30);
    }

    /// Returns true iff the region is empty.
    bool empty() pure const nothrow
    {
        return _intervals.length == 0;
    }

    ///
    unittest
    {
        alias R = Region!(int, int);
        alias TI = R.TaggedInterval;

        R emptyRegion1;
        auto emptyRegion2 = R([TI(0, 0, 0), TI(0, 10, 10)]);
        auto emptyRegion3 = R([TI(0, 0, 0), TI(1, 0, 0)]);
        auto region1 = R(0, 0, 10);

        assert(emptyRegion1.empty);
        assert(emptyRegion2.empty);
        assert(emptyRegion3.empty);
        assert(!region1.empty);
    }

    protected void normalize()
    {
        if (_intervals.length == 0)
        {
            return;
        }

        _intervals.sort();
        TaggedInterval accInterval = _intervals[0];
        size_t insertIdx = 0;

        foreach (i, intervalB; _intervals[1 .. $])
        {
            if (intervalB.empty)
            {
                continue;
            }
            else if (accInterval.intersects(intervalB))
            {
                // If two intervals intersect their union is the same as the convex hull of both.
                accInterval = TaggedInterval.convexHull(accInterval, intervalB);
            }
            else
            {
                if (!accInterval.empty)
                {
                    _intervals[insertIdx++] = accInterval;
                }
                accInterval = intervalB;
            }
        }

        if (!accInterval.empty)
        {
            _intervals[insertIdx++] = accInterval;
        }
        _intervals.length = insertIdx;
    }

    unittest
    {
        alias R = Region!(int, int);
        alias TI = R.TaggedInterval;

        auto region = R([
            TI(3, 8349, 8600),
            TI(3, 8349, 8349),
            TI(3, 8349, 8850),
            TI(3, 8349, 9100),
            TI(3, 8349, 8349),
            TI(3, 8349, 9350),
            TI(3, 8349, 9600),
            TI(3, 8349, 8349),
            TI(3, 8349, 9900),
            TI(3, 8349, 10150),
            TI(3, 8349, 8349),
            TI(3, 8349, 10400),
            TI(3, 8349, 10650),
            TI(3, 8349, 10800),
            TI(3, 8499, 10800),
            TI(3, 8749, 10800),
            TI(3, 8749, 8749),
            TI(3, 8999, 10800),
            TI(3, 9249, 10800),
            TI(3, 9549, 10800),
            TI(3, 9799, 10800),
            TI(3, 10049, 10800),
            TI(3, 10299, 10800),
            TI(3, 10549, 10800),
        ]);
        auto normalizedRegion = R([
            TI(3, 8349, 10800),
        ]);

        assert(region == normalizedRegion);
    }

    /// Computes the union of all tagged intervals.
    Region opBinary(string op)(in Region other) const if (op == "|")
    {
        return Region(this._intervals.dup ~ other._intervals.dup);
    }

    ///
    unittest
    {
        alias R = Region!(int, int);
        alias TI = R.TaggedInterval;

        assert((R(0, 10, 20) | R(0, 0, 5)) == R([TI(0, 10, 20), TI(0, 0, 5)]));
        assert((R(0, 10, 20) | R(0, 5, 15)) == R(0, 5, 20));
        assert((R(0, 10, 20) | R(0, 12, 18)) == R(0, 10, 20));
        assert((R(0, 10, 20) | R(0, 10, 20)) == R(0, 10, 20));
        assert((R(0, 10, 20) | R(0, 15, 25)) == R(0, 10, 25));
        assert((R(0, 10, 20) | R(0, 25, 30)) == R([TI(0, 10, 20), TI(0, 25, 30)]));
        assert((R(0, 10, 20) | R(1, 25, 30)) == R([TI(0, 10, 20), TI(1, 25, 30)]));
    }

    /// Computes the intersection of the two regions.
    Region opBinary(string op)(in Region other) const if (op == "&")
    {
        if (this.empty || other.empty)
        {
            return Region();
        }

        auto intersectionAcc = appender!(TaggedInterval[]);
        intersectionAcc.reserve(this._intervals.length + other._intervals.length);

        foreach (lhsInterval; this._intervals)
        {
            // dfmt off
            auto intersectingRhs = other._intervals
                .assumeSorted!"a.isStrictlyBefore(b)"
                .equalRange(lhsInterval);
            // dfmt on

            if (intersectingRhs.empty)
            {
                continue;
            }

            TaggedInterval tmpIntersection = lhsInterval;

            foreach (rhsInterval; intersectingRhs)
            {
                tmpIntersection &= rhsInterval;

                if (tmpIntersection.empty)
                {
                    // Intersection is empty; stop trying rhsIntervals.
                    break;
                }
            }

            if (!tmpIntersection.empty)
            {
                intersectionAcc ~= tmpIntersection;
            }
        }

        return Region(intersectionAcc.data);
    }

    ///
    unittest
    {
        alias R = Region!(int, int);
        alias TI = R.TaggedInterval;

        assert((R(0, 10, 20) & R(0, 0, 5)) == R([]));
        assert((R(0, 10, 20) & R(0, 5, 15)) == R(0, 10, 15));
        assert((R(0, 10, 20) & R(0, 12, 18)) == R(0, 12, 18));
        assert((R(0, 10, 20) & R(0, 10, 20)) == R(0, 10, 20));
        assert((R(0, 10, 20) & R(0, 15, 25)) == R(0, 15, 20));
        assert((R(0, 10, 20) & R(0, 25, 30)) == R([]));
        assert((R(0, 10, 20) & R(1, 25, 30)) == R([]));
    }

    Region opBinary(string op)(in TaggedInterval interval) const if (op == "-")
    {
        if (interval.empty)
        {
            return Region(this._intervals.dup);
        }

        auto differenceAcc = appender!(TaggedInterval[]);
        differenceAcc.reserve(this._intervals.length + 1);

        foreach (lhsInterval; this._intervals)
        {
            if (lhsInterval.intersects(interval))
            {
                auto tmpDifference = lhsInterval - interval;

                differenceAcc ~= tmpDifference._intervals;
            }
            else
            {
                differenceAcc ~= lhsInterval;
            }
        }

        return Region(differenceAcc.data);
    }

    Region opBinary(string op)(in Region other) const if (op == "-")
    {
        if (other.empty)
        {
            return Region(this._intervals.dup);
        }

        auto differenceAcc = appender!(TaggedInterval[]);
        differenceAcc.reserve(this._intervals.length + other._intervals.length);

        foreach (lhsInterval; this._intervals)
        {
            // dfmt off
            auto intersectingRhs = other._intervals
                .assumeSorted!"a.isStrictlyBefore(b)"
                .equalRange(lhsInterval);
            // dfmt on

            auto tmpDifference = Region(lhsInterval);

            foreach (rhsInterval; intersectingRhs)
            {
                if (tmpDifference.empty
                        || tmpDifference._intervals[$ - 1].isStrictlyBefore(rhsInterval))
                {
                    // Remaining rhsItervals will not intersect anymore.
                    break;
                }

                tmpDifference -= rhsInterval;
            }

            if (!tmpDifference.empty)
            {
                differenceAcc ~= tmpDifference._intervals;
            }
        }

        return Region(differenceAcc.data);
    }

    ///
    unittest
    {
        alias R = Region!(int, int);
        alias TI = R.TaggedInterval;

        assert((R(0, 10, 20) - R(0, 0, 5)) == R(0, 10, 20));
        assert((R(0, 10, 20) - R(0, 5, 15)) == R(0, 15, 20));
        assert((R(0, 10, 20) - R(0, 12, 18)) == R([TI(0, 10, 12), TI(0, 18, 20)]));
        assert((R(0, 10, 20) - R(0, 10, 20)).empty);
        assert((R(0, 10, 20) - R(0, 15, 25)) == R(0, 10, 15));
        assert((R(0, 10, 20) - R(0, 25, 30)) == R(0, 10, 20));
        assert((R(0, 10, 20) - R(1, 25, 30)) == R(0, 10, 20));
    }

    Region opOpAssign(string op, T)(in T other)
            if (is(T : Region) || is(T : TaggedInterval))
    {
        static if (op == "|" || op == "&" || op == "-")
        {
            mixin("auto tmp = this " ~ op ~ " other;");

            this._intervals = tmp._intervals;

            return this;
        }
        else
        {
            static assert(0, "unsupported operator: " ~ op);
        }
    }

    unittest
    {
        alias R = Region!(int, int);
        alias TI = R.TaggedInterval;

        R accRegion;
        auto inputRegion1 = R([
            TI(3, 8349, 8600),
            TI(3, 8349, 8850),
            TI(3, 8349, 9100),
            TI(3, 8349, 9350),
            TI(3, 8349, 9600),
            TI(3, 8349, 9900),
            TI(3, 8349, 10150),
            TI(3, 8349, 10400),
            TI(3, 8349, 10650),
            TI(3, 8349, 10800),
            TI(3, 8499, 10800),
            TI(3, 8749, 10800),
            TI(3, 8999, 10800),
            TI(3, 9249, 10800),
            TI(3, 9549, 10800),
            TI(3, 9799, 10800),
            TI(3, 10049, 10800),
            TI(3, 10299, 10800),
            TI(3, 10549, 10800),
        ]);
        auto inputRegion2 = R([
            TI(3, 2297, 11371),
        ]);
        auto expectedResult = R([
            TI(3, 2297, 11371),
        ]);

        accRegion |= inputRegion1;

        assert(accRegion == inputRegion1);

        accRegion |= inputRegion2;

        assert(accRegion == expectedResult);
        assert((inputRegion1 | inputRegion2) == expectedResult);
    }
}

unittest
{
    Region!(int, int) r;
}

/**
    Returns true iff `thing` is empty

    See_Also: Region.empty, Region.TaggedInterval.empty
*/
bool empty(T)(in T thing) pure nothrow
        if (is(T : Region!Args, Args...) || is(T : Region!Args.TaggedInterval, Args...))
{
    return thing.empty;
}

///
unittest
{
    alias R = Region!(int, int);
    alias TI = R.TaggedInterval;

    R emptyRegion;
    TI emptyTI;
    auto region = R(0, 0, 10);
    auto ti = TI(0, 0, 10);

    assert(empty(emptyRegion));
    assert(empty(emptyTI));
    assert(!empty(region));
    assert(!empty(ti));
}

/**
    Returns the union of all elements.

    See_Also: Region.opBinary!"|", Region.TaggedInterval.opBinary!"|"
*/
auto union_(Range)(Range regions)
        if (isInputRange!Range && is(ElementType!Range : Region!Args, Args...))
{
    alias Region = Unqual!(ElementType!Range);

    // dfmt off
    return Region(regions
        .map!"a._intervals.dup"
        .join);
    // dfmt on
}

///
unittest
{
    alias R = Region!(int, int);
    alias TI = R.TaggedInterval;

    R emptyRegion;
    TI emptyTI;
    auto region = R(0, 0, 10);
    auto ti = TI(0, 0, 10);

    assert(empty(emptyRegion));
    assert(empty(emptyTI));
    assert(!empty(region));
    assert(!empty(ti));
}
