/**
    Some additional mathematical functions.

    Copyright: © 2018 Arne Ludwig <arne.ludwig@posteo.de>
    License: Subject to the terms of the MIT license, as written in the
             included LICENSE file.
    Authors: Arne Ludwig <arne.ludwig@posteo.de>
*/
module djunctor.util.math;

import std.algorithm : sort, sum;
import std.conv : to;
import std.range : ElementType, isForwardRange, walkLength;

/// Caluclate the mean of range.
ElementType!Range mean(Range)(Range values) if (isForwardRange!Range)
{
    auto sum = values.sum;
    auto length = values.length.to!(ElementType!Range);

    return sum / length;
}

unittest
{
    {
        auto values = [2, 4, 8];
        assert(values.mean == 4);
    }
    {
        auto values = [1.0, 2.0, 3.0, 4.0];
        assert(values.mean == 2.5);
    }
}

/// Caluclate the median of range.
ElementType!Range median(Range)(Range values) if (__traits(compiles, sort(values)))
{
    assert(values.length > 0, "median is undefined for empty set");
    auto middleIdx = values.length / 2;
    auto useAverage = values.length > 1 && values.length % 2 == 0;
    auto sortedValues = values.sort;

    if (useAverage)
        return (sortedValues[middleIdx] + sortedValues[middleIdx + 1]) / 2;
    else
        return values[middleIdx];
}

unittest
{
    {
        auto values = [4, 2, 8];
        assert(values.median == 4);
    }
    {
        auto values = [2, 1, 3, 0, 4, 9, 8, 5, 6, 3, 9];
        assert(values.median == 4);
    }
    {
        auto values = [2.0, 1.0, 4.0, 3.0, 5.0];
        assert(values.median == 3.0);
    }
}