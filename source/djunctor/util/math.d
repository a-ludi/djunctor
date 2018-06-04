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
import std.traits : isIntegral, isNumeric;

/// Calculate the mean of range.
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

/// Calculate the median of range.
ElementType!Range median(Range)(Range values) if (__traits(compiles, sort(values)))
{
    assert(values.length > 0, "median is undefined for empty set");
    auto middleIdx = values.length / 2;
    auto useAverage = values.length > 1 && values.length % 2 == 0;
    auto sortedValues = values.sort;

    if (useAverage)
        return (sortedValues[middleIdx - 1] + sortedValues[middleIdx]) / 2;
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
        auto values = [4, 3, 2, 8];
        assert(values.median == 3);
    }
    {
        auto values = [4, 6, 2, 8];
        assert(values.median == 5);
    }
    {
        auto values = [2, 1, 3, 0, 4, 9, 8, 5, 6, 3, 9];
        assert(values.median == 4);
    }
    {
        auto values = [2.0, 1.0, 4.0, 3.0];
        assert(values.median == 2.5);
    }
    {
        auto values = [2.0, 1.0, 4.0, 3.0, 5.0];
        assert(values.median == 3.0);
    }
}

/**
    Round x upward according to base, ie. returns the next integer larger or
    equal to x which is divisible by base.

    Returns: x rounded upward according to base.
*/
Integer ceil(Integer)(in Integer x, in Integer base) pure nothrow 
        if (isIntegral!Integer)
{
    // dfmt off
    return x % base == 0
        ? x
        : (x / base + 1) * base;
    // dfmt on
}

///
unittest
{
    assert(ceil(8, 10) == 10);
    assert(ceil(32, 16) == 32);
    assert(ceil(101, 100) == 200);
}

/**
    Round x downward according to base, ie. returns the next integer smaller or
    equal to x which is divisible by base.

    Returns: x rounded downward according to base.
*/
Integer floor(Integer)(in Integer x, in Integer base) pure nothrow 
        if (isIntegral!Integer)
{
    return (x / base) * base;
}

///
unittest
{
    assert(floor(8, 10) == 0);
    assert(floor(32, 16) == 32);
    assert(floor(101, 100) == 100);
}

/// Returns the absolute difference between two numbers.
Num absdiff(Num)(in Num a, in Num b) pure nothrow if (isNumeric!Num)
{
    // dfmt off
    return a > b
        ? a - b
        : b - a;
    // dfmt on
}

///
unittest
{
    assert(absdiff(2UL, 3UL) == 1UL);
    assert(absdiff(-42, 13) == 55);
    assert(absdiff(2.5, 5) == 2.5);
}
