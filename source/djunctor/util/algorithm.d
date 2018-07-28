/**
    Some additional alogorithm functions.

    Copyright: © 2018 Arne Ludwig <arne.ludwig@posteo.de>
    License: Subject to the terms of the MIT license, as written in the
             included LICENSE file.
    Authors: Arne Ludwig <arne.ludwig@posteo.de>
*/
module djunctor.util.algorithm;

import std.conv : to;
import std.functional : binaryFun, unaryFun;
import std.traits : Unqual;

/**
    Order `a` and `b` lexicographically by applying each `fun` to them. For
    unary functions compares `fun(a) < fun(b)`.
*/
bool orderLexicographically(T, fun...)(T a, T b) pure nothrow
{
    static foreach (i, getFieldValue; fun)
    {
        {
            auto aValue = unaryFun!getFieldValue(a);
            auto bValue = unaryFun!getFieldValue(b);

            if (aValue != bValue)
            {
                return aValue < bValue;
            }
        }
    }

    return false;
}

/**
    Compare `a` and `b` lexicographically by applying each `fun` to them. For
    unary functions compares `fun(a) < fun(b)`.
*/
int cmpLexicographically(T, fun...)(T a, T b) pure nothrow
{
    static foreach (i, getFieldValue; fun)
    {
        {
            auto aValue = unaryFun!getFieldValue(a);
            auto bValue = unaryFun!getFieldValue(b);

            if (aValue < bValue)
            {
                return -1;
            }
            else if (aValue > bValue)
            {
                return 1;
            }
        }
    }

    return 0;
}

/**
    Slices an input array into slices of equivalent adjacent elements.
    In other languages this is often called `partitionBy`, `groupBy`
    or `sliceWhen`.

    Equivalence is defined by the predicate `pred`, which can be binary,
    which is passed to `std.functional.binaryFun`. Two range elements
    `a` and `b` are considered equivalent if `pred(a,b)` is true.

    This predicate must be an equivalence relation, that is, it must be
    reflexive (`pred(x,x)` is always true), symmetric
    (`pred(x,y) == pred(y,x)`), and transitive (`pred(x,y) && pred(y,z)`
    implies `pred(x,z)`). If this is not the case, the range returned by
    chunkBy may assert at runtime or behave erratically.

    Params:
     pred = Predicate for determining equivalence.
     r = An array to be sliced.

    Returns: With a binary predicate, a range of slices is returned in which
    all elements in a given slice are equivalent under the given predicate.

    Notes:

    Equivalent elements separated by an intervening non-equivalent element will
    appear in separate subranges; this function only considers adjacent
    equivalence. Elements in the subranges will always appear in the same order
    they appear in the original range.
*/
auto sliceBy(alias pred, Array)(ref Array array) pure nothrow
        if (is(Unqual!Array == T[], T))
{
    return SliceByImpl!(pred, Array)(array);
}

///
unittest
{
    import std.algorithm.comparison : equal;

    // Grouping by particular attribute of each element:
    auto data = [
        [1, 1],
        [1, 2],
        [2, 2],
        [2, 3]
    ];

    auto r1 = data.sliceBy!((a,b) => a[0] == b[0]);
    assert(r1.equal([
        data[0 .. 2],
        data[2 .. 4]
    ]));

    auto r2 = data.sliceBy!((a,b) => a[1] == b[1]);
    assert(r2.equal([
        data[0 .. 1],
        data[1 .. 3],
        data[3 .. 4],
    ]));
}

private struct SliceByImpl(alias pred, Array)
        if (is(Unqual!Array == T[], T))
{
    private alias equivalent = binaryFun!pred;

    private Array array;
    private size_t sliceStart;
    private size_t sliceEnd;

    this(ref Array array)
    {
        this.array = array;

        if (!empty)
        {
            popFront();
        }
    }

    void popFront()
    {
        assert(!empty, "Attempting to popFront an empty SliceByImpl");

        sliceStart = sliceEnd++;

        if (empty)
        {
            return;
        }

        auto refElement = array[sliceStart];
        while (sliceEnd < array.length && equivalent(refElement, array[sliceEnd]))
        {
            ++sliceEnd;
        }
    }

    @property bool empty() const pure nothrow
    {
        return sliceStart >= array.length;
    }

    @property auto front()
    {
        assert(!empty, "Attempting to fetch the front of an empty SliceByImpl");

        return array[sliceStart .. sliceEnd];
    }
}
