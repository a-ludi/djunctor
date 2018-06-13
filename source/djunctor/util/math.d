/**
    Some additional mathematical functions.

    Copyright: © 2018 Arne Ludwig <arne.ludwig@posteo.de>
    License: Subject to the terms of the MIT license, as written in the
             included LICENSE file.
    Authors: Arne Ludwig <arne.ludwig@posteo.de>
*/
module djunctor.util.math;

import std.algorithm : filter, sort, sum, swap;
import std.array : Appender;
import std.conv : to;
import std.exception : assertThrown;
import std.functional : binaryFun;
import std.range : assumeSorted, ElementType, isForwardRange, walkLength;
import std.traits : isIntegral, isNumeric;
import std.typecons : Flag, No, Yes;

debug import std.stdio : writeln;

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

class EdgeExistsException : Exception
{
    this()
    {
        super("edge cannot be inserted: edge already exists");
    }
}

class MissingEdgeException : Exception
{
    this()
    {
        super("edge not found");
    }
}

class MissingNodeException : Exception
{
    this()
    {
        super("node not found");
    }
}

/// This structure represents a graph with optional edge
/// payloads. The graph is represented as a list of edges which is
/// particularly suited for sparse graphs. While the set of nodes is fixed the
/// set of edges is mutable.
struct Graph(Node, Weight=void, Flag!"isDirected" isDirected = No.isDirected)
{
    static immutable isWeighted = is(typeof(Weight()));

    static struct Edge
    {
        Node start;
        Node end;

        static if (isWeighted)
            Weight weight;

        this(Node start, Node end)
        {
            this.start = start;
            this.end = end;

            static if (!isDirected)
            {
                if (end < start)
                {
                    swap(this.start, this.end);
                }
            }
        }

        static if (isWeighted)
        {

            this(Node start, Node end, Weight weight)
            {
                this(start, end);
                this.weight = weight;
            }
        }

        bool opEquals(in Edge other) const pure nothrow
        {
            static if (isWeighted)
            {
                return this.start == other.start && this.end == other.end && this.weight == other.weight;
            }
            else
            {
                return this.start == other.start && this.end == other.end;
            }
        }

        int opCmp(in Edge other) const pure nothrow
        {
            // dfmt off
            immutable cmpFields = isWeighted
                ? ["start", "end", "weight"]
                : ["start", "end"];
            // dfmt on

            static foreach(field; cmpFields)
            {
                if (mixin("this."~field~" < other."~field))
                {
                    return -1;
                }
                else if (mixin("this."~field~" > other."~field))
                {
                    return 1;
                }
            }

            return 0;
        }

        private int compareNodes(in Edge other) const pure nothrow
        {
            immutable cmpFields = ["start", "end"];

            static foreach(field; cmpFields)
            {
                if (mixin("this."~field~" < other."~field))
                {
                    return -1;
                }
                else if (mixin("this."~field~" > other."~field))
                {
                    return 1;
                }
            }

            return 0;
        }
    }

    protected static bool orderByNodes(in Edge a, in Edge b) nothrow pure
    {
        return a.compareNodes(b) < 0;
    }

    static Edge edge(T...)(T args)
    {
        return Edge(args);
    }

    const(Node[]) nodes;
    protected Appender!(Edge[]) _edges;

    @property const(Edge[]) edges() const nothrow pure
    {
        return _edges.data;
    }

    this(in Node[] nodes)
    {
        Node[] sortedNodes = nodes.dup;

        sortedNodes.sort();

        this.nodes = sortedNodes;
    }

    /// Add an edge to this graph and handle existing edges with `handleConflict`.
    /// The handler must have this signature `Edge handleConflict(Edge, Edge)`.
    void add(alias handleConflict = ConflictStrategy.error)(Edge edge)
    {
        if (!has(edge.start) || !has(edge.end))
        {
            throw new MissingNodeException();
        }

        auto sortedEdges = assumeSorted!orderByNodes(_edges.data);
        auto trisectedEdges = sortedEdges.trisect(edge);
        auto existingEdges = trisectedEdges[1];
        auto existingEdgeIdx = trisectedEdges[0].length;

        if (existingEdges.empty)
        {
            forceAdd(edge);
        }
        else
        {
            auto newEdge = binaryFun!handleConflict(existingEdges.front, edge);

            _edges.data[existingEdgeIdx] = newEdge;
        }
    }

    /// ditto
    void opOpAssign(string op)(Edge edge) if (op == "~")
    {
        add(edge);
    }

    ///
    unittest
    {
        auto g1 = Graph!(int, int)([1, 2]);

        auto e1 = g1.edge(1, 2, 1);
        auto e2 = g1.edge(1, 2, 2);

        g1 ~= e1;

        assertThrown!EdgeExistsException(g1.add(e2));

        with (g1.ConflictStrategy)
        {
            g1.add!replace(e2);

            assert(g1.get(g1.edge(1, 2)) == e2);

            g1.add!keep(e1);

            assert(g1.get(g1.edge(1, 2)) == e2);

            g1.add!sumWeights(e2);

            assert(g1.get(g1.edge(1, 2)).weight == 2 * e2.weight);
        }
    }

    /// Some pre-defined conflict handlers for `add`.
    static struct ConflictStrategy
    {
        static if (isWeighted)
        {
            /// Return an edge with sum of both weights. If given payload will be
            /// kept from existingEdge .
            static Edge sumWeights(Edge existingEdge, Edge newEdge)
            {
                existingEdge.weight += newEdge.weight;

                return existingEdge;
            }
        }

        /// Throw `EdgeExistsException`.
        static Edge error(Edge existingEdge, Edge newEdge)
        {
            throw new EdgeExistsException();
        }

        /// Replace the existingEdge by newEdge.
        static Edge replace(Edge existingEdge, Edge newEdge)
        {
            return newEdge;
        }

        /// Keep existingEdge – discard newEdge.
        static Edge keep(Edge existingEdge, Edge newEdge)
        {
            return existingEdge;
        }
    }

    /// Forcibly add an edge to this graph.
    protected void forceAdd(Edge edge)
    {
        _edges ~= edge;
        _edges.data.sort;
    }

    /// Check if edge/node exists in this graph. Ignores the weight if weighted.
    bool opBinaryRight(string op)(Node node) const pure nothrow if (op == "in")
    {
        auto sortedNodes = assumeSorted(nodes);

        return sortedNodes.contains(node);
    }

    /// ditto
    bool has(Node node) const pure nothrow
    {
        return node in this;
    }

    /// Check if edge exists in this graph. Only the `start` and `end` node
    /// will be compared.
    bool opBinaryRight(string op)(Edge edge) const pure nothrow if (op == "in")
    {
        auto sortedEdges = assumeSorted!orderByNodes(edges);

        return sortedEdges.contains(edge);
    }

    /// ditto
    bool has(Edge edge) const pure nothrow
    {
        return edge in this;
    }

    /// Get the designated edge from this graph. Only the `start` and `end`
    /// node will be compared.
    Edge get(Edge edge)
    {
        auto sortedEdges = assumeSorted!orderByNodes(edges);
        auto existingEdges = sortedEdges.equalRange(edge);

        if (existingEdges.empty)
        {
            throw new MissingEdgeException();
        }
        else
        {
            return existingEdges.front;
        }
    }

    ///
    unittest
    {
        auto g1 = Graph!(int, int)([1, 2]);

        auto e1 = g1.edge(1, 2, 1);

        g1 ~= e1;

        assert(g1.get(g1.edge(1, 2)) == e1);
        assertThrown!MissingEdgeException(g1.get(g1.edge(1, 1)));
    }

    static if (isDirected)
    {
        /// Returns a range of in/outgoing edges of node `n`.
        auto inEdges(Node n) const nothrow pure
        {
            return edges[].filter!(e => e.end == n);
        }

        /// ditto
        auto outEdges(Node n) const nothrow pure
        {
            return edges[].filter!(e => e.start == n);
        }

        ///
        unittest
        {
            import std.algorithm : equal;

            auto g1 = Graph!(int, void, Yes.isDirected)([1, 2, 3]);

            g1 ~= g1.edge(1, 1);
            g1 ~= g1.edge(1, 2);
            g1 ~= g1.edge(2, 2);
            g1 ~= g1.edge(2, 3);

            // dfmt off
            assert(g1.inEdges(1).equal([
                g1.edge(1, 1),
            ]));
            assert(g1.outEdges(1).equal([
                g1.edge(1, 1),
                g1.edge(1, 2),
            ]));
            assert(g1.inEdges(2).equal([
                g1.edge(1, 2),
                g1.edge(2, 2),
            ]));
            assert(g1.outEdges(2).equal([
                g1.edge(2, 2),
                g1.edge(2, 3),
            ]));
            assert(g1.inEdges(3).equal([
                g1.edge(2, 3),
            ]));
            assert(g1.outEdges(3).empty);
            // dfmt on
        }

        /// Get the in/out degree of node `n`.
        size_t inDegree(Node n) const nothrow pure
        {
            return inEdges(n).walkLength;
        }

        /// ditto
        size_t outDegree(Node n) const nothrow pure
        {
            return outEdges(n).walkLength;
        }

        ///
        unittest
        {
            auto g1 = Graph!(int, void, Yes.isDirected)([1, 2, 3]);

            g1 ~= g1.edge(1, 1);
            g1 ~= g1.edge(1, 2);
            g1 ~= g1.edge(2, 2);
            g1 ~= g1.edge(2, 3);

            assert(g1.inDegree(1) == 1);
            assert(g1.outDegree(1) == 2);
            assert(g1.inDegree(2) == 2);
            assert(g1.outDegree(2) == 2);
            assert(g1.inDegree(3) == 1);
            assert(g1.outDegree(3) == 0);
        }
    }
    else
    {
        /// Returns a range of all edges incident to node `n`.
        auto incidentEdges(Node n) const nothrow pure
        {
            return edges[].filter!(e => e.start == n || e.end == n);
        }

        /// ditto
        alias inEdges = incidentEdges;

        /// ditto
        alias outEdges = incidentEdges;

        ///
        unittest
        {
            import std.algorithm : equal;

            auto g1 = Graph!int([1, 2, 3]);

            g1 ~= g1.edge(1, 1);
            g1 ~= g1.edge(1, 2);
            g1 ~= g1.edge(2, 2);
            g1 ~= g1.edge(2, 3);

            // dfmt off
            assert(g1.incidentEdges(1).equal([
                g1.edge(1, 1),
                g1.edge(1, 2),
            ]));
            assert(g1.incidentEdges(2).equal([
                g1.edge(1, 2),
                g1.edge(2, 2),
                g1.edge(2, 3),
            ]));
            assert(g1.incidentEdges(3).equal([
                g1.edge(2, 3),
            ]));
            // dfmt on
        }

        /// Get the degree of node `n`.
        size_t degree(Node n) const nothrow pure
        {
            return incidentEdges(n).walkLength;
        }

        /// ditto
        alias inDegree = degree;

        /// ditto
        alias outDegree = degree;
    }
}

///
unittest
{
    //   +-+  +-+
    //   \ /  \ /
    //   (1)--(2)
    auto g1 = Graph!int([1, 2]);

    g1 ~= g1.edge(1, 1);
    g1 ~= g1.edge(1, 2);
    g1.add(g1.edge(2, 2));

    assert(g1.edge(1, 1) in g1);
    assert(g1.edge(1, 2) in g1);
    assert(g1.edge(2, 1) in g1);
    assert(g1.has(g1.edge(2, 2)));


    //   0.5     0.5
    //   +-+     +-+
    //   \ /     \ /
    //   (1)-----(2)
    //       1.0
    auto g2 = Graph!(int, double)([1, 2]);

    g2 ~= g2.edge(1, 1, 0.5);
    g2 ~= g2.edge(1, 2, 1.0);
    g2.add(g2.edge(2, 2, 0.5));

    assert(g2.edge(1, 1) in g2);
    assert(g2.edge(1, 2) in g2);
    assert(g2.edge(2, 1) in g2);
    assert(g2.has(g2.edge(2, 2)));


    //   0.5     0.5
    //   +-+     +-+
    //   \ v     v /
    //   (1)---->(2)
    //       1.0
    auto g3 = Graph!(int, double, Yes.isDirected)([1, 2]);

    g3 ~= g3.edge(1, 1, 0.5);
    g3 ~= g3.edge(1, 2, 1.0);
    g3.add(g3.edge(2, 2, 0.5));

    assert(g3.edge(1, 1) in g3);
    assert(g3.edge(1, 2) in g3);
    assert(!(g3.edge(2, 1) in g3));
    assert(g3.has(g3.edge(2, 2)));


    //   +-+   +-+
    //   \ v   v /
    //   (1)-->(2)
    auto g4 = Graph!(int, void, Yes.isDirected)([1, 2]);

    g4 ~= g4.edge(1, 1);
    g4 ~= g4.edge(1, 2);
    g4.add(g4.edge(2, 2));

    assert(g4.edge(1, 1) in g4);
    assert(g4.edge(1, 2) in g4);
    assert(!(g4.edge(2, 1) in g4));
    assert(g4.has(g4.edge(2, 2)));
}

///
unittest
{
    //     -1     1         1
    // (1)----(2)---(3) (4)---(5) (6)
    size_t[] contigs = [1, 2, 3, 4, 5, 6];
    auto contigGraph = Graph!(size_t, int)([1, 2, 3, 4, 5, 6]);

    contigGraph.add(contigGraph.edge(1, 2, -1));
    contigGraph.add(contigGraph.edge(2, 3, 1));
    contigGraph.add(contigGraph.edge(4, 5, 1));

    foreach (contig; contigs)
    {
        assert(contigGraph.degree(contig) <= 2);
    }
}
