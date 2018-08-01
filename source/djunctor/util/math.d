/**
    Some additional mathematical functions.

    Copyright: © 2018 Arne Ludwig <arne.ludwig@posteo.de>
    License: Subject to the terms of the MIT license, as written in the
             included LICENSE file.
    Authors: Arne Ludwig <arne.ludwig@posteo.de>
*/
module djunctor.util.math;

import djunctor.util.algorithm : cmpLexicographically, sliceBy;
import std.algorithm : all, copy, filter, map, max, sort, sum, swap, uniq;
import std.array : Appender, array;
import std.conv : to;
import std.exception : assertThrown;
import std.functional : binaryFun, unaryFun;
import std.range : assumeSorted, chain, ElementType, enumerate, isForwardRange,
    retro, walkLength;
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
struct Graph(Node, Weight = void, Flag!"isDirected" isDirected = No.isDirected, EdgePayload = void)
{
    static immutable isWeighted = !is(Weight == void);
    static immutable hasEdgePayload = !is(EdgePayload == void);

    static struct Edge
    {
        protected Node _start;
        protected Node _end;

        static if (isWeighted)
            Weight weight;

        static if (hasEdgePayload)
            EdgePayload payload;

        /// Construct an edge.
        this(Node start, Node end)
        {
            this._start = start;
            this._end = end;

            static if (!isDirected)
            {
                if (end < start)
                {
                    swap(this._start, this._end);
                }
            }
        }

        static if (isWeighted)
        {
            /// ditto
            this(Node start, Node end, Weight weight)
            {
                this(start, end);
                this.weight = weight;
            }
        }

        static if (hasEdgePayload && !is(EdgePayload : Weight))
        {
            /// ditto
            this(Node start, Node end, EdgePayload payload)
            {
                this(start, end);
                this.payload = payload;
            }
        }

        static if (isWeighted && hasEdgePayload)
        {
            /// ditto
            this(Node start, Node end, Weight weight, EdgePayload payload)
            {
                this(start, end);
                this.weight = weight;
                this.payload = payload;
            }
        }

        /// Get the start of this edge. For undirected graphs this is the
        /// smaller of both incident nodes.
        @property Node start() const pure nothrow
        {
            return _start;
        }

        /// Get the end of this edge. For undirected graphs this is the
        /// larger of both incident nodes.
        @property Node end() const pure nothrow
        {
            return _end;
        }

        /**
            Get target of this edge beginning at node `from`. For undirected
            graphs returns the other node of this edge.

            Throws: MissingNodeException if this edge does not start in node `from`.
        */
        Node target(Node from) const
        {
            static if (isDirected)
            {
                if (start == from)
                {
                    return end;
                }
                else
                {
                    throw new MissingNodeException();
                }
            }
            else
            {
                if (start == from)
                {
                    return end;
                }
                else if (end == from)
                {
                    return start;
                }
                else
                {
                    throw new MissingNodeException();
                }
            }
        }

        /**
            Get source of this edge beginning at node `from`. For undirected
            graphs returns the other node of this edge.

            Throws: MissingNodeException if this edge does not end in node `from`.
        */
        static if (isDirected)
        {
            Node source(Node from) const
            {
                if (end == from)
                {
                    return start;
                }
                else
                {
                    throw new MissingNodeException();
                }
            }
        }
        else
        {
            alias source = target;
        }

        /// Two edges are equal iff their incident nodes (and weight) are the
        /// same.
        bool opEquals(in Edge other) const pure nothrow
        {
            static if (isWeighted)
            {
                return this.start == other.start && this.end == other.end
                    && this.weight == other.weight;
            }
            else
            {
                return this.start == other.start && this.end == other.end;
            }
        }

        /// Orders edge lexicographically by `start`, `end`(, `weight`).
        int opCmp(in Edge other) const pure nothrow
        {
            static if (isWeighted)
            {
                // dfmt off
                return cmpLexicographically!(
                    typeof(this),
                    "a.start",
                    "a.end",
                    "a.weight",
                )(this, other);
                // dfmt on
            }
            else
            {
                // dfmt off
                return cmpLexicographically!(
                    typeof(this),
                    "a.start",
                    "a.end",
                )(this, other);
                // dfmt on
            }
        }

        private int compareNodes(in Edge other) const pure nothrow
        {
            // dfmt off
            return cmpLexicographically!(
                typeof(this),
                "a.start",
                "a.end",
            )(this, other);
            // dfmt on
        }

        /**
            Returns the node that is connects this edge with other edge. In
            case of undirected graphs this is just the common node of both
            edges; in directed case this is the end node of this edge if it
            matches the start node of other edge.

            Throws: MissingNodeException if the connecting node is undefined.
        */
        Node getConnectingNode(in Edge other) const
        {
            static if (isDirected)
            {
                if (this.end == other.start)
                {
                    return this.end;
                }
            }
            else
            {
                if (this.end == other.start || this.end == other.end)
                {
                    return this.end;
                }
                else if (this.start == other.start || this.start == other.end)
                {
                    return this.start;
                }
            }

            throw new MissingNodeException();
        }
    }

    static bool orderByNodes(in Edge a, in Edge b) nothrow pure
    {
        return a.compareNodes(b) < 0;
    }

    static bool groupByNodes(in Edge a, in Edge b) nothrow pure
    {
        return a.compareNodes(b) == 0;
    }

    /// Construct an edge for this graph.
    static Edge edge(T...)(T args)
    {
        return Edge(args);
    }

    protected Node[] _nodes;
    protected Appender!(Edge[]) _edges;

    /// The set (ordered list) of nodes.
    @property const(Node[]) nodes() const nothrow pure
    {
        return _nodes;
    }

    private @property void nodes(Node[] nodes)
    {
        nodes.sort();

        this._nodes = nodes.uniq.array;
    }

    /// Get the set (ordered list) of edges in this graph.
    @property auto edges() nothrow pure
    {
        return chain(_edges.data);
    }

    /// ditto
    @property auto edges() const nothrow pure
    {
        return chain(_edges.data);
    }

    /**
        Construct a graph from a set of nodes (and edges). Modifies `nodes`
        while sorting but releases it after construction.

        Throws: MissingNodeException if an edge has a node that is not present
                in this graph .
        Throws: EdgeExistsException if an edge already exists when trying
                inserting it.
    */
    this(Node[] nodes)
    {
        this.nodes = nodes;
    }

    /// ditto
    this(Node[] nodes, Edge[] edges)
    {
        this(nodes);

        _edges.reserve(edges.length);
        foreach (edge; edges)
        {
            this.add(edge);
        }
    }

    this(this)
    {
        _nodes = _nodes.dup;
    }

    /// Add a set of edges to this graph without any checks.
    void bulkAdd(Edge[] edges)
    {
        this._edges ~= edges;
        _edges.data.sort;
    }

    /// Add a set of edges to this graph and merge mutli-edges using `merge`.
    void addAndMerge(alias merge)(Edge[] edges)
    {
        //static assert(is(typeof(merge(edges)) == Edge), "expected `Edge merge(Edge[] multiEdge)`");
        bulkAdd(edges);

        // dfmt off
        auto bufferRest = _edges
            .data
            .sliceBy!groupByNodes
            .map!(unaryFun!merge)
            .copy(_edges.data);
        // dfmt on
        _edges.shrinkTo(_edges.data.length - bufferRest.length);
    }

    ///
    unittest
    {
        auto g1 = Graph!(int, int)([1, 2]);

        static g1.Edge sumWeights(g1.Edge[] multiEdge)
        {
            auto sumOfWeights = multiEdge.map!"a.weight".sum;
            auto mergedEdge = multiEdge[0];
            mergedEdge.weight = sumOfWeights;

            return mergedEdge;
        }

        // dfmt off
        auto edges = [
            g1.edge(1, 2, 1),
            g1.edge(1, 2, 1),
            g1.edge(1, 2, 1),
            g1.edge(2, 3, 2),
            g1.edge(2, 3, 2),
            g1.edge(3, 4, 3),
        ];
        g1.addAndMerge!sumWeights(edges);
        assert(g1.edges == [
            g1.edge(1, 2, 3),
            g1.edge(2, 3, 4),
            g1.edge(3, 4, 3),
        ]);
        // dfmt on
    }

    /// Add an edge to this graph and handle existing edges with `handleConflict`.
    /// The handler must have this signature `Edge handleConflict(Edge, Edge)`.
    Edge add(alias handleConflict = ConflictStrategy.error)(Edge edge)
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
            return forceAdd(edge);
        }
        else
        {
            auto newEdge = binaryFun!handleConflict(existingEdges.front, edge);

            return replaceEdge(existingEdgeIdx, newEdge);
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

            ///
            unittest
            {
                auto g1 = Graph!(int, int)([1, 2]);
                alias CS = g1.ConflictStrategy;

                g1 ~= g1.edge(1, 2, 1);

                auto addedEdge = g1.add!(CS.sumWeights)(g1.edge(1, 2, 1));

                assert(addedEdge.weight == 2);
            }
        }

        /// Throw `EdgeExistsException`.
        static inout(Edge) error(inout(Edge) existingEdge, inout(Edge) newEdge)
        {
            throw new EdgeExistsException();
        }

        ///
        unittest
        {
            auto g1 = Graph!int([1, 2]);
            alias CS = g1.ConflictStrategy;

            g1 ~= g1.edge(1, 2);

            assertThrown!EdgeExistsException(g1.add!(CS.error)(g1.edge(1, 2)));
        }

        /// Replace the existingEdge by newEdge.
        static inout(Edge) replace(inout(Edge) existingEdge, inout(Edge) newEdge)
        {
            return newEdge;
        }

        ///
        unittest
        {
            auto g1 = Graph!(int, int)([1, 2]);
            alias CS = g1.ConflictStrategy;

            g1 ~= g1.edge(1, 2, 1);

            auto addedEdge = g1.add!(CS.replace)(g1.edge(1, 2, 2));

            assert(addedEdge.weight == 2);
        }

        /// Keep existingEdge – discard newEdge.
        static inout(Edge) keep(inout(Edge) existingEdge, inout(Edge) newEdge)
        {
            return existingEdge;
        }

        ///
        unittest
        {
            auto g1 = Graph!(int, int)([1, 2]);
            alias CS = g1.ConflictStrategy;

            g1 ~= g1.edge(1, 2, 1);

            auto addedEdge = g1.add!(CS.keep)(g1.edge(1, 2, 2));

            assert(addedEdge.weight == 1);
        }
    }

    /// Forcibly add an edge to this graph.
    protected Edge forceAdd(Edge edge)
    {
        _edges ~= edge;
        _edges.data.sort;

        return edge;
    }

    /// Replace an edge in this graph.
    protected Edge replaceEdge(in size_t edgeIdx, Edge newEdge)
    {
        auto shouldSort = _edges.data[edgeIdx] != newEdge;

        _edges.data[edgeIdx] = newEdge;

        if (shouldSort)
        {
            _edges.data.sort;
        }

        return newEdge;
    }

    void filterEdges(alias pred)()
    {
        // dfmt off
        auto bufferRest = _edges
            .data
            .filter!pred
            .copy(_edges.data);
        // dfmt on
        _edges.shrinkTo(_edges.data.length - bufferRest.length);
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
    auto get(Edge edge)
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

    /// Returns the ndex of node `n` in the list of nodes.
    size_t indexOf(Node n) const
    {
        auto sortedNodes = assumeSorted(nodes);
        auto tristectedNodes = sortedNodes.trisect(n);

        if (tristectedNodes[1].empty)
        {
            throw new MissingNodeException();
        }

        return tristectedNodes[0].length;
    }

    ///
    unittest
    {

        auto g1 = Graph!(int, int)([1, 2]);

        assert(g1.indexOf(1) == 0);
        assert(g1.indexOf(2) == 1);
        assertThrown!MissingNodeException(g1.indexOf(3));
    }

    static if (isDirected)
    {
        /// Returns a range of in/outgoing edges of node `n`.
        auto inEdges(Node n) nothrow pure
        {
            return _edges.data[].filter!(e => e.end == n);
        }

        /// ditto
        auto inEdges(Node n) const nothrow pure
        {
            return edges[].filter!(e => e.end == n);
        }

        /// ditto
        auto outEdges(Node n) nothrow pure
        {
            return _edges.data[].filter!(e => e.start == n);
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
        auto incidentEdges(Node n) nothrow pure
        {
            return _edges.data[].filter!(e => e.start == n || e.end == n);
        }

        /// ditto
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

    //   +-+  +-+
    //   \ /  \ /
    //   (1)--(2)
    //
    // payload(1, 1) = [1];
    // payload(1, 2) = [2];
    // payload(2, 2) = [3];
    auto g5 = Graph!(int, void, No.isDirected, int[])([1, 2]);

    g5 ~= g5.edge(1, 1, [1]);
    g5 ~= g5.edge(1, 2, [2]);
    g5.add(g5.edge(2, 2, [3]));

    assert(g5.edge(1, 1) in g5);
    assert(g5.get(g5.edge(1, 1)).payload == [1]);
    assert(g5.edge(1, 2) in g5);
    assert(g5.get(g5.edge(1, 2)).payload == [2]);
    assert(g5.edge(2, 1) in g5);
    assert(g5.get(g5.edge(2, 1)).payload == [2]);
    assert(g5.has(g5.edge(2, 2)));
    assert(g5.get(g5.edge(2, 2)).payload == [3]);
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

class EmptySetException : Exception
{
    this(string msg)
    {
        super(msg);
    }
}

struct NaturalNumberSet
{
    private static immutable partSize = 8 * size_t.sizeof;
    private static immutable size_t firstBit = 1;
    private static immutable size_t lastBit = firstBit << (partSize - 1);
    private static immutable size_t emptyPart = 0;
    private static immutable size_t fullPart = ~emptyPart;

    private size_t[] parts;
    private size_t nMax;

    this(this)
    {
        parts = parts.dup;
    }

    private this(size_t[] parts)
    {
        this.parts = parts;
    }

    private bool inBounds(in size_t n) const pure nothrow
    {
        return n < nMax;
    }

    void reserveFor(in size_t n)
    {
        if (parts.length == 0)
        {
            parts.length = max(1, ceil(n, partSize) / partSize);
        }

        while (!inBounds(n))
        {
            parts.length *= 2;
            nMax = parts.length * partSize;
        }
    }

    @property size_t capacity() pure const nothrow
    {
        return nMax;
    }

    private size_t partIdx(in size_t n) const pure nothrow
    {
        return n / partSize;
    }

    private size_t idxInPart(in size_t n) const pure nothrow
    {
        return n % partSize;
    }

    private size_t itemMask(in size_t n) const pure nothrow
    {
        return firstBit << idxInPart(n);
    }

    static size_t inverse(in size_t n) pure nothrow
    {
        return n ^ fullPart;
    }

    void add(in size_t n)
    {
        reserveFor(n);

        parts[partIdx(n)] |= itemMask(n);
    }

    void remove(in size_t n)
    {
        if (!inBounds(n))
        {
            return;
        }

        parts[partIdx(n)] &= inverse(itemMask(n));
    }

    bool has(in size_t n) const pure nothrow
    {
        if (!inBounds(n))
        {
            return false;
        }

        return (parts[partIdx(n)] & itemMask(n)) != emptyPart;
    }

    bool empty() const pure nothrow
    {
        return parts.all!(part => part == emptyPart);
    }

    size_t minElement() const
    {
        foreach (i, part; parts)
        {
            if (part != emptyPart)
            {
                size_t j = 0;

                while (((part >> j) & firstBit) != firstBit)
                {
                    ++j;
                }

                return i * partSize + j;
            }
        }

        throw new EmptySetException("empty set has no minElement");
    }

    size_t maxElement() const
    {
        foreach (i, part; parts.retro.enumerate)
        {
            if (part != emptyPart)
            {
                size_t j = 0;

                while (((part << j) & lastBit) != lastBit)
                {
                    ++j;
                }

                return (parts.length - i - 1) * partSize + (partSize - j - 1);
            }
        }

        throw new EmptySetException("empty set has no maxElement");
    }

    unittest
    {
        foreach (i; 0 .. 2 * NaturalNumberSet.partSize)
        {
            NaturalNumberSet set;

            set.add(i + 5);
            set.add(i + 7);

            assert(set.minElement() == i + 5);
            assert(set.maxElement() == i + 7);
        }
    }

    /// Returns a range of the elements in this set. The elements are ordered
    /// ascending.
    @property auto elements() const pure nothrow
    {
        static struct ElementsRange
        {
            const NaturalNumberSet* set;
            bool _empty = false;
            size_t i = 0;
            size_t part;
            size_t j = 0;

            this(const NaturalNumberSet* set) pure nothrow
            {
                this.set = set;
                this._empty = set.empty;

                if (!this.empty)
                {
                    this.part = set.parts[i];
                    if (!set.has(front))
                    {
                        popFront();
                    }
                }
            }

            void popFront() pure nothrow
            {
                assert(!empty, "Attempting to popFront an empty elements range");
                ++j;

                while (shiftedPartEmpty)
                {
                    nextPart();

                    if (empty)
                    {
                        return;
                    }
                }

                while (((part >> j) & firstBit) != firstBit && !shiftedPartEmpty)
                {
                    ++j;
                }

                if (shiftedPartEmpty)
                {
                    popFront();
                }
            }

            @property size_t front() const pure nothrow
            {
                assert(!empty, "Attempting to fetch the front of an empty elements range");
                return i * partSize + j;
            }

            @property bool empty() const pure nothrow
            {
                return _empty;
            }

            private @property bool shiftedPartEmpty() const pure nothrow
            {
                return (part >> j) == emptyPart || j >= partSize;
            }

            private void nextPart() pure nothrow
            {
                // move to start of next part
                ++i;
                j = 0;

                if (i < set.parts.length)
                {
                    part = set.parts[i];
                }
                else
                {
                    _empty = true;
                }
            }
        }

        return ElementsRange(&this);
    }

    ///
    unittest
    {
        import std.algorithm : equal;
        import std.range : iota;

        NaturalNumberSet set;
        auto someNumbers = iota(set.partSize).filter!"a % 3 == 0";

        foreach (i; someNumbers)
        {
            set.add(i);
        }

        assert(equal(someNumbers, set.elements));
    }
}

unittest
{
    NaturalNumberSet set;

    // add some numbers
    foreach (i; 0 .. set.partSize)
    {
        if (i % 2 == 0)
        {
            set.add(i);
        }
    }

    // force extension of set
    foreach (i; set.partSize .. 2 * set.partSize)
    {
        if (i % 3 == 0)
        {
            set.add(i);
        }
    }

    // validate presence
    foreach (i; 0 .. 2 * set.partSize)
    {
        if (i / set.partSize == 0 && i % 2 == 0)
        {
            assert(set.has(i));
        }
        else if (i / set.partSize == 1 && i % 3 == 0)
        {
            assert(set.has(i));
        }
        else
        {
            assert(!set.has(i));
        }
    }
}
