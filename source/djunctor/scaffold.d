/**
    Work with scaffold graphs.

    Copyright: © 2018 Arne Ludwig <arne.ludwig@posteo.de>
    License: Subject to the terms of the MIT license, as written in the
             included LICENSE file.
    Authors: Arne Ludwig <arne.ludwig@posteo.de>
*/
module djunctor.util.scaffold;

import djunctor.util.log;
import djunctor.util.math : Graph, MissingNodeException, NaturalNumberSet;
import std.algorithm : canFind, countUntil, equal, filter, joiner, map, sort, sum, swap;
import std.array : Appender, array;
import std.exception : assertThrown;
import std.functional : binaryFun;
import std.range : assumeSorted, iota, only, refRange, retro, walkLength;
import std.typecons : Flag, No, Tuple, tuple, Yes;
import vibe.data.json : toJson = serializeToJson;

debug import std.conv : to;
debug import std.stdio : writeln;

///
unittest
{
    //   contig 1      contig 2
    //
    //  o        o     o        o
    //          / e1 e2 \      / e4
    //    o -- o ------- o -- o
    //     \        e3         \
    //      \                   \ e5
    //   e10 \   ____________   /
    //        \ /   e6       \ /
    //    o -- o         o -- o
    //          \ e7 e8 /      \ e9
    //  o        o     o        o
    //
    //   contig 4      contig 3
    //
    alias J = Join!int;
    alias S = Scaffold!int;
    alias CN = ContigNode;
    alias CP = ContigPart;
    auto scaffold = buildScaffold!(sumPayloads!int, int)(5, [
        J(CN(1, CP.end), CN(1, CP.post ), 1), // e1
        J(CN(1, CP.end), CN(1, CP.post ), 1), // e1
        J(CN(2, CP.pre), CN(2, CP.begin), 1), // e2
        J(CN(1, CP.end), CN(2, CP.begin), 1), // e3
        J(CN(2, CP.end), CN(2, CP.post ), 1), // e4
        J(CN(2, CP.end), CN(3, CP.end  ), 1), // e5
        J(CN(4, CP.end), CN(3, CP.end  ), 1), // e6
        J(CN(4, CP.end), CN(4, CP.post ), 1), // e7
        J(CN(3, CP.pre), CN(3, CP.begin), 1), // e8
        J(CN(3, CP.end), CN(3, CP.post ), 1), // e9
        J(CN(4, CP.end), CN(1, CP.begin), 1), // e10
    ]).discardAmbiguousJoins!int;
    //
    //   contig 1      contig 2
    //
    //  o        o     o        o
    //          / e1 e2 \      / e4
    //    o -- o ------- o -- o
    //              e3
    //
    //    o -- o         o -- o
    //          \ e7 e8 /      \ e9
    //  o        o     o        o
    //
    //   contig 4      contig 3

    assert(J(CN(1, CP.end), CN(1, CP.post ))  in scaffold); // e1
    assert(J(CN(2, CP.pre), CN(2, CP.begin))  in scaffold); // e2
    assert(J(CN(1, CP.end), CN(2, CP.begin))  in scaffold); // e3
    assert(J(CN(2, CP.end), CN(2, CP.post ))  in scaffold); // e4
    assert(J(CN(2, CP.end), CN(3, CP.end  )) !in scaffold); // e5
    assert(J(CN(4, CP.end), CN(3, CP.end  )) !in scaffold); // e6
    assert(J(CN(4, CP.end), CN(4, CP.post ))  in scaffold); // e7
    assert(J(CN(3, CP.pre), CN(3, CP.begin))  in scaffold); // e8
    assert(J(CN(3, CP.end), CN(3, CP.post ))  in scaffold); // e9
    assert(J(CN(4, CP.end), CN(1, CP.begin)) !in scaffold); // e10

    assert(scaffold.get(J(CN(1, CP.end), CN(1, CP.post ))).payload == 2); // e1
}

/// Each contig has four designated parts where joins can start or end.
static enum ContigPart
{
    /// Designates a transcendent point *before* the contig where
    /// front extensions end.
    pre,
    /// Designates the begin of the contig.
    begin,
    /// Designates the end of the contig.
    end,
    /// Designates a transcendent point *after* the contig where
    /// back extensions end.
    post,
}

bool isReal(in ContigPart contigPart) pure nothrow
{
    return contigPart == ContigPart.begin || contigPart == ContigPart.end;
}
bool isTranscendent(in ContigPart contigPart) pure nothrow
{
    return contigPart == ContigPart.pre || contigPart == ContigPart.post;
}


/**
    A contig is represented by four `ContigNodes` in the scaffold graph.

    See_Also: ContigPart
*/
alias ContigNode = Tuple!(size_t, "contigId", ContigPart, "contigPart");
/// Represents a set of joins which conclude possibly several scaffolding
/// variants.
alias Scaffold(T) = Graph!(ContigNode, void, No.isDirected, T);
alias Join(T) = Scaffold!T.Edge;

Join!T opOpAssignPayloads(string op, T)(Join!T j1, Join!T j2) pure nothrow
{
    mixin("j1.payload " ~ op ~ " j2.payload;");

    return j1;
}

alias sumPayloads(T) = opOpAssignPayloads!("+=", T);
alias concatenatePayloads(T) = opOpAssignPayloads!("~=", T);

J dontJoin(J)(J j1, J j2) pure nothrow
{
    j1.payload = T.init;

    return j1;
}

/// Returns true iff join is a default edge of the scaffold graph.
bool isDefault(J)(in J join) pure nothrow
{
    // dfmt off
    return join.start.contigPart == ContigPart.begin &&
           join.end.contigPart == ContigPart.end &&
           join.start.contigId == join.end.contigId;
    // dfmt on
}

/// Returns true iff join is a gap edge of the scaffold graph.
bool isGap(J)(in J join) pure nothrow
{
    // dfmt off
    return join.start.contigId != join.end.contigId &&
        (join.start.contigPart.isReal) &&
        (join.end.contigPart.isReal);
    // dfmt on
}

/// Returns true iff join is an extension edge of the scaffold graph.
bool isExtension(J)(in J join) pure nothrow
{
    return isFrontExtension(join) ^ isBackExtension(join);
}

/// Returns true iff join is a front extension edge of the scaffold graph.
bool isFrontExtension(J)(in J join) pure nothrow
{
    // dfmt off
    return join.start.contigId == join.end.contigId &&
        join.start.contigPart == ContigPart.pre &&
        join.end.contigPart == ContigPart.begin;
    // dfmt on
}

/// Returns true iff join is a back extension edge of the scaffold graph.
bool isBackExtension(J)(in J join) pure nothrow
{
    // dfmt off
    return join.start.contigId == join.end.contigId &&
        join.start.contigPart == ContigPart.end &&
        join.end.contigPart == ContigPart.post;
    // dfmt on
}

/// Returns true iff join is a valid edge of the scaffold graph.
bool isValid(J)(in J join) pure nothrow
{
    return isDefault(join) ^ isGap(join) ^ isExtension(join);
}

/// Build a scaffold graph using `rawJoins`. This creates default edges and
/// inserts the rawJoins.
Scaffold!T buildScaffold(alias handleConflict, T)(in size_t numReferenceContigs, Join!T[] rawJoins)
{
    // dfmt off
    auto scaffold = initScaffold!T(numReferenceContigs)
        .addJoins!(handleConflict, T)(rawJoins);
    // dfmt on

    return scaffold;
}

private Scaffold!T initScaffold(T)(in size_t numReferenceContigs)
{
    // dfmt off
    auto contigIds = iota(1, numReferenceContigs + 1);
    auto contigNodes = contigIds
        .map!(contigId => only(
            ContigNode(contigId, ContigPart.pre),
            ContigNode(contigId, ContigPart.begin),
            ContigNode(contigId, ContigPart.end),
            ContigNode(contigId, ContigPart.post),
        ))
        .joiner
        .array;
    // dfmt on
    auto initialScaffold = Scaffold!T(contigNodes);

    foreach (contigId; contigIds)
    {
        initialScaffold ~= getDefaultJoin!T(contigId);
    }

    return initialScaffold;
}

Join!T getDefaultJoin(T)(size_t contigId) pure nothrow
{
    // dfmt off
    return Join!T(
        ContigNode(contigId, ContigPart.begin),
        ContigNode(contigId, ContigPart.end),
    );
    // dfmt on
}

private Scaffold!T addJoins(alias handleConflict, T)(Scaffold!T scaffold, Join!T[] rawJoins)
{
    foreach (join; rawJoins)
    {
        assert(join.isValid && !join.isDefault);

        scaffold.add!handleConflict(join);
    }

    return scaffold;
}

/// This marks ambiguous gap insertion for removal.
Scaffold!T discardAmbiguousJoins(T)(Scaffold!T scaffold)
{
    foreach (contigNode; scaffold.nodes)
    {
        assert(!contigNode.contigPart.isTranscendent || scaffold.degree(contigNode) <= 1);

        if (contigNode.contigPart.isReal && scaffold.degree(contigNode) > 2)
        {
            auto incidentGapJoins = scaffold.incidentEdges(contigNode).filter!isGap.array;

            if (incidentGapJoins.length > 1)
            {
                // dfmt off
                logJsonDiagnostic(
                    "info", "skipping ambiguous gap joins",
                    "sourceContigNode", contigNode.toJson,
                    "joins", incidentGapJoins.toJson,
                );
                // dfmt on

                foreach (join; incidentGapJoins)
                {
                    join.payload = T.init;
                    scaffold.add!(scaffold.ConflictStrategy.replace)(join);
                }
            }
        }
    }

    return removeNoneJoins!T(scaffold);
}

/// Remove marked edges from the graph. This always keeps the default edges.
Scaffold!T removeNoneJoins(T)(Scaffold!T scaffold)
{
    // dfmt off
    return Scaffold!T(scaffold.nodes.dup, scaffold
        .edges
        .filter!(join => isDefault(join) || join.payload != T.init)
        .map!(e => cast(Join!T) e)
        .array
    );
    // dfmt on
}

/// Remove extension edges were they coincide with a gap edge combining their
/// payload. This is intended to build pile ups with all reads that contribute
/// to each gap.
Scaffold!T mergeExtensionsWithGaps(alias mergePayloads, T)(Scaffold!T scaffold)
{
    foreach (contigNode; scaffold.nodes)
    {
        assert(!contigNode.contigPart.isTranscendent || scaffold.degree(contigNode) <= 1);
        assert(scaffold.degree(contigNode) <= 3);

        if (contigNode.contigPart.isReal && scaffold.degree(contigNode) == 3)
        {
            auto incidentJoins = scaffold.incidentEdges(contigNode).filter!(j => !isDefault(j)).array;
            assert(incidentJoins.length == 2);
            // The gap join has real `contigPart`s on both ends.
            int gapJoinIdx = incidentJoins[0].target(contigNode).contigPart.isReal ? 0 : 1;
            auto gapJoin = incidentJoins[gapJoinIdx];
            auto extensionJoin = incidentJoins[$ - gapJoinIdx - 1];

            gapJoin.payload = binaryFun!mergePayloads(gapJoin.payload, extensionJoin.payload);
            extensionJoin.payload = T.init;

            scaffold.add!(scaffold.ConflictStrategy.replace)(gapJoin);
            scaffold.add!(scaffold.ConflictStrategy.replace)(extensionJoin);
        }
    }

    return removeNoneJoins!T(scaffold);
}

///
unittest
{
    alias J = Join!int;
    alias S = Scaffold!int;
    alias CN = ContigNode;
    alias CP = ContigPart;
    //   contig 1      contig 2
    //
    //  o        o     o        o
    //          / e1 e2 \      / e4
    //    o -- o ------- o -- o
    //              e3
    //
    //    o -- o         o -- o
    //          \ e7 e8 /      \ e9
    //  o        o     o        o
    //
    //   contig 4      contig 3
    //
    auto scaffold = buildScaffold!(sumPayloads!int, int)(5, [
        J(CN(1, CP.end), CN(1, CP.post ), 1), // e1
        J(CN(1, CP.end), CN(1, CP.post ), 1), // e1
        J(CN(2, CP.pre), CN(2, CP.begin), 1), // e2
        J(CN(1, CP.end), CN(2, CP.begin), 1), // e3
        J(CN(2, CP.end), CN(2, CP.post ), 1), // e4
        J(CN(4, CP.end), CN(4, CP.post ), 1), // e7
        J(CN(3, CP.pre), CN(3, CP.begin), 1), // e8
        J(CN(3, CP.end), CN(3, CP.post ), 1), // e9
    ]).mergeExtensionsWithGaps!("a + b", int);

    assert(J(CN(1, CP.end), CN(1, CP.post )) !in scaffold); // e1
    assert(J(CN(2, CP.pre), CN(2, CP.begin)) !in scaffold); // e2
    assert(J(CN(1, CP.end), CN(2, CP.begin))  in scaffold); // e3
    assert(J(CN(2, CP.end), CN(2, CP.post ))  in scaffold); // e4
    assert(J(CN(4, CP.end), CN(4, CP.post ))  in scaffold); // e7
    assert(J(CN(3, CP.pre), CN(3, CP.begin))  in scaffold); // e8
    assert(J(CN(3, CP.end), CN(3, CP.post ))  in scaffold); // e9

    assert(scaffold.get(J(CN(1, CP.end), CN(2, CP.begin))).payload == 4); // merged 2 * e1 + e2 + e3
}

/**
    Performs a linear walk through a scaffold graph starting in startNode.
    A linear walk is a sequence of adjacent joins where no node is visited
    twice unless the graph is cyclic in which case the first node will appear
    twice. The implementation requires the graph to have linear components,
    ie. for every node the degree must be at most two. If the component of
    `startNode` is cyclic then the walk will in `startNode` and the `isCyclic`
    flag will be set.

    The direction of the walk can be influenced by giving `firstJoin`.

    **Note:** if one wants to read the `isCyclic` flag it is required to use
    `std.range.refRange` in most cases.

    Returns: range of joins in the scaffold graph.
    Throws: MissingNodeException if any node is encountered that is not part
            of the graph.
*/
LinearWalk!T linearWalk(T)(Scaffold!T scaffold, ContigNode startNode)
{
    return LinearWalk!T(scaffold, startNode);
}

/// ditto
LinearWalk!T linearWalk(T)(Scaffold!T scaffold, ContigNode startNode, Join!T firstJoin)
{
    return LinearWalk!T(scaffold, startNode, firstJoin);
}

///
unittest
{
    alias Payload = int;
    alias J = Join!Payload;
    alias S = Scaffold!Payload;
    alias CN = ContigNode;
    alias CP = ContigPart;
    //   contig 1      contig 2
    //
    //  o        o     o        o
    //                         / e4
    //    o -- o ------- o -- o
    //              e3
    //
    //    o -- o         o -- o
    //          \ e7 e8 /      \ e9
    //  o        o     o        o
    //
    //   contig 4      contig 3
    //
    auto joins1 = [
        J(CN(1, CP.end), CN(2, CP.begin), 1), // e3
        J(CN(2, CP.end), CN(2, CP.post ), 1), // e4
        J(CN(4, CP.end), CN(4, CP.post ), 1), // e7
        J(CN(3, CP.pre), CN(3, CP.begin), 1), // e8
        J(CN(3, CP.end), CN(3, CP.post ), 1), // e9
    ];
    auto scaffold1 = buildScaffold!(sumPayloads!Payload, Payload)(5, joins1);
    auto walks1 = [
        [
            getDefaultJoin!Payload(1),
            joins1[0],
            getDefaultJoin!Payload(2),
            joins1[1],
        ],
        [
            getDefaultJoin!Payload(4),
            joins1[2],
        ],
        [
            joins1[3],
            getDefaultJoin!Payload(3),
            joins1[4],
        ],
    ];

    alias getWalkStart = (walk) => walk[0].source(walk[0].getConnectingNode(walk[1]));

    foreach (walk; walks1)
    {
        auto reverseWalk = walk.retro.array;
        auto computedWalk = linearWalk!Payload(scaffold1, getWalkStart(walk));
        auto computedReverseWalk = linearWalk!Payload(scaffold1, getWalkStart(reverseWalk));

        assert(equal(walk[], refRange(&computedWalk)));
        assert(!computedWalk.isCyclic);
        assert(equal(reverseWalk[], refRange(&computedReverseWalk)));
        assert(!computedReverseWalk.isCyclic);
    }

    //   contig 1      contig 2
    //
    //  o        o     o        o
    //
    //              e1
    //    o -- o ------- o -- o
    //     \_________________/
    //              e2
    //
    auto joins2 = [
        J(CN(1, CP.end), CN(2, CP.begin), 1), // e1
        J(CN(2, CP.end), CN(1, CP.begin ), 1), // e2
    ];
    auto scaffold2 = buildScaffold!(sumPayloads!Payload, Payload)(5, joins2);
    auto walk2 = [
        getDefaultJoin!Payload(1),
        joins2[0],
        getDefaultJoin!Payload(2),
        joins2[1],
    ];

    {
        auto computedWalk = linearWalk!Payload(scaffold2, getWalkStart(walk2), walk2[0]);
        auto reverseWalk2 = walk2.retro.array;
        auto computedReverseWalk2 = linearWalk!Payload(scaffold2, getWalkStart(reverseWalk2), reverseWalk2[0]);

        assert(equal(walk2[], refRange(&computedWalk)));
        assert(computedWalk.isCyclic);
        assert(equal(reverseWalk2[], refRange(&computedReverseWalk2)));
        assert(computedWalk.isCyclic);
    }
}

struct LinearWalk(T)
{
    private Scaffold!T scaffold;
    private size_t currentNodeIdx;
    private Join!T currentJoin;
    private bool isEmpty = false;
    private Flag!"isCyclic" isCyclic = No.isCyclic;
    private NaturalNumberSet visitedNodes;

    /// Start linear walk through a scaffold graph in startNode.
    this(Scaffold!T scaffold, ContigNode startNode)
    {
        this.scaffold = scaffold;
        this.currentNode = startNode;
        this.visitedNodes.reserveFor(this.scaffold.nodes.length - 1);
        this.markVisited(this.currentNodeIdx);
        this.popFront();
    }

    /// Start linear walk through a scaffold graph in startNode.
    this(Scaffold!T scaffold, ContigNode startNode, Join!T firstJoin)
    {
        this.scaffold = scaffold;
        this.currentNode = startNode;
        this.visitedNodes.reserveFor(this.scaffold.nodes.length - 1);
        this.markVisited(this.currentNodeIdx);
        this.currentJoin = firstJoin;
        currentNode = currentJoin.target(currentNode);
        this.markVisited(this.currentNodeIdx);
    }

    void popFront()
    {
        assert(scaffold.degree(currentNode) <= 2, "fork in linear walk");

        // If isCyclic is set then the last edge of the cycle is being popped.
        // Thus we stop.
        if (isCyclic)
        {
            return endOfWalk();
        }

        auto candidateEdges = scaffold
            .incidentEdges(currentNode)
            .filter!(join => !visitedNodes.has(scaffold.indexOf(join.target(currentNode))));
        bool noSuccessorNodes = candidateEdges.empty;

        if (noSuccessorNodes)
        {
            // Iff the current node has more than one neighbor the graph has
            // a cycle.
            if (scaffold.degree(currentNode) > 1)
            {
                assert(scaffold.degree(currentNode) == 2);

                return lastEdgeOfCycle();
            }
            else
            {
                return endOfWalk();
            }
        }

        currentJoin = candidateEdges.front;
        currentNode = currentJoin.target(currentNode);
        markVisited(currentNodeIdx);
    }

    @property Join!T front()
    {
        return currentJoin;
    }

    @property bool empty()
    {
        return isEmpty;
    }

    private void lastEdgeOfCycle()
    {
        isCyclic = Yes.isCyclic;
        // Find the missing edge.
        currentJoin = scaffold
            .incidentEdges(currentNode)
            .filter!(join => join != currentJoin)
            .front;
    }

    private void endOfWalk()
    {
        currentNodeIdx = scaffold.nodes.length;
        isEmpty = true;
    }

    private @property ContigNode currentNode()
    {
        return scaffold.nodes[currentNodeIdx];
    }

    private @property void currentNode(ContigNode node)
    {
        this.currentNodeIdx = this.scaffold.indexOf(node);
    }

    private void markVisited(size_t nodeIdx)
    {
        this.visitedNodes.add(nodeIdx);
    }
}
