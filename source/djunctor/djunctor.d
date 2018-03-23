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
import std.stdio : writeln;

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

        static auto getOrdered(string order, T)(T thingA, T thingB) pure nothrow
                if (order == "a2b" || order == "b2a")
        {
            static if (order == "a2b")
                return tuple(thingA, thingB);
            else
                return tuple(thingB, thingA);
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
        yes,
        no,
    }

    Contig contigA;
    Contig contigB;
    Complement complement;
    LocalAlignment[] localAlignments;

    invariant
    {
        assert(localAlignments.length >= 1);
    }

    int opCmp(ref const AlignmentChain other) const pure nothrow
    {
        import std.math : sgn;

        size_t idCompare = other.contigA.id - this.contigA.id;
        if (idCompare != 0)
            return cast(int) sgn(idCompare);

        idCompare = other.contigB.id - this.contigB.id;
        if (idCompare != 0)
            return cast(int) sgn(idCompare);

        auto thisFirstLA = this.localAlignments[0];
        auto otherFirstLA = other.localAlignments[0];
        size_t locusCompare = otherFirstLA.contigA.begin - thisFirstLA.contigA.begin;
        if (locusCompare != 0)
            return cast(int) sgn(locusCompare);

        locusCompare = otherFirstLA.contigB.begin - thisFirstLA.contigB.begin;
        if (locusCompare != 0)
            return cast(int) sgn(locusCompare);

        auto thisLastLA = this.localAlignments[$ - 1];
        auto otherLastLA = other.localAlignments[$ - 1];
        locusCompare = otherLastLA.contigA.end - thisLastLA.contigA.end;
        if (locusCompare != 0)
            return cast(int) sgn(locusCompare);

        locusCompare = otherLastLA.contigB.end - thisLastLA.contigB.end;
        if (locusCompare != 0)
            return cast(int) sgn(locusCompare);

        return 0;
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
