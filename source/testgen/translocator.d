/**
    This package is used to translocate assembly gaps from assembly to another.

    Copyright: © 2018 Arne Ludwig <arne.ludwig@posteo.de>
    License: Subject to the terms of the MIT license, as written in the
             included LICENSE file.
    Authors: Arne Ludwig <arne.ludwig@posteo.de>
*/
module testgen.translocator;

import djunctor.alignments : AlignmentChain;
import djunctor.dazzler : getAlignments, getNumContigs, getFastaEntries;
import djunctor.djunctor : getRegion;
import djunctor.util.fasta : parseFastaRecord;
import djunctor.util.log;
import djunctor.util.region : Region;
import std.algorithm : cache, copy, filter, joiner, map;
import std.array : array;
import std.conv : to;
import std.format : format;
import std.range : assumeSorted, chain, chunks, enumerate, only, repeat, slide, takeExactly;
import std.stdio : stdout;
import testgen.commandline : Options;
import vibe.data.json : toJson = serializeToJson;

/// Start the translocator algorithm with preprocessed options.
void runWithOptions(in ref Options options)
{
    Translocator(options).run();
}

alias AssemblyMask = Region!(size_t, size_t, "contigId");
alias AssemblyInterval = AssemblyMask.TaggedInterval;
alias AssemblyPoint = AssemblyMask.TaggedPoint;

private struct Translocator
{
    const(Options) options;
    size_t numContigsAssembly2;
    AlignmentChain[] alignments;
    AssemblyMask mappedRegions;

    this(in Options options)
    {
        this.options = options;
    }

    void run()
    {
        logJsonDiagnostic("state", "enter", "function", "run");
        init();
        writeOutputAssembly();
        logJsonDiagnostic("state", "exit", "function", "run");
    }

    protected void init()
    {
        logJsonDiagnostic("state", "enter", "function", "init");
        // dfmt off
        alignments = getAlignments(
            options.assembly2Db,
            options.assembly1Db,
            options.alignmentFile,
            options
        );
        mappedRegions = AssemblyMask(alignments
            .filter!"a.isProper"
            .map!(getRegion!AssemblyMask)
            .map!"a.intervals.dup"
            .joiner
            .array);
        // dfmt on

        logJsonDebug("mappedRegions", mappedRegions.intervals.toJson);
        logJsonDiagnostic("state", "exit", "function", "init");
    }

    protected void writeOutputAssembly()
    {
        logJsonDiagnostic("state", "enter", "function", "writeOutputAssembly");
        immutable dchar unknownBase = 'n';
        auto assembly2Contigs = getFastaEntries(options.assembly2Db, cast(size_t[]) [], options).map!parseFastaRecord;
        auto mappedRegions = mappedRegions.intervals.assumeSorted!"a.contigId < b.contigId";
        AssemblyInterval needle;

        // dfmt off
        assembly2Contigs
            .enumerate(1)
            .map!(
                (enumContigs) => getScaffoldHeader(enumContigs[0]),
                (enumContigs) => {
                    needle.contigId = enumContigs[0];
                    auto contigSequence = enumContigs[1][].array;
                    auto contigMappedRegions = chain(
                        only(needle),
                        mappedRegions.equalRange(needle)
                    );

                    return contigMappedRegions
                        .slide(2)
                        .map!((keepRegions) => {
                            auto numNs = keepRegions[0].end == 0
                                ? 0
                                : keepRegions[1].begin - keepRegions[0].end;
                            auto sequencePiece = chain(
                                repeat(unknownBase).takeExactly(numNs),
                                contigSequence[keepRegions[1].begin .. keepRegions[1].end],
                            );
                            debug logJsonDebug(
                                "keepRegions", keepRegions.array.toJson,
                                "numNs", numNs,
                                "numBases", keepRegions[1].end - keepRegions[1].begin,
                            );

                            return sequencePiece;
                        })
                        .map!"a()"
                        .cache
                        .joiner
                        .chunks(options.fastaLineWidth)
                        .joiner("\n");
                }
            )
            .map!(scaffoldParts => chain(scaffoldParts[0], "\n", scaffoldParts[1]()))
            .joiner("\n")
            .copy(stdout.lockingTextWriter);
        // dfmt on
        logJsonDiagnostic("state", "exit", "function", "writeOutputAssembly");
    }

    static protected string getScaffoldHeader(in size_t scaffoldId)
    {
        return format!">translocated_gaps_%d"(scaffoldId);
    }
}

