/**
    Defines bindings and utilities to/for the dazzler commands.

    Copyright: © 2018 Arne Ludwig <arne.ludwig@posteo.de>
    License: Subject to the terms of the MIT license, as written in the
             included LICENSE file.
    Authors: Arne Ludwig <arne.ludwig@posteo.de>
*/
module djunctor.dazzler;

import djunctor.commandline : hasOption, isOptionsList;
import djunctor.alignments : AlignmentChain;
import djunctor.util.fasta : parseFastaRecord;
import djunctor.util.log;
import djunctor.util.range : arrayChunks, takeExactly;
import djunctor.util.tempfile : mkstemp;
import std.algorithm : all, cache, canFind, endsWith, equal, filter, isSorted, joiner,
    map, min, sort, splitter, startsWith, SwapStrategy, uniq;
import std.array : appender, Appender, array, uninitializedArray;
import std.conv : to;
import std.exception : enforce;
import std.file : exists, remove;
import std.format : format, formattedRead;
import std.meta : Instantiate;
import std.path : absolutePath, baseName, buildPath, dirName, relativePath,
    stripExtension, withExtension;
import std.process : Config, escapeShellCommand, kill, pipeProcess,
    ProcessPipes, Redirect, wait;
import std.range : chain, chunks, drop, only, take, slide;
import std.range.primitives : ElementType, empty, isForwardRange, isInputRange;
import std.stdio : File, writeln;
import std.string : lineSplitter, outdent;
import std.traits : hasMember, isIntegral, isSomeChar, isSomeString, Unqual;
import std.typecons : Flag, No, tuple, Tuple, Yes;
import std.variant : Algebraic;
import vibe.data.json : Json, toJson = serializeToJson;

/// File suffixes of hidden DB files.
private immutable hiddenDbFileSuffixes = [".bps", ".hdr", ".idx"];

/// Constant holding the .dam file extension.
immutable damFileExtension = ".dam";

/// The Dazzler tools require sequence of a least minSequenceLength base pairs.
immutable minSequenceLength = 14;

/**
    Return a list of hidden files associated to every `.dam`/`.db` file. These
    files contain the actual data used in all the computation. Thus, we
    carefully check for their existence.
*/
auto getHiddenDbFiles(string dbFile)
{
    import std.algorithm : map;

    return hiddenDbFileSuffixes.map!(suffix => buildPath(dbFile.dirName,
            "." ~ dbFile.baseName.withExtension(suffix).to!string));
}

class DazzlerCommandException : Exception
{
    pure nothrow @nogc @safe this(string msg, string file = __FILE__,
            size_t line = __LINE__, Throwable nextInChain = null)
    {
        super(msg, file, line, nextInChain);
    }
}

enum ProvideMethod
{
    copy,
    symlink,
}

/**
    Provide dbFile in `workdir`.

    Returns: Path of the dbFile in `workdir`.
*/
string provideDamFileInWorkdir(in string dbFile, ProvideMethod provideMethod, in string workdir)
{
    import std.file : copy, symlink;
    import std.range : chain, only;

    alias inWorkdir = anyDbFile => buildPath(workdir, anyDbFile.baseName);
    auto allDbFiles = chain(only(dbFile), getHiddenDbFiles(dbFile));

    foreach (anyDbFile; allDbFiles)
    {
        switch (provideMethod)
        {
        case ProvideMethod.copy:
            copy(anyDbFile, inWorkdir(anyDbFile));
            break;
        case ProvideMethod.symlink:
            symlink(anyDbFile, inWorkdir(anyDbFile));
            break;
        default:
            assert(0);
        }
    }

    return inWorkdir(dbFile);
}

/// Build a new .dam file by using the given subset of reads in inDbFile.
string dbSubset(Options, R)(in string inDbFile, R readIds, in Options options)
        if (hasOption!(Options, "workdir", isSomeString)
            && hasOption!(Options, "dbsplitOptions", isOptionsList))
{
    immutable outDbNameTemplate = "subset-XXXXXX";

    auto outDbTemplate = buildPath(options.workdir, outDbNameTemplate);
    auto outDb = mkstemp(outDbTemplate, damFileExtension);

    outDb.file.close();
    remove(outDb.name);
    buildSubsetDb(inDbFile, outDb.name, readIds, options.workdir);
    dbsplit(outDb.name, options.dbsplitOptions, options.workdir);

    return outDb.name;
}

AlignmentChain[] getLocalAlignments(Options)(in string dbA, in Options options)
        if (hasOption!(Options, "dalignerOptions", isOptionsList) && hasOption!(Options,
            "ladumpOptions", isOptionsList) && hasOption!(Options, "workdir", isSomeString))
{
    if (!lasFilesGenerated(dbA, options.workdir))
    {
        dalign(dbA, options.dalignerOptions, options.workdir);
    }

    return processGeneratedLasFiles(dbA, null, options);
}

AlignmentChain[] getLocalAlignments(Options)(in string dbA, in string dbB, in Options options)
        if (hasOption!(Options, "dalignerOptions", isOptionsList) && hasOption!(Options,
            "ladumpOptions", isOptionsList) && hasOption!(Options, "workdir", isSomeString))
{
    if (!lasFilesGenerated(dbA, dbB, options.workdir))
    {
        dalign(dbA, dbB, options.dalignerOptions, options.workdir);
    }

    return processGeneratedLasFiles(dbA, dbB, options);
}

void computeLocalAlignments(Options)(in string[] dbList, in Options options)
        if (hasOption!(Options, "dalignerOptions", isOptionsList)
            && hasOption!(Options, "workdir", isSomeString))
{
    dalign(dbList, options.dalignerOptions, options.workdir);
}

AlignmentChain[] getMappings(Options)(in string dbA, in string dbB, in Options options)
        if (hasOption!(Options, "damapperOptions", isOptionsList) && hasOption!(Options,
            "ladumpOptions", isOptionsList) && hasOption!(Options, "workdir", isSomeString))
{
    if (!lasFilesGenerated(dbA, dbB, options.workdir))
    {
        damapper(dbA, dbB, options.damapperOptions, options.workdir);
    }

    return processGeneratedLasFiles(dbA, dbB, options);
}

void computeMappings(Options)(in string[] dbList, in Options options)
        if (hasOption!(Options, "damapperOptions", isOptionsList)
            && hasOption!(Options, "workdir", isSomeString))
{
    damapper(dbList, options.damapperOptions, options.workdir);
}

private AlignmentChain[] processGeneratedLasFiles(Options)(in string dbA,
        in string dbB, in Options options)
        if (hasOption!(Options, "ladumpOptions", isOptionsList)
            && hasOption!(Options, "workdir", isSomeString))
{
    auto lasFileList = getLasFiles(dbA, dbB, options.workdir);
    auto dbFiles = tuple(dbA, dbB);

    // dfmt off
    auto alignmentChains =
        lasFileList
            .map!(lasFile => ladump(lasFile, dbFiles.expand,
                    options.ladumpOptions, options.workdir))
            .map!(lasDump => readAlignmentList(lasDump))
            .joiner
            .array;
    // dfmt on
    alignmentChains.sort!("a < b", SwapStrategy.stable);

    return alignmentChains;
}

private enum ChainPartType
{
    start = '>',
    continuation = '-',
    alternateStart = '+',
    noChainInFile = '.',
}

private auto readAlignmentList(S)(in S lasDump) if (isSomeString!S)
{
    import std.algorithm : chunkBy, count, filter;
    import std.range : chunks, drop, enumerate;

    immutable recordSeparator = ';';
    immutable chainPartFormat = "P %d %d %c %c;L %d %d;C %d %d %d %d;D %d";
    immutable numChainPartLines = chainPartFormat.count(recordSeparator) + 1;
    // dfmt off
    alias RawChainPart = Tuple!(
        AlignmentChain.Contig, "contigA",
        AlignmentChain.Contig, "contigB",
        AlignmentChain.Complement, "complement",
        ChainPartType, "chainPartType",
        AlignmentChain.LocalAlignment, "localAlignment",
    );
    // dfmt on

    /// Build chunks of numChainPartLines lines.
    alias byChainPartSplitter = lasDump => lasDump.lineSplitter.drop(2).chunks(numChainPartLines);
    /// Parse chunks of numChainPartLines lines into RawChainPart s.
    alias parseChainPart = chainPartLines => {
        AlignmentChain.Contig contigA, contigB;
        char complement;
        char chainPartType;
        AlignmentChain.LocalAlignment localAlignment;

        auto joinedLines = chainPartLines.joiner(only(recordSeparator));
        // dfmt off
        int numMatches = joinedLines
            .array
            .formattedRead!chainPartFormat(
                contigA.id,
                contigB.id,
                complement,
                chainPartType,
                contigA.length,
                contigB.length,
                localAlignment.contigA.begin,
                localAlignment.contigA.end,
                localAlignment.contigB.begin,
                localAlignment.contigB.end,
                localAlignment.numDiffs,
            );
        // dfmt on
        assert(numMatches == 11, format!"%d matches in chunk: `%s`"(numMatches, joinedLines.array));
        AlignmentChain.Complement flagComplement = complement == 'c'
            ? AlignmentChain.Complement.yes : AlignmentChain.Complement.no;

        return RawChainPart(contigA, contigB, flagComplement,
                chainPartType.to!ChainPartType, localAlignment);
    };
    // dfmt off
    /// Returns true iff both parts have the same contigs.
    alias sameContigsInvolved = (part1, part2) =>
        part1.contigA.id == part2.contigA.id &&
        part1.contigB.id == part2.contigB.id;
    /// Returns true iff `part1` and `part2` belong to the same chain.
    alias isChainContinuation = (part1, part2) =>
        (
            part1.chainPartType == ChainPartType.start &&
            part2.chainPartType == ChainPartType.continuation
        )
        ||
        (
            part1.chainPartType == ChainPartType.alternateStart &&
            part2.chainPartType == ChainPartType.continuation
        )
        ||
        (
            part1.chainPartType == ChainPartType.continuation &&
            part2.chainPartType == ChainPartType.continuation
        );
    /// Make sure the complement flag is the same for all parts of a chain.
    alias assertSameComplementValue = (chainPart1, chainPart2) => {
        assert(chainPart1.complement == chainPart2.complement);

        return true;
    };
    alias belongToSameChain = (chainPart1, chainPart2) pure =>
        chainPart1() == chainPart2() || (
            sameContigsInvolved(chainPart1(), chainPart2()) &&
            isChainContinuation(chainPart1(), chainPart2()) &&
            assertSameComplementValue(chainPart1(), chainPart2()));
    /// Convert the chunk of chain parts into an AlignmentChain.
    alias buildAlignmentChain = chainPartChunk => AlignmentChain(
        chainPartChunk[0],
        chainPartChunk[1].front()().contigA,
        chainPartChunk[1].front()().contigB,
        chainPartChunk[1].front()().complement,
        chainPartChunk[1]
            .map!(chainPart => chainPart().localAlignment)
            .array);

    return byChainPartSplitter(lasDump)
        .map!parseChainPart
        .chunkBy!belongToSameChain
        .enumerate
        .map!buildAlignmentChain;
    // dfmt on
}

unittest
{
    immutable testLasDump = q"EOF
        + P 9
        % P 1337
        P 1 2 n >
        L 8 9
        C 3 4 5 6
        D 7
        P 1 2 n -
        L 17 18
        C 12 13 14 15
        D 16
        P 19 20 c +
        L 26 27
        C 21 22 23 24
        D 25
        P 19 20 c -
        L 35 36
        C 30 31 32 33
        D 34
        P 37 38 n .
        L 35 36
        C 39 40 41 42
        D 43
        P 46 47 c .
        L 53 54
        C 48 49 50 51
        D 52
        P 46 47 n .
        L 53 54
        C 57 58 59 60
        D 61
        P 64 65 c >
        L 71 72
        C 66 67 68 69
        D 70
        P 55 56 c -
        L 80 81
        C 75 76 77 78
        D 79
EOF".outdent;

    auto alignmentChains = readAlignmentList(testLasDump).array;
    with (AlignmentChain) with (Complement) with (LocalAlignment)
            {
                // dfmt off
                assert(alignmentChains ==
                    [
                        AlignmentChain(
                            0,
                            Contig(1, 8),
                            Contig(2, 9),
                            no,
                            [
                                LocalAlignment(
                                    Locus(3, 4),
                                    Locus(5, 6),
                                    7
                                ),
                                LocalAlignment(
                                    Locus(12, 13),
                                    Locus(14, 15),
                                    16
                                )
                            ]
                        ),
                        AlignmentChain(
                            1,
                            Contig(19, 26),
                            Contig(20, 27),
                            yes,
                            [
                                LocalAlignment(
                                    Locus(21, 22),
                                    Locus(23, 24),
                                    25
                                ),
                                LocalAlignment(
                                    Locus(30, 31),
                                    Locus(32, 33),
                                    34
                                )
                            ]
                        ),
                        AlignmentChain(
                            2,
                            Contig(37, 35),
                            Contig(38, 36),
                            no,
                            [
                                LocalAlignment(
                                    Locus(39, 40),
                                    Locus(41, 42),
                                    43
                                )
                            ]
                        ),
                        AlignmentChain(
                            3,
                            Contig(46, 53),
                            Contig(47, 54),
                            yes,
                            [
                                LocalAlignment(
                                    Locus(48, 49),
                                    Locus(50, 51),
                                    52
                                ),
                            ]
                        ),
                        AlignmentChain(
                            4,
                            Contig(46, 53),
                            Contig(47, 54),
                            no,
                            [
                                LocalAlignment(
                                    Locus(57, 58),
                                    Locus(59, 60),
                                    61
                                )
                            ]
                        ),
                        AlignmentChain(
                            5,
                            Contig(64, 71),
                            Contig(65, 72),
                            yes,
                            [
                                LocalAlignment(
                                    Locus(66, 67),
                                    Locus(68, 69),
                                    70
                                )
                            ]
                        ),
                        AlignmentChain(
                            6,
                            Contig(55, 80),
                            Contig(56, 81),
                            yes,
                            [
                                LocalAlignment(
                                    Locus(75, 76),
                                    Locus(77, 78),
                                    79
                                )
                            ]
                        )
                    ],
                    format!"error parsing test input; result was: %s"(alignmentChains)
                );
                // dfmt on
            }
}

void attachTracePoints(Options)(AlignmentChain*[] alignmentChains, in string dbA,
        in string dbB, in string[] dazzlerOptions, in Options options)
        if (hasOption!(Options, "workdir", isSomeString)
            && hasOption!(Options, "ladumpTraceOptions", isOptionsList))
{
    auto lasFileList = getLasFiles(dbA, dbB, options.workdir);
    // NOTE: dump only for matching A reads; better would be for matching
    //       B reads but that is not possible with `LAdump`.
    auto aReadIds = alignmentChains.map!"a.contigA.id".uniq.array;

    // dfmt off
    auto tracePointDumps =
        lasFileList
            .map!(lasFile => ladump(lasFile, dbA, dbB, aReadIds,
                    options.ladumpTraceOptions, options.workdir))
            .map!(lasDump => readTracePointList(lasDump))
            .joiner
            .array;
    // dfmt on
    tracePointDumps.sort!("a < b", SwapStrategy.stable);
    assert(isSorted!"*a < *b"(alignmentChains), "alignmentChains must be sorted");

    auto numAttached = alignmentChains.attachTracePointDumps(tracePointDumps,
            getTracePointDistance(dazzlerOptions));
    assert(numAttached == alignmentChains.length,
            "missing trace point lists for some alignment chains");

    foreach (idx, alignmentChain; alignmentChains)
    {
        assert(alignmentChain.tracePointDistance > 0);

        foreach (localAlignment; alignmentChain.localAlignments)
        {
            assert(localAlignment.tracePoints.length > 0);
        }
    }
}

// dfmt off
private alias TracePointDump = Tuple!(
    size_t, "contigAId",
    size_t, "contigBId",
    size_t, "contigABegin",
    size_t, "contigBBegin",
    size_t, "contigAEnd",
    size_t, "contigBEnd",
    AlignmentChain.LocalAlignment.TracePoint[][], "tracePointLists",
);
// dfmt on

private TracePointDump[] readTracePointList(S)(in S lasDump) if (isSomeString!S)
{
    alias TracePoint = AlignmentChain.LocalAlignment.TracePoint;

    debug auto nextLineTypes = "+%@";
    size_t contigAId;
    size_t contigBId;
    size_t contigAFirstBegin;
    size_t contigBFirstBegin;
    size_t contigALastEnd;
    size_t contigBLastEnd;
    bool collectingChainParts = false;
    bool chainContinuation = false;
    bool collectingTracePoints = false;
    Appender!(TracePoint[]) tracePointList;
    Appender!(TracePoint[][]) tracePointLists;
    auto tracePointDumps = appender!(TracePointDump[]);

    void startTracePointList(size_t numTracePoints)
    {
        tracePointList = appender!(TracePoint[]);
        tracePointList.reserve(numTracePoints);
        collectingTracePoints = true;
    }

    void finishTracePointList()
    {
        if (collectingTracePoints)
        {
            tracePointLists ~= tracePointList.data;
            collectingTracePoints = false;
        }
    }

    void startChain()
    {
        tracePointLists = appender!(TracePoint[][]);
        collectingChainParts = true;
    }

    void nextChainPart()
    {
        finishTracePointList();
        chainContinuation = true;
    }

    void finishChain()
    {
        if (collectingChainParts)
        {
            finishTracePointList();
            // dfmt off
            tracePointDumps ~= TracePointDump(
                contigAId,
                contigBId,
                contigAFirstBegin,
                contigBFirstBegin,
                contigALastEnd,
                contigBLastEnd,
                tracePointLists.data,
            );
            // dfmt on
            collectingChainParts = false;
        }
        chainContinuation = false;
    }

    foreach (dumpLine; lasDump.lineSplitter)
    {
        if (dumpLine.length == 0)
            continue;

        debug assert(nextLineTypes.canFind(dumpLine[0 .. 1]));
        switch (dumpLine[0])
        {
        case '+':
            if (dumpLine[2] == 'T')
            {
                size_t numTotalTracePoints;

                auto numMatches = dumpLine.formattedRead!"+ T %d"(numTotalTracePoints);
                assert(numMatches == 1);

                tracePointDumps.reserve(numTotalTracePoints);
                debug nextLineTypes = "P+%@";
            }
            break;
        case 'P':
            size_t nextContigAId;
            size_t nextContigBId;
            char complement;
            char chainPartType;

            // dfmt off
            auto numMatches = dumpLine.formattedRead!"P %d %d %c %c"(
                nextContigAId,
                nextContigBId,
                complement,
                chainPartType,
            );
            // dfmt on
            assert(numMatches == 4);

            final switch (chainPartType.to!ChainPartType)
            {
            case ChainPartType.start:
            case ChainPartType.alternateStart:
            case ChainPartType.noChainInFile:
                finishChain();
                contigAId = nextContigAId;
                contigBId = nextContigBId;
                startChain();
                break;
            case ChainPartType.continuation:
                if (contigAId == nextContigAId
                        && contigBId == nextContigBId)
                {
                    nextChainPart();
                }
                else
                {
                    debug logJsonDebug("info", "missing chain start");
                    // NOTE: This is a work-around for a bug in LAdump: if the
                    //       beginning of a chain is not a "proper overlap" but
                    //       a continuation is then LAdump will print the
                    //       continuation only. Thus, a part marked as
                    //       continuation still may be a chain start.
                    finishChain();
                    contigAId = nextContigAId;
                    contigBId = nextContigBId;
                    startChain();
                }
                break;
            }
            debug nextLineTypes = "C";
            break;
        case 'C':
            size_t contigABegin;
            size_t contigBBegin;
            size_t contigAEnd;
            size_t contigBEnd;

            // dfmt off
                auto numMatches = dumpLine.formattedRead!"C %d %d %d %d"(
                    contigABegin,
                    contigAEnd,
                    contigBBegin,
                    contigBEnd,
                );
                // dfmt on
            assert(numMatches == 4);

            if (chainContinuation)
            {
                contigALastEnd = contigAEnd;
                contigBLastEnd = contigBEnd;
            }
            else
            {
                contigAFirstBegin = contigABegin;
                contigBFirstBegin = contigBBegin;
                contigALastEnd = contigAEnd;
                contigBLastEnd = contigBEnd;
            }

            debug nextLineTypes = "T";
            break;
        case 'T':
            size_t numTracePoints;

            auto numMatches = dumpLine.formattedRead!"T %d"(numTracePoints);
            assert(numMatches == 1);

            startTracePointList(numTracePoints);
            debug nextLineTypes = " ";
            break;
        case ' ':
            if (collectingTracePoints)
            {
                TracePoint tracePoint;

                auto numMatches = dumpLine.formattedRead!"   %d %d"(tracePoint.numDiffs,
                        tracePoint.numBasePairs);
                assert(numMatches == 2);

                tracePointList ~= tracePoint;
            }
            debug nextLineTypes = " P";
            break;
        default:
            break; // ignore
        }
    }

    finishChain();

    return tracePointDumps.data;
}

unittest
{
    immutable testLasDump = q"EOF
        + P 3
        % P 42
        + T 9
        % T 42
        @ T 42
        P 1 1 n >
        C 0 1337 0 42
        T 3
           2 102
           3 101
           4 104
        P 1 1 c +
        C 0 1338 0 43
        T 3
           3 101
           4 104
           2 102
        P 1 2 n >
        C 0 1339 0 44
        T 4
           6 105
           1 101
           2 100
           3  97
        P 1 2 n -
        C 0 1340 0 45
        T 2
           0   2
           2 102
EOF".outdent;

    with (AlignmentChain.LocalAlignment)
    {
        auto tracePointDumps = readTracePointList(testLasDump).array;
        // dfmt off
        auto expectedResult = [
            TracePointDump(1, 1, 0, 0, 1337, 42, [[
                TracePoint(2, 102),
                TracePoint(3, 101),
                TracePoint(4, 104),
            ]]),
            TracePointDump(1, 1, 0, 0, 1338, 43, [[
                TracePoint(3, 101),
                TracePoint(4, 104),
                TracePoint(2, 102),
            ]]),
            TracePointDump(1, 2, 0, 0, 1340, 45, [
                [
                    TracePoint(6, 105),
                    TracePoint(1, 101),
                    TracePoint(2, 100),
                    TracePoint(3, 97),
                ],
                [
                    TracePoint(0, 2),
                    TracePoint(2, 102),
                ]
            ]),
        ];
        // dfmt on

        assert(tracePointDumps == expectedResult);
    }
}

auto fingerprint(AlignmentChain* alignmentChain) pure nothrow
{
    // dfmt off
    return tuple(
        alignmentChain.contigA.id,
        alignmentChain.contigB.id,
        alignmentChain.first.contigA.begin + 0,
        alignmentChain.first.contigB.begin + 0,
        alignmentChain.last.contigA.end + 0,
        alignmentChain.last.contigB.end + 0,
    );
    // dfmt on
}

auto fingerprint(TracePointDump* tracePointDump) pure nothrow
{
    // dfmt off
    return tuple(
        tracePointDump.contigAId,
        tracePointDump.contigBId,
        tracePointDump.contigABegin,
        tracePointDump.contigBBegin,
        tracePointDump.contigAEnd,
        tracePointDump.contigBEnd,
    );
    // dfmt on
}

private size_t attachTracePointDumps(AlignmentChain*[] alignmentChains,
        TracePointDump[] tracePointDumps, in size_t tracePointDistance) pure
{
    assert(tracePointDistance > 0);
    assert(isSorted!"*a < *b"(alignmentChains));

    size_t numAlignmentChainsAffected = 0;
    size_t numLoops = 0;
    size_t i = 0;
    size_t j = 0;

    while (i < alignmentChains.length && j < tracePointDumps.length
            && numLoops < alignmentChains.length + tracePointDumps.length)
    {
        auto alignmentChain = alignmentChains[i];
        auto tracePointDump = &tracePointDumps[j];
        auto acFingerprint = alignmentChain.fingerprint;
        auto tpdFingerprint = tracePointDump.fingerprint;

        // dfmt off
        debug logJsonDebug(
            "acFingerprint", acFingerprint.toJson,
            "tpdFingerprint", tpdFingerprint.toJson,
        );
        // dfmt on

        if (acFingerprint == tpdFingerprint)
        {
            foreach (k, ref localAlignment; alignmentChain.localAlignments)
            {
                localAlignment.tracePoints = tracePointDump.tracePointLists[k];
            }
            alignmentChain.tracePointDistance = tracePointDistance;
            numAlignmentChainsAffected = i + 1;
            ++i;
        }
        else if (tpdFingerprint < acFingerprint)
        {
            ++j;
        }
        else
        {
            assert(tpdFingerprint > acFingerprint);
            throw new Exception(
                    format!"missing trace point data for alignment chain: %s"(*alignmentChain));
        }

        ++numLoops;
    }

    return numAlignmentChainsAffected;
}

unittest
{
    with (AlignmentChain) with (LocalAlignment)
        {
            // dfmt off
        auto alignmentChains = [
            new AlignmentChain(
                1,
                Contig(1, 1337),
                Contig(1, 42),
                Complement.no,
                [
                    LocalAlignment(Locus(0, 1337), Locus(0, 42), 0),
                ],
            ),
            new AlignmentChain(
                2,
                Contig(1, 1338),
                Contig(1, 43),
                Complement.no,
                [
                    LocalAlignment(Locus(0, 1338), Locus(0, 43), 0),
                ],
            ),
            new AlignmentChain(
                3,
                Contig(1, 1340),
                Contig(3, 45),
                Complement.no,
                [
                    LocalAlignment(Locus(0, 1339), Locus(0, 44), 0),
                    LocalAlignment(Locus(0, 1340), Locus(0, 45), 0),
                ],
            ),
        ];
        auto tracePointDumps = [
            TracePointDump(1, 1, 0, 0, 1337, 42, [[
                TracePoint(1, 102),
                TracePoint(2, 101),
                TracePoint(3, 104),
            ]]),
            TracePointDump(1, 1, 0, 0, 1338, 43, [[
                TracePoint(4, 101),
                TracePoint(5, 104),
                TracePoint(6, 102),
            ]]),
            TracePointDump(1, 2, 0, 0, 1, 2, [[]]),
            TracePointDump(1, 2, 0, 0, 3, 4, [[]]),
            TracePointDump(1, 3, 0, 0, 1340, 45, [
                [
                    TracePoint(7, 105),
                    TracePoint(8, 101),
                    TracePoint(9, 100),
                    TracePoint(10, 97),
                ],
                [
                    TracePoint(11, 2),
                    TracePoint(12, 102),
                ],
            ]),
        ];
        auto expectedAlignmentChains = [
            AlignmentChain(
                1,
                Contig(1, 1337),
                Contig(1, 42),
                Complement.no,
                [
                    LocalAlignment(Locus(0, 1337), Locus(0, 42), 0, [
                        TracePoint(1, 102),
                        TracePoint(2, 101),
                        TracePoint(3, 104),
                    ]),
                ],
                101
            ),
            AlignmentChain(
                2,
                Contig(1, 1338),
                Contig(1, 43),
                Complement.no,
                [
                    LocalAlignment(Locus(0, 1338), Locus(0, 43), 0, [
                        TracePoint(4, 101),
                        TracePoint(5, 104),
                        TracePoint(6, 102),
                    ]),
                ],
                101
            ),
            AlignmentChain(
                3,
                Contig(1, 1340),
                Contig(3, 45),
                Complement.no,
                [
                    LocalAlignment(Locus(0, 1339), Locus(0, 44), 0, [
                        TracePoint(7, 105),
                        TracePoint(8, 101),
                        TracePoint(9, 100),
                        TracePoint(10, 97),
                    ]),
                    LocalAlignment(Locus(0, 1340), Locus(0, 45), 0, [
                        TracePoint(11, 2),
                        TracePoint(12, 102),
                    ]),
                ],
                101
            ),
        ];
        // dfmt on

            assert(alignmentChains.attachTracePointDumps(tracePointDumps, 101) == 3);
            assert(alignmentChains.map!"*a".array == expectedAlignmentChains);
        }
}

size_t getTracePointDistance(in string[] dazzlerOptions) pure
{
    immutable defaultTracePointDistance = 100;

    foreach (option; dazzlerOptions)
    {
        if (option.startsWith(cast(const(string)) DalignerOptions.tracePointDistance))
        {
            return option[2 .. $].to!size_t;
        }
    }

    return defaultTracePointDistance;
}

unittest
{
    with (DalignerOptions)
    {
        assert(getTracePointDistance([]) == 100);
        assert(getTracePointDistance([identity]) == 100);
        // dfmt off
        assert(getTracePointDistance([
            tracePointDistance ~ "42",
            averageCorrelationRate ~ ".8",
            identity,
        ]) == 42);
        assert(getTracePointDistance([
            identity,
            tracePointDistance ~ "42",
            averageCorrelationRate ~ ".8",
        ]) == 42);
        assert(getTracePointDistance([
            averageCorrelationRate ~ ".8",
            identity,
            tracePointDistance ~ "42",
        ]) == 42);
        // dfmt on
    }
}

/**
    Get the designated set of records in FASTA format. If recordNumbers is
    empty the whole DB will be converted.
*/
auto getFastaEntries(Options, Range)(in string dbFile, Range recordNumbers, in Options options)
        if (hasOption!(Options, "fastaLineWidth", isIntegral) && hasOption!(Options,
            "workdir", isSomeString) && isInputRange!Range && is(ElementType!Range : size_t))
{
    // dfmt off
    string[] dbdumpOptions = [
        DBdumpOptions.readNumber,
        DBdumpOptions.originalHeader,
        DBdumpOptions.sequenceString,
    ];
    // dfmt on

    return readDbDump(dbdump(dbFile, recordNumbers, dbdumpOptions,
            options.workdir), recordNumbers, options.fastaLineWidth);
}

private auto readDbDump(S, Range)(S dbDump, Range recordNumbers, in size_t lineLength)
        if (isInputRange!S && isSomeString!(ElementType!S)
            && isInputRange!Range && is(ElementType!Range : size_t))
{
    import std.algorithm : count, filter, sort;
    import std.array : appender;
    import std.range : chunks, drop;

    immutable lineSeparator = '\n';
    immutable subrecordSeparator = ';';
    immutable recordFormat = "R %d;H %d %s;L %d %d %d;S %d %s";
    immutable numRecordLines = recordFormat.count(subrecordSeparator) + 1;

    /// Build chunks of numRecordLines lines.
    alias byRecordSplitter = dbDump => dbDump.drop(6).arrayChunks(numRecordLines);
    /// Parse chunks of numRecordLines lines into FASTA format.
    alias parseRecord = recordLines => {
        size_t recordNumber;
        size_t headerLineLength;
        string headerLine;
        size_t locationWell;
        size_t locationPulseStart;
        size_t locationPulseEnd;
        size_t sequenceLength;
        string sequence;

        auto joinedLines = recordLines.joiner(only(subrecordSeparator)).array;

        // dfmt off
        int numMatches = joinedLines
            .formattedRead!recordFormat(
                recordNumber,
                headerLineLength,
                headerLine,
                locationWell,
                locationPulseStart,
                locationPulseEnd,
                sequenceLength,
                sequence,
            );
        // dfmt on
        assert(numMatches == 8, format!"%d matches in chunk: `%s`"(numMatches, joinedLines.array));

        bool isSkipping = recordNumbers.length > 0 && !recordNumbers.canFind(recordNumber);
        // dfmt off
        debug logJsonDebug(
            "isSkipping", isSkipping,
            "wantedRecordNumbers", recordNumbers.toJson,
            "recordNumber", recordNumber,
            "headerLine", headerLine,
        );
        // dfmt on

        // skip unwanted records
        if (isSkipping)
            return null;

        auto fastaData = appender!string;
        fastaData.reserve(headerLine.length + sequence.length + sequence.length / lineLength + 1);

        fastaData ~= headerLine ~ lineSeparator;
        fastaData ~= sequence.chunks(lineLength).joiner(only(lineSeparator));

        return fastaData.data;
    };

    return byRecordSplitter(dbDump).map!parseRecord
        .map!"a()"
        .cache;
}

unittest
{
    immutable testDbDump = q"EOF
        + R 4
        + M 0
        + H 104539
        @ H 26
        + S 38
        @ S 19574
        R 1
        H 22 >Sim/1/0_14 RQ=0.975
        L 0 0 14
        S 14 ggcccaggcagccc
        R 2
        H 22 >Sim/2/0_9 RQ=0.975
        L 0 0 9
        S 9 cacattgtg
        R 3
        H 23 >Sim/3/0_11 RQ=0.975
        L 0 0 11
        S 11 gagtgcagtgg
        R 4
        H 23 >Sim/4/0_4 RQ=0.975
        L 0 0 4
        S 4 gagc
        R 5
        H 24 >Sim/5/0_60 RQ=0.975
        L 0 0 60
        S 60 gagcgagcgagcgagcgagcgagcgagcgagcgagcgagcgagcgagcgagcgagcgagc
EOF".outdent;

    {
        size_t[] recordIds = [];
        auto fastaEntries = readDbDump(testDbDump.lineSplitter, recordIds, 50).array;
        // dfmt off
        assert(fastaEntries == [
            ">Sim/1/0_14 RQ=0.975\nggcccaggcagccc",
            ">Sim/2/0_9 RQ=0.975\ncacattgtg",
            ">Sim/3/0_11 RQ=0.975\ngagtgcagtgg",
            ">Sim/4/0_4 RQ=0.975\ngagc",
            ">Sim/5/0_60 RQ=0.975\ngagcgagcgagcgagcgagcgagcgagcgagcgagcgagcgagcgagcga\ngcgagcgagc",
        ], fastaEntries.to!string);
        // dfmt on
    }
    {
        size_t[] recordIds = [1, 3];
        auto fastaEntries = readDbDump(testDbDump.lineSplitter, recordIds, 50).array;
        // dfmt off
        assert(fastaEntries == [
            ">Sim/1/0_14 RQ=0.975\nggcccaggcagccc",
            ">Sim/3/0_11 RQ=0.975\ngagtgcagtgg",
        ], fastaEntries.to!string);
        // dfmt on
    }
}

/// Build a .dam file with the given set of FASTA records.
string buildDamFile(Range, Options)(Range fastaRecords, Options options)
        if (hasOption!(Options, "dbsplitOptions", isOptionsList) && hasOption!(Options, "workdir",
            isSomeString) && isInputRange!Range && isSomeString!(ElementType!Range))
{
    immutable tempDbNameTemplate = "auxiliary-XXXXXX";

    auto tempDbTemplate = buildPath(options.workdir, tempDbNameTemplate);
    auto tempDb = mkstemp(tempDbTemplate, damFileExtension);

    tempDb.file.close();
    remove(tempDb.name);
    fasta2dam(tempDb.name, fastaRecords, options.workdir);
    dbsplit(tempDb.name, options.dbsplitOptions, options.workdir);

    return tempDb.name;
}

unittest
{
    import djunctor.util.tempfile : mkdtemp;
    import std.file : rmdirRecurse, isFile;

    // dfmt off
    auto fastaRecords = [
        ">Sim/1/0_14 RQ=0.975\nggcccacccaggcagccc",
        ">Sim/3/0_11 RQ=0.975\ngagtgcgtgcagtgg",
    ];
    // dfmt on
    struct Options
    {
        string[] dbsplitOptions;
        string workdir;
    }

    auto tmpDir = mkdtemp("./.unittest-XXXXXX");
    auto options = Options([], tmpDir);
    scope (exit)
        rmdirRecurse(tmpDir);

    string dbName = buildDamFile(fastaRecords[], options);

    assert(dbName.isFile);
    foreach (hiddenDbFile; getHiddenDbFiles(dbName))
    {
        assert(hiddenDbFile.isFile);
    }
}

/**
    Self-dalign dbFile and build consensus using daccord.

    Returns: list of consensus DBs.
*/
string getConsensus(Options)(in string dbFile, in size_t readId, in Options options)
        if (hasOption!(Options, "daccordOptions", isOptionsList)
            && hasOption!(Options, "dalignerOptions",
            isOptionsList) && hasOption!(Options, "dbsplitOptions",
            isOptionsList) && hasOption!(Options, "workdir", isSomeString))
{
    static struct ModifiedOptions
    {
        string[] daccordOptions;
        string[] dalignerOptions;
        string[] dbsplitOptions;
        string workdir;
    }

    auto readIdx = readId - 1;
    // dfmt off
    auto consensusDbs = getConsensus(dbFile, const(ModifiedOptions)(
        options.daccordOptions ~ format!"%s%d,%d"(cast(string) DaccordOptions.readInterval, readIdx, readIdx),
        options.dalignerOptions,
        options.dbsplitOptions,
        options.workdir,
    ));
    // dfmt on

    if (consensusDbs.length == 0)
    {
        throw new Exception("empty consensus");
    }

    assert(consensusDbs.length == 1, "too many consensus DBs");

    return consensusDbs[0];
}

///
string[] getConsensus(Options)(in string dbFile, in Options options)
        if (hasOption!(Options, "daccordOptions", isOptionsList)
            && hasOption!(Options, "dalignerOptions",
            isOptionsList) && hasOption!(Options, "dbsplitOptions",
            isOptionsList) && hasOption!(Options, "workdir", isSomeString))
{
    dalign(dbFile, options.dalignerOptions, options.workdir);
    computeErrorProfile(dbFile, options);

    // dfmt off
    auto consensusDbs = getLasFiles(dbFile, options.workdir)
        .filter!(lasFile => !lasEmpty(lasFile, dbFile, null, options.workdir))
        .map!(lasFile => daccord(dbFile, lasFile, options.daccordOptions, options.workdir))
        .array;
    // dfmt on
    foreach (consensusDb; consensusDbs)
        dbsplit(consensusDb, options.dbsplitOptions, options.workdir);

    return consensusDbs;
}

private void computeErrorProfile(Options)(in string dbFile, in Options options)
        if (hasOption!(Options, "daccordOptions", isOptionsList)
            && hasOption!(Options, "workdir", isSomeString))
{
    // dfmt off
    auto eProfOptions = options
        .daccordOptions
        .filter!(option => !option.startsWith(
            cast(string) DaccordOptions.produceFullSequences,
            //cast(string) DaccordOptions.readInterval,
            cast(string) DaccordOptions.readsPart,
            cast(string) DaccordOptions.errorProfileFileName,
        ))
        .chain(only(DaccordOptions.computeErrorProfileOnly))
        .array;
    auto lasFiles = getLasFiles(dbFile, options.workdir)
        .filter!(lasFile => !lasEmpty(lasFile, dbFile, null, options.workdir));
    // dfmt on

    foreach (lasFile; lasFiles)
    {
        // Produce error profile
        silentDaccord(dbFile, lasFile, eProfOptions, options.workdir);
    }
}

unittest
{
    import djunctor.util.tempfile : mkdtemp;
    import std.file : rmdirRecurse, isFile;

    // dfmt off
    auto fastaRecords = [
        ">Sim/1/0_1050 RQ=0.975\nattTgaggcatcagccactgcacccagccttgtgccctttctgagagccgggaagatgctcccaggagccctcg\nggaggcttccctccggtcgtcgtggccagaattgtctcctgctcgtgtggagtcggtggcgggccaggcgaatg\nggagctaccggggctgccgctttggactgctcggcatttgccccatggggctgcacaggggcccaggctggctg\nagaatgtccctgggtccaggaggcagacggaggtacagcccagcagccaggaggtgttcaggatgttccccagt\ncagcacccgtggaggggagggaggaggcagggtgggcgaggaaggtccaacagtggacggcctgcccacaagag\nagctctgagctgggagctggcagagttgctgcaagtgggtgtgggccaggactgactgggcctgtgcacctgcc\ntggatgcatcagtggtcgtggtgctgcccgggaagggcgtgaagctccctgcagccaaggatcctggaggtgca\ngacatcacccagcccaccggacaacagcctgccctacttcgaggagctctgggcagcccagccccatgtccccc\ntcacgccccaccccacactgacaaaaagaccacaggattccaacagtccaaccagggggaggccgttgaattcg\nggggacaaccagaaacgcctgaaacagagataaagagactgatatggaaaagactgggctggcatggtggctcc\ncaactgggatcccagtgcttgtgagaggccgaggcgggaggatcacttgagcccagaagttcaagaccagcgtg\nggcaacatagtgagaccccgtctcttttaaaaatccttttttaattaggcaggcataggtagttgcgtgcctgc\nttttcccagctgctagggaggtagaggcaggagaatcacgggagtttcgaagtccaaggtcacagtgagctgtg\nattgcaccactgcactccagcctgggcaacatggcaagaccccatctctaaaagaaagaaacaagaagacatgg\nagagaaatatccaa",
        ">Sim/2/0_1050 RQ=0.975\nattagagCcatcagccactgcacccagccttgtgccctttctgagagccgggaagatgctcccaggagccctcg\nggaggcttccctccggtcgtcgtggccagaattgtctcctgctcgtgtggagtcggtggcgggccaggcgaatg\nggagctaccggggctgccgctttggactgctcggcatttgccccatggggctgcacaggggcccaggctggctg\nagaatgtccctgggtccaggaggcagacggaggtacagcccagcagccaggaggtgttcaggatgttccccagt\ncagcacccgtggaggggagggaggaggcagggtgggcgaggaaggtccaacagtggacggcctgcccacaagag\nagctctgagctgggagctggcagagttgctgcaagtgggtgtgggccaggactgactgggcctgtgcacctgcc\ntggatgcatcagtggtcgtggtgctgcccgggaagggcgtgaagctccctgcagccaaggatcctggaggtgca\ngacatcacccagcccaccggacaacagcctgccctacttcgaggagctctgggcagcccagccccatgtccccc\ntcacgccccaccccacactgacaaaaagaccacaggattccaacagtccaaccagggggaggccgttgaattcg\nggggacaaccagaaacgcctgaaacagagataaagagactgatatggaaaagactgggctggcatggtggctcc\ncaactgggatcccagtgcttgtgagaggccgaggcgggaggatcacttgagcccagaagttcaagaccagcgtg\nggcaacatagtgagaccccgtctcttttaaaaatccttttttaattaggcaggcataggtagttgcgtgcctgc\nttttcccagctgctagggaggtagaggcaggagaatcacgggagtttcgaagtccaaggtcacagtgagctgtg\nattgcaccactgcactccagcctgggcaacatggcaagaccccatctctaaaagaaagaaacaagaagacatgg\nagagaaatatccaa",
        ">Sim/3/0_1050 RQ=0.975\nattagaggcatcagccactgcacccagccttgtgccctttctgagagccgggaagatgctcccaggagccctcg\nggaggcttccctccggtcgtcgtggccagaattgtctcctgctcgtgtggagtcggtggcgggccaggcgaatg\nggagctaccggggctgccgctttggactgctcggcatttgccccatggggctgcacaggggcccaggctggctg\nagaatgtccctgggtccaggaggcagacggaggtacagcccagcagccaggaggtgttcaggatgttccccagt\ncagcacccgtggaggggagggaggaggcagggtgggcgaggaaggtccaacagtggacggcctgcccacaagag\nagctctgagctgggagctggcagagttgctgcaagtgggtgtgggccaggactgactgggcctgtgcacctgcc\ntggatgcatcagtggtcgtggtgctgcccgggaagggcgtgaagctccctgcagccaaggatcctggaggtgca\ngacatcacccagcccaccggacaacagcctgccctacttcgaggagctctgggcagcccagccccatgtccccc\ntcacgccccaccccacactgacaaaaagaccacaggattccaacagtccaaccagggggaggccgttgaattcg\nggggacaaccagaaacgcctgaaacagagataaagagactgatatggaaaagactgggctggcatggtggctcc\ncaactgggatcccagtgcttgtgagaggccgaggcgggaggatcacttgagcccagaagttcaagaccagcgtg\nggcaacatagtgagaccccgtctcttttaaaaatccttttttaattaggcaggcataggtagttgcgtgcctgc\nttttcccagctgctagggaggtagaggcaggagaatcacgggagtttcgaagtccaaggtcacagtgagctgtg\nattgcaccactgcactccagcctgggcaacatggcaagaccccatctctaaaagaaagaaacaagaagacatgg\nagagaaatatccaa",
    ];

    // dfmt on
    struct Options
    {
        string[] dbsplitOptions;
        string[] dalignerOptions;
        string[] daccordOptions;
        string[] dbdumpOptions;
        size_t fastaLineWidth;
        string workdir;
    }

    auto tmpDir = mkdtemp("./.unittest-XXXXXX");
    // dfmt off
    auto options = Options(
        [],
        [DalignerOptions.minAlignmentLength ~ "15"],
        [],
        [
            DBdumpOptions.readNumber,
            DBdumpOptions.originalHeader,
            DBdumpOptions.sequenceString,
        ],
        74,
        tmpDir,
    );
    // dfmt on
    scope (exit)
        rmdirRecurse(tmpDir);

    string dbName = buildDamFile(fastaRecords[], options);
    string[] consensusDbs = getConsensus(dbName, options);
    assert(consensusDbs.length >= 1);
    auto consensusFasta = getFastaEntries(consensusDbs[0], cast(size_t[])[], options);
    auto expectedSequence = fastaRecords[$ - 1].lineSplitter.drop(1).joiner.array;
    auto consensusSequence = consensusFasta.front.lineSplitter.drop(1).joiner.array;

    assert(expectedSequence == consensusSequence,
            format!"expected %s but got %s"(expectedSequence, consensusSequence));
}

string[] getLasFiles(in string dbA, in string baseDirectory)
{
    return getLasFiles(dbA, null, baseDirectory);
}

string[] getLasFiles(in string dbA, in string dbB, in string baseDirectory)
{
    import std.algorithm : max;

    static struct InferredParams
    {
        size_t numBlocks;
        string directory;
        string fileStem;

        this(string dbFile)
        {
            if (dbFile != null)
            {
                this.numBlocks = getNumBlocks(dbFile);
                this.directory = dbFile.dirName;
                this.fileStem = dbFile.baseName.stripExtension;
            }
        }

        string fileNamePart(in size_t blockIdx = 0)
        {
            assert(blockIdx <= this.numBlocks);
            if (this.numBlocks == 1)
                return this.fileStem;
            else
                return format!"%s.%d"(this.fileStem, blockIdx + 1);
        }
    }

    immutable fileTemplate = "%s/%s.%s.las";
    auto aOptions = InferredParams(dbA);
    auto bOptions = dbB == null ? aOptions : InferredParams(dbB);
    size_t numLasFiles = aOptions.numBlocks * max(1, bOptions.numBlocks);
    string[] fileList;

    fileList.length = numLasFiles;
    foreach (size_t i; 0 .. aOptions.numBlocks)
    {
        foreach (size_t j; 0 .. bOptions.numBlocks)
        {
            fileList[i * bOptions.numBlocks + j] = format!fileTemplate(baseDirectory,
                    aOptions.fileNamePart(i), bOptions.fileNamePart(j));
        }
    }

    return fileList;
}

bool lasFilesGenerated(in string dbA, in string baseDirectory)
{
    return lasFilesGenerated(dbA, null, baseDirectory);
}

bool lasFilesGenerated(in string dbA, in string dbB, in string baseDirectory)
{
    // dfmt off
    return getLasFiles(dbA, dbB, baseDirectory)
        .map!(lasFile => lasFile.exists)
        .all;
    // dfmt on
}

size_t getNumContigs(Options)(in string damFile, in Options options)
        if (hasOption!(Options, "workdir", isSomeString))
{
    immutable contigNumFormat = "+ R %d";
    immutable contigNumFormatStart = contigNumFormat[0 .. 4];
    size_t numContigs;
    size_t[] empty;
    // dfmt off
    auto matchingLine = dbdump(damFile, empty, [], options.workdir)
        .filter!(line => line.startsWith(contigNumFormatStart))
        .front;
    // dfmt on

    if (!matchingLine)
    {
        auto errorMessage = format!"could not read the contig count in `%s`"(damFile);
        throw new DazzlerCommandException(errorMessage);
    }

    if (formattedRead!contigNumFormat(matchingLine, numContigs) != 1)
    {
        auto errorMessage = format!"could not read the contig count in `%s`"(damFile);
        throw new DazzlerCommandException(errorMessage);
    }

    return numContigs;
}

auto getScaffoldStructure(Options)(in string damFile, in Options options)
        if (hasOption!(Options, "workdir", isSomeString))
{
    immutable string[] dbshowOptions = [DBshowOptions.noSequence];

    auto rawScaffoldInfo = dbshow(damFile, dbshowOptions, options.workdir);

    return ScaffoldStructureReader(rawScaffoldInfo);
}

alias ScaffoldPart = Algebraic!(ContigPart, GapPart);

struct ContigPart
{
    size_t globalContigId;
    size_t scaffoldId;
    size_t contigId;
    size_t begin;
    size_t end;
    string header;

    invariant
    {
        assert(begin < end);
    }

    @property size_t length() const pure nothrow
    {
        return end - begin;
    }
}

struct GapPart
{
    size_t beginGlobalContigId;
    size_t endGlobalContigId;
    size_t scaffoldId;
    size_t beginContigId;
    size_t endContigId;
    size_t begin;
    size_t end;

    invariant
    {
        assert(begin < end);
    }

    @property size_t length() const pure nothrow
    {
        return end - begin;
    }
}

private struct ScaffoldStructureReader
{
    static immutable scaffoldInfoLineFormat = "%s:: Contig %d[%d,%d]";
    alias RawScaffoldInfo = typeof("".lineSplitter);

    private RawScaffoldInfo rawScaffoldInfo;
    private ContigPart lastContigPart;
    private ScaffoldPart currentPart;
    private bool _empty;

    this(string rawScaffoldInfo)
    {
        this.rawScaffoldInfo = rawScaffoldInfo.lineSplitter;
        // Force the first element to be a contigPart.
        this.currentPart = GapPart();
        // Make `scaffoldId`s start at 0.
        this.lastContigPart.scaffoldId = -1UL;

        if (!empty)
        {
            popFront();
        }
    }

    void popFront()
    {
        assert(!empty, "Attempting to popFront an empty ScaffoldStructureReader");

        if (rawScaffoldInfo.empty)
        {
            _empty = true;

            return;
        }

        auto nextContigPart = ContigPart(lastContigPart.globalContigId + 1);

        // dfmt off
        rawScaffoldInfo.front.formattedRead!scaffoldInfoLineFormat(
            nextContigPart.header,
            nextContigPart.contigId,
            nextContigPart.begin,
            nextContigPart.end,
        );
        // dfmt on
        if (lastContigPart.contigId >= nextContigPart.contigId)
        {
            nextContigPart.scaffoldId = lastContigPart.scaffoldId + 1;
        }

        if (currentPart.peek!GapPart !is null
                || lastContigPart.scaffoldId != nextContigPart.scaffoldId)
        {
            assert(nextContigPart.header[$ - 1] == ' ');
            // Remove the trailing space
            nextContigPart.header = nextContigPart.header[0 .. $ - 1];
            lastContigPart = nextContigPart;
            currentPart = nextContigPart;
            rawScaffoldInfo.popFront();
        }
        else
        {
            // dfmt off
            currentPart = GapPart(
                lastContigPart.globalContigId,
                nextContigPart.globalContigId,
                lastContigPart.scaffoldId,
                lastContigPart.contigId,
                nextContigPart.contigId,
                lastContigPart.end,
                nextContigPart.begin,
            );
            // dfmt on
        }
    }

    @property ScaffoldPart front() const
    {
        assert(!empty, "Attempting to fetch the front of an empty ScaffoldStructureReader");
        return currentPart;
    }

    @property bool empty() const
    {
        return _empty;
    }

    ScaffoldStructureReader save()
    {
        ScaffoldStructureReader copy;

        copy.rawScaffoldInfo = this.rawScaffoldInfo.save;
        copy.lastContigPart = this.lastContigPart;
        copy.currentPart = this.currentPart;
        copy._empty = this._empty;

        return copy;
    }
}

unittest
{
    auto exampleDump = q"EOS
>reference_mod/1/0_837550 RQ=0.850 :: Contig 0[0,8300]
>reference_mod/1/0_837550 RQ=0.850 :: Contig 1[12400,20750]
>reference_mod/1/0_837550 RQ=0.850 :: Contig 2[29200,154900]
>reference_mod/1/0_837550 RQ=0.850 :: Contig 3[159900,169900]
>reference_mod/1/0_837550 RQ=0.850 :: Contig 4[174900,200650]
>reference_mod/1/0_837550 RQ=0.850 :: Contig 5[203650,216400]
>reference_mod/1/0_837550 RQ=0.850 :: Contig 6[218900,235150]
>reference_mod/1/0_837550 RQ=0.850 :: Contig 7[238750,260150]
>reference_mod/1/0_837550 RQ=0.850 :: Contig 8[263650,837500]
>reference_mod/2/0_1450 RQ=0.850 :: Contig 0[0,1450]
EOS";

    auto reader = ScaffoldStructureReader(exampleDump);
    auto scaffoldStructure = reader.array;

    // dfmt off
    assert(scaffoldStructure == [
        ScaffoldPart(ContigPart(
            1, 0, 0, 0, 8300,
            ">reference_mod/1/0_837550 RQ=0.850",
        )),
        ScaffoldPart(GapPart(1, 2, 0, 0, 1, 8300, 12400)),
        ScaffoldPart(ContigPart(
            2, 0, 1, 12400, 20750,
            ">reference_mod/1/0_837550 RQ=0.850",
        )),
        ScaffoldPart(GapPart(2, 3, 0, 1, 2, 20750, 29200)),
        ScaffoldPart(ContigPart(
            3, 0, 2, 29200, 154900,
            ">reference_mod/1/0_837550 RQ=0.850",
        )),
        ScaffoldPart(GapPart(3, 4, 0, 2, 3, 154900, 159900)),
        ScaffoldPart(ContigPart(
            4, 0, 3, 159900, 169900,
            ">reference_mod/1/0_837550 RQ=0.850",
        )),
        ScaffoldPart(GapPart(4, 5, 0, 3, 4, 169900, 174900)),
        ScaffoldPart(ContigPart(
            5, 0, 4, 174900, 200650,
            ">reference_mod/1/0_837550 RQ=0.850",
        )),
        ScaffoldPart(GapPart(5, 6, 0, 4, 5, 200650, 203650)),
        ScaffoldPart(ContigPart(
            6, 0, 5, 203650, 216400,
            ">reference_mod/1/0_837550 RQ=0.850",
        )),
        ScaffoldPart(GapPart(6, 7, 0, 5, 6, 216400, 218900)),
        ScaffoldPart(ContigPart(
            7, 0, 6, 218900, 235150,
            ">reference_mod/1/0_837550 RQ=0.850",
        )),
        ScaffoldPart(GapPart(7, 8, 0, 6, 7, 235150, 238750)),
        ScaffoldPart(ContigPart(
            8, 0, 7, 238750, 260150,
            ">reference_mod/1/0_837550 RQ=0.850",
        )),
        ScaffoldPart(GapPart(8, 9, 0, 7, 8, 260150, 263650)),
        ScaffoldPart(ContigPart(
            9, 0, 8, 263650, 837500,
            ">reference_mod/1/0_837550 RQ=0.850",
        )),
        ScaffoldPart(ContigPart(
            10, 1, 0, 0, 1450,
            ">reference_mod/2/0_1450 RQ=0.850",
        )),
    ]);
    // dfmt on
}

/**
    Get the hidden files comprising the designated mask.
*/
auto getMaskFiles(in string dbFile, in string maskDestination)
{
    auto destinationDir = maskDestination.dirName;
    auto maskName = maskDestination.baseName;
    auto dbName = dbFile.baseName.stripExtension;
    auto maskHeader = format!"%s/.%s.%s.anno"(destinationDir, dbName, maskName);
    auto maskData = format!"%s/.%s.%s.data"(destinationDir, dbName, maskName);

    return tuple!("header", "data")(maskHeader, maskData);
}

/// Thrown on failure while reading a Dazzler mask.
///
/// See_Also: `readMask`
class MaskReaderException : Exception
{
    pure nothrow @nogc @safe this(string msg, string file = __FILE__,
            size_t line = __LINE__, Throwable nextInChain = null)
    {
        super(msg, file, line, nextInChain);
    }
}

private
{
    alias MaskHeaderEntry = int;
    alias MaskDataPointer = long;
    alias MaskDataEntry = int;
}

/**
    Read the `Region`s of a Dazzler mask for `dbFile`.

    Throws: MaskReaderException
    See_Also: `writeMask`, `getMaskFiles`
*/
Region[] readMask(Region, Options)(in string dbFile, in string maskDestination, in Options options)
        if (hasOption!(Options, "workdir", isSomeString))
{
    alias _enforce = enforce!MaskReaderException;

    auto maskFileNames = getMaskFiles(dbFile, maskDestination);
    auto maskHeader = readMaskHeader(maskFileNames.header);
    auto maskData = getBinaryFile!MaskDataEntry(maskFileNames.data);

    auto maskRegions = appender!(Region[]);
    alias RegionContigId = typeof(maskRegions.data[0].tag);
    alias RegionBegin = typeof(maskRegions.data[0].begin);
    alias RegionEnd = typeof(maskRegions.data[0].end);
    auto numReads = getNumContigs(dbFile, options).to!int;

    size_t currentContig = 1;

    _enforce(maskHeader.numReads == numReads, "mask does not match DB");
    _enforce(maskHeader.size == 0, "corrupted mask: expected 0");
    _enforce(maskHeader.dataPointers.length == numReads + 1,
            "corrupted mask: unexpected number of data pointers");

    foreach (dataPtrRange; maskHeader.dataPointers[].slide!(No.withPartial)(2))
    {
        auto dataPtrs = dataPtrRange.map!(ptr => ptr / MaskDataEntry.sizeof)
            .takeExactly!2;
        _enforce(0 <= dataPtrs[0] && dataPtrs[0] <= dataPtrs[1]
                && dataPtrs[1] <= maskData.length, "corrupted mask: data pointer out of bounds");
        _enforce(dataPtrs[0] % 2 == 0 && dataPtrs[1] % 2 == 0,
                "corrupted mask: non-sense data pointers");

        foreach (interval; maskData[dataPtrs[0] .. dataPtrs[1]].chunks(2))
        {
            enforce!MaskReaderException(interval.length == 2 && 0 <= interval[0]
                    && interval[0] <= interval[1], "corrupted mask: invalid interval");

            Region newRegion;
            newRegion.tag = currentContig.to!RegionContigId;
            newRegion.begin = interval[0].to!RegionBegin;
            newRegion.end = interval[1].to!RegionEnd;

            maskRegions ~= newRegion;
        }

        ++currentContig;
    }

    return maskRegions.data;
}

private auto readMaskHeader(in string fileName)
{
    auto headerFile = File(fileName, "rb");
    MaskHeaderEntry[2] headerBuffer;
    auto numPointers = (headerFile.size - headerBuffer.sizeof) / MaskDataPointer.sizeof;
    auto pointerBuffer = uninitializedArray!(MaskDataPointer[])(numPointers);

    enforce!MaskReaderException(headerFile.rawRead(headerBuffer).length == headerBuffer.length,
            format!"error while reading mask header `%s`: file too short"(fileName));
    enforce!MaskReaderException(headerFile.rawRead(pointerBuffer).length == numPointers,
            format!"error while reading mask header `%s`: file too short"(fileName));

    return tuple!("numReads", "size", "dataPointers")(headerBuffer[0],
            headerBuffer[1], pointerBuffer);
}

private T[] getBinaryFile(T)(in string fileName)
{
    auto file = File(fileName, "rb");
    auto bufferLength = file.size / T.sizeof;
    auto dataBuffer = file.rawRead(uninitializedArray!(T[])(bufferLength));

    enforce!MaskReaderException(dataBuffer.length == bufferLength,
            format!"error while reading binary file `%s`: expected %d bytes of data but read only %d"(
                fileName, bufferLength * T.sizeof, dataBuffer.length * T.sizeof));

    return dataBuffer;
}

/**
    Write the list of regions to a Dazzler mask for `dbFile`.

    See_Also: `readMask`, `getMaskFiles`
*/
void writeMask(Region, Options)(in string dbFile, in string maskDestination,
        in Region[] regions, in Options options)
        if (hasOption!(Options, "workdir", isSomeString))
{
    // dfmt off
    alias MaskRegion = Tuple!(
        MaskHeaderEntry, "tag",
        MaskDataEntry, "begin",
        MaskDataEntry, "end",
    );
    // dfmt on

    if (regions.length == 0)
    {
        // dfmt off
        logJsonDiagnostic(
            "notice", "skipping empty mask",
            "dbFile", dbFile,
            "maskDestination", maskDestination,
        );
        // dfmt on

        return;
    }

    auto maskFileNames = getMaskFiles(dbFile, maskDestination);
    auto maskHeader = File(maskFileNames.header, "wb");
    auto maskData = File(maskFileNames.data, "wb");

    // dfmt off
    auto maskRegions = regions
        .map!(region => MaskRegion(
            region.tag.to!MaskHeaderEntry,
            region.begin.to!MaskDataEntry,
            region.end.to!MaskDataEntry,
        ))
        .array;
    // dfmt on
    maskRegions.sort();

    auto numReads = getNumContigs(dbFile, options).to!MaskHeaderEntry;
    MaskHeaderEntry size = 0; // this seems to be zero always (see DAMASKER/TANmask.c:422)
    MaskHeaderEntry currentContig = 1;
    MaskDataPointer dataPointer = 0;

    maskHeader.rawWrite([numReads, size]);
    maskHeader.rawWrite([dataPointer]);
    foreach (maskRegion; maskRegions)
    {
        assert(maskRegion.tag >= currentContig);

        while (maskRegion.tag > currentContig)
        {
            maskHeader.rawWrite([dataPointer]);
            ++currentContig;
        }

        if (maskRegion.tag == currentContig)
        {
            maskData.rawWrite([maskRegion.begin, maskRegion.end]);
            dataPointer += typeof(maskRegion.begin).sizeof + typeof(maskRegion.end).sizeof;
        }
    }

    foreach (emptyContig; currentContig .. numReads + 1)
    {
        maskHeader.rawWrite([dataPointer]);
    }
}

/// Options for `daccord`.
enum DaccordOptions : string
{
    /// number of threads (default 4)
    numberOfThreads = "-t",
    /// window size (default 40)
    windowSize = "-w",
    /// advance size (default 10)
    advanceSize = "-a",
    /// max depth (default 18446744073709551615)
    maxDepth = "-d",
    /// produce full sequences (default 0)
    produceFullSequences = "-f",
    /// verbosity (default 18446744073709551615)
    verbosity = "-V",
    /// read interval (default 0,18446744073709551615)
    readInterval = "-I",
    /// reads part (default 0,1)
    readsPart = "-J",
    /// error profile file name (default input.las.eprof)
    errorProfileFileName = "-E",
    /// minimum window coverage (default 3)
    minWindowCoverage = "-m",
    /// maximum window error (default 18446744073709551615)
    maxWindowError = "-e",
    /// minimum length of output (default 0)
    minLengthOfOutput = "-l",
    /// minimum k-mer filter frequency (default 0)
    minKMerFilterFrequency = "--minfilterfreq",
    /// maximum k-mer filter frequency (default 2)
    maxKMerFilterFrequency = "--maxfilterfreq",
    /// temporary file prefix (default daccord_ozelot_4500_1529654843)
    temporaryFilePrefix = "-T",
    /// maximum number of alignments considered per read (default 5000)
    maxAlignmentsPerRead = "-D",
    /// maximum number of alignments considered per read (default 0)
    maxAlignmentsPerReadVard = "--vard",
    /// compute error profile only (default disable)
    computeErrorProfileOnly = "--eprofonly",
    /// compute error distribution estimate (default disable)
    computeErrorDistributionEstimate = "--deepprofileonly",
    /// kmer size (default 8)
    kmerSize = "-k",
}

/// Options for `daligner`.
enum DalignerOptions : string
{
    verbose = "-v",
    /// If the -b option is set, then the daligner assumes the data has a
    /// strong compositional bias (e.g. >65% AT rich).
    strongCompositionalBias = "-b",
    /// If the -A option is set (“A” for “asymmetric”) then just overlaps
    /// where the a-read is in block X and the b-read is in block Y are
    /// reported, and if X = Y then it further reports only those overlaps
    /// where the a-read index is less than the b-read index.
    asymmetric = "-A",
    /// If the -I option is set (“I” for “identity”) then when X = Y, overlaps
    /// between different portions of the same read will also be found and
    /// reported.
    identity = "-I",
    /// Search code looks for a pair of diagonal bands of width 2^^w
    /// (default 26 = 64) that contain a collection of exact matching k-mers
    /// (default 14) between the two reads, such that the total number of
    /// bases covered by the k-mer hits is h (default 35).
    kMerSize = "-k",
    /// ditto
    bandWidth = "-w",
    /// ditto
    hitBaseCoverage = "-h",
    /// Suppresses the use of any k-mer that occurs more than t times in
    /// either the subject or target block.
    maxKmerOccurence = "-t",
    /// Let the program automatically select a value of t that meets a given
    /// memory usage limit specified (in Gb) by the -M parameter.
    maxKmerMemory = "-M",
    tempDir = "-P",
    /// Searching for local alignments involving at least -l base pairs
    /// (default 1000) or more, that have an average correlation rate of
    /// -e (default 70%).
    minAlignmentLength = "-l",
    /// ditto
    averageCorrelationRate = "-e",
    /// The local alignments found will be output in a sparse encoding where
    /// a trace point on the alignment is recorded every -s base pairs of
    /// the a-read (default 100bp).
    tracePointDistance = "-s",
    /// By setting the -H parameter to say N, one alters daligner so that it
    /// only reports overlaps where the a-read is over N base-pairs long.
    minAReadLength = "-H",
    /// The program runs with 4 threads by default, but this may be set to
    /// any power of 2 with the -T option.
    numThreads = "-T",
    /// If there are one or more interval tracks specified with the -m option
    /// (m for mask), then the reads of the DB or DB’s to which the track
    /// applies are soft masked with the union of the intervals of all the
    /// interval tracks that apply, that is any k-mers that contain any bases
    /// in any of the masked intervals are ignored for the purposes of seeding
    /// a match.
    masks = "-m",
}

/// Options for `damapper`.
enum DamapperOptions : string
{
    verbose = "-v",
    /// If the -b option is set, then the daligner assumes the data has a
    /// strong compositional bias (e.g. >65% AT rich).
    strongCompositionalBias = "-b",
    /// Search code looks for a pair of diagonal bands of width 2^^w
    /// (default 26 = 64) that contain a collection of exact matching k-mers
    /// (default 14) between the two reads, such that the total number of
    /// bases covered by the k-mer hits is h (default 35).
    kMerSize = "-k",
    /// Suppresses the use of any k-mer that occurs more than t times in
    /// either the subject or target block.
    maxKmerOccurence = "-t",
    /// Let the program automatically select a value of t that meets a given
    /// memory usage limit specified (in Gb) by the -M parameter.
    maxKmerMemory = "-M",
    /// ditto
    averageCorrelationRate = "-e",
    /// The local alignments found will be output in a sparse encoding where
    /// a trace point on the alignment is recorded every -s base pairs of
    /// the a-read (default 100bp).
    tracePointDistance = "-s",
    /// The program runs with 4 threads by default, but this may be set to
    /// any power of 2 with the -T option.
    numThreads = "-T",
    /// If there are one or more interval tracks specified with the -m option
    /// (m for mask), then the reads of the DB or DB’s to which the track
    /// applies are soft masked with the union of the intervals of all the
    /// interval tracks that apply, that is any k-mers that contain any bases
    /// in any of the masked intervals are ignored for the purposes of seeding
    /// a match.
    masks = "-m",
    /// If the -n option is given then all chains that are within the given
    /// fraction of the best are also reported, e.g. -n.95 reports all
    /// matches within 95% of the top match.
    bestMatches = "-n",
    /// The -p option requests that damapper produce a repeat profile track
    /// for each read.
    repeatProfileTrack = "-p",
    /// The parameter -z asks that LAs are sorted in pile order as opposed to
    /// map order (see the -a option of daligner for which this is the
    /// negation).
    sortPileOrder = "-z",
    /// If the -C option is set, then damapper also outputs a file Y.X.las
    /// for a given block pair that contains all the same matches as in
    /// X.Y.las but where the A-read is a contig of the reference and the
    /// B-read is a mapped read. And if the -N options is set, then the file
    /// Y.X.las is not produced.
    symmetric = "-C",
    /// ditto
    oneDirection = "-N",
}

/// Options for `DBdump`.
enum DBdumpOptions
{
    readNumber = "-r",
    originalHeader = "-h",
    sequenceString = "-s",
    sNROfACGTChannels = "-a",
    intrinsicQualityVector = "-i",
    quivaValues = "-q",
    repeatProfileVector = "-p",
    masks = "-m",
    untrimmedDatabase = "-u",
    upperCase = "-U",
}

/// Options for `DBshow`.
enum DBshowOptions
{
    untrimmedDatabase = "-u",
    showQuiva = "-q",
    showArrowPulseSequence = "-a",
    noSequence = "-n",
    masks = "-m",
    produceQuivaFile = "-Q",
    produceArrowFile = "-A",
    upperCase = "-U",
    fastaLineWidth = "-w",
}

/// Options for `fasta2DAM` and `fasta2DB`.
enum Fasta2DazzlerOptions
{
    verbose = "-v",
    /// Import files listed 1/line in given file.
    fromFile = "-f",
    /// Import data from stdin, use optional name as data source.
    fromStdin = "-i",
}

/// Options for `LAdump`.
enum LAdumpOptions
{
    coordinates = "-c",
    numDiffs = "-d",
    tracePoints = "-t",
    lengths = "-l",
    properOverlapsOnly = "-o",
}

private
{
    bool lasEmpty(in string lasFile, in string dbA, in string dbB, in string workdir)
    {
        auto dumpHeader = ladump(lasFile, dbA, dbB, [], workdir);
        size_t numParts;

        dumpHeader.formattedRead!"+ P %d"(numParts);

        return numParts == 0;
    }

    void dalign(in string refDam, in string[] dalignerOpts, in string workdir)
    {
        dalign([refDam], dalignerOpts, workdir);
    }

    void dalign(in string refDam, in string readsDam, in string[] dalignerOpts, in string workdir)
    {
        dalign([refDam, readsDam], dalignerOpts, workdir);
    }

    void dalign(in string[] dbList, in string[] dalignerOpts, in string workdir)
    {
        assert(dbList.length >= 1);
        auto isSelfAlignment = dbList.length == 1;
        auto additionalOptions = only(isSelfAlignment ? DalignerOptions.identity : null);
        auto inputFiles = isSelfAlignment ? [dbList[0], dbList[0]] : dbList;
        const(string[]) inputFilesRelativeToWorkDir = inputFiles.map!(
                f => f.relativeToWorkdir(workdir)).array;

        executeCommand(chain(only("daligner"), additionalOptions, dalignerOpts,
                inputFilesRelativeToWorkDir), workdir);
    }

    void damapper(in string refDam, in string readsDam, in string[] damapperOpts, in string workdir)
    {
        damapper([refDam, readsDam], damapperOpts, workdir);
    }

    void damapper(in string[] dbList, in string[] damapperOpts, in string workdir)
    {
        const(string[]) dbListRelativeToWorkDir = dbList.map!(
                f => f.relativeToWorkdir(workdir)).array;

        executeCommand(chain(only("damapper", DamapperOptions.symmetric),
                damapperOpts, dbListRelativeToWorkDir), workdir);
    }

    string daccord(in string dbFile, in string lasFile, in string[] daccordOpts, in string workdir)
    {
        alias esc = escapeShellCommand;
        string daccordedDb = dbFile.stripExtension.to!string ~ "-daccord.dam";

        // dfmt off
        executeShell(chain(
            only("daccord"),
            only(esc(daccordOpts)),
            only(esc(lasFile.relativeToWorkdir(workdir))),
            only(esc(dbFile.relativeToWorkdir(workdir))),
            only("|"),
            only("fasta2DAM", Fasta2DazzlerOptions.fromStdin),
            only(esc(daccordedDb.relativeToWorkdir(workdir))),
        ), workdir);
        // dfmt on

        return daccordedDb;
    }

    void silentDaccord(in string dbFile, in string lasFile, in string[] daccordOpts,
            in string workdir)
    {
        // dfmt off
        executeCommand(chain(
            only("daccord"),
            daccordOpts,
            only(lasFile.relativeToWorkdir(workdir)),
            only(dbFile.relativeToWorkdir(workdir)),
        ), workdir);
        // dfmt on
    }

    void buildSubsetDb(R)(in string inDbFile, in string outDbFile, R readIds, in string workdir)
    {
        alias esc = escapeShellCommand;
        // dfmt off
        auto escapedReadIds = readIds
            .map!(to!size_t)
            .map!(to!string)
            .map!esc;
        // dfmt on

        // dfmt off
        executeShell(chain(
            only("DBshow"),
            only(esc(inDbFile.relativeToWorkdir(workdir))),
            escapedReadIds,
            only("|"),
            only("fasta2DAM", Fasta2DazzlerOptions.fromStdin),
            only(esc(outDbFile.relativeToWorkdir(workdir))),
        ), workdir);
        // dfmt on
    }

    void fasta2dam(Range)(in string outFile, Range fastaRecords, in string workdir)
            if (isInputRange!(Unqual!Range) && isSomeString!(ElementType!(Unqual!Range)))
    {
        import std.algorithm : each, joiner;
        import std.file : getcwd;
        import std.process : Config, pipeProcess, Redirect, wait;
        import std.range : chunks;

        immutable writeChunkSize = 1024 * 1024;
        auto outFileArg = outFile.relativeToWorkdir(workdir);
        auto command = ["fasta2DAM", Fasta2DazzlerOptions.fromStdin, outFileArg];
        //dfmt off
        logJsonDebug(
            "action", "execute",
            "type", "pipe",
            "command", command.map!Json.array,
            "input", fastaRecords.map!Json.array,
            "state", "pre",
        );
        //dfmt on
        auto process = pipeProcess(["fasta2DAM", Fasta2DazzlerOptions.fromStdin,
                outFileArg], Redirect.stdin, null, // env
                Config.none, workdir);
        //dfmt off
        fastaRecords
            .filter!(fastaRecord => parseFastaRecord(fastaRecord).length >= minSequenceLength)
            .joiner(only('\n'))
            .chain("\n")
            .chunks(writeChunkSize)
            .each!(chunk => process.stdin.write(chunk.array));
        //dfmt on
        process.stdin.close();
        auto exitStatus = wait(process.pid);
        if (exitStatus != 0)
        {
            throw new DazzlerCommandException(
                    format!"command `fasta2dam` failed with exit code %d"(exitStatus));
        }

        return;
    }

    void fasta2dam(in string inFile, in string outFile, in string workdir)
    {
        executeCommand(only("fasta2DAM", outFile.relativeToWorkdir(workdir), inFile), workdir);
    }

    void dbsplit(in string dbFile, in string[] dbsplitOptions, in string workdir)
    {
        executeCommand(chain(only("DBsplit"), dbsplitOptions,
                only(dbFile.relativeToWorkdir(workdir))), workdir);
    }

    void lasort(in string[] lasFiles, in string workdir)
    {
        import std.file : rename;
        import std.path : setExtension;

        alias sortedName = (lasFile) => lasFile.setExtension(".S.las");

        executeCommand(chain(only("LAsort"), lasFiles), workdir);
        foreach (lasFile; lasFiles)
        {
            rename(sortedName(lasFile), lasFile);
        }
    }

    string ladump(in string lasFile, in string dbA, in string dbB,
            in string[] ladumpOpts, in string workdir)
    {
        return ladump(lasFile, dbA, dbB, [], ladumpOpts, workdir);
    }

    string ladump(in string lasFile, in string dbA, in string dbB, in size_t[] readIds,
            in string[] ladumpOpts, in string workdir)
    {
        // dfmt off
        return executeCommand(chain(
            only("LAdump"),
            ladumpOpts,
            only(
                dbA.relativeToWorkdir(workdir),
                dbB.relativeToWorkdir(workdir),
                lasFile.relativeToWorkdir(workdir)
            ),
            readIds.map!(to!string),
        ), workdir);
        // dfmt on
    }

    auto dbdump(Range)(in string dbFile, Range recordNumbers,
            in string[] dbdumpOptions, in string workdir)
            if (isForwardRange!Range && is(ElementType!Range : size_t))
    {
        static struct DBDump
        {
            static immutable lineTerminator = "\n";

            const string dbFile;
            const Range recordNumbers;
            const string[] dbdumpOptions;
            const string workdir;
            ProcessPipes dbdump;
            string currentLine;

            ~this()
            {
                if (!(dbdump.pid is null))
                    releaseProcess();
            }

            void releaseProcess()
            {
                if (!dbdump.stdout.isOpen)
                    return;

                dbdump.stdout.close();

                version (Posix)
                {
                    import core.sys.posix.signal : SIGKILL;

                    dbdump.pid.kill(SIGKILL);
                }
                else
                {
                    static assert(0, "Only intended for use on POSIX compliant OS.");
                }
                dbdump.pid.wait();
            }

            private void assertInitialized()
            {
                if (!(dbdump.pid is null))
                    return;

                auto command = chain(only("DBdump"), dbdumpOptions,
                        only(dbFile.relativeToWorkdir(workdir)), recordNumbers.map!(to!string));
                // dfmt off
                logJsonDebug(
                    "action", "execute",
                    "type", "pipe",
                    "command", command.map!Json.array,
                    "state", "pre",
                );
                // dfmt on
                dbdump = pipeProcess(command.array, Redirect.stdout, null, Config.none, workdir);

                popFront();
            }

            void popFront()
            {
                import std.string : stripRight;

                assertInitialized();
                currentLine = dbdump.stdout.readln();

                if (currentLine.empty)
                {
                    currentLine = null;
                    releaseProcess();
                }

                if (currentLine.endsWith(lineTerminator))
                    currentLine = currentLine[0 .. $ - lineTerminator.length];
            }

            @property string front()
            {
                assertInitialized();

                return currentLine;
            }

            @property bool empty()
            {
                assertInitialized();

                if (currentLine is null)
                {
                    releaseProcess();

                    return true;
                }
                else
                {
                    return false;
                }
            }
        }

        return new DBDump(dbFile, recordNumbers, dbdumpOptions, workdir);
    }

    string dbshow(in string dbFile, in string contigId, in string workdir)
    {
        return executeCommand(only("DBshow", dbFile.relativeToWorkdir(workdir), contigId), workdir);
    }

    string dbshow(in string dbFile, in string[] dbshowOptions, in string workdir)
    {
        return executeCommand(chain(only("DBshow"), dbshowOptions,
                only(dbFile.relativeToWorkdir(workdir))), workdir);
    }

    size_t getNumBlocks(in string damFile)
    {
        // see also in dazzler's DB.h:394
        //     #define DB_NBLOCK "blocks = %9d\n"  //  number of blocks
        immutable blockNumFormat = "blocks = %d";
        immutable blockNumFormatStart = blockNumFormat[0 .. 6];
        size_t numBlocks;
        auto matchingLine = File(damFile).byLine.filter!(
                line => line.startsWith(blockNumFormatStart)).front;

        if (!matchingLine)
        {
            auto errorMessage = format!"could not read the block count in `%s`"(damFile);
            throw new DazzlerCommandException(errorMessage);
        }

        if (formattedRead!blockNumFormat(matchingLine, numBlocks) != 1)
        {
            auto errorMessage = format!"could not read the block count in `%s`"(damFile);
            throw new DazzlerCommandException(errorMessage);
        }

        return numBlocks;
    }

    string executeCommand(Range)(in Range command, in string workdir = null)
            if (isInputRange!(Unqual!Range) && isSomeString!(ElementType!(Unqual!Range)))
    {
        import std.process : Config, execute;

        string output = command.executeWrapper!("command",
                sCmd => execute(sCmd, null, // env
                    Config.none, size_t.max, workdir));
        return output;
    }

    void executeShell(Range)(in Range command, in string workdir = null)
            if (isInputRange!(Unqual!Range) && isSomeString!(ElementType!(Unqual!Range)))
    {
        import std.algorithm : joiner;
        import std.process : Config, executeShell;

        string output = command.executeWrapper!("shell",
                sCmd => executeShell(sCmd.joiner(" ").array.to!string, null, // env
                    Config.none, size_t.max, workdir));
    }

    void executeScript(Range)(in Range command, in string workdir = null)
            if (isInputRange!(Unqual!Range) && isSomeString!(ElementType!(Unqual!Range)))
    {
        import std.process : Config, executeShell;

        string output = command.executeWrapper!("script",
                sCmd => executeShell(sCmd.buildScriptLine, null, // env
                    Config.none, size_t.max, workdir));
    }

    string executeWrapper(string type, alias execCall, Range)(in Range command)
            if (isInputRange!(Unqual!Range) && isSomeString!(ElementType!(Unqual!Range)))
    {
        import std.array : array;
        import std.algorithm : filter;
        import std.string : lineSplitter;

        auto sanitizedCommand = command.filter!"a != null".array;

        // dfmt off
        logJsonDebug(
            "action", "execute",
            "type", type,
            "command", sanitizedCommand.map!Json.array,
            "state", "pre",
        );
        // dfmt on
        auto result = execCall(sanitizedCommand);
        // dfmt off
        logJsonDebug(
            "action", "execute",
            "type", type,
            "command", sanitizedCommand.map!Json.array,
            "output", result
                .output[0 .. min(1024, $)]
                .lineSplitter
                .map!Json
                .array,
            "exitStatus", result.status,
            "state", "post",
        );
        // dfmt on
        if (result.status > 0)
        {
            throw new DazzlerCommandException(
                    format("process %s returned with non-zero exit code %d: %s",
                    sanitizedCommand[0], result.status, result.output));
        }

        return result.output;
    }

    string buildScriptLine(in string[] command)
    {
        return escapeShellCommand(command) ~ " | sh -sve";
    }

    string relativeToWorkdir(in string fileName, in string workdir)
    {
        return relativePath(absolutePath(fileName), absolutePath(workdir));
    }
}
