/**
    Defines bindinds and utilities to/for the dazzler commands.

    Copyright: © 2018 Arne Ludwig <arne.ludwig@posteo.de>
    License: Subject to the terms of the MIT license, as written in the
             included LICENSE file.
    Authors: Arne Ludwig <arne.ludwig@posteo.de>
*/
module djunctor.dazzler;

import djunctor.commandline : hasOption, isOptionsList;
import djunctor.djunctor : AlignmentChain, AlignmentContainer;
import djunctor.util.log;
import djunctor.util.tempfile : mkstemp;
import std.algorithm : equal, map, joiner, sort, splitter, SwapStrategy;
import std.array : array;
import std.conv : to;
import std.file : remove;
import std.format : format, formattedRead;
import std.meta : Instantiate;
import std.path : absolutePath, baseName, buildPath, dirName, relativePath,
    stripExtension, withExtension;
import std.process : escapeShellCommand;
import std.range : chain, drop, only;
import std.range.primitives : ElementType, isInputRange;
import std.stdio : File, writeln;
import std.string : lineSplitter, outdent;
import std.traits : hasMember, isIntegral, isSomeString, Unqual;
import std.typecons : Flag, No, tuple, Tuple, Yes;

/// File suffixes of hidden DB files.
private immutable hiddenDbFileSuffixes = [".bps", ".hdr", ".idx"];

/// Constant holding the .dam file extension.
immutable damFileExtension = ".dam";

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
    Provide dbFile in workdir.

    Returns: Workdir location of the dbFile.
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

AlignmentContainer!(AlignmentChain[]) getLocalAlignments(Options)(in string dbA, in Options options)
        if (hasOption!(Options, "dalignerOptions", isOptionsList) && hasOption!(Options,
            "ladumpOptions", isOptionsList) && hasOption!(Options, "workdir", isSomeString))
{
    return getLocalAlignments(dbA, null, options);
}

AlignmentContainer!(AlignmentChain[]) getLocalAlignments(Options)(in string dbA,
        in string dbB, in Options options)
        if (hasOption!(Options, "dalignerOptions", isOptionsList) && hasOption!(Options,
            "ladumpOptions", isOptionsList) && hasOption!(Options, "workdir", isSomeString))
{
    dalign(dbA, dbB, options.dalignerOptions, options.workdir);

    return processGeneratedLasFiles(dbA, dbB, options);
}

AlignmentContainer!(AlignmentChain[]) getMappings(Options)(in string dbA,
        in string dbB, in Options options)
        if (hasOption!(Options, "damapperOptions", isOptionsList) && hasOption!(Options,
            "ladumpOptions", isOptionsList) && hasOption!(Options, "workdir", isSomeString))
{
    damapper(dbA, dbB, options.damapperOptions, options.workdir);

    return processGeneratedLasFiles(dbA, dbB, options);
}

private auto processGeneratedLasFiles(Options)(in string dbA, in string dbB, in Options options)
        if (hasOption!(Options, "ladumpOptions", isOptionsList)
            && hasOption!(Options, "workdir", isSomeString))
{
    auto lasFileLists = getLasFiles(dbA, dbB);
    AlignmentContainer!(AlignmentChain[]) results;

    // TODO prevent double execution if dbB == null
    static foreach (order; typeof(lasFileLists).orders)
    {
        {
            auto lasFileList = __traits(getMember, lasFileLists, order);
            auto dbFiles = typeof(lasFileLists).getOrdered!order(dbA, dbB);

            //lasort(lasFileList);
            // dfmt off
            auto alignmentChains =
                lasFileList
                    .map!(lasFile => ladump(lasFile, dbFiles.expand, options.ladumpOptions, options.workdir))
                    .map!(lasDump => readAlignmentList(lasDump))
                    .joiner
                    .array;
            alignmentChains.sort!("a < b", SwapStrategy.stable);
            __traits(getMember, results, order) = alignmentChains;
            // dfmt on
        }
    }

    return results;
}

private auto readAlignmentList(S)(in S lasDump) if (isSomeString!S)
{
    import std.algorithm : chunkBy, count, filter;
    import std.range : chunks, drop, enumerate;

    immutable recordSeparator = ';';
    immutable chainPartFormat = "P %d %d %c %c;L %d %d;C %d %d %d %d;D %d";
    immutable numChainPartLines = chainPartFormat.count(recordSeparator) + 1;
    enum ChainPartType
    {
        start = '>',
        continuation = '-',
        alternateStart = '+',
        noChainInFile = '.',
    }
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
    /// Returns true iff part1 and part2 belong to the same chain.
    alias isChainContinuation = (part1, part2) =>
        part1.chainPartType == ChainPartType.start && part2.chainPartType == ChainPartType.continuation ||
        part1.chainPartType == ChainPartType.alternateStart && part2.chainPartType == ChainPartType.continuation ||
        part1.chainPartType == ChainPartType.continuation && part2.chainPartType == ChainPartType.continuation;
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

/**
    Get the designated set of records in FASTA format. If recordNumbers is
    empty the whole DB will be converted.
*/
auto getFastaEntries(Options, Range)(in string dbFile, Range recordNumbers, in Options options)
        if (hasOption!(Options, "dbdumpOptions", isOptionsList)
            && hasOption!(Options, "fastaLineWidth",
            isIntegral) && hasOption!(Options, "workdir", isSomeString)
            && isInputRange!Range && is(ElementType!Range == size_t))
{
    return readDbDump(dbdump(dbFile, options.dbdumpOptions, options.workdir),
            recordNumbers, options.fastaLineWidth);
}

private auto readDbDump(S, Range)(in S dbDump, Range recordNumbers, in size_t lineLength)
        if (isSomeString!S && isInputRange!Range && is(ElementType!Range == size_t))
{
    import std.algorithm : count, filter, sort;
    import std.array : appender;
    import std.range : chunks, drop;

    immutable lineSeparator = '\n';
    immutable subrecordSeparator = ';';
    immutable recordFormat = "R %d;H %d %s;L %d %d %d;S %d %s";
    immutable numRecordLines = recordFormat.count(subrecordSeparator) + 1;

    auto sortedRecordNumbers = recordNumbers.sort;
    /// Build chunks of numRecordLines lines.
    alias byRecordSplitter = dbDump => dbDump.lineSplitter.drop(6).chunks(numRecordLines);
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

        auto joinedLines = recordLines.joiner(only(subrecordSeparator));
        // dfmt off
        int numMatches = joinedLines
            .array
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

        // skip unwanted records
        if (sortedRecordNumbers.length > 0 && !sortedRecordNumbers.contains(recordNumber))
            return "";

        auto fastaData = appender!string;
        fastaData.reserve(headerLine.length + sequence.length + sequence.length / lineLength + 1);

        fastaData ~= headerLine ~ lineSeparator;
        fastaData ~= sequence.chunks(lineLength).joiner(only(lineSeparator));

        return fastaData.data;
    };

    return byRecordSplitter(dbDump).map!parseRecord.filter!`a() != ""`;
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
        auto fastaEntries = readDbDump(testDbDump, recordIds, 50).map!"a()".array;
        // dfmt off
        assert(fastaEntries == [
            ">Sim/1/0_14 RQ=0.975\nggcccaggcagccc",
            ">Sim/2/0_9 RQ=0.975\ncacattgtg",
            ">Sim/3/0_11 RQ=0.975\ngagtgcagtgg",
            ">Sim/4/0_4 RQ=0.975\ngagc",
            ">Sim/5/0_60 RQ=0.975\ngagcgagcgagcgagcgagcgagcgagcgagcgagcgagcgagcgagcga\ngcgagcgagc",
        ]);
        // dfmt on
    }
    {
        size_t[] recordIds = [1, 3];
        auto fastaEntries = readDbDump(testDbDump, recordIds, 50).map!"a()".array;
        // dfmt off
        assert(fastaEntries == [
            ">Sim/1/0_14 RQ=0.975\nggcccaggcagccc",
            ">Sim/3/0_11 RQ=0.975\ngagtgcagtgg",
        ]);
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
        ">Sim/1/0_14 RQ=0.975\nggcccaggcagccc",
        ">Sim/3/0_11 RQ=0.975\ngagtgcagtgg",
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

string[] getConsensus(Options)(in string dbFile, in Options options)
        if (hasOption!(Options, "daccordOptions", isOptionsList) && hasOption!(Options,
            "dalignerOptions", isOptionsList) && hasOption!(Options, "workdir", isSomeString))
{
    dalign(dbFile, options.dalignerOptions, options.workdir);

    // dfmt off
    return getLasFiles(dbFile)
        .a2b
        .map!(lasFile => daccord(dbFile, lasFile, options.daccordOptions, options.workdir))
        .array;
    // dfmt on
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
    auto options = Options([], ["-l15"], [], ["-r", "-h", "-s"], 74, tmpDir);
    scope (exit)
        rmdirRecurse(tmpDir);

    string dbName = buildDamFile(fastaRecords[], options);
    string[] consensusDbs = getConsensus(dbName, options);
    assert(consensusDbs.length >= 1);
    auto consensusFasta = getFastaEntries(consensusDbs[0], cast(size_t[])[], options);
    auto expectedSequence = fastaRecords[$ - 1].lineSplitter.drop(1).joiner.array;
    auto consensusSequence = consensusFasta.front()().lineSplitter.drop(1).joiner.array;

    assert(expectedSequence == consensusSequence);
}

AlignmentContainer!(string[]) getLasFiles(in string dbA)
{
    return getLasFiles(dbA, null);
}

AlignmentContainer!(string[]) getLasFiles(in string dbA, in string dbB)
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

    immutable singleFileTemplate = "%s/%s.las";
    immutable twoFileTemplate = "%s/%s.%s.las";
    auto aOptions = InferredParams(dbA);
    auto bOptions = InferredParams(dbB);
    size_t numLasFilesPerList = aOptions.numBlocks * max(1, bOptions.numBlocks);
    AlignmentContainer!(string[]) fileLists;

    fileLists.a2b.length = numLasFilesPerList;
    fileLists.b2a.length = numLasFilesPerList;
    foreach (size_t i; 0 .. aOptions.numBlocks)
    {
        if (bOptions.numBlocks > 0)
        {
            foreach (size_t j; 0 .. bOptions.numBlocks)
            {
                fileLists.a2b[i * bOptions.numBlocks + j] = format!twoFileTemplate(
                        aOptions.directory, aOptions.fileNamePart(i), bOptions.fileNamePart(j));
                fileLists.b2a[i * bOptions.numBlocks + j] = format!twoFileTemplate(
                        bOptions.directory, bOptions.fileNamePart(j), aOptions.fileNamePart(i));
            }
        }
        else
        {
            fileLists.a2b[i] = fileLists.b2a[i] = format!singleFileTemplate(aOptions.directory,
                    aOptions.fileNamePart(i));
        }
    }

    return fileLists;
}

private
{
    void dalign(in string refDam, in string[] dalignerOpts, in string workdir)
    {
        dalign(refDam, null, dalignerOpts, workdir);
    }

    void dalign(in string refDam, in string readsDam, in string[] dalignerOpts, in string workdir)
    {
        executeScript(chain(only("HPC.daligner"), dalignerOpts,
                only(refDam.relativeToWorkdir(workdir)),
                only(readsDam.relativeToWorkdir(workdir))), workdir);
    }

    void damapper(in string refDam, in string readsDam, in string[] damapperOpts, in string workdir)
    {
        executeScript(chain(only("HPC.damapper", "-C"), damapperOpts,
                only(refDam.relativeToWorkdir(workdir)), only(readsDam.relativeToWorkdir(workdir))), workdir);
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
            only("fasta2DAM", "-i"),
            only(esc(daccordedDb.relativeToWorkdir(workdir))),
        ), workdir);
        // dfmt on

        return daccordedDb;
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
        auto process = pipeProcess(["fasta2DAM", "-i", outFileArg],
                Redirect.stdin, null, // env
                Config.none, workdir);//dfmt off
        fastaRecords
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
        return executeCommand(chain(only("LAdump"), ladumpOpts,
                only(dbA.relativeToWorkdir(workdir)),
                only(dbB.relativeToWorkdir(workdir)), only(lasFile)), workdir);
    }

    string dbdump(in string dbFile, in string[] dbdumpOptions, in string workdir)
    {
        return executeCommand(chain(only("DBdump"), dbdumpOptions,
                only(dbFile.relativeToWorkdir(workdir))), workdir);
    }

    string dbshow(in string dbFile, in string contigId, in string workdir)
    {
        return executeCommand(only("DBshow", dbFile.relativeToWorkdir(workdir), contigId), workdir);
    }

    size_t getNumBlocks(in string damFile)
    {
        import std.algorithm : filter, startsWith;

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
        logDiagnostic(output);
    }

    void executeScript(Range)(in Range command, in string workdir = null)
            if (isInputRange!(Unqual!Range) && isSomeString!(ElementType!(Unqual!Range)))
    {
        import std.process : Config, executeShell;

        string output = command.executeWrapper!("script",
                sCmd => executeShell(sCmd.buildScriptLine, null, // env
                    Config.none, size_t.max, workdir));
        logDiagnostic(output);
    }

    string executeWrapper(string type, alias execCall, Range)(in Range command)
            if (isInputRange!(Unqual!Range) && isSomeString!(ElementType!(Unqual!Range)))
    {
        import std.array : array;
        import std.algorithm : filter;

        auto sanitizedCommand = command.filter!"a != null".array;

        logDiagnostic("executing " ~ type ~ ": %s", sanitizedCommand);
        auto result = execCall(sanitizedCommand);
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
        return escapeShellCommand(command) ~ " | sh -sv";
    }

    string relativeToWorkdir(in string fileName, in string workdir)
    {
        return relativePath(absolutePath(fileName), absolutePath(workdir));
    }
}
