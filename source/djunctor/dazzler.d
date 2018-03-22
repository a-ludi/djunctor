/**
    Defines bindinds and utilities to/for the dazzler commands.

    Copyright: © 2018 Arne Ludwig <arne.ludwig@posteo.de>
    License: Subject to the terms of the MIT license, as written in the
             included LICENSE file.
    Authors: Arne Ludwig <arne.ludwig@posteo.de>
*/
module djunctor.dazzler;

import djunctor.djunctor : AlignmentContainer;
import djunctor.log;
import std.algorithm : equal, splitter;
import std.array : array;
import std.conv : to;
import std.format : format, formattedRead;
import std.range : chain, only;
import std.range.primitives : ElementType, isInputRange;
import std.stdio : File, writeln;
import std.string : outdent;
import std.traits : Unqual, isSomeString;
import std.typecons : Flag, No, Yes;

/// File suffixes of hidden DB files.
private immutable hiddenDbFileSuffixes = [".bps", ".hdr", ".idx"];

/**
    Return a list of hidden files associated to every `.dam`/`.db` file. These
    files contain the actual data used in all the computation. Thus, we
    carefully check for their existence.
*/
auto getHiddenDbFiles(string dbFile)
{
    import std.algorithm : map;
    import std.path : baseName, buildPath, dirName, withExtension;

    return hiddenDbFileSuffixes.map!(suffix => buildPath(dbFile.dirName,
            "." ~ dbFile.baseName.withExtension(suffix).to!string));
}

private string workdir;

void setWorkdir(string workdir_)
{
    workdir = workdir_;
    logDebug(format!"using workdir `%s`"(workdir));
}

string getWorkdir()
{
    assert(workdir != null);

    return workdir;
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

void provideDamFileInWorkdir(in string dbFile, ProvideMethod provideMethod)
{
    import std.file : copy, symlink;
    import std.path : baseName, buildPath;
    import std.range : chain, only;

    auto allDbFiles = chain(only(dbFile), getHiddenDbFiles(dbFile));

    foreach (anyDbFile; allDbFiles)
    {
        switch (provideMethod)
        {
        case ProvideMethod.copy:
            copy(anyDbFile, buildPath(getWorkdir(), anyDbFile.baseName));
            break;
        case ProvideMethod.symlink:
            symlink(anyDbFile, buildPath(getWorkdir(), anyDbFile.baseName));
            break;
        default:
            assert(0);
        }
    }
}

/*
def generateDamFile(fasta_file, workdir, *, dbsplit_opts, logger):
def getLocalAlignments(ref_file, reads_file=None, *, daligner_opts, ladump_opts,
def getMappings(ref_file, reads_file, *, damapper_opts, ladump_opts, logger):
*/

/**
    Holds a chain of local alignments that form a compound alignment. An AlignmentChain should
    contain at least one element.
*/
struct AlignmentChain
{
    struct LocalAlignment
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

    struct Contig
    {
        size_t id;
        size_t length;
    }

    enum Complement
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
}

private auto readAlignmentList(S)(in S lasDump) if (isSomeString!S)
{
    import std.algorithm : chunkBy, count, filter, joiner, map;
    import std.range : chunks;
    import std.string : lineSplitter;
    import std.typecons : Tuple;

    immutable recordSeparator = ';';
    immutable chainPartFormat = "P %d %d %c %c;C %d %d %d %d;D %d;L %d %d";
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
    alias byChainPartSplitter = lasDump => lasDump.lineSplitter.chunks(numChainPartLines);
    /// Parse chunks of numChainPartLines lines into RawChainPart s.
    alias parseChainPart = chainPartLines => {
        AlignmentChain.Contig contigA, contigB;
        char complement;
        char chainPartType;
        AlignmentChain.LocalAlignment localAlignment;

        auto joinedLines = chainPartLines.joiner(only(recordSeparator));
        // dfmt off
        int numMatches = joinedLines
            .save
            .array
            .formattedRead!chainPartFormat(
                contigA.id,
                contigB.id,
                complement,
                chainPartType,
                localAlignment.contigA.begin,
                localAlignment.contigA.end,
                localAlignment.contigB.begin,
                localAlignment.contigB.end,
                localAlignment.numDiffs,
                contigA.length,
                contigB.length,
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
        chainPartChunk.front()().contigA,
        chainPartChunk.front()().contigB,
        chainPartChunk.front()().complement,
        chainPartChunk
            .map!(chainPart => chainPart().localAlignment)
            .array);

    return byChainPartSplitter(lasDump)
        .map!parseChainPart
        .chunkBy!belongToSameChain
        .map!buildAlignmentChain;
    // dfmt on
}

unittest
{
    immutable testLasDump = q"EOF
        P 1 2 n >
        C 3 4 5 6
        D 7
        L 8 9
        P 1 2 n -
        C 12 13 14 15
        D 16
        L 17 18
        P 19 20 c +
        C 21 22 23 24
        D 25
        L 26 27
        P 19 20 c -
        C 30 31 32 33
        D 34
        L 35 36
        P 37 38 n .
        C 39 40 41 42
        D 43
        L 35 36
        P 46 47 c .
        C 48 49 50 51
        D 52
        L 53 54
        P 46 47 n .
        C 57 58 59 60
        D 61
        L 53 54
        P 64 65 c >
        C 66 67 68 69
        D 70
        L 71 72
        P 55 56 c -
        C 75 76 77 78
        D 79
        L 80 81
EOF".outdent;

    auto alignmentChains = readAlignmentList(testLasDump).array;
    with (AlignmentChain) with (Complement) with (LocalAlignment)
            {
                // dfmt off
            assert(alignmentChains == [
                AlignmentChain(
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

AlignmentContainer!(string[]) getLasFiles(string dbA)
{
    return getLasFiles(dbA, null);
}

AlignmentContainer!(string[]) getLasFiles(string dbA, string dbB)
{
    immutable fileTemplate = "%s.%s.las";
    auto aOptions = LasFileListOptions(dbA);
    auto bOptions = LasFileListOptions(dbB);
    size_t numLasFilesPerList = aOptions.numBlocks * bOptions.numBlocks;
    AlignmentContainer!(string[]) fileLists;

    fileLists.a2b.length = numLasFilesPerList;
    fileLists.b2a.length = numLasFilesPerList;
    foreach (size_t i; 0 .. aOptions.numBlocks)
    {
        foreach (size_t j; 0 .. bOptions.numBlocks)
        {
            fileLists.a2b[i * bOptions.numBlocks + j] = format!fileTemplate(
                    aOptions.fileNamePart(i), bOptions.fileNamePart(j),);
            fileLists.b2a[i * bOptions.numBlocks + j] = format!fileTemplate(
                    bOptions.fileNamePart(j), aOptions.fileNamePart(i),);
        }
    }

    return fileLists;
}

//string getDamFileName(fastaFile, workdir);
//void processLasFiles(refFile, readsFile, ladumpOpts, logger);

private
{
    struct LasFileListOptions
    {
        size_t numBlocks;
        string directory;
        string fileStem;

        this(string dbFile)
        {
            import std.path : baseName, dirName, stripExtension;

            this.numBlocks = getNumBlocks(dbFile);
            this.directory = dbFile.dirName;
            this.fileStem = dbFile.baseName.stripExtension;
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

    void dalign(string refDam, string readsDam, string[] dalignerOpts)
    {
        executeScript(chain(only("HPC.daligner"), dalignerOpts, only(refDam, readsDam),));
    }

    void damapper(string refDam, string readsDam, string[] damapperOpts)
    {
        executeScript(chain(only("HPC.damapper", "-C"), damapperOpts, only(refDam, readsDam),));
    }

    void fasta2dam(string inFile, string outFile)
    {
        executeCommand(only("fasta2DAM", outFile, inFile));
    }

    void dbsplit(string dbFile, string[] dbsplitOptions)
    {
        executeCommand(chain(only("DBsplit"), dbsplitOptions, only(dbFile),));
    }

    void lasort(string[] lasFiles)
    {
        executeCommand(chain(only("LAsort"), lasFiles,));
    }

    string ladump(string lasFile, string dbA, string dbB, string[] ladumpOpts)
    {
        return executeCommand(chain(only("LAdump"), ladumpOpts, only(dbA, dbB, lasFile),));
    }

    string dbshow(string dbFile, string contigId)
    {
        return executeCommand(only("DBshow", dbFile, contigId));
    }

    size_t getNumBlocks(string damFile)
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

    string executeCommand(Range)(Range command)
            if (isInputRange!(Unqual!Range) && isSomeString!(ElementType!(Unqual!Range)))
    {
        import std.process : Config, execute;

        string output = command.executeWrapper!(sCmd => execute(sCmd, null, // env
                Config.none, size_t.max, getWorkdir()));

        return output;
    }

    void executeScript(Range)(Range command)
            if (isInputRange!(Unqual!Range) && isSomeString!(ElementType!(Unqual!Range)))
    {
        import std.process : Config, executeShell;

        string output = command.executeWrapper!(sCmd => executeShell(sCmd.buildScriptLine,
                null, // env
                Config.none, size_t.max, getWorkdir()));

        logDiagnostic(output);
    }

    string executeWrapper(alias execCall, Range)(Range command)
            if (isInputRange!(Unqual!Range) && isSomeString!(ElementType!(Unqual!Range)))
    {
        import std.array : array;
        import std.algorithm : filter;

        auto sanitizedCommand = command.filter!"a != null".array;

        logDiagnostic("executing command: %s", sanitizedCommand);
        auto result = execCall(sanitizedCommand);

        if (result.status > 0)
        {
            throw new DazzlerCommandException(
                    format("process %s returned with non-zero exit code %d: %s",
                    sanitizedCommand[0], result.status, result.output));
        }

        return result.output;
    }

    string buildScriptLine(string[] command)
    {
        import std.process : escapeShellCommand;

        return escapeShellCommand(command) ~ " | sh -sv";
    }
}
