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
import std.algorithm : equal, map, joiner, sort, splitter, SwapStrategy;
import std.array : array;
import std.conv : to;
import std.format : format, formattedRead;
import std.meta : Instantiate;
import std.range : chain, only;
import std.range.primitives : ElementType, isInputRange;
import std.stdio : File, writeln;
import std.string : outdent;
import std.traits : hasMember, isSomeString, Unqual;
import std.typecons : Flag, No, tuple, Tuple, Yes;

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

/**
    Provide dbFile in workdir.

    Returns: Workdir location of the dbFile.
*/
string provideDamFileInWorkdir(in string dbFile, ProvideMethod provideMethod)
{
    import std.file : copy, symlink;
    import std.path : baseName, buildPath;
    import std.range : chain, only;

    alias inWorkdir = anyDbFile => buildPath(getWorkdir(), anyDbFile.baseName);
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
        if (hasOption!(Options, "dalignerOptions", isOptionsList)
            && hasOption!(Options, "ladumpOptions", isOptionsList))
{
    return getLocalAlignments(dbA, null, options);
}

AlignmentContainer!(AlignmentChain[]) getLocalAlignments(Options)(in string dbA,
        in string dbB, in Options options)
        if (hasOption!(Options, "dalignerOptions", isOptionsList)
            && hasOption!(Options, "ladumpOptions", isOptionsList))
{
    dalign(dbA, dbB, options.dalignerOptions);

    return processGeneratedLasFiles(dbA, dbB, options);
}

AlignmentContainer!(AlignmentChain[]) getMappings(Options)(in string dbA,
        in string dbB, in Options options)
        if (hasOption!(Options, "damapperOptions", isOptionsList)
            && hasOption!(Options, "ladumpOptions", isOptionsList))
{
    damapper(dbA, dbB, options.damapperOptions);

    return processGeneratedLasFiles(dbA, dbB, options);
}

private auto processGeneratedLasFiles(Options)(in string dbA, in string dbB, in Options options)
        if (hasOption!(Options, "ladumpOptions", isOptionsList))
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
                    .map!(lasFile => ladump(lasFile, dbFiles.expand, options.ladumpOptions))
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
    import std.range : chunks, drop;
    import std.string : lineSplitter;

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
            .save
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
    import std.algorithm : max;

    static struct InferredParams
    {
        size_t numBlocks;
        string directory;
        string fileStem;

        this(string dbFile)
        {
            import std.path : baseName, dirName, stripExtension;

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

    immutable singleFileTemplate = "%s.las";
    immutable twoFileTemplate = "%s.%s.las";
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
                        aOptions.fileNamePart(i), bOptions.fileNamePart(j));
                fileLists.b2a[i * bOptions.numBlocks + j] = format!twoFileTemplate(
                        bOptions.fileNamePart(j), aOptions.fileNamePart(i));
            }
        }
        else
        {
            fileLists.a2b[i] = fileLists.b2a[i] = format!singleFileTemplate(
                    aOptions.fileNamePart(i));
        }
    }

    return fileLists;
}

private
{
    void dalign(in string refDam, in string[] dalignerOpts)
    {
        dalign(refDam, null, dalignerOpts);
    }

    void dalign(in string refDam, in string readsDam, in string[] dalignerOpts)
    {
        executeScript(chain(only("HPC.daligner"), dalignerOpts, only(refDam), only(readsDam)));
    }

    void damapper(in string refDam, in string readsDam, in string[] damapperOpts)
    {
        executeScript(chain(only("HPC.damapper", "-C"), damapperOpts,
                only(refDam), only(readsDam)));
    }

    void fasta2dam(in string inFile, in string outFile)
    {
        executeCommand(only("fasta2DAM", outFile, inFile));
    }

    void dbsplit(in string dbFile, in string[] dbsplitOptions)
    {
        executeCommand(chain(only("DBsplit"), dbsplitOptions, only(dbFile)));
    }

    void lasort(in string[] lasFiles)
    {
        import std.file : rename;
        import std.path : setExtension;

        alias sortedName = (lasFile) => lasFile.setExtension(".S.las");

        executeCommand(chain(only("LAsort"), lasFiles));
        foreach (lasFile; lasFiles)
        {
            rename(sortedName(lasFile), lasFile);
        }
    }

    string ladump(in string lasFile, in string dbA, in string dbB, in string[] ladumpOpts)
    {
        return executeCommand(chain(only("LAdump"), ladumpOpts, only(dbA),
                only(dbB), only(lasFile)));
    }

    string dbshow(in string dbFile, in string contigId)
    {
        return executeCommand(only("DBshow", dbFile, contigId));
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

    string executeCommand(Range)(in Range command)
            if (isInputRange!(Unqual!Range) && isSomeString!(ElementType!(Unqual!Range)))
    {
        import std.process : Config, execute;

        string output = command.executeWrapper!("command",
                sCmd => execute(sCmd, null, // env
                    Config.none, size_t.max, getWorkdir()));

        return output;
    }

    void executeScript(Range)(in Range command)
            if (isInputRange!(Unqual!Range) && isSomeString!(ElementType!(Unqual!Range)))
    {
        import std.process : Config, executeShell;

        string output = command.executeWrapper!("script",
                sCmd => executeShell(sCmd.buildScriptLine, null, // env
                    Config.none,
                    size_t.max, getWorkdir()));

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
        import std.process : escapeShellCommand;

        return escapeShellCommand(command) ~ " | sh -sv";
    }
}
