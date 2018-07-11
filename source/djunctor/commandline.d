/**
    Defines the behavior of the djunctor command line client.

    Copyright: © 2018 Arne Ludwig <arne.ludwig@posteo.de>
    License: Subject to the terms of the MIT license, as written in the
             included LICENSE file.
    Authors: Arne Ludwig <arne.ludwig@posteo.de>
*/
module djunctor.commandline;

import darg : ArgParseHelp, Argument, Help, helpString, Option, OptionFlag,
    parseArgs, usageString;
import djunctor.dazzler : provideDamFileInWorkdir, provideLasFileInWorkdir,
    ProvideMethod, DaccordOptions, DalignerOptions, DamapperOptions,
    LAdumpOptions;
import djunctor.scaffold : JoinPolicy;
import djunctor.util.log;
import std.conv;
import std.stdio;
import std.meta : Instantiate;
import std.range.primitives : ElementType, isForwardRange;
import std.traits : hasMember, isSomeString;
import vibe.data.json : serializeToJsonString;

/// Possible returns codes of the command line execution.
enum ReturnCode
{
    ok,
    commandlineError,
    djunctorError,
}

/// Start `djunctor` with the given set of arguments.
ReturnCode runDjunctorCommandline(string[] args)
{
    Options options;

    try
    {
        options = processOptions(args);
    }
    catch (ArgParseHelp e)
    {
        // Help was requested
        writeln(usage);
        write(help);

        return ReturnCode.ok;
    }
    catch (Exception e)
    {
        logError(e.msg);
        writeln(usage);

        return ReturnCode.commandlineError;
    }

    const Options finalOptions = options;
    logInfo(finalOptions.serializeToJsonString());

    scope (exit)
    {
        if (!(finalOptions.keepTemp || options.generateCLIOptions))
        {
            cleanWorkDir(finalOptions);
        }
    }

    try
    {
        if (options.generateCLIOptions)
        {
            printCLIOptionsForAlignment(finalOptions);
        }
        else
        {
            import djunctor.djunctor : runWithOptions;

            runWithOptions(finalOptions);
        }

        return ReturnCode.ok;
    }
    catch (Exception e)
    {
        logError(e.to!string);

        return ReturnCode.djunctorError;
    }
}

Options processOptions(string[] args)
{
    Options options = parseArgs!Options(args[1 .. $]);

    initLogger(options);
    normalizePaths(options);
    addDefaultOptions(options);
    verifyOptions(options);

    if (options.generateCLIOptions)
    {
        // Do not do any "real" work.
        return options;
    }

    createWorkDir(options);
    options.refDb = getRefDb(options);
    options.readsDb = getReadsDb(options);
    options.selfAlignmentFile = provideLasFileInWorkdir(options.selfAlignmentInputFile, options.provideMethod, options.workdir);
    options.refVsReadsAlignmentFile = provideLasFileInWorkdir(options.refVsReadsAlignmentInputFile, options.provideMethod, options.workdir);

    return options;
}

void printCLIOptionsForAlignment(in ref Options options)
{
    import std.string : join;

    // dfmt off
    writeln("# self alignment options (consider using `HPC.daligner`)");
    writeln((
        ["daligner"] ~
        options.selfAlignmentOptions ~
        [options.refFile, options.refFile]
    ).join(" "));
    writeln("# ref vs reads alignment options (consider using `HPC.damapper`)");
    writeln((
        ["damapper"] ~
        options.refVsReadsAlignmentOptions ~
        [options.refFile, options.readsFile]
    ).join(" "));
    // dfmt on
}

/**
    The set of options used by `djunctor`. For details see the source code or
    run `./djunctor -h`.
*/
struct Options
{
    @Option("help", "h")
    @Help("Prints this help.")
    OptionFlag help;

    @Argument("REFERENCE")
    @Help("reference assembly in .dam format")
    string refFile;
    string refDb;

    @Argument("READS")
    @Help("set of PacBio reads in .dam format")
    string readsFile;
    string readsDb;

    @Argument("SELF-ALIGNMENT")
    @Help(q"{
        local alignments of the reference against itself in form of a .las
        file as produced by `daligner`
    }")
    string selfAlignmentInputFile;
    @Option()
    string selfAlignmentFile;

    @Argument("REF-VS-READS-ALIGNMENT")
    @Help(q"{
        alignments chains of the reads against the reference in form of a .las
        file as produced by `damapper`
    }")
    string refVsReadsAlignmentInputFile;
    @Option()
    string refVsReadsAlignmentFile;

    @Option("generate-options")
    @Help("print recommended CLI options for the alignments (daligner/damapper) and exit")
    OptionFlag generateCLIOptions;

    @Option("block")
    @Help("if given only block <ulong> of the reference DB will be processed")
    size_t refBlockNum;

    @Option("reads")
    @Help(q"{
        restrict the set of reads to this list of IDs given in a file as
        produced by --unused-reads
    }")
    string readsListFile = null;
    @Option()
    size_t[] readsList;

    @Option("reads-error")
    @Help("estimated error rate in reads")
    double readsErrorRate = .15;

    @Option("reference-error")
    @Help("estimated error rate in reference")
    double referenceErrorRate = .01;

    @Option("min-anchor-length")
    @Help("alignment need to have at least this length of unique anchoring sequence")
    size_t minAnchorLength = 200;

    @Option("good-anchor-length")
    @Help("alignment anchors with at least this length will get no penalty")
    size_t goodAnchorLength = 1000;

    @Option("min-reads-per-pile-up")
    @Help("alignment anchors with at least this length will get no penalty")
    size_t minReadsPerPileUp = 5;

    @Option("min-extension-length")
    @Help("extensions must have at least <ulong> bps of consensus to be inserted")
    size_t minExtensionLength = 100;

    @Option("out-mask")
    @Help(q"{
        write inferred repeat mask into a Dazzler mask. Given a path-like
        string without extension: the `dirname` designates the directory to
        write the mask to. The mask comprises two hidden files
        `.[REFERENCE].[MASK].{anno,data}`.
    }")
    string outMask = null;

    @Option("in-mask")
    @Help("use Dazzler mask to initialize the repeat mask (see also `--out-mask`)")
    string inMask = null;

    @Option("unused-reads")
    @Help("if given write unused read IDs to the designated file as JSON array")
    string unusedReadsList = null;

    @Option("join-policy")
    @Help(q"{
        allow only joins (gap filling) in the given mode:
        `scaffoldGaps` (only join gaps inside of scaffolds –
        marked by `n`s in FASTA),
        `scaffolds` (join gaps inside of scaffolds and try to join scaffolds),
        `contigs` (break input into contigs and re-scaffold everything;
        maintains scaffold gaps where new scaffolds are consistent)
    }")
    JoinPolicy joinPolicy = JoinPolicy.scaffoldGaps;

    @Option("extend-contigs")
    @Help("if given extend contigs even if no spanning reads can be found")
    OptionFlag shouldExtendContigs;

    /// List of options to pass to `daligner`.
    @Option()
    string[] dalignerOptions = [];
    string[] selfAlignmentOptions;
    string[] pileUpAlignmentOptions;

    /// List of options to pass to `damapper`.
    @Option()
    // dfmt off
    string[] damapperOptions = [
        DamapperOptions.bestMatches ~ ".7",
    ];
    string[] refVsReadsAlignmentOptions;
    // dfmt on

    /// List of options to pass to `daccord`.
    @Option()
    // dfmt off
    string[] daccordOptions = [
        DaccordOptions.produceFullSequences,
    ];
    // dfmt on

    /// List of options to pass to `DBsplit`.
    @Option()
    string[] dbsplitOptions = [];

    @Option("fasta-line-width", "w")
    @Help("list of options to pass to `DBsplit`")
    size_t fastaLineWidth = 50;

    /// List of options to pass to `LAdump` for simple dump.
    @Option()
    // dfmt off
    string[] ladumpOptions = [
        LAdumpOptions.coordinates,
        LAdumpOptions.numDiffs,
        LAdumpOptions.lengths,
    ];
    // dfmt on

    /// List of options to pass to `LAdump` for trace point dump.
    @Option()
    // dfmt off
    string[] ladumpTraceOptions = [
        LAdumpOptions.coordinates,
        LAdumpOptions.tracePoints,
    ];
    // dfmt on

    @Option("confidence", "c")
    @Help("mask region were coverage is out of confidence interval with <double> confidence")
    double confidence = .95;

    @Option("input-provide-method", "p")
    @Help(q"{
        use <ProvideMethod> to provide the input files in the working directory;
        either `symlink` or `copy` (default: `symlink`)
    }")
    ProvideMethod provideMethod = ProvideMethod.symlink;

    @Option("keep-temp", "k")
    @Help("keep the temporary files; outputs the exact location")
    OptionFlag keepTemp;

    /**
        Last part of the working directory name. A directory in the temp
        directory as returned by `std.file.tmpDir` with the naming scheme will
        be created to hold all data for the computation.
    */
    string workdirTemplate = "djunctor-XXXXXX";

    /// This is a temporary directory to store all working data.
    @Option()
    string workdir;

    @Option()
    size_t verbose = 0;
    @Option("verbose", "v")
    @Help("increase output to help identify problems; use up to three times")
    void increaseVerbosity() pure
    {
        ++verbose;
    }

    /// Return a table that lists all `@Option`s in this struct.
    string toString() const
    {
        import std.algorithm : joiner, map, maxElement;
        import std.array : array;
        import std.format : format;
        import std.range : zip;
        import std.traits : hasUDA;

        const(string)[] options = [];
        const(string)[] values = [];

        static foreach (member; __traits(allMembers, Options))
        {
            static if (!is(typeof(__traits(getMember, this, member)) == function))
            {
                static if (hasUDA!(__traits(getMember, this, member), Option)
                        || hasUDA!(__traits(getMember, this, member), Argument))
                {
                    options ~= member;
                    values ~= to!string(__traits(getMember, this, member));
                }
            }
        }

        const size_t width = options.map!"a.length".maxElement;

        // dfmt off
        return zip(options, values)
            .map!(entry => format!"%*s: %s\n"(width, entry[0], entry[1]))
            .joiner(" ")
            .array
            .to!string;
        // dfmt on
    }
}

private
{
    immutable usage = usageString!Options("djunctor");
    immutable help = helpString!Options;

    void initLogger(ref Options options) nothrow
    {
        switch (options.verbose)
        {
        case 3:
            setLogLevel(LogLevel.debug_);
            break;
        case 2:
            setLogLevel(LogLevel.diagnostic);
            break;
        case 1:
            setLogLevel(LogLevel.info);
            break;
        case 0:
        default:
            setLogLevel(LogLevel.error);
            break;
        }
    }

    void normalizePaths(ref Options options)
    {
        import std.path : absolutePath;

        // dfmt off
        immutable fileOptions = [
            "refFile",
            "readsFile",
            "selfAlignmentFile",
            "refVsReadsAlignmentFile",
            "readsListFile",
            "outMask",
            "inMask",
            "unusedReadsList",
        ];
        // dfmt on

        static foreach (fileOption; fileOptions)
        {
            mixin("options." ~ fileOption ~ " = absolutePath(options." ~ fileOption ~ ");");
        }
    }

    void addDefaultOptions(ref Options options)
    {
        import std.format : format;

        // dfmt off
        options.selfAlignmentOptions = options.dalignerOptions ~ [
            DalignerOptions.identity,
            format!(DalignerOptions.minAlignmentLength ~ "%d")(options.minAnchorLength),
            format!(DalignerOptions.averageCorrelationRate ~ "%f")((1 - options.referenceErrorRate)^^2),
        ];

        options.refVsReadsAlignmentOptions = options.damapperOptions ~ [
            DamapperOptions.symmetric,
            format!(DamapperOptions.averageCorrelationRate ~ "%f")((1 - options.referenceErrorRate) * (1 - options.readsErrorRate)),
        ];

        options.pileUpAlignmentOptions = options.dalignerOptions ~ [
            DalignerOptions.identity,
            format!(DalignerOptions.minAlignmentLength ~ "%d")(options.minAnchorLength),
            format!(DalignerOptions.averageCorrelationRate ~ "%f")((1 - options.readsErrorRate)^^2),
        ];
        // dfmt on
    }

    void verifyOptions(ref Options options)
    {
        verifyInputFiles(options);
        verifyOutputFiles(options);
    }

    void verifyInputFiles(ref Options options)
    {
        verifyDamFile(options.refFile);
        verifyDamFile(options.readsFile);
        options.readsList = verifyReadsListFile(options.readsListFile);
        verifyInMask(options.refFile, options.inMask);
    }

    void verifyDamFile(in string damFile, in size_t blockNum = 0)
    {
        import djunctor.dazzler : getHiddenDbFiles, getNumBlocks;
        import std.algorithm : endsWith;
        import std.exception : enforce;
        import std.file : exists;
        import std.format : format;

        enforce!Exception(damFile.endsWith(".dam"), format!"expected .dam file, got `%s`"(damFile));
        enforce!Exception(damFile.exists, format!"cannot open file `%s`"(damFile));

        foreach (hiddenDbFile; getHiddenDbFiles(damFile))
        {
            enforce!Exception(hiddenDbFile.exists,
                    format!"cannot open hidden database file `%s`"(hiddenDbFile));
        }

        size_t numBlocks = getNumBlocks(damFile);
        enforce!Exception(blockNum == 0 || blockNum <= numBlocks,
                format!"cannot select block %d; databse has only %d blocks"(blockNum, numBlocks));
    }

    size_t[] verifyReadsListFile(ref string readsListFile)
    {
        import std.algorithm : endsWith, joiner, map;
        import std.array : array;
        import std.exception : enforce;
        import std.file : exists, readText;
        import std.format : format;
        import std.range : tee;
        import vibe.data.json : parseJsonString;

        if (readsListFile is null)
        {
            return [];
        }

        enforce!Exception(readsListFile.exists, format!"cannot open file `%s`"(readsListFile));
        auto inputData = readText(readsListFile).parseJsonString(readsListFile);
        // dfmt off
        return inputData[]
            .map!"a.get!long"
            .tee!(n => enforce!Exception(n >= 1, format!"unexpected ID %d must be >= 1"(n)))
            .array
            .to!(size_t[]);
        // dfmt on
    }

    void verifyInMask(in string refFile, in string maskDestination)
    {
        import djunctor.dazzler : getMaskFiles;
        import std.algorithm : endsWith;
        import std.exception : enforce;
        import std.file : exists;
        import std.format : format;

        if (maskDestination is null)
        {
            return;
        }

        foreach (maskFile; getMaskFiles(refFile, maskDestination))
        {
            enforce!Exception(maskFile.exists, format!"cannot open file `%s`"(maskFile));
        }
    }

    void verifyOutputFiles(ref Options options)
    {
        import djunctor.dazzler : getMaskFiles;

        enforceCanWriteIfPresent(options.unusedReadsList);

        if (options.outMask !is null)
        {
            foreach (maskFile; getMaskFiles(options.refFile, options.outMask))
            {
                enforceCanWriteIfPresent(maskFile);
            }
        }
    }

    void createWorkDir(ref Options options)
    {
        import std.file : tempDir;
        import std.path : buildPath;
        import djunctor.util.tempfile : mkdtemp;

        auto workdirTemplate = buildPath(tempDir(), options.workdirTemplate);

        options.workdir = mkdtemp(workdirTemplate);
    }

    void cleanWorkDir(in ref Options options)
    {
        import std.file : rmdirRecurse;

        try
        {
            rmdirRecurse(options.workdir);
        }
        catch (Exception e)
        {
            logWarn(to!string(e));
        }
    }

    string getRefDb(in Options options)
    {
        import std.path : extension, setExtension;

        string dbFile = provideDamFileInWorkdir(options.refFile,
                options.provideMethod, options.workdir);

        if (options.refBlockNum > 0)
        {
            dbFile = dbFile.setExtension(options.refBlockNum.to!string ~ dbFile.extension);
        }

        return dbFile;
    }

    string getReadsDb(in Options options)
    {
        import djunctor.dazzler : dbSubset;

        if (options.readsList.length == 0)
        {
            return provideDamFileInWorkdir(options.readsFile,
                    options.provideMethod, options.workdir);
        }
        else
        {
            return dbSubset(options.readsFile, options.readsList[], options);
        }
    }

    void enforceCanWriteIfPresent(string fileName)
    {
        if (fileName is null)
        {
            return;
        }

        import std.exception : ErrnoException;
        import std.file : exists, FileException, remove;
        import std.format : format;
        import std.stdio : File;

        auto deleteAfterwards = !fileName.exists;

        try
        {
            auto unusedReadsList = File(fileName, "a");
        }
        catch (ErrnoException e)
        {
            throw new Exception(format!"cannot open file `%s` for writing: %s"(fileName, e));
        }

        if (deleteAfterwards)
        {
            try
            {
                remove(fileName);
            }
            catch (FileException e)
            {
                // dfmt off
                logJsonDebug(
                    "info", "failed to delete file after testing",
                    "error", e.toString(),
                    "file", fileName,
                );
                // dfmt on
            }
        }
    }
}

template hasOption(T, string name, alias S)
{
    static if (hasMember!(T, name))
    {
        enum hasOption = Instantiate!(S, typeof(__traits(getMember, T, name)));
    }
    else
    {
        enum hasOption = false;
    }
}

template hasOption(T, string name, S)
{
    static if (hasMember!(T, name))
    {
        enum hasOption = is(typeof(__traits(getMember, T, name)) == S);
    }
    else
    {
        enum hasOption = false;
    }
}

unittest
{
    struct Mock
    {
        string a;
        string b;
        int c;
    }

    static assert(hasOption!(Mock, "a", string));
    static assert(hasOption!(Mock, "a", isSomeString));
    static assert(!hasOption!(Mock, "c", string));
    static assert(!hasOption!(Mock, "c", isSomeString));
    static assert(!hasOption!(Mock, "d", string));
    static assert(!hasOption!(Mock, "d", isSomeString));
}

template isOptionsList(T)
{
    enum isOptionsList = isForwardRange!T && isSomeString!(ElementType!T);
}

unittest
{
    static assert(isOptionsList!(string[]));
    static assert(!isOptionsList!(string));
    static assert(!isOptionsList!(int[]));
    static assert(!isOptionsList!(int));
}
