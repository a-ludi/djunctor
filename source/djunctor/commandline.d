/**
    Defines the behavior of the djunctor command line client.

    Copyright: © 2018 Arne Ludwig <arne.ludwig@posteo.de>
    License: Subject to the terms of the MIT license, as written in the
             included LICENSE file.
    Authors: Arne Ludwig <arne.ludwig@posteo.de>
*/
module djunctor.commandline;

import darg : ArgParseError, ArgParseHelp, Argument, Help, helpString, Option,
    OptionFlag, parseArgs, usageString;
import std.conv;
import std.stdio;
import djunctor.dazzler : provideDamFileInWorkdir, ProvideMethod;
import djunctor.util.log;
import std.meta : Instantiate;
import std.range.primitives : ElementType, isForwardRange;
import std.traits : hasMember, isSomeString;
import vibe.data.json : serializeToJsonString;

/// Possible returns codes of the commandline execution.
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
        if (!finalOptions.keepTemp)
        {
            cleanWorkDir(finalOptions);
        }
    }

    try
    {
        import djunctor.djunctor : runWithOptions;

        runWithOptions(finalOptions);

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
    createWorkDir(options);
    options.refDb = provideDamFileInWorkdir(options.refFile, options.provideMethod, options.workdir);
    options.readsDb = provideDamFileInWorkdir(options.readsFile,
            options.provideMethod, options.workdir);

    return options;
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

    @Option("coord-transform", "T")
    @Help("write a Python 2.7 compatible script which transform input coordinates to output coordinates")
    string coordTransform = null;

    @Option("reference-error")
    @Help("estimated error rate in reference")
    double referenceErrorRate = .95;

    @Option("min-anchor-length")
    @Help("alignment need to have at least this length of unique anchoring sequence")
    size_t minAnchorLength = 200;

    @Option("good-anchor-length")
    @Help("alignment anchors with at least this length will get no penalty")
    size_t goodAnchorLength = 1000;

    @Option("min-absolute-pile-up-size")
    @Help("alignment anchors with at least this length will get no penalty")
    size_t minAbsolutePileUpSize = 5;

    @Option("daligner-options")
    @Help("list of options to pass to `daligner`")
    string[] dalignerOptions = [];

    @Option("dmapper-options")
    @Help("list of options to pass to `damapper`")
    string[] damapperOptions = [];

    @Option("daccord-options")
    @Help("list of options to pass to `daccord`")
    string[] daccordOptions = [];

    @Option("dbsplit-options")
    @Help("list of options to pass to `DBsplit`")
    string[] dbsplitOptions = [];

    @Option("fasta-line-width", "w")
    @Help("list of options to pass to `DBsplit`")
    size_t fastaLineWidth = 50;

    /// List of options to pass to `LAdump` for simple dump.
    @Option()
    // dfmt off
    string[] ladumpOptions = [
        "-c", // output alignment coordinates
        "-d", // output number of differences for each local alignment
        "-l", // output lengths of the contigs
        "-o", // output proper overlaps only
    ];
    // dfmt on

    /// List of options to pass to `LAdump` for trace point dump.
    @Option()
    // dfmt off
    string[] ladumpTraceOptions = [
        "-c", // output alignment coordinates
        "-t", // output number of differences for each local alignment
        "-o", // output proper overlaps only
    ];
    // dfmt on

    /// List of options to pass to `DBdump`
    @Option()
    // dfmt off
    string[] dbdumpOptions = [
        "-r", // read number
        "-h", // original file name string (header) and location: well, pulse start, pulse end
        "-s", // sequence string
    ];
    // dfmt on

    @Option("confidence", "c")
    @Help("discard pile ups with <ulong>% confidence if too large/small")
    size_t confidence = 95;

    @Option("input-provide-method", "p")
    @Help("use this method to provide the input files in the working directory")
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
    string toString() pure const
    {
        import std.algorithm : joiner, map, maxElement, sort;
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

        options.refFile = absolutePath(options.refFile);
        options.readsFile = absolutePath(options.readsFile);
    }

    void addDefaultOptions(ref Options options) nothrow
    {
        import std.format : format;

        if (!options.dalignerOptions || options.dalignerOptions.length == 0)
        {
            // dfmt off
            options.dalignerOptions = [
            ];
            // dfmt on
        }

        if (!options.damapperOptions || options.damapperOptions.length == 0)
        {
            // dfmt off
            options.damapperOptions = options.dalignerOptions ~ [
            ];
            // dfmt on
        }
    }

    void verifyOptions(ref Options options)
    {
        verifyInputFiles(options);
        verifyOutputFiles(options);
    }

    void verifyInputFiles(ref Options options)
    {
        import std.algorithm : endsWith;
        import std.file : exists;
        import std.format : format;
        import djunctor.dazzler : getHiddenDbFiles;

        foreach (inputFile; [options.refFile, options.readsFile])
        {
            if (!inputFile.endsWith(".dam"))
            {
                throw new Exception(format!"expected .dam file, got `%s`"(inputFile));
            }

            if (!inputFile.exists)
            {
                throw new Exception(format!"cannot open file `%s`"(inputFile));
            }

            foreach (hiddenDbFile; getHiddenDbFiles(inputFile))
            {
                if (!hiddenDbFile.exists)
                {
                    throw new Exception(
                            format!"cannot open hidden database file `%s`"(hiddenDbFile));
                }
            }
        }
    }

    void verifyOutputFiles(ref Options options)
    {
        import std.exception : ErrnoException;
        import std.format : format;
        import std.stdio : File;

        if (!(options.coordTransform is null))
        {
            try
            {
                auto coordTransformScript = File(options.coordTransform, "a");
            }
            catch (ErrnoException e)
            {
                throw new Exception(format!"cannot write coord transform file `%s`: %s"(
                        options.coordTransform, e));
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
