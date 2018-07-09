/**
    Defines the behavior of the djunctor command line client.

    Copyright: © 2018 Arne Ludwig <arne.ludwig@posteo.de>
    License: Subject to the terms of the MIT license, as written in the
             included LICENSE file.
    Authors: Arne Ludwig <arne.ludwig@posteo.de>
*/
module testgen.commandline;

import darg : ArgParseHelp, Argument, Help, helpString, Multiplicity, Option, OptionFlag,
    parseArgs, usageString;
import djunctor.dazzler : provideDamFileInWorkdir, provideLasFileInWorkdir,
    ProvideMethod, DaccordOptions, DalignerOptions, DamapperOptions,
    LAdumpOptions;
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
    testgenError,
}

/// Start `testgen` with the given set of arguments.
ReturnCode runTestGenCommandline(string[] args)
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
        final switch (options.action)
        {
        case Action.translocate:
            import testgen.translocator : runWithOptions;

            runWithOptions(finalOptions);
            break;
        }

        return ReturnCode.ok;
    }
    catch (Exception e)
    {
        logError(e.to!string);

        return ReturnCode.testgenError;
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
    options.assembly1Db = getDb(options.assembly1File, options);
    options.assembly2Db = getDb(options.assembly2File, options);
    provideLasFileInWorkdir(options.alignmentFile, options.provideMethod, options.workdir);

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

    @Argument("ACTION")
    @Help(q"{
        execute the given ACTION: `translocate` assembly gaps from one
        assembly to another printing the result in FASTA on stdout.
    }")
    Action action;

    @Argument("ASSEMBLY1")
    @Help("assembly1 in .dam format")
    string assembly1File;
    string assembly1Db;

    @Argument("ASSEMBLY2")
    @Help("assembly2 in .dam format")
    string assembly2File;
    string assembly2Db;

    @Argument("ALIGNMENT")
    @Help("alignment of ASSEMBLY1 vs. ASSEMBLY2")
    string alignmentFile;

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

enum Action
{
    translocate,
}

private
{
    immutable usage = usageString!Options("testgen");
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
            "assembly1File",
            "assembly2File",
            "alignmentFile",
        ];
        // dfmt on

        static foreach (fileOption; fileOptions)
        {
            mixin("options." ~ fileOption ~ " = absolutePath(options." ~ fileOption ~ ");");
        }
    }

    void addDefaultOptions(ref Options options)
    {
    }

    void verifyOptions(ref Options options)
    {
        verifyInputFiles(options);
    }

    void verifyInputFiles(ref Options options)
    {
        verifyDamFile(options.assembly1File);
        verifyDamFile(options.assembly2File);
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

    string getDb(in string dbFile, in Options options, in size_t blockNum = 0)
    {
        import std.path : extension, setExtension;

        string workdirDbFile = provideDamFileInWorkdir(dbFile,
                options.provideMethod, options.workdir);

        if (blockNum > 0)
        {
            workdirDbFile = workdirDbFile.setExtension(blockNum.to!string ~ workdirDbFile.extension);
        }

        return workdirDbFile;
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
