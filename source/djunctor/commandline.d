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
import djunctor.log;

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
        options = parseArgs!Options(args[1 .. $]);
    }
    catch (ArgParseError e)
    {
        logError(e.msg);
        writeln(usage);

        return ReturnCode.commandlineError;
    }
    catch (ArgParseHelp e)
    {
        // Help was requested
        writeln(usage);
        write(help);

        return ReturnCode.ok;
    }

    initLogger(options);
    addDefaultOptions(options);
    if (!verifyOptions(options))
    {
        // some options are invalid
        return ReturnCode.commandlineError;
    }
    createWorkDir(options);

    const Options finalOptions = options;

    logInfo(finalOptions.to!string);

    try
    {
        import djunctor.djunctor : runWithOptions;

        runWithOptions(finalOptions);

        return ReturnCode.ok;
    }
    catch (Exception e)
    {
        writeln(e);

        return ReturnCode.djunctorError;
    }
    finally
    {
        if (!finalOptions.keepTemp)
        {
            cleanWorkDir(finalOptions);
        }
    }
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

    @Argument("READS")
    @Help("set of PacBio reads in .dam format")
    string readsFile;

    @Option("daligner-options", "dao")
    @Help("list of options to pass to `daligner`")
    string[] dalignerOptions = [];

    @Option("dmapper-options", "dmo")
    @Help("list of options to pass to `damapper`")
    string[] damapperOptions = [];

    @Option("dbsplit-options", "dmo")
    @Help("list of options to pass to `DBsplit`")
    string[] dbsplitOptions = [];

    /// List of options to pass to `LAdump`
    @Option()
    string[] ladumpOptions = ["-c", // output alignment coordinates
        "-d", // output number of differences for each local alignment
        "-l", // output lengths of the contigs
        ];

    @Option("confidence", "c")
    @Help("discard pile ups with <ulong>% confidence if too large/small")
    size_t confidence = 95;

    @Option("keep-temp")
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
        import std.traits : getUDAs;

        const(string)[] options = [];
        const(string)[] values = [];

        static foreach (member; __traits(allMembers, Options))
        {
            static if (!is(typeof(__traits(getMember, this, member)) == function))
            {
                static if (getUDAs!(__traits(getMember, this, member), Option)
                        .length || getUDAs!(__traits(getMember, this, member), Argument).length)
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

    void addDefaultOptions(ref Options options) nothrow
    {
        import std.format : format;

        if (!options.dalignerOptions || options.dalignerOptions.length == 0)
        {
            options.dalignerOptions = ["-s126", // prevent integer overflow in trace point procedure
                "-t20", // ignore k-mers occuring more than this number
                ];
        }

        if (!options.damapperOptions || options.damapperOptions.length == 0)
        {
            options.damapperOptions = options.dalignerOptions;
        }
    }

    bool verifyOptions(ref Options options)
    {
        if (!verifyInputFiles(options))
            return false;

        return true;
    }

    bool verifyInputFiles(ref Options options)
    {
        import std.algorithm : endsWith;
        import std.file : exists;
        import djunctor.dazzler : getHiddenDbFiles;

        foreach (inputFile; [options.refFile, options.readsFile])
        {
            if (!inputFile.endsWith(".dam"))
            {
                logError("expected .dam file, got `%s`", inputFile);

                return false;
            }

            if (!inputFile.exists)
            {
                logError("cannot open file `%s`", inputFile);

                return false;
            }

            foreach (hiddenDbFile; getHiddenDbFiles(inputFile))
            {
                if (!hiddenDbFile.exists)
                {
                    logError("cannot open hidden database file `%s`", hiddenDbFile);

                    return false;
                }
            }
        }

        return true;
    }

    void createWorkDir(ref Options options)
    {
        import std.algorithm : endsWith;
        import std.exception : ErrnoException;
        import std.file : tempDir;
        import std.path : buildPath;
        import std.string : fromStringz, toStringz;

        version (Posix)
        {
            import core.sys.posix.stdlib : mkdtemp;

            char[255] workdirNameBuffer;
            auto workdirTemplate = buildPath(tempDir(), options.workdirTemplate);
            auto len = workdirTemplate.length;
            assert(len < workdirNameBuffer.length);
            assert(workdirTemplate.endsWith("XXXXXX"));

            workdirNameBuffer[0 .. len] = workdirTemplate[];
            workdirNameBuffer[len] = 0;

            if (null == mkdtemp(workdirNameBuffer.ptr))
            {
                throw new ErrnoException("cannot create workdir", __FILE__, __LINE__);
            }

            options.workdir = to!string(fromStringz(workdirNameBuffer.ptr));
        }
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
