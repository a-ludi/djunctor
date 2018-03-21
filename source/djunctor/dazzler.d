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
import std.conv;
import std.range : chain, only;
import std.range.primitives : ElementType, isInputRange;
import std.stdio : File, writeln;
import std.traits : Unqual, isSomeString;

/// File suffixes of hidden DB files.
private immutable hiddenDbFileSuffixes = [".bps", ".hdr", ".idx"];

/**
    Return a list o hidden files associated to every `.dam`/`.db` file. These
    files contain the actual data used in all the computation. Thus, we
    carefully check for their existence.
*/
auto getHiddenDbFiles(string dbFile)
{
    import std.algorithm : map;
    import std.path : baseName, chainPath, dirName, withExtension;

    return hiddenDbFileSuffixes.map!(delegate(suffix) {
        return chainPath(dbFile.dirName, "." ~ dbFile.baseName.withExtension(suffix).to!string,);
    });
}

private string workdir;

void setWorkdir(string workdir_)
{
    workdir = workdir_;
}

string getWorkdir()
{
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

/*
def generateDamFile(fasta_file, workdir, *, dbsplit_opts, logger):
def getLocalAlignments(ref_file, reads_file=None, *, daligner_opts, ladump_opts,
def getMappings(ref_file, reads_file, *, damapper_opts, ladump_opts, logger):
class AlignmentListReader:
*/

AlignmentContainer!(string[]) getLasFiles(string dbA)
{
    return getLasFiles(dbA, null);
}

AlignmentContainer!(string[]) getLasFiles(string dbA, string dbB)
{
    import std.format : format;

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
            import std.format : format;

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
        import std.format : format, formattedRead;

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
        import std.format : format;

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
