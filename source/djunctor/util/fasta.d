/**
    Some functions to work with FASTA data.

    Copyright: © 2018 Arne Ludwig <arne.ludwig@posteo.de>
    License: Subject to the terms of the MIT license, as written in the
             included LICENSE file.
    Authors: Arne Ludwig <arne.ludwig@posteo.de>
*/
module djunctor.util.fasta;

import std.algorithm : count, equal, joiner, startsWith;
import std.array : appender, array;
import std.conv : to;
import std.format : format, formattedRead;
import std.range : chain, chunks, drop, ElementType, isBidirectionalRange, only,
    take, walkLength;
import std.string : indexOf, lineSplitter, outdent;
import std.traits : isSomeChar, isSomeString;
import std.typecons : tuple, Tuple;

/**
    Gives access to FASTA data. Does not copy the input sequence.
*/
template Fasta(T) if (isSomeString!T)
{
    struct Fasta
    {
        static immutable headerIndicator = '>';
        const T data;
        private size_t[] recordIndex;

        alias data this;

        /**
            Build an index in order to give fast access to individual records.
            This is called implicitly when accessing individual records using
            `opIndex` or `length`.
        */
        void buildIndex()
        {
            if ((data.length > 0 && recordIndex.length > 0) || (data.length == 0
                    && recordIndex.length == 0))
                return;

            recordIndex.reserve(data.count(headerIndicator) + 1);
            long currIdx = data.indexOf(headerIndicator);

            while (currIdx >= 0)
            {
                recordIndex ~= currIdx.to!size_t;
                currIdx = data.indexOf(headerIndicator, currIdx + 1);
            }

            recordIndex ~= data.length;
        }

        /**
            Get the FASTA record at idx (zero-based).

            Returns: `FastaRecord!T` at index idx.
        */
        FastaRecord!T opIndex(size_t idx)
        {
            assert(0 <= idx && idx < length, "index out of bounds");
            buildIndex();
            auto recordBegin = recordIndex[idx];
            auto recordEnd = recordIndex[idx + 1];

            return data[recordBegin .. recordEnd].parseFastaRecord();
        }

        /// Get the number of FASTA records.
        @property size_t length()
        {
            buildIndex();

            return recordIndex.length - 1;
        }

        /// Returns true iff line starts with '>'.
        static bool isHeaderLine(in T line) pure
        {
            return line.startsWith(only(headerIndicator));
        }
    }
}

///
unittest
{
    auto fasta1 = Fasta!string(q"EOF
        >sequence1
        CTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCT
        AACCCTAACCCTAACCCTAACCCTAACCCTAACAACCCTAACCCTAACCC
        >sequence2
        AAGCTGAGCAGGGCTTTAAAGCTATCTTATTAATAATTATTTCTGTATTG
        TAACCCTAACCCTAAACCTAACCCTAACCCTAACCCTAACAACCCTAACC
EOF".outdent);
    // dfmt off
    auto fasta1Records = [
        q"EOF
            >sequence1
            CTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCT
            AACCCTAACCCTAACCCTAACCCTAACCCTAACAACCCTAACCCTAACCC
EOF".outdent.parseFastaRecord,
        q"EOF
            >sequence2
            AAGCTGAGCAGGGCTTTAAAGCTATCTTATTAATAATTATTTCTGTATTG
            TAACCCTAACCCTAAACCTAACCCTAACCCTAACCCTAACAACCCTAACC
EOF".outdent.parseFastaRecord,
    ];
    // dfmt on

    assert(fasta1.length == 2, fasta1.length.to!string);
    assert(fasta1[0] == fasta1Records[0]);
    assert(fasta1[1] == fasta1Records[1]);
}

/// Convenience wrapper around `Fasta!T(T data)`.
Fasta!T parseFasta(T)(T data)
{
    return typeof(return)(data);
}

///
unittest
{
    string fastaData = q"EOF
        >sequence1
        CTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCT
        AACCCTAACCCTAACCCTAACCCTAACCCTAACAACCCTAACCCTAACCC
        >sequence2
        AAGCTGAGCAGGGCTTTAAAGCTATCTTATTAATAATTATTTCTGTATTG
        TAACCCTAACCCTAAACCTAACCCTAACCCTAACCCTAACAACCCTAACC
EOF".outdent;
    auto fasta = fastaData.parseFasta();
    // dfmt off
    auto fastaRecords = [
        FastaRecord!string(q"EOF
            >sequence1
            CTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCT
            AACCCTAACCCTAACCCTAACCCTAACCCTAACAACCCTAACCCTAACCC
EOF".outdent),
        FastaRecord!string(q"EOF
            >sequence2
            AAGCTGAGCAGGGCTTTAAAGCTATCTTATTAATAATTATTTCTGTATTG
            TAACCCTAACCCTAAACCTAACCCTAACCCTAACCCTAACAACCCTAACC
EOF".outdent),
    ];
    // dfmt on

    assert(fasta.length == 2, fasta.length.to!string);
    assert(fasta[0] == fastaRecords[0]);
    assert(fasta[1] == fastaRecords[1]);
}

/**
    Gives access to a single FASTA record. Does not copy the input sequence.
*/
template FastaRecord(T) if (isSomeString!T)
{
    struct FastaRecord
    {
        static immutable lineSep = "\n";
        alias Slice = Tuple!(int, int);

        const T data;

        alias data this;

        /// Get this record in FASTA format.
        auto toFasta(int lineWidth = 50) pure const
        {
            auto formattedBody = this[].chunks(lineWidth).joiner(lineSep);

            return chain(header, lineSep, formattedBody, lineSep);
        }

        /// Get the complete header line.
        @property auto header() pure const
        {
            return data.lineSplitter.front;
        }

        /// Get the length of the sequence (in characters).
        @property size_t length() pure const
        {
            return this[].walkLength;
        }

        @property size_t opDollar(size_t dim : 0)()
        {
            return length;
        }

        /// Get the sequence of this FASTA record without newlines.
        auto opIndex() pure const
        {
            return data.lineSplitter.drop(1).joiner;
        }

        /// Get the sequence character at index `i` of this FASTA record.
        auto opIndex(int i) pure const
        {
            i = normalizeIndex(i);
            assert(0 <= i && i < length, "index out of bounds");

            return this[].drop(i).front;
        }

        auto opIndex(in Slice slice) pure const
        {
            auto i = normalizeIndex(slice[0]);
            auto j = normalizeIndex(slice[1]);
            assert(0 <= i && i <= j && j < length, format!"index out of bounds: %s"(slice));

            return this[].drop(i).take(j - i);
        }

        /// Get sub-sequence from `i` to `j` (exclusive) of this FASTA record.
        auto opSlice(size_t dim : 0)(int i, int j)
        {
            return tuple(i, j);
        }

        /// ditto
        auto opSlice(size_t dim : 0)(size_t i, size_t j)
        {
            return tuple(i.to!int, j.to!int);
        }

        private int normalizeIndex(int i) const
        {
            auto length = this.length;

            while (i < 0)
                i += length;

            return i;
        }
    }
}

///
unittest
{
    // dfmt off
    auto fastaRecord1 = q"EOF
        >sequence1
        CTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCT
        AACCCTAACCCTAACCCTAACCCTAACCCTAACAACCCTAACCCTAACCC
EOF".outdent.parseFastaRecord;
    // dfmt on

    assert(fastaRecord1.header == ">sequence1");
    assert(fastaRecord1[].equal("CTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACAACCCTAACCCTAACCC"));
    assert(fastaRecord1[0 .. 5].equal("CTAAC"));
    assert(fastaRecord1.toFasta(13).equal(q"EOF
        >sequence1
        CTAACCCTAACCC
        TAACCCTAACCCT
        AACCCTAACCCTA
        ACCCTAACCCTAA
        CCCTAACCCTAAC
        CCTAACCCTAACC
        CTAACAACCCTAA
        CCCTAACCC
EOF".outdent));

    // dfmt off
    auto fastaRecord2 = q"EOF
        >sequence2
        AAGCTGAGCAGGGCTTTAAAGCTATCTTATTAATAATTATTTCTGTATTG
        TAACCCTAACCCTAAACCTAACCCTAACCCTAACCCTAACAACCCTAACC
EOF".outdent.parseFastaRecord;
    // dfmt on

    assert(fastaRecord2.header == ">sequence2");
    assert(fastaRecord2[].equal("AAGCTGAGCAGGGCTTTAAAGCTATCTTATTAATAATTATTTCTGTATTGTAACCCTAACCCTAAACCTAACCCTAACCCTAACCCTAACAACCCTAACC"));
    assert(fastaRecord2[5 .. 10].equal("GAGCA"));
    assert(fastaRecord2.toFasta(45).equal(q"EOF
        >sequence2
        AAGCTGAGCAGGGCTTTAAAGCTATCTTATTAATAATTATTTCTG
        TATTGTAACCCTAACCCTAAACCTAACCCTAACCCTAACCCTAAC
        AACCCTAACC
EOF".outdent));
}

/// Convenience wrapper around `FastaRecord!T(T data)`.
FastaRecord!T parseFastaRecord(T)(T data)
{
    return typeof(return)(data);
}

///
unittest
{
    string fastaRecordData = q"EOF
        >sequence1
        CTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCT
        AACCCTAACCCTAACCCTAACCCTAACCCTAACAACCCTAACCCTAACCC
EOF".outdent;
    auto fastaRecord = fastaRecordData.parseFastaRecord();

    assert(fastaRecord.header == ">sequence1");
    assert(fastaRecord[].equal("CTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACAACCCTAACCCTAACCC"));
    assert(fastaRecord[0 .. 5].equal("CTAAC"));
}

/// Build a `FastaRecord!T` from header and sequence.
FastaRecord!T buildFastaRecord(T)(in T header, in T sequence) if (isSomeString!T)
{
    immutable lineSep = FastaRecord!T.lineSep;
    auto builder = appender!T;

    builder.reserve(header.length + sequence.length + 2 * lineSep.length);

    builder ~= header;
    builder ~= lineSep;
    builder ~= sequence;
    builder ~= lineSep;

    return typeof(return)(builder.data);
}

// TODO test

template PacBioHeader(T) if (isSomeString!T)
{
    struct PacBioHeader
    {
        static immutable headerFormat = ">%s/%d/%d_%d RQ=%f";

        T name;
        size_t id;
        size_t well;
        size_t sequenceLength;
        float readQuality;

        /// Construct a `PacBioHeader!T` from `header`.
        this(T header)
        {
            this.parse(header);
        }

        /// Assign new `header` data.
        void opAssign(T header)
        {
            this.parse(header);
        }

        /// Builds the header string.
        S to(S : T)() const
        {
            return buildHeader();
        }

        private T buildHeader() const
        {
            // dfmt off
            return format!headerFormat(
                name,
                id,
                well,
                sequenceLength,
                readQuality,
            );
            // dfmt on
        }

        private void parse(in T header)
        {
            // dfmt off
            auto numMatches = header[].formattedRead!headerFormat(
                name,
                id,
                well,
                sequenceLength,
                readQuality,
            );
            // dfmt on

            assert(numMatches == 5);
        }
    }
}

///
unittest
{
    string header = ">name/1/0_1337 RQ=0.75";
    auto pbHeader1 = PacBioHeader!string(header);

    assert(pbHeader1.to!string == ">name/1/0_1337 RQ=0.750000");
    assert(pbHeader1.name == "name");
    assert(pbHeader1.id == 1);
    assert(pbHeader1.well == 0);
    assert(pbHeader1.sequenceLength == 1337);
    assert(pbHeader1.readQuality == 0.75);

    PacBioHeader!string pbHeader2 = header;

    assert(pbHeader2 == pbHeader1);
}

/// Convenience wrapper around `PacBioHeader!T(T header)`.
PacBioHeader!T parsePacBioHeader(T)(T header)
{
    return typeof(return)(header);
}

///
unittest
{
    string header = ">name/1/0_1337 RQ=0.75";
    auto pbHeader1 = header.parsePacBioHeader();

    assert(pbHeader1.to!string == ">name/1/0_1337 RQ=0.750000");
    assert(pbHeader1.name == "name");
    assert(pbHeader1.id == 1);
    assert(pbHeader1.well == 0);
    assert(pbHeader1.sequenceLength == 1337);
    assert(pbHeader1.readQuality == 0.75);
}

/**
    Compute the reverse complement of a DNA sequence. Only bases A, T, C, G
    will be translated; all other characters are left as is. Replacement
    preserves casing of the characters.
*/
auto reverseComplementer(Range)(Range sequence)
        if (isBidirectionalRange!Range && isSomeChar!(ElementType!Range))
{
    import std.algorithm : substitute;
    import std.range : retro;

    // dfmt off
    return sequence
        .retro
        .substitute!(
            "A", "T",
            "G", "C",
            "T", "A",
            "C", "G",
            "a", "t",
            "g", "c",
            "t", "a",
            "c", "g",
        );
    // dfmt on
}

/// ditto
T reverseComplement(T)(in T sequence) if (isSomeString!T)
{
    import std.array : array;

    return sequence[].reverseComplementer.array.to!T;
}

FastaRecord!T reverseComplement(T)(in FastaRecord!T fastaRecord) if (isSomeString!T)
{
    immutable lineSep = FastaRecord!T.lineSep;
    auto header = fastaRecord.header;
    auto sequence = fastaRecord[].array.reverseComplement;
    auto builder = appender!T;

    builder.reserve(header.length + sequence.length + 2 * lineSep.length);

    builder ~= header;
    builder ~= lineSep;
    builder ~= sequence;
    builder ~= lineSep;

    return typeof(return)(builder.data);
}

///
unittest
{
    auto seq = "GGTTGTAAATTGACTGTTGTCTGCT\ngccaatctactggtgggggagagat";
    auto revComp = "atctctcccccaccagtagattggc\nAGCAGACAACAGTCAATTTACAACC";

    assert(seq.reverseComplement == revComp);
    assert(seq.reverseComplementer.equal(revComp));

    // dfmt off
    auto fastaRecord1 = q"EOF
        >sequence1
        CTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCT
        AACCCTAACCCTAACCCTAACCCTAACCCTAACAACCCTAACCCTAACCC
EOF".outdent.parseFastaRecord;
    // dfmt on
    auto fastaRecord1RevComp = fastaRecord1.reverseComplement;

    assert(fastaRecord1RevComp.header == ">sequence1");
    assert(fastaRecord1RevComp[].equal("GGGTTAGGGTTAGGGTTGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAG"));
    assert(fastaRecord1RevComp[0 .. 5].equal("GGGTT"));
    assert(fastaRecord1RevComp.toFasta(13).equal(q"EOF
        >sequence1
        GGGTTAGGGTTAG
        GGTTGTTAGGGTT
        AGGGTTAGGGTTA
        GGGTTAGGGTTAG
        GGTTAGGGTTAGG
        GTTAGGGTTAGGG
        TTAGGGTTAGGGT
        TAGGGTTAG
EOF".outdent));
}
