/**
    This is the main algorithm of this package.

    Copyright: © 2018 Arne Ludwig <arne.ludwig@posteo.de>
    License: Subject to the terms of the MIT license, as written in the
             included LICENSE file.
    Authors: Arne Ludwig <arne.ludwig@posteo.de>
*/
module djunctor.djunctor;

import djunctor.alignments : AlignmentChain, alignmentCoverage,
    buildPileUpsFromReadAlignments = buildPileUps, getAlignmentChainRefs,
    getType, haveEqualIds, isExtension, isGap, isValid, makeJoin, PileUp,
    ReadAlignment, ReadAlignmentType;
import djunctor.commandline : Options;
import djunctor.dazzler : attachTracePoints, buildDamFile, DalignerOptions,
    DamapperOptions, GapPart, getConsensus, getFastaEntries, getLocalAlignments,
    getMappings, getNumContigs, getScaffoldStructure, LAdumpOptions, readMask,
    ScaffoldPart, writeMask;
import djunctor.util.fasta : parseFastaRecord, parsePacBioHeader,
    reverseComplement;
import djunctor.util.log;
import djunctor.util.math : ceil, floor, NaturalNumberSet;
import djunctor.util.region : empty, min, Region, sup;
import djunctor.util.scaffold : ContigNode, contigStarts, enforceJoinPolicy,
    getDefaultJoin, getUnkownJoin, initScaffold, isAntiParallel, isDefault,
    isExtension, isFrontExtension, isGap, isUnkown, isValid, linearWalk,
    normalizeUnkownJoins, Scaffold;
import dstats.distrib : invPoissonCDF;
import std.algorithm : all, chunkBy, each, equal, filter, find, fold, isSorted,
    joiner, map, max, maxIndex, min, setDifference, sort, sum, swap,
    SwapStrategy, uniq;
import std.array : appender, array, join, minimallyInitializedArray;
import std.conv;
import std.exception : enforce, ErrnoException;
import std.format : format;
import std.math : abs, floor;
import std.range : assumeSorted, chain, chunks, drop, ElementType, InputRange,
    inputRangeObject, iota, isInputRange, only, repeat, retro, take, walkLength,
    zip;
import std.stdio : File, write, writeln;
import std.typecons : Flag, No, tuple, Tuple, Yes;
import vibe.data.json : Json, toJson = serializeToJson;

/// Start the `djunctor` algorithm with preprocessed options.
void runWithOptions(in ref Options options)
{
    new DJunctor(options).run();
}

alias ReferenceMask = Region!(size_t, size_t, "contigId");
alias ReferenceInterval = ReferenceMask.TaggedInterval;
alias ReferencePoint = ReferenceMask.TaggedPoint;
alias ReadMask = Region!(size_t, size_t, "readId");
alias ReadInterval = ReadMask.TaggedInterval;
alias ReadPoint = ReadMask.TaggedPoint;

/// Returns the alignment region of alignmentChain.
Mask getRegion(Mask = ReferenceMask)(in AlignmentChain alignmentChain) pure
        if (is(Mask : Region!Args, Args...))
{
    auto contigAId = alignmentChain.contigA.id;
    // dfmt off
    return Mask(alignmentChain
        .localAlignments
        .map!(la => Mask.TaggedInterval(
            contigAId,
            la.contigA.begin,
            la.contigA.end,
        ))
        .array
    );
    // dfmt on
}

/**
    Get points on the reference where the pileUp should be cropped. Returns
    one common trace point for each involved contig.

    See_Also: `getCommonTracePoint`
*/
private ReferencePoint[] getCroppingRefPositions(PileUp pileUp,
        in ReferenceMask mask, in size_t tracePointDistance)
{
    auto alignmentsByContig = pileUp.splitAlignmentsByContigA();
    auto commonTracePoints = alignmentsByContig.map!(acs => getCommonTracePoint(acs,
            mask, tracePointDistance)).array;
    auto contigAIds = alignmentsByContig.map!"a[0].contigA.id".array;

    auto cropPos = zip(contigAIds, commonTracePoints).map!(
            zipped => ReferencePoint(zipped.expand)).array;

    // dfmt off
    debug logJsonDebug(
        "alignmentsByContig", alignmentsByContig.toJson,
        "pileUp", pileUp.toJson,
        "cropPos", cropPos.toJson,
    );
    // dfmt on

    return cropPos;
}

/// Sort alignment chains of pileUp into groups with the same `contigA.id`.
private const(AlignmentChain)[][] splitAlignmentsByContigA(PileUp pileUp)
{
    if (pileUp.length == 0)
    {
        return [];
    }

    // dfmt off
    auto orderedAlignments = pileUp
        .map!"a[]"
        .joiner
        .array
        .sort;
    // dfmt on
    auto alignmentsByContig = orderedAlignments.chunkBy!"a.contigA.id == b.contigA.id";

    // dfmt off
    return alignmentsByContig
        .filter!"!a.empty"
        .map!(acsByContig => acsByContig
            .map!(ac => cast(const(typeof(ac))) ac))
        .map!array
        .array;
    // dfmt on
}

/// Returns a common trace points wrt. contigA that is not in mask and is on
/// the relevant half of the contig.
private long getCommonTracePoint(in AlignmentChain[] alignmentChains,
        in ReferenceMask mask, in size_t tracePointDistance)
{
    static long getCommonTracePoint(R)(R tracePointCandidates, ReferenceMask allowedTracePointRegion) pure
    {
        auto commonTracePoints = tracePointCandidates.filter!(c => c in allowedTracePointRegion);

        return commonTracePoints.empty ? -1 : commonTracePoints.front.value.to!long;
    }

    auto contigA = alignmentChains[0].contigA;
    auto alignmentType = ReadAlignment(alignmentChains[0]).type;
    assert(alignmentType != ReadAlignmentType.gap);
    // dfmt off
    auto commonAlignmentRegion = alignmentChains
        .map!getRegion
        .fold!"a & b";
    auto relevantContigHalf = alignmentType == ReadAlignmentType.front
        ? ReferenceMask(contigA.id, 0, contigA.length / 2)
        : ReferenceMask(contigA.id, contigA.length / 2, contigA.length);
    // dfmt on
    auto allowedTracePointRegion = (relevantContigHalf & commonAlignmentRegion) - mask;

    if (allowedTracePointRegion.empty)
    {
        return -1;
    }

    auto tracePointMin = ceil(min(allowedTracePointRegion), tracePointDistance);
    auto tracePointSup = ceil(sup(allowedTracePointRegion), tracePointDistance);
    // dfmt off
    auto tracePointCandidates = iota(tracePointMin, tracePointSup, tracePointDistance)
        .map!(tracePointCandidate => ReferencePoint(contigA.id, tracePointCandidate));
    // dfmt on

    // dfmt off
    return alignmentType == ReadAlignmentType.front
        ? getCommonTracePoint(tracePointCandidates.retro, allowedTracePointRegion)
        : getCommonTracePoint(tracePointCandidates, allowedTracePointRegion);
    // dfmt on
}

///
unittest
{
    immutable tracePointDistance = 5UL;

    size_t readId = 0;
    AlignmentChain getDummyRead(size_t begin, size_t end)
    {
        immutable referenceLength = 25;
        immutable extensionLength = 5;
        auto readLength = end - begin + extensionLength;

        with (AlignmentChain) with (LocalAlignment)
            {
                // dfmt off
                return AlignmentChain(
                    readId,
                    Contig(1, 25),
                    Contig(1, readLength),
                    Complement.no,
                    [LocalAlignment(
                        Locus(begin, end),
                        begin == 0
                            ? Locus(extensionLength, readLength)
                            : Locus(0, readLength - extensionLength),
                    )],
                );
                // dfmt on
            }
    }

    //        0    5   10   15   20   25   30
    //   ref: |------------------------|.....
    //        :    .    .    .    .    :
    // reads: :|----------------------------| #1
    //        :    |------------------------| #2
    //        :    .  |---------------------| #3
    //        :    .|-----------------------| #4
    //        :    .   |--------------------| #5
    //        :    .    .    .    .    :    .
    // mask1: :    [=======) .    .    :    .
    // mask2: :    .    .  [=======)   :    .
    // mask3: [========================)    .
    // dfmt off
    auto backExtensions = [
        getDummyRead(1, 25),
        getDummyRead(5, 25),
        getDummyRead(8, 25),
        getDummyRead(6, 25),
        getDummyRead(9, 25),
    ];
    // dfmt on

    //       -5    0    5   10   15   20   25
    //   ref: .....|------------------------|
    //        .    :    .    .    .    .    :
    // reads: |----------------------------|: #1
    //        |------------------------|    : #2
    //        |---------------------|  .    : #3
    //        |-----------------------|.    : #4
    //        |--------------------|   .    : #5
    //        .    :    .    .    .    .    :
    // mask1: .    :    [=======) .    .    :
    // mask2: .    :    .    .  [=======)   :
    // mask3: .    [========================)
    // dfmt off
    auto frontExtensions = [
        getDummyRead(0, 24),
        getDummyRead(0, 20),
        getDummyRead(0, 17),
        getDummyRead(0, 19),
        getDummyRead(0, 16),
    ];
    // dfmt on
    auto mask1 = ReferenceMask(1, 5, 13);
    auto mask2 = ReferenceMask(1, 13, 21);
    auto mask3 = ReferenceMask(1, 0, 25);

    assert(getCommonTracePoint(backExtensions, mask1, tracePointDistance) == 15);
    assert(getCommonTracePoint(backExtensions, mask2, tracePointDistance) == -1);
    assert(getCommonTracePoint(backExtensions, mask3, tracePointDistance) == -1);
    assert(getCommonTracePoint(frontExtensions, mask1, tracePointDistance) == 0);
    assert(getCommonTracePoint(frontExtensions, mask2, tracePointDistance) == 10);
    assert(getCommonTracePoint(frontExtensions, mask3, tracePointDistance) == -1);
}

ReadInterval getCroppingSlice(const AlignmentChain alignmentChain,
        in ReferencePoint[] croppingRefPoints)
{
    auto alignmentType = ReadAlignment(alignmentChain).type;
    auto read = alignmentChain.contigB;
    auto tracePointDistance = alignmentChain.tracePointDistance;
    // dfmt off
    auto croppingRefPos = croppingRefPoints
        .filter!(p => p.contigId == alignmentChain.contigA.id)
        .front
        .value;
    // dfmt on

    size_t croppingTracePointIdx(in AlignmentChain.LocalAlignment localAlignment)
    {
        auto firstTracePointRefPos = ceil(localAlignment.contigA.begin, tracePointDistance);
        auto openInterval = alignmentType == ReadAlignmentType.front;
        assert(croppingRefPos >= firstTracePointRefPos);

        return (croppingRefPos - firstTracePointRefPos) / tracePointDistance + (openInterval ? 0 : 1);
    }

    bool coversCroppingRefPos(in AlignmentChain.LocalAlignment localAlignment)
    {
        return localAlignment.contigA.begin <= croppingRefPos
            && croppingRefPos < localAlignment.contigA.end;
    }

    // dfmt off
    auto coveringLocalAlignment = alignmentChain
        .localAlignments
        .filter!coversCroppingRefPos
        .front;
    auto readCroppingPos =
        coveringLocalAlignment.contigB.begin +
        coveringLocalAlignment
            .tracePoints[0 .. croppingTracePointIdx(coveringLocalAlignment)]
            .map!"a.numBasePairs"
            .sum;
    // dfmt on
    size_t readBeginIdx;
    size_t readEndIdx;

    final switch (alignmentType)
    {
    case ReadAlignmentType.front:
        readBeginIdx = 0;
        readEndIdx = readCroppingPos;
        break;
    case ReadAlignmentType.back:
        readBeginIdx = readCroppingPos;
        readEndIdx = read.length;
        break;
    case ReadAlignmentType.gap:
        assert(0, "unreachable");
    }

    if (alignmentChain.complement)
    {
        swap(readBeginIdx, readEndIdx);
        readBeginIdx = read.length - readBeginIdx;
        readEndIdx = read.length - readEndIdx;
    }

    return ReadInterval(read.id, readBeginIdx, readEndIdx);
}

auto cropPileUp(PileUp pileUp, in ReferenceMask mask, in Options options)
{
    auto cropper = PileUpCropper(pileUp, mask, options);
    cropper.buildDb();

    return cropper.result;
}

private struct PileUpCropper
{
    PileUp pileUp;
    const(ReferenceMask) repeatMask;
    const(Options) options;
    private ReferencePoint[] croppingRefPositions;
    private string croppedDb;

    @property auto result()
    {
        return tuple!("db", "referencePositions")(croppedDb, croppingRefPositions);
    }

    void buildDb()
    {
        preparePileUp();
        fetchCroppingRefPositions();
        auto croppedFastEntries = getCroppedFastaEntries();
        croppedDb = buildDamFile(croppedFastEntries, options);
    }

    private void preparePileUp()
    {
        fetchTracePoints(pileUp, options);
        pileUp.sort!"a[0].contigB.id < b[0].contigB.id"();
    }

    private void fetchCroppingRefPositions()
    {
        auto tracePointDistance = pileUp[0][0].tracePointDistance;
        croppingRefPositions = pileUp.getCroppingRefPositions(repeatMask, tracePointDistance);
        logJsonDebug("croppingRefPositions", croppingRefPositions.toJson);

        foreach (refPos; croppingRefPositions)
        {
            enforce!Exception(refPos.value != -1, "could not find a common trace point");
        }
    }

    private auto getCroppedFastaEntries()
    {
        auto croppedFastaEntries = minimallyInitializedArray!(string[])(pileUp.length);

        foreach (croppedDbIdx, readAlignment, rawFastaEntry; pileUpWithSequence())
        {
            auto readCroppingSlice = getReadCroppingSlice(readAlignment);
            // dfmt off
            debug logJsonDebug(
                "croppedDbIdx", croppedDbIdx,
                "readCroppingSlice", readCroppingSlice.toJson,
            );
            // dfmt on

            croppedFastaEntries[croppedDbIdx] = cropFastaEntry(rawFastaEntry, readCroppingSlice);
        }

        return croppedFastaEntries;
    }

    private auto pileUpWithSequence()
    {
        size_t[] readIds = pileUp.map!"a[0].contigB.id + 0".array;

        // dfmt off
        return zip(
            iota(pileUp.length),
            pileUp,
            getFastaEntries(options.readsDb, readIds, options).map!parseFastaRecord,
        );
        // dfmt on
    }

    private auto getReadCroppingSlice(ReadAlignment readAlignment)
    {
        // dfmt off
        auto readCroppingSlice = readAlignment[]
            .map!(alignmentChain => alignmentChain.getCroppingSlice(croppingRefPositions))
            .fold!"a & b";
        // dfmt on
        assert(!readCroppingSlice.empty, "invalid/empty read cropping slice");

        return readCroppingSlice;
    }

    private string cropFastaEntry(T)(T rawFastaEntry, in ReadInterval readCroppingSlice)
    {
        immutable lineSep = typeof(pileUpWithSequence().front[2]).lineSep;

        auto croppedLength = readCroppingSlice.size;
        auto newHeader = buildNewHeader(rawFastaEntry, readCroppingSlice);
        auto expectedFastaLength = newHeader.length + croppedLength
            + croppedLength / options.fastaLineWidth + 2;
        auto croppedFastaEntry = appender!string;

        croppedFastaEntry.reserve(expectedFastaLength);
        croppedFastaEntry ~= newHeader;
        croppedFastaEntry ~= lineSep;
        // dfmt off
        croppedFastaEntry ~= rawFastaEntry[readCroppingSlice.begin .. readCroppingSlice.end]
            .chunks(options.fastaLineWidth)
            .joiner(lineSep);
        // dfmt on
        croppedFastaEntry ~= lineSep;

        return croppedFastaEntry.data;
    }

    private string buildNewHeader(T)(T rawFastaEntry, in ReadInterval readCroppingSlice)
    {
        auto headerBuilder = parsePacBioHeader(rawFastaEntry.header);
        headerBuilder.qualityRegionBegin = 0;
        headerBuilder.qualityRegionBegin = readCroppingSlice.size;

        return headerBuilder.to!string;
    }
}

ref PileUp fetchTracePoints(ref PileUp pileUp, in Options options)
{
    auto allAlignmentChains = pileUp.getAlignmentChainRefs();
    allAlignmentChains.sort!("*a < *b", SwapStrategy.stable);
    allAlignmentChains.attachTracePoints(options.refDb, options.readsDb,
            options.damapperOptions, options);

    return pileUp;
}

// dfmt off
/// Information about the point where the two sequences should be spliced.
alias SpliceSite = Tuple!(
    ReferencePoint, "croppingRefPosition",
    AlignmentChain.Complement, "complement",
);
/// This characterizes an insertion.
alias InsertionInfo = Tuple!(
    string, "sequenceDb",
    SpliceSite[], "spliceSites",
);
// dfmt on
alias ResultScaffold = Scaffold!InsertionInfo;
/// This characterizes an insertion.
alias Insertion = ResultScaffold.Edge;

Insertion concatenateSpliceSites(Insertion existingJoin, Insertion newJoin)
{
    if (existingJoin.payload.sequenceDb is null)
    {
        existingJoin.payload.sequenceDb = newJoin.payload.sequenceDb;
    }
    existingJoin.payload.spliceSites ~= newJoin.payload.spliceSites;

    return existingJoin;
}

auto getFastaRecord(in string sequenceDb, in size_t sequenceId, in Options options)
{
    if (sequenceDb is null)
    {
        return typeof(parseFastaRecord("")).init;
    }

    auto records = getFastaEntries(sequenceDb, [sequenceId + 0], options).map!parseFastaRecord;

    if (records.empty)
    {
        throw new Exception(format!"could not fetch %d from %s"(sequenceId, sequenceDb));
    }

    return records.front;
}

class DJunctor
{
    /// Stop after maxLoops `mainLoop`s at the latest.
    static immutable maxLoops = 1;
    /**
        Two scores are considered similar if the relative "error" of them is
        smaller than defaulMaxRelativeDiff.

        The relative error is defined as:

            relError(a, b) = abs(a - b)/max(a - b)

        **Implementation note:** all computations are done in integer
        arithmetic hence AlignmentChain.maxScore corresponds to 1 in the above
        equation.
    */
    static immutable maxRelativeDiff = AlignmentChain.maxScore / 20; // 5% of larger value
    /**
        Two scores are considered similar if the absolute "error" of them is
        smaller than defaulMaxAbsoluteDiff.

        The relative error is defined as:

            absError(a, b) = abs(a - b)

        **Implementation note:** all computations are done in integer
        arithmetic hence all scores are lower than or equal to
        AlignmentChain.maxScore .
    */
    static immutable maxAbsoluteDiff = AlignmentChain.maxScore / 100; // 1% wrt. score

    size_t numReferenceContigs;
    size_t numReads;
    AlignmentChain[] selfAlignment;
    ReferenceMask repetitiveRegions;
    AlignmentChain[] readsAlignment;
    const Options options;
    /// Set of "good" alignments to be considered in further processing.
    AlignmentChain[] catCandidates;
    /// A scaffold graph accumulating all planned insertions.
    ResultScaffold catHits;
    /// Scaffold structure of the reference.
    const(ScaffoldPart)[] scaffoldStructure;
    /// Keep track of unused reads.
    NaturalNumberSet unusedReads;

    this(in ref Options options)
    {
        this.options = options;
    }

    void run()
    {
        logJsonDiagnostic("state", "enter", "function", "run");
        // dfmt off
        init();
        assessRepeatStructure();
        filterAlignments();

        PileUp[] pileUps = buildPileUps();
        foreach (ref pileUp; pileUps)
        {
            processPileUp(pileUp);
        }

        // dfmt off
        logJsonDebug("catHits", [
            "nodes": catHits.nodes.toJson,
            "joins": catHits.edges.map!(join => [
                "start": join.start.toJson,
                "end": join.end.toJson,
                "payload": join.payload.toJson,
            ]).array.toJson,
        ]);
        // dfmt on

        writeNewAssembly();
        writeUnusedReadsList();
        // dfmt on
        logJsonDiagnostic("state", "exit", "function", "run");
    }

    protected void init()
    {
        static struct SelfAlignmentOptions
        {
            string[] dalignerOptions;
            string[] ladumpOptions;
            string workdir;
        }

        static struct ReadsToRefAlignmentOptions
        {
            string[] damapperOptions;
            string[] ladumpOptions;
            string workdir;
        }

        logJsonDiagnostic("state", "enter", "function", "djunctor.init");
        numReferenceContigs = getNumContigs(options.refDb, options);
        numReads = getNumContigs(options.readsDb, options);
        scaffoldStructure = getScaffoldStructure(options.refDb, options).array;

        if (options.inMask !is null)
        {
            this.repetitiveRegions = ReferenceMask(readMask!ReferenceInterval(options.refDb,
                    options.inMask, options));
        }

        // dfmt off
        selfAlignment = getLocalAlignments(options.refDb, const(SelfAlignmentOptions)(
            options.dalignerOptions ~ [
                format!(DalignerOptions.minAlignmentLength ~ "%d")(options.minAnchorLength),
                format!(DalignerOptions.averageCorrelationRate ~ "%f")((1 - options.referenceErrorRate)^^2),
            ],
            options.ladumpOptions[].filter!(o => o != LAdumpOptions.properOverlapsOnly).array,
            options.workdir,
        ));
        readsAlignment = getMappings(options.refDb, options.readsDb, const(ReadsToRefAlignmentOptions)(
            options.damapperOptions ~ [
                format!(DamapperOptions.averageCorrelationRate ~ "%f")((1 - options.referenceErrorRate) * (1 - options.readsErrorRate)),
            ],
            options.ladumpOptions[].array,
            options.workdir,
        ));
        // dfmt on
        catCandidates = readsAlignment.dup;

        // dfmt off
        catHits = initScaffold!(
            (contigId) => InsertionInfo(options.refDb, []),
            InsertionInfo
        )(numReferenceContigs);
        // dfmt on
        insertUnkownJoins();

        unusedReads.reserveFor(numReads);
        foreach (readId; iota(1, numReads + 1))
        {
            unusedReads.add(readId);
        }
        logJsonDiagnostic("state", "exit", "function", "djunctor.init");
    }

    protected void assessRepeatStructure()
    {
        logJsonDiagnostic("state", "enter", "function", "djunctor.assessRepeatStructure");

        auto alphaHalf = (1 - options.confidence) / 2;
        auto selfCoverage = alignmentCoverage(selfAlignment);
        auto readsCoverage = alignmentCoverage(readsAlignment);
        // dfmt off
        auto selfCoverageConfidenceInterval = tuple(
            invPoissonCDF(alphaHalf, selfCoverage),
            invPoissonCDF(1 - alphaHalf, selfCoverage),
        );
        auto readsCoverageConfidenceInterval = tuple(
            invPoissonCDF(alphaHalf, readsCoverage),
            invPoissonCDF(1 - alphaHalf, readsCoverage),
        );
        logJsonDebug(
            "selfCoverage", selfCoverage,
            "selfCoverageConfidenceInterval", selfCoverageConfidenceInterval.toJson,
            "readsCoverage", readsCoverage,
            "readsCoverageConfidenceInterval", readsCoverageConfidenceInterval.toJson,
        );
        alias AssessmentStage(Assessor : RepeatAssessor) = Tuple!(
            string, "name",
            Assessor, "assessor",
            const(AlignmentChain[]), "input",
        );
        // dfmt on

        AssessmentStage!Assessor assessmentStage(Assessor : RepeatAssessor)(string name,
                Assessor assessor, const(AlignmentChain[]) input)
        {
            return typeof(return)(name, assessor, input);
        }

        // dfmt off
        auto assessmentStages = tuple(
            assessmentStage("self-alignment", new BadAlignmentCoverageAssessor(selfCoverageConfidenceInterval.expand), selfAlignment),
            assessmentStage("reads-alignment", new BadAlignmentCoverageAssessor(readsCoverageConfidenceInterval.expand), readsAlignment),
        );
        // dfmt on

        foreach (stage; assessmentStages)
        {
            auto repetitiveRegions = stage.assessor(stage.input);

            // dfmt off
            logJsonDebug(
                "assessor", stage.name,
                "repetitiveRegions", repetitiveRegions.intervals.toJson,
            );
            // dfmt on

            if (shouldLog(LogLevel.diagnostic) && options.outMask != null)
            {
                auto maskName = format!"%s-%s"(options.outMask, stage.name);

                writeMask(options.refDb, maskName, repetitiveRegions.intervals, options);
            }

            this.repetitiveRegions |= repetitiveRegions;
        }

        // dfmt off
        logJsonDebug(
            "assessor", "finalResult",
            "repetitiveRegions", this.repetitiveRegions.intervals.toJson,
        );
        // dfmt on

        if (options.outMask != null)
        {
            writeMask(options.refDb, options.outMask, this.repetitiveRegions.intervals, options);
        }

        logJsonDiagnostic("state", "exit", "function", "djunctor.assessRepeatStructure");
    }

    protected void filterAlignments()
    {
        logJsonDiagnostic("state", "enter", "function", "djunctor.filterReads");

        // dfmt off
        auto filters = tuple(
            new WeaklyAnchoredAlignmentChainsFilter(repetitiveRegions, options.minAnchorLength),
            new ImproperAlignmentChainsFilter(),
            new AmbiguousAlignmentChainsFilter(&unusedReads),
            new RedundantAlignmentChainsFilter(&unusedReads),
        );
        // dfmt on
        AlignmentChain[] filterInput = catCandidates[];
        // dfmt off
        logJsonDiagnostic(
            "filterStage", "Input",
            "discardedAlignmentChains", Json.emptyArray,
            "keptAlignmentChains", filterInput.toJson,
        );
        // dfmt on
        foreach (i, filter; filters)
        {
            auto filterOutput = filter(filterInput[]);

            assert(isSorted(filterOutput));
            if (shouldLog(LogLevel.diagnostic))
            {
                // dfmt off
                auto discardedAlignmentChains = setDifference(
                    filterInput,
                    filterOutput,
                ).array;
                logJsonDiagnostic(
                    "filterStage", typeof(filter).stringof,
                    "discardedAlignmentChains", discardedAlignmentChains.toJson,
                    "keptAlignmentChains", filterOutput.toJson,
                );
                // dfmt on
            }

            filterInput = filterOutput;
        }
        auto filterPipelineOutput = filterInput;

        this.catCandidates = filterPipelineOutput;

        logJsonDiagnostic("state", "exit", "function", "djunctor.filterReads");
    }

    protected static bool similarScore(size_t a, size_t b) pure
    {
        auto diff = a > b ? a - b : b - a;
        auto magnitude = a > b ? a : b;

        return diff < maxAbsoluteDiff || (diff * AlignmentChain.maxScore / magnitude) < maxRelativeDiff;
    }

    unittest
    {
        immutable eps = 1;
        immutable refValueSmall = 2 * maxAbsoluteDiff;
        immutable refValueLarge = AlignmentChain.maxScore / 2;

        // test absolute part
        assert(similarScore(refValueSmall, refValueSmall));
        assert(similarScore(refValueSmall, refValueSmall + (maxAbsoluteDiff - 1)));
        assert(similarScore(refValueSmall, refValueSmall - (maxAbsoluteDiff - 1)));
        assert(!similarScore(refValueSmall, refValueSmall + (maxAbsoluteDiff + 1)));
        assert(!similarScore(refValueSmall, refValueSmall - (maxAbsoluteDiff + 1)));
        // test relative part
        assert(similarScore(refValueLarge, refValueLarge));
        assert(similarScore(refValueLarge,
                refValueLarge + refValueLarge * maxRelativeDiff / AlignmentChain.maxScore));
        assert(similarScore(refValueLarge,
                refValueLarge - refValueLarge * maxRelativeDiff / AlignmentChain.maxScore) + eps);
        assert(!similarScore(refValueLarge,
                refValueLarge + refValueLarge * 2 * maxRelativeDiff / AlignmentChain.maxScore));
        assert(!similarScore(refValueLarge,
                refValueLarge - refValueLarge * 2 * maxRelativeDiff / AlignmentChain.maxScore));
    }

    protected void insertUnkownJoins()
    {
        // dfmt off
        scaffoldStructure[]
            .filter!(part => part.peek!GapPart !is null)
            .map!(gapPart => gapPart.get!GapPart)
            .map!(gapPart => getUnkownJoin(
                gapPart.beginGlobalContigId,
                gapPart.endGlobalContigId,
                InsertionInfo(
                    null,
                    [
                        SpliceSite(ReferencePoint(
                            gapPart.beginGlobalContigId,
                            gapPart.begin,
                        ), AlignmentChain.Complement.no),
                        SpliceSite(ReferencePoint(
                            gapPart.beginGlobalContigId,
                            gapPart.end,
                        ), AlignmentChain.Complement.no),
                    ]
                ),
            ))
            .each!(unkownJoin => catHits.add(unkownJoin));
        // dfmt on
    }

    protected PileUp[] buildPileUps()
    {
        logJsonDiagnostic("state", "enter", "function", "djunctor.buildPileUps");

        // dfmt off
        auto pileUps = buildPileUpsFromReadAlignments(numReferenceContigs, catCandidates)
            .filter!(pileUp => pileUp.length >= options.minReadsPerPileUp)
            .array;
        // dfmt on

        // dfmt off
        logJsonDebug("pileUps", pileUps
            .map!(pileUp => Json([
                "type": Json(pileUp.getType.to!string),
                "readAlignments": pileUp.map!"a[]".array.toJson,
            ]))
            .array
            .toJson);
        // dfmt on

        logJsonDiagnostic("state", "exit", "function", "djunctor.buildPileUps");

        return pileUps;
    }

    protected void processPileUp(ref PileUp pileUp)
    {
        auto croppingResult = cropPileUp(pileUp, repetitiveRegions, options);
        auto referenceReadIdx = bestReadAlignmentIndex(pileUp, Yes.preferSpanning, options);
        auto referenceRead = pileUp[referenceReadIdx];
        auto consensusDb = buildConsensus(croppingResult.db, referenceReadIdx + 1);

        if (pileUp.isExtension && shouldSkipShortExtension(croppingResult,
                consensusDb, referenceRead))
        {
            return;
        }

        addInsertionToScaffold(referenceRead, consensusDb, croppingResult.referencePositions);
        addFlankingContigSlicesToScaffold(croppingResult.referencePositions);
        markReadsAsUsed(pileUp);
    }

    protected size_t bestReadAlignmentIndex(in PileUp pileUp,
            Flag!"preferSpanning" preferSpanning, in Options options) pure
    {
        auto pileUpType = pileUp.getType;

        // dfmt off
        return pileUp
            .map!((readAlignment) => insertionScore(readAlignment, pileUpType, preferSpanning))
            .maxIndex;
        // dfmt on
    }

    protected bool shouldSkipShortExtension(T)(T croppingResult,
            in string consensusDb, in ReadAlignment referenceRead)
    {
        assert(croppingResult.referencePositions.length == 1);

        auto consensusSequence = getFastaRecord(consensusDb, 1, options);
        auto refPos = croppingResult.referencePositions[0].value;
        auto refLength = referenceRead[0].contigA.length;
        ulong extensionLength = referenceRead.isFrontExtension
            ? consensusSequence.length.to!ulong - refPos.to!ulong
            : consensusSequence.length.to!ulong - (refLength - refPos).to!ulong;

        return extensionLength < options.minExtensionLength;
    }

    protected void addInsertionToScaffold(ref ReadAlignment referenceRead,
            in string consensusDb, ref ReferencePoint[] referencePositions)
    {
        auto insertion = makeJoin!Insertion(referenceRead);
        // dfmt off
        insertion.payload = InsertionInfo(
            consensusDb,
            zip(referencePositions, referenceRead[].map!"a.complement")
                .map!(spliceSite => cast(SpliceSite) spliceSite)
                .array,
        );
        // dfmt on
        catHits.add(insertion);
        logJsonDebug("insertion", insertion.toJson);
    }

    protected void addFlankingContigSlicesToScaffold(ref ReferencePoint[] referencePositions)
    {
        foreach (ref referencePosition; referencePositions)
        {
            auto contigEdge = getDefaultJoin!InsertionInfo(referencePosition.contigId);
            // dfmt off
            contigEdge.payload = InsertionInfo(
                options.refDb,
                [SpliceSite(referencePosition, AlignmentChain.Complement.no)],
            );
            // dfmt on
            auto insertedEdge = catHits.add!concatenateSpliceSites(contigEdge);
            // dfmt off
            logJsonDebug(
                "info", "contigEdgeUpdate",
                "insertion", insertedEdge.toJson,
            );
            // dfmt on
        }
    }

    protected void markReadsAsUsed(ref in PileUp pileUp)
    {
        foreach (readAlignment; pileUp)
        {
            foreach (alignmentChain; readAlignment)
            {
                unusedReads.remove(alignmentChain.contigB.id);
            }
        }
    }

    protected size_t insertionScore(in ReadAlignment readAlignment,
            in ReadAlignmentType pileUpType,
            in Flag!"preferSpanning" preferSpanning = No.preferSpanning) pure
    {
        immutable shortAnchorPenaltyMagnitude = AlignmentChain.maxScore / 512;
        immutable notSpanningPenaltyMagnitude = AlignmentChain.maxScore / 2;
        immutable improperAlignmentPenaltyMagnitude = AlignmentChain.maxScore / 8;

        long numAlignments = readAlignment.length;
        long expectedAlignmentCount = pileUpType == ReadAlignmentType.gap ? 2 : 1;
        auto alignmentAnchor = readAlignment[].map!getRegion
            .fold!"a | b" - repetitiveRegions;
        long avgAnchorSize = alignmentAnchor.size / numAlignments;
        debug long avgAlignmentLength = readAlignment[].map!"a.totalLength".sum / numAlignments;
        debug assert(avgAnchorSize <= avgAlignmentLength);
        long avgAlignmentScore = readAlignment[].map!"a.score".sum / numAlignments;
        long shortAnchorPenalty = floor(shortAnchorPenaltyMagnitude * (
                (options.goodAnchorLength + 1) / avgAnchorSize.to!float) ^^ 2).to!size_t;
        // dfmt off
        long notSpanningPenalty = preferSpanning
            ? (expectedAlignmentCount - numAlignments) * notSpanningPenaltyMagnitude
            : 0;
        // dfmt on
        long improperAlignmentPenalty = readAlignment[].map!"a.isProper ? 0 : 1".sum
            * improperAlignmentPenaltyMagnitude / numAlignments;
        // dfmt off
        size_t score = max(0, (
              avgAlignmentScore
            - shortAnchorPenalty
            - notSpanningPenalty
            - improperAlignmentPenalty
        ));
        // dfmt on

        debug
        {
            size_t readId = readAlignment[0].contigB.id;
            auto contigIds = readAlignment[].map!"a.contigA.id".array;
            // dfmt off
            debug logJsonDebug(
                "readId", readId,
                "contigIds", contigIds.toJson,
                "expectedAlignmentCount", expectedAlignmentCount,
                "avgAnchorSize", avgAnchorSize,
                "avgAlignmentLength", avgAlignmentLength,
                "avgAlignmentScore", avgAlignmentScore,
                "shortAnchorPenalty", shortAnchorPenalty,
                "score", score,
            );
            // dfmt on
        }

        return score;
    }

    protected string buildConsensus(string croppedDb, size_t referenceReadId)
    {
        static struct ConsensusAlignmentOptions
        {
            string[] daccordOptions;
            string[] dalignerOptions;
            string[] dbsplitOptions;
            string workdir;
        }

        // dfmt off
        string consensusDb = getConsensus(croppedDb, referenceReadId, const(ConsensusAlignmentOptions)(
            options.daccordOptions,
            options.dalignerOptions ~ [
                format!(DalignerOptions.minAlignmentLength ~ "%d")(options.minAnchorLength),
                format!(DalignerOptions.averageCorrelationRate ~ "%f")((1 - options.readsErrorRate)^^2),
            ],
            options.dbsplitOptions,
            options.workdir,
        ));
        // dfmt on

        debug logJsonDebug("consensusDb", consensusDb);
        assert(getNumContigs(consensusDb, options) == 1);

        return consensusDb;
    }

    protected void writeNewAssembly()
    {
        logJsonDiagnostic("state", "enter", "function", "djunctor.writeNewAssembly");

        // dfmt off
        catHits = catHits
            .enforceJoinPolicy!InsertionInfo(options.joinPolicy)
            .normalizeUnkownJoins!InsertionInfo();
        // dfmt on

        foreach (startNode; contigStarts!InsertionInfo(catHits))
        {
            writeNewContig(startNode);
        }

        logJsonDiagnostic("state", "exit", "function", "djunctor.writeNewAssembly");
    }

    protected void writeNewContig(ContigNode startNode)
    {
        auto globalComplement = AlignmentChain.Complement.no;
        auto insertionBegin = startNode;

        writeHeader(startNode);
        size_t currentOutputColumn = 0;
        foreach (currentInsertion; linearWalk!InsertionInfo(catHits, startNode))
        {
            currentOutputColumn = writeInsertion(insertionBegin,
                    currentInsertion, globalComplement, currentOutputColumn);

            insertionBegin = currentInsertion.target(insertionBegin);
            if (currentInsertion.isAntiParallel)
            {
                globalComplement = cast(AlignmentChain.Complement) !globalComplement;
            }
        }
        writeln();
    }

    protected void writeHeader(in ContigNode begin) const
    {
        writeln(format!"> contig %d"(begin.contigId));
    }

    protected size_t writeInsertion(in ContigNode begin, in Insertion insertion,
            in AlignmentChain.Complement globalComplement, in size_t currentOutputColumn) const
    {
        assert(insertion.isValid, "invalid insertion");

        auto isContig = insertion.isDefault;
        auto isUnkownSequence = insertion.isUnkown;
        auto sequenceDb = insertion.payload.sequenceDb;
        auto sequenceId = isContig ? begin.contigId : 1;
        auto fastaRecord = getFastaRecord(sequenceDb, sequenceId, options);
        SpliceSite[] spliceSites = insertion.payload.spliceSites.dup;

        bool localComplement = false;
        alias Sequence = typeof(fastaRecord[].array);
        Sequence sequence;
        size_t contigSpliceStart;
        size_t contigSpliceEnd;

        if (isContig)
        {
            auto contigLength = fastaRecord.length;

            switch (spliceSites.length)
            {
            case 0:
                contigSpliceStart = 0;
                contigSpliceEnd = contigLength;
                break;
            case 1:
                auto splicePosition = spliceSites[0].croppingRefPosition.value;

                if (splicePosition < contigLength / 2)
                {
                    contigSpliceStart = splicePosition;
                    contigSpliceEnd = contigLength;
                }
                else
                {
                    contigSpliceStart = 0;
                    contigSpliceEnd = splicePosition;
                }
                break;
            case 2:
                assert(spliceSites.length == 2);
                assert(spliceSites[0].croppingRefPosition.contigId
                        == spliceSites[1].croppingRefPosition.contigId);

                contigSpliceStart = spliceSites[0].croppingRefPosition.value;
                contigSpliceEnd = spliceSites[1].croppingRefPosition.value;

                if (contigSpliceEnd < contigSpliceStart)
                {
                    swap(contigSpliceStart, contigSpliceEnd);
                }
                break;
            default:
                assert(0, "too many spliceSites");
            }

            sequence = fastaRecord[contigSpliceStart .. contigSpliceEnd].array;
            assert(globalComplement != (begin < insertion.target(begin)));
        }
        else if (isUnkownSequence)
        {
            assert(spliceSites.length == 2);
            auto sequenceLength = spliceSites[1].croppingRefPosition.value
                - spliceSites[0].croppingRefPosition.value;

            sequence = (cast(dchar) 'n').repeat.take(sequenceLength).array;
        }
        else
        {
            if (insertion.isExtension)
            {
                assert(spliceSites.length == 1);
            }
            else
            {
                assert(insertion.isGap);
                assert(spliceSites.length == 2);

                if (begin.contigId != spliceSites[0].croppingRefPosition.contigId)
                {
                    swap(spliceSites[0], spliceSites[1]);
                }
            }

            localComplement = spliceSites[0].complement;
            sequence = fastaRecord[].array;
        }

        bool effectiveComplement = globalComplement ^ localComplement;

        // dfmt off
        logJsonDebug(
            "info", "writing insertion",
            "isContig", isContig,
            "sequenceDb", sequenceDb,
            "sequenceId", sequenceId,
            "spliceSites", spliceSites.toJson,
            "localComplement", localComplement,
            "globalComplement", globalComplement,
            "effectiveComplement", effectiveComplement,
        );
        // dfmt on

        if (effectiveComplement)
        {
            sequence = reverseComplement(sequence);
        }

        // BEGIN wrapped output
        if (sequence.length == 0)
        {
            return currentOutputColumn;
        }

        writeln(sequence.take(options.fastaLineWidth - currentOutputColumn));

        foreach (sequenceChunk; sequence.drop(options.fastaLineWidth - currentOutputColumn)
                .chunks(options.fastaLineWidth))
        {
            if (sequenceChunk.walkLength == options.fastaLineWidth)
            {
                writeln(sequenceChunk);
            }
            else
            {
                write(sequenceChunk);
            }
        }

        return (currentOutputColumn + sequence.length) % options.fastaLineWidth;
        // END wrapped output
    }

    protected void writeUnusedReadsList()
    {
        auto unusedReads = this.unusedReads.elements;
        logJsonDebug("unusedReads", unusedReads.array.toJson);

        if (options.unusedReadsList is null)
        {
            return;
        }

        try
        {
            auto unusedReadsList = File(options.unusedReadsList, "w");

            unusedReadsList.write(unusedReads.array.toJson);
        }
        catch (ErrnoException e)
        {
            // dfmt off
            logJsonError(
                "info", "cannot write unused reads list",
                "error", e.to!string,
                "file", options.unusedReadsList,
                "unusedReads", unusedReads.toJson,
            );
            // dfmt on
        }
    }
}

interface AlignmentChainFilter
{
    AlignmentChain[] opCall(AlignmentChain[] alignmentChains);
}

abstract class ReadFilter : AlignmentChainFilter
{
    NaturalNumberSet* unusedReads;

    this(NaturalNumberSet* unusedReads)
    {
        this.unusedReads = unusedReads;
    }

    override AlignmentChain[] opCall(AlignmentChain[] alignmentChains)
    {
        // dfmt off
        auto discardedReadIds = getDiscardedReadIds(alignmentChains)
            .map!(ac => ac.contigB.id)
            .array
            .sort
            .uniq
            .array
            .assumeSorted;
        // dfmt on
        markReadsAsUsed(discardedReadIds);

        return alignmentChains.filter!(ac => !discardedReadIds.contains(ac.contigB.id)).array;
    }

    private void markReadsAsUsed(ReadIds)(ref ReadIds discardedReadIds)
    {
        foreach (discardedReadId; discardedReadIds)
        {
            unusedReads.remove(discardedReadId);
        }
    }

    InputRange!(AlignmentChain) getDiscardedReadIds(AlignmentChain[] alignmentChains);
}

/// Discard improper alignments.
class ImproperAlignmentChainsFilter : AlignmentChainFilter
{
    override AlignmentChain[] opCall(AlignmentChain[] alignmentChains)
    {
        return alignmentChains.filter!"a.isProper".array;
    }
}

/// Discard read if it has an alignment that – extended with the
/// exceeding read sequence on either end – is fully contained in
/// a single contig.
class RedundantAlignmentChainsFilter : ReadFilter
{
    this(NaturalNumberSet* unusedReads)
    {
        super(unusedReads);
    }

    override InputRange!(AlignmentChain) getDiscardedReadIds(AlignmentChain[] alignmentChains)
    {
        assert(alignmentChains.map!"a.isProper".all);

        return inputRangeObject(alignmentChains.filter!(ac => ac.isFullyContained));
    }
}

/// Discard read if it has an alignment that aligns in multiple
/// locations of one contig with approximately equal quality.
class AmbiguousAlignmentChainsFilter : ReadFilter
{
    this(NaturalNumberSet* unusedReads)
    {
        super(unusedReads);
    }

    override InputRange!(AlignmentChain) getDiscardedReadIds(AlignmentChain[] alignmentChains)
    {
        alias AlignmentsChunk = typeof(alignmentChains.chunkBy!haveEqualIds.front);

        bool isAmgiguouslyAlignedRead(AlignmentsChunk alignmentsChunk)
        {
            auto readAlignments = alignmentsChunk.array;
            readAlignments.sort!"a.score > b.score";

            // dfmt off
            return (
                readAlignments.length > 1 &&
                DJunctor.similarScore(readAlignments[0].score, readAlignments[1].score)
            );
            // dfmt on
        }

        assert(alignmentChains.map!"a.isProper".all);

        // dfmt off
        return alignmentChains
            .chunkBy!haveEqualIds
            .filter!isAmgiguouslyAlignedRead
            .joiner
            .inputRangeObject;
        // dfmt on
    }
}

class WeaklyAnchoredAlignmentChainsFilter : AlignmentChainFilter
{
    size_t minAnchorLength;
    const(ReferenceMask) repetitiveRegions;

    this(const(ReferenceMask) repetitiveRegions, size_t minAnchorLength)
    {
        this.repetitiveRegions = repetitiveRegions;
        this.minAnchorLength = minAnchorLength;
    }

    AlignmentChain[] opCall(AlignmentChain[] alignmentChains)
    {
        bool isStronglyAnchored(AlignmentChain alignment)
        {
            // Mark reads as weakly anchored if they are mostly anchored in a
            // repetitive region
            auto uniqueAlignmentRegion = getRegion(alignment) - repetitiveRegions;
            bool isWeaklyAnchored = uniqueAlignmentRegion.size <= minAnchorLength;

            return !isWeaklyAnchored;
        }

        // dfmt off
        return alignmentChains
            .filter!isStronglyAnchored
            .array;
        // dfmt on
    }
}

interface MaskAssessor(T)
{
    ReferenceMask opCall(T input);
}

alias RepeatAssessor = MaskAssessor!(const(AlignmentChain[]));

/**
    Mask reference regions where the alignment coverage is not within set limits. This helps to
    identify repetitive or bad quality regions.
*/
class BadAlignmentCoverageAssessor : RepeatAssessor
{
    double lowerLimit;
    double upperLimit;

    /// Create an assessor with these limits.
    this(double lowerLimit, double upperLimit)
    {
        this.lowerLimit = lowerLimit;
        this.upperLimit = upperLimit;
    }

    version (unittest)
    {
        /**
            Generate a list of alignment chains for testing:

            ```
                                 c1                             c2                    c3
                   0    5   10   15   20   25   30       0    5   10   15      0    5   10   15
            ref:   [-----------------------------)       [--------------)      [--------------)
            reads: .    .    .    .    .    .    .       .    .    .    .      .    .    .    .
                   . #1 [------------) .    .    .   #12 [--) .    .    .  #23 .[--).    .    .
                   . #2 [------------) .    .    .   #13 [--) .    .    .  #24 . [--)    .    .
                   . #3 [--------------)    .    .   #14 [----)    .    .  #25 .  [--)   .    .
                   .    . #4 [---------)    .    .   #15 [----)    .    .  #26 .   [--)  .    .
                   .    . #5 [-------------------)   #16 [--------------)  #27 .    [--) .    .
                   .    . #6 [-------------------)   #17 [--------------)  #28 .    .[--).    .
                   .    .    #7 [----------------)   #18 [--------------)  #29 .    . [--)    .
                   .    .    .    . #8 [---------)       .#19 [---------)  #30 .    .  [--)   .
                   .    .    .    . #9 [---------)       .#20 [---------)  #31 .    .   [--)  .
                   .    .    .    .#10 [---------)       .#21 [---------)  #32 .    .    [--) .
                   .    .    .    .    #11 [-----)       .    #22 [-----)  #33 .    .    .[--).
                   .    .    .    .    .    .    .       .    .    .    .      .    .    .    .
            mask:  [====)    [=======) [=========)       [==) [=========)      [==) .    . [==)
                   .    .    .    .    .    .    .       .    .    .    .      .    .    .    .
            cov:   ^    .    .    .    .    .    .       .    .    .    .      .    .    .    .
                   |    .    .    .    .    .    .       .    .    .    .      .    .    .    .
                   |    .    .  +----+ .   +----+.       +-+  .   +----+.      .    .    .    .
                   |    .    +--+####| +---+####|.       |#|  +---+####|.      .    .    .    .
                 5 |.........|#######+-+########|.       |#+--+########|.      .    .    .    .
                   |    .    |##################|.       |#############|.      .    .    .    .
                   |    +----+##################|.       |#############|.      .  +-------+   .
                   |    |#######################|.       |#############|.      . ++#######++  .
                   |    |#######################|.       |#############|.      .++#########++ .
                 0 +----+-----------------------+--------+-------------+--------+-----------+---->
            ```
        */
        private static AlignmentChain[] getTestAlignments()
        {
            with (AlignmentChain) with (LocalAlignment)
                {
                    size_t alignmentChainId = 0;
                    size_t contReadId = 0;
                    AlignmentChain getDummyAlignment(size_t contigId,
                            size_t contigLength, size_t beginIdx, size_t endIdx)
                    {
                        // dfmt off
                        return AlignmentChain(
                            ++alignmentChainId,
                            Contig(contigId, contigLength),
                            Contig(alignmentChainId, endIdx),
                            Complement.no,
                            [
                                LocalAlignment(
                                    Locus(beginIdx, beginIdx + 1),
                                    Locus(0, 1),
                                    0,
                                ),
                                LocalAlignment(
                                    Locus(endIdx - 1, endIdx),
                                    Locus(2, 3),
                                    0,
                                ),
                            ],
                        );
                        // dfmt on
                    }

                    // dfmt off
                    return [
                        getDummyAlignment(1, 30,  5, 18), //  #1
                        getDummyAlignment(1, 30,  5, 18), //  #2
                        getDummyAlignment(1, 30,  5, 20), //  #3
                        getDummyAlignment(1, 30, 10, 20), //  #4
                        getDummyAlignment(1, 30, 10, 30), //  #5
                        getDummyAlignment(1, 30, 10, 30), //  #6
                        getDummyAlignment(1, 30, 13, 30), //  #7
                        getDummyAlignment(1, 30, 20, 30), //  #8
                        getDummyAlignment(1, 30, 20, 30), //  #9
                        getDummyAlignment(1, 30, 20, 30), // #10
                        getDummyAlignment(1, 30, 24, 30), // #11
                        getDummyAlignment(2, 15,  0,  3), // #12
                        getDummyAlignment(2, 15,  0,  3), // #13
                        getDummyAlignment(2, 15,  0,  5), // #14
                        getDummyAlignment(2, 15,  0,  5), // #15
                        getDummyAlignment(2, 15,  0, 15), // #16
                        getDummyAlignment(2, 15,  0, 15), // #17
                        getDummyAlignment(2, 15,  0, 15), // #18
                        getDummyAlignment(2, 15,  5, 15), // #19
                        getDummyAlignment(2, 15,  5, 15), // #20
                        getDummyAlignment(2, 15,  5, 15), // #21
                        getDummyAlignment(2, 15,  9, 15), // #22
                        getDummyAlignment(3, 15,  1,  4), // #23
                        getDummyAlignment(3, 15,  2,  5), // #24
                        getDummyAlignment(3, 15,  3,  6), // #25
                        getDummyAlignment(3, 15,  4,  7), // #26
                        getDummyAlignment(3, 15,  5,  8), // #27
                        getDummyAlignment(3, 15,  6,  9), // #28
                        getDummyAlignment(3, 15,  7, 10), // #29
                        getDummyAlignment(3, 15,  8, 11), // #30
                        getDummyAlignment(3, 15,  9, 12), // #31
                        getDummyAlignment(3, 15, 10, 13), // #32
                        getDummyAlignment(3, 15, 11, 14), // #33
                    ];
                    // dfmt on
                }
        }
    }

    /// Apply the assessor to the given set of alignment.
    override ReferenceMask opCall(const(AlignmentChain[]) alignments)
    {
        if (alignments.length == 0)
        {
            return ReferenceMask();
        }

        static immutable OK = CoverageZone.ok;
        auto maskAcc = appender!(ReferenceInterval[]);
        auto masker = Masker();
        auto changeEvents = coverageChanges(alignments);
        auto lastEvent = changeEvents.front;

        foreach (event; changeEvents)
        {
            auto currentZone = coverageZone(event.currentCoverage);
            auto newZone = coverageZone(event.newCoverage);

            if (masker.isMasking && event.contigId != lastEvent.contigId)
            {
                maskAcc ~= masker.finish(lastEvent.position);
            }

            // dfmt off
            if (
                !masker.isMasking &&
                (newZone != OK || (currentZone == OK && currentZone != newZone))
            )
            // dfmt on
            {
                masker.start(event.contigId, event.position);
            }
            else if (masker.isMasking && currentZone != OK && newZone == OK)
            {
                maskAcc ~= masker.finish(event.position);
            }

            lastEvent = event;
        }

        if (masker.isMasking)
        {
            maskAcc ~= masker.finish(lastEvent.position);
        }

        return ReferenceMask(maskAcc.data);
    }

    ///
    unittest
    {
        auto alignments = getTestAlignments();
        alias CoverageChange = CoverageChangeRange.CoverageChange;

        auto assessor = new BadAlignmentCoverageAssessor(3, 5);
        // dfmt off
        assert(assessor(alignments) == ReferenceMask([
            ReferenceInterval(1,  0,  5),
            ReferenceInterval(1, 10, 18),
            ReferenceInterval(1, 20, 30),
            ReferenceInterval(2,  0,  3),
            ReferenceInterval(2,  5, 15),
            ReferenceInterval(3,  0,  3),
            ReferenceInterval(3, 12, 15),
        ]));
        // dfmt on
    }

    private static enum CoverageZone
    {
        low,
        ok,
        high,
    }

    private CoverageZone coverageZone(in int coverage) const pure nothrow
    {
        // dfmt off
        return coverage < lowerLimit
            ? CoverageZone.low
            : coverage > upperLimit
                ? CoverageZone.high
                : CoverageZone.ok;
        // dfmt on
    }

    private static struct Masker
    {
        private bool _isMasking = false;
        private size_t contigId;
        private size_t maskStart;

        void start(in size_t contigId, in size_t maskStart) pure nothrow
        {
            this._isMasking = true;
            this.contigId = contigId;
            this.maskStart = maskStart;
        }

        ReferenceInterval finish(in size_t maskEnd) pure nothrow
        {
            this._isMasking = false;

            // dfmt off
            return ReferenceInterval(
                this.contigId,
                this.maskStart,
                maskEnd,
            );
            // dfmt on
        }

        @property bool isMasking() const pure nothrow
        {
            return this._isMasking;
        }
    }

    /// Transforms a range of alignment chains into a range of coverage
    /// change events.
    private struct CoverageChangeRange
    {
        // dfmt off
        static alias AlignmentEvent = Tuple!
        (
            size_t, "contigId",
            size_t, "position",
            int, "diff",
        );
        // dfmt on

        static struct CoverageChange
        {
            size_t contigId;
            size_t position;
            int currentCoverage;
            int newCoverage;

            bool hasChanged() const pure nothrow
            {
                return currentCoverage != newCoverage;
            }
        }

        AlignmentEvent[] alignmentEvents;
        size_t currentEventIdx = 0;
        CoverageChange _front;

        private this(AlignmentEvent[] alignmentEvents)
        {
            this.alignmentEvents = alignmentEvents;
            if (!empty)
            {
                // Advance to first event.
                popFront();
            }
        }

        private static CoverageChangeRange create(Range)(Range alignments)
                if (isInputRange!Range && is(ElementType!Range : const(AlignmentChain)))
        {
            if (alignments.length == 0)
            {
                return CoverageChangeRange();
            }

            // dfmt off
            auto alignmentEvents = alignments
                .map!(alignment => only(
                    AlignmentEvent(alignment.contigA.id, alignment.first.contigA.begin, 1),
                    AlignmentEvent(alignment.contigA.id, alignment.last.contigA.end, -1),
                ))
                .joiner;
            auto contigBoundaryEvents = alignments
                .map!"a.contigA"
                .uniq
                .map!(contig => only(
                    AlignmentEvent(contig.id, 0, 0),
                    AlignmentEvent(contig.id, contig.length, 0),
                ))
                .joiner;
            // dfmt on
            auto changeEvents = chain(alignmentEvents, contigBoundaryEvents).array;
            changeEvents.sort();

            return CoverageChangeRange(changeEvents);
        }

        void popFront()
        {
            _front.currentCoverage = _front.newCoverage;

            if (currentEventIdx == alignmentEvents.length)
            {
                // There is no more data; just pop the front and return;
                assert(_front.currentCoverage == 0, "coverage should drop to zero in the end");
                ++currentEventIdx;

                return;
            }

            auto currentEvent = alignmentEvents[currentEventIdx];
            auto eventAcc = currentEvent;
            eventAcc.diff = 0;

            // Collect all events for one position.
            for (; currentEventIdx < alignmentEvents.length; ++currentEventIdx)
            {
                currentEvent = alignmentEvents[currentEventIdx];

                // dfmt off
                if (
                    eventAcc.contigId != currentEvent.contigId ||
                    eventAcc.position != currentEvent.position
                )
                // dfmt on
                {
                    break;
                }

                eventAcc.diff += currentEvent.diff;
            }

            assert(eventAcc.contigId == _front.contigId || _front.currentCoverage == 0,
                    "coverage should drop to zero between contigs");
            _front.contigId = eventAcc.contigId;
            _front.position = eventAcc.position;
            _front.newCoverage = _front.currentCoverage + eventAcc.diff;
        }

        @property bool empty()
        {
            return currentEventIdx > alignmentEvents.length;
        }

        @property CoverageChange front()
        {
            return _front;
        }
    }

    static CoverageChangeRange coverageChanges(Range)(Range alignments)
    {
        return CoverageChangeRange.create(alignments);
    }

    unittest
    {
        auto alignments = getTestAlignments();

        alias CoverageChange = CoverageChangeRange.CoverageChange;

        // dfmt off
        assert(coverageChanges(alignments).equal([
            CoverageChange(1,  0, 0, 0),
            CoverageChange(1,  5, 0, 3),
            CoverageChange(1, 10, 3, 6),
            CoverageChange(1, 13, 6, 7),
            CoverageChange(1, 18, 7, 5),
            CoverageChange(1, 20, 5, 6),
            CoverageChange(1, 24, 6, 7),
            CoverageChange(1, 30, 7, 0),
            CoverageChange(2,  0, 0, 7),
            CoverageChange(2,  3, 7, 5),
            CoverageChange(2,  5, 5, 6),
            CoverageChange(2,  9, 6, 7),
            CoverageChange(2, 15, 7, 0),
            CoverageChange(3,  0, 0, 0),
            CoverageChange(3,  1, 0, 1),
            CoverageChange(3,  2, 1, 2),
            CoverageChange(3,  3, 2, 3),
            CoverageChange(3,  4, 3, 3),
            CoverageChange(3,  5, 3, 3),
            CoverageChange(3,  6, 3, 3),
            CoverageChange(3,  7, 3, 3),
            CoverageChange(3,  8, 3, 3),
            CoverageChange(3,  9, 3, 3),
            CoverageChange(3, 10, 3, 3),
            CoverageChange(3, 11, 3, 3),
            CoverageChange(3, 12, 3, 2),
            CoverageChange(3, 13, 2, 1),
            CoverageChange(3, 14, 1, 0),
            CoverageChange(3, 15, 0, 0),
        ]));
        // dfmt on
    }
}
