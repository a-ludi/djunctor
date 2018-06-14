/**
    This is the main algorithm of this package.

    Copyright: © 2018 Arne Ludwig <arne.ludwig@posteo.de>
    License: Subject to the terms of the MIT license, as written in the
             included LICENSE file.
    Authors: Arne Ludwig <arne.ludwig@posteo.de>
*/
module djunctor.djunctor;

import djunctor.alignments : AlignmentChain, alignmentCoverage, buildPileUps, getAlignmentChainRefs,
    getType, haveEqualIds, isGap, isValid, PileUp, ReadAlignment, ReadAlignmentType;
import djunctor.commandline : Options;
import djunctor.dazzler : buildDamFile, getConsensus, getFastaEntries,
    getLocalAlignments, getMappings, getNumContigs, getTracePointDistance,
    attachTracePoints, writeMask;
import djunctor.util.fasta : parseFasta, parseFastaRecord, parsePacBioHeader,
    reverseComplement;
import djunctor.util.log;
import djunctor.util.math : absdiff, ceil, floor, mean, median;
import djunctor.util.range : Comparator;
import djunctor.util.region : empty, Region, union_;
import djunctor.util.string : indent;
import dstats.distrib : invPoissonCDF;
import core.exception : AssertError;
import std.algorithm : all, any, cache, canFind, chunkBy, each, equal, filter,
    find, fold, group, isSorted, joiner, map, max, maxElement, min,
    multiwayUnion, setDifference, sort, sum, swap, SwapStrategy, uniq;
import std.array : appender, Appender, array, join;
import std.container : BinaryHeap, heapify, make;
import std.conv;
import std.exception : assertNotThrown, assertThrown;
import std.format : format, formattedWrite;
import std.math : abs, floor, sgn;
import std.range : assumeSorted, chain, chunks, ElementType, enumerate,
    ForwardRange, InputRange, inputRangeObject, iota, isForwardRange,
    isInputRange, only, retro, slide, SortedRange, tail, take, walkLength, zip;
import std.stdio : File, write, writeln;
import std.string : outdent;
import std.typecons : Flag, No, tuple, Tuple, Yes;
import vibe.data.json : Json, toJson = serializeToJson;

/// Start the `djunctor` algorithm with preprocessed options.
void runWithOptions(in ref Options options)
{
    new DJunctor(options).run();
}

alias ReferenceMask = Region!(size_t, size_t, "contigId");
alias ReferenceInterval = ReferenceMask.TaggedInterval;

/// Returns the alignment region of alignmentChain.
ReferenceMask getRegion(in AlignmentChain alignmentChain) pure nothrow
{
    // dfmt off
    return ReferenceMask(
        alignmentChain.contigA.id,
        alignmentChain.first.contigA.begin,
        alignmentChain.last.contigA.end,
    );
    // dfmt on
}

// dfmt off
alias CroppingSlice = Tuple!(
    size_t, "begin",
    size_t, "end",
);
alias CroppingRefPosition = Tuple!(
    size_t, "frontIdx",
    size_t, "frontContigId",
    size_t, "backIdx",
    size_t, "backContigId",
);
// dfmt on

CroppingSlice intersection(in CroppingSlice aSlice, in CroppingSlice bSlice) pure nothrow
{
    return CroppingSlice(max(aSlice[0], bSlice[0]), min(aSlice[1], bSlice[1]));
}

CroppingRefPosition getCroppingRefPositions(const PileUp pileUp, in size_t tracePointDistance) pure nothrow
{
    // dfmt off
    auto frontAlignments = pileUp
        .filter!(readAlignment => readAlignment.isFrontExtension || readAlignment.isGap)
        .map!"a[$ - 1]"
        .array;
    auto backAlignments = pileUp
        .filter!(readAlignment => readAlignment.isBackExtension || readAlignment.isGap)
        .map!"a[0]"
        .array;
    // dfmt on
    auto commonFrontTracePoint = frontAlignments.getCommonFrontTracePoint(tracePointDistance);
    auto commonBackTracePoint = backAlignments.getCommonBackTracePoint(tracePointDistance);
    size_t frontContig;
    size_t backContig;

    if (frontAlignments.length > 0 && backAlignments.length > 0)
    {
        frontContig = frontAlignments[0].contigA.id;
        backContig = backAlignments[0].contigA.id;
    }
    else if (frontAlignments.length > 0)
    {
        frontContig = backContig = frontAlignments[0].contigA.id;
    }
    else if (backAlignments.length > 0)
    {
        backContig = frontContig = backAlignments[0].contigA.id;
    }
    else
    {
        assert(0, "empty pile up");
    }

    // dfmt off
    return CroppingRefPosition(
        commonFrontTracePoint,
        frontContig,
        commonBackTracePoint,
        backContig,
    );
    // dfmt on
}

CroppingSlice getCroppingSlice(const AlignmentChain alignmentChain,
        in CroppingRefPosition croppingRefPosition)
{
    auto alignment = const(ReadAlignment)(alignmentChain);
    auto alignmentType = alignment.type;
    auto tracePointDistance = alignmentChain.tracePointDistance;
    // dfmt off
    auto refPos = alignmentType == ReadAlignmentType.front
        ? croppingRefPosition.frontIdx
        : croppingRefPosition.backIdx;
    // dfmt on
    auto openInterval = alignmentType == ReadAlignmentType.front;
    size_t readCroppingPos;

    foreach (localAlignment; alignmentChain.localAlignments)
    {
        auto firstTracePointRefPos = ceil(localAlignment.contigA.begin, tracePointDistance);
        auto numTracePoints = localAlignment.tracePoints.length;
        assert(refPos >= firstTracePointRefPos);
        auto tracePointIdx = (refPos - firstTracePointRefPos) / tracePointDistance;

        if (tracePointIdx < numTracePoints)
        {
            auto endIdx = tracePointIdx + (openInterval ? 0 : 1);
            // dfmt off
            readCroppingPos =
                localAlignment.contigB.begin +
                localAlignment
                    .tracePoints[0 .. endIdx]
                    .map!"a.numBasePairs"
                    .sum;
            // dfmt on
            break;
        }
    }

    // dfmt off
    return alignmentType == ReadAlignmentType.front
        ? CroppingSlice(0, readCroppingPos)
        : CroppingSlice(readCroppingPos, alignmentChain.contigB.length);
    // dfmt on
}

/**
    Returns a back "common" trace points wrt. contigB. The returned
    trace point is not necessarily common to *all* alignment chains, ie. the
    alignment intervals are incrementally intersected reject those that would
    cause an empty intersection. The intervals are ordered descending by
    end point and starting point. Example:

    ---
    :....:....:....:....:....:....:....:....:
    :    :    :    :    |----+--------------|
    :    :    :    |---------+--------------|
    :    :    :    :    :  |-+---------|    :
    :    :    :    :    :   |+--------|:    :
    :    :    :    :    :|---+-------| :    :
    :    : |-----------|:    :    :    :    :
    |---------|    :    :    :    :    :    :

     :  = trace point
     +  = "common" trace point
    |-| = alignment interval
    ---
*/
size_t getCommonBackTracePoint(const AlignmentChain[] alignmentChains, in size_t tracePointDistance) pure nothrow
{
    // dfmt off
    auto alignmentSlices = alignmentChains
        .map!(ac => tuple(
            ac.last.contigA.end + 0,
            ac.first.contigA.begin + 0,
        ))
        .array;
    // dfmt on

    foreach (alignmentSlice; alignmentSlices)
    {
        assert(alignmentSlice[0] > alignmentSlice[1]);
    }
    // dfmt off
    debug logJsonDebug(
        "alignmentSlices", alignmentSlices.map!"[a[0], a[1]]".array.toJson,
        "contigAIds", alignmentChains.map!"a.contigA.id".array.toJson,
        "type", "back",
    );
    // dfmt on

    return getCommonBackTracePoint(alignmentSlices, tracePointDistance);
}

private size_t getCommonBackTracePoint(Tuple!(size_t,
        size_t)[] alignmentSlices, in size_t tracePointDistance) pure nothrow
{
    if (alignmentSlices.length == 0)
        return 0;

    alignmentSlices.sort!"a > b";
    auto commonSlice = tuple(0UL, size_t.max);
    auto lastCommonSlice = commonSlice;

    foreach (alignmentSlice; alignmentSlices)
    {
        commonSlice[0] = ceil(max(commonSlice[0], alignmentSlice[1]), tracePointDistance);
        commonSlice[1] = floor(min(commonSlice[1], alignmentSlice[0]), tracePointDistance);

        if (commonSlice[1] < commonSlice[0])
        {
            break;
        }

        lastCommonSlice = commonSlice;
    }

    return lastCommonSlice[0];
}

unittest
{
    immutable tracePointDistance = 5UL;
    // dfmt off
    // 0    5   10   15   20   25
    // :....:....:....:....:....:
    // |---------+--------------|
    //      |----+--------------|
    //         |-+---------|
    //       |---+-------|
    //          |+--------|
    auto alignmentSlices1 = [
        tuple(25UL,  0UL),
        tuple(25UL,  5UL),
        tuple(20UL,  8UL),
        tuple(18UL,  6UL),
        tuple(19UL,  9UL),
    ];
    // dfmt on
    auto commonBackTracePoint1 = getCommonBackTracePoint(alignmentSlices1, tracePointDistance);

    assert(commonBackTracePoint1 == 10);

    // dfmt off
    // 0    5   10   15   20   25
    // :....:....:....:....:....:
    // |-------------|
    //               |+----|
    //            |---+-------|
    //              |-+---------|
    auto alignmentSlices2 = [
        tuple(14UL,  0UL),
        tuple(20UL, 14UL),
        tuple(23UL, 11UL),
        tuple(25UL, 13UL),
    ];
    // dfmt on
    auto commonBackTracePoint2 = getCommonBackTracePoint(alignmentSlices2, tracePointDistance);

    assert(commonBackTracePoint2 == 15);

    Tuple!(size_t, size_t)[] alignmentSlices3;
    auto commonBackTracePoint3 = getCommonBackTracePoint(alignmentSlices3, tracePointDistance);

    assert(commonBackTracePoint3 == 0);
}

/**
    Returns a back "common" trace points wrt. contigB. The returned
    trace point is not necessarily common to *all* alignment chains, ie. the
    alignment intervals are incrementally intersected reject those that would
    cause an empty intersection. The intervals are ordered ascending by
    starting point and end point. Example:

    ---
    :....:....:....:....:....:....:....:....:
    |--------------+----|    :    :    :    :
    |--------------+---------|    :    :    :
    :    |---------+-|  :    :    :    :    :
    :    :|--------+|   :    :    :    :    :
    :    : |-------+---|:    :    :    :    :
    :    :    :    :    :|-----------| :    :
    :    :    :    :    :    :    |---------|

     :  = trace point
     +  = "common" trace point
    |-| = alignment interval
    ---
*/
size_t getCommonFrontTracePoint(const AlignmentChain[] alignmentChains, in size_t tracePointDistance) pure nothrow
{
    // dfmt off
    auto alignmentSlices = alignmentChains
        .map!(ac => tuple(
            ac.first.contigA.begin + 0,
            ac.last.contigA.end + 0,
        ))
        .array;
    // dfmt on

    foreach (alignmentSlice; alignmentSlices)
    {
        assert(alignmentSlice[0] < alignmentSlice[1]);
    }
    // dfmt off
    debug logJsonDebug(
        "alignmentSlices", alignmentSlices.map!"[a[0], a[1]]".array.toJson,
        "contigAIds", alignmentChains.map!"a.contigA.id".array.toJson,
        "type", "front",
    );
    // dfmt on

    return getCommonFrontTracePoint(alignmentSlices, tracePointDistance);
}

private size_t getCommonFrontTracePoint(Tuple!(size_t,
        size_t)[] alignmentSlices, in size_t tracePointDistance) pure nothrow
{
    if (alignmentSlices.length == 0)
        return size_t.max;

    alignmentSlices.sort!"a < b";
    auto commonSlice = tuple(0UL, size_t.max);
    auto lastCommonSlice = commonSlice;

    foreach (alignmentSlice; alignmentSlices)
    {
        commonSlice[0] = max(commonSlice[0], alignmentSlice[0]);
        commonSlice[0] = ceil(commonSlice[0], tracePointDistance);
        commonSlice[1] = min(commonSlice[1], alignmentSlice[1]);
        commonSlice[1] = floor(commonSlice[1], tracePointDistance);

        if (commonSlice[1] < commonSlice[0])
        {
            break;
        }

        lastCommonSlice = commonSlice;
    }

    return lastCommonSlice[1];
}

///
unittest
{
    immutable tracePointDistance = 5UL;
    // dfmt off
    // 0    5   10   15   20   25
    // :....:....:....:....:....:
    // |--------------+---------|
    // |--------------+----|
    //      |---------+-|
    //        |-------+---|
    //       |--------+|
    auto alignmentSlices1 = [
        tuple( 0UL, 25UL),
        tuple( 0UL, 20UL),
        tuple( 5UL, 17UL),
        tuple( 7UL, 19UL),
        tuple( 6UL, 16UL),
    ];
    // dfmt on
    auto commonFrontTracePoint1 = getCommonFrontTracePoint(alignmentSlices1, tracePointDistance);

    assert(commonFrontTracePoint1 == 15);

    // dfmt off
    // 0    5   10   15   20   25
    // :....:....:....:....:....:
    //            |-------------|
    //      |----+|
    //   |-------+---|
    // |---------+-|
    auto alignmentSlices2 = [
        tuple(11UL, 25UL),
        tuple( 5UL, 11UL),
        tuple( 2UL, 14UL),
        tuple( 0UL, 12UL),
    ];
    // dfmt on
    auto commonFrontTracePoint2 = getCommonFrontTracePoint(alignmentSlices2, tracePointDistance);

    assert(commonFrontTracePoint2 == 10);

    Tuple!(size_t, size_t)[] alignmentSlices3;
    auto commonFrontTracePoint3 = getCommonFrontTracePoint(alignmentSlices3, tracePointDistance);

    assert(commonFrontTracePoint3 == size_t.max);
}

/// This characterizes an insertion.
// dfmt off
alias Hit = Tuple!(
    CoordinateTransform.Insertion, "insertion",
    size_t, "readId",
    AlignmentChain.Complement, "complement",
    string, "dbFile",
);
// dfmt on

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
    alias numEstimateUseless = (numCatUnsure, iteration) => (numCatUnsure - numCatUnsure / 20) / (
            iteration + 1);
    /// Do not insert extensions that are improbable short after consensus.
    static immutable minExtensionLength = 100;

    size_t numReferenceContigs;
    size_t numReads;
    AlignmentChain[] selfAlignment;
    ReferenceMask repetitiveRegions;
    AlignmentChain[] readsAlignment;
    const Options options;
    /// Set of read ids not to be considered in further processing.
    size_t[] catUseless;
    /// Set of alignments to be considered in further processing.
    AlignmentChain[] catCandidates;
    /// Set of alignments to be used for gap filling.
    Hit[] catHits;
    /// Access contigs distributed over several DBs by ID.
    DBUnion referenceDb;
    /// Transform coordinates before insertions to coordinates after them.
    CoordinateTransform coordTransform;
    alias T = coordTransform;
    size_t iteration;

    this(in ref Options options)
    {
        this.options = options;
        this.catUseless = [];
        this.catHits = [];
        this.iteration = 0;
        this.referenceDb.baseDb = options.refDb;
    }

    DJunctor run()
    {
        logJsonDiagnostic("state", "enter", "function", "run");
        // dfmt off
        this
            .init
            .mainLoop
            .finish;
        // dfmt on
        logJsonDiagnostic("state", "exit", "function", "run");

        return this;
    }

    protected DJunctor init()
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
        // dfmt off
        selfAlignment = getLocalAlignments(options.refDb, const(SelfAlignmentOptions)(
            options.dalignerOptions ~ [
                format!"-l%d"(options.minAnchorLength),
                format!"-e%f"((1 - options.referenceErrorRate)^^2),
            ],
            options.ladumpOptions[].filter!"a != \"-o\"".array,
            options.workdir,
        ));
        readsAlignment = getMappings(options.refDb, options.readsDb, const(ReadsToRefAlignmentOptions)(
            options.damapperOptions ~ [
                format!"-e%f"((1 - options.referenceErrorRate) * (1 - options.readsErrorRate)),
            ],
            options.ladumpOptions[].array,
            options.workdir,
        ));
        // dfmt on
        catCandidates = readsAlignment.dup;
        logJsonDiagnostic("state", "exit", "function", "djunctor.init");

        return this;
    }

    protected DJunctor mainLoop()
    {
        logJsonDiagnostic("state", "enter", "function", "djunctor.mainLoop");
        do
        {
            assessRepeatStructure();
            filterAlignments();
            findHits();

            // dfmt off
            logJsonDiagnostic(
                "iteration", iteration,
                "uselessReads", catUseless.toJson,
                "numCandiates", catCandidates.length,
                "numHits", catHits.length
            );
            // dfmt on

            if (catHits.length > 0)
            {
                insertHits();
            }

            ++iteration;
        }
        while (catHits.length > 0 && iteration < maxLoops);
        logJsonDiagnostic("state", "exit", "function", "djunctor.mainLoop");

        return this;
    }

    protected DJunctor assessRepeatStructure()
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

            if (shouldLog(LogLevel.diagnostic) && options.repeatMask != null)
            {
                auto maskName = format!"%s-%s"(options.repeatMask, stage.name);

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

        if (options.repeatMask != null)
        {
            writeMask(options.refDb, options.repeatMask,
                    this.repetitiveRegions.intervals, options);
        }

        logJsonDiagnostic("state", "exit", "function", "djunctor.assessRepeatStructure");

        return this;
    }

    protected DJunctor filterAlignments()
    {
        logJsonDiagnostic("state", "enter", "function", "djunctor.filterReads");

        // dfmt off
        auto filters = tuple(
            new WeaklyAnchoredAlignmentChainsFilter(repetitiveRegions, options.minAnchorLength),
            new ImproperAlignmentChainsFilter(),
            new AmbiguousAlignmentChainsFilter(),
            new RedundantAlignmentChainsFilter(),
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

        return this;
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

    protected DJunctor findHits()
    {
        logJsonDiagnostic("state", "enter", "function", "djunctor.findHits");

        // dfmt off
        auto pileUps = buildPileUps(catCandidates)
            .filter!(pileUp => pileUp.length >= options.minReadsPerPileUp)
            .array;
        // dfmt on

        logFillingInfo!("findHits", "raw")(pileUps);
        // dfmt off
        logJsonDebug("pileUps", pileUps
            .map!(pileUp => Json([
                "type": Json(pileUp.getType.to!string),
                "readAlignments": pileUp.map!"a[]".array.toJson,
            ]))
            .array
            .toJson);
        // dfmt on

        foreach (pileUp; pileUps)
        {
            assert(pileUp.isValid, "ambiguous pile up");

            auto croppedDbResult = getCroppedPileUpDb(pileUp);
            auto hit = buildConsensus(croppedDbResult);

            // Mark best scoring read for later sequence insertion.
            catHits ~= hit;
            logJsonDebug("hit", hit.toJson);
        }

        logJsonDiagnostic("state", "exit", "function", "djunctor.findHits");

        return this;
    }

    protected Hit getCroppedPileUpDb(PileUp pileUp)
    {
        auto pileUpType = pileUp.getType;
        fetchTracePoints(pileUp);
        auto tracePointDistance = pileUp[0][0].tracePointDistance;
        auto croppingRefPositions = pileUp.getCroppingRefPositions(tracePointDistance);
        logJsonDebug("croppingRefPositions", croppingRefPositions.toJson);
        size_t[] readIds = pileUp.map!"a[0].contigB.id + 0".array;
        // dfmt off
        auto rawFastaEntries = zip(
            readIds,
            getFastaEntries(options.readsDb, readIds, options).map!parseFastaRecord,
        ).array.assumeSorted!"a[0] < b[0]";
        // dfmt on
        auto croppedFastaEntries = appender!(string[]);
        croppedFastaEntries.reserve(rawFastaEntries.length);
        size_t bestInsertionScore;
        size_t bestInsertionLength;
        size_t bestInsertionIdx;
        AlignmentChain.Complement bestInsertionComplement;

        foreach (croppedDbIdx, readAlignment; pileUp)
        {
            immutable dummyFastaEntry = typeof(rawFastaEntries.front[1])();
            auto readId = readAlignment[0].contigB.id;
            auto rawFastaEntry = rawFastaEntries.equalRange(tuple(readId,
                    dummyFastaEntry)).front[1];
            immutable lineSep = typeof(rawFastaEntry).lineSep;
            auto readCroppingSlice = readAlignment[].map!(
                    alignmentChain => alignmentChain.getCroppingSlice(croppingRefPositions))
                .fold!((aSlice, bSlice) => intersection(aSlice, bSlice));
            auto insertionLength = readCroppingSlice[1] - readCroppingSlice[0];

            // dfmt off
            debug logJsonDebug(
                "readId", readId,
                "complement", readAlignment[0].complement,
                "rawReadCroppingSlices", readAlignment[].map!(alignmentChain => alignmentChain.getCroppingSlice(croppingRefPositions)).array.toJson,
                "readCroppingSlice", readCroppingSlice.toJson,
                "croppedReadLength", readCroppingSlice[1] - readCroppingSlice[0],
            );
            // dfmt on

            if (readAlignment[0].complement)
            {
                auto readLength = readAlignment[0].contigB.length;

                foreach (ref endpoint; readCroppingSlice)
                    endpoint = readLength - endpoint;
                swap(readCroppingSlice.expand);
            }
            assert(readCroppingSlice[0] < readCroppingSlice[1], "invalid/empty read cropping slice");

            size_t insertionScore = getInsertionScore(readAlignment,
                    pileUpType, Yes.preferSpanning);

            if (insertionScore > bestInsertionScore)
            {
                bestInsertionScore = insertionScore;
                bestInsertionLength = insertionLength;
                bestInsertionIdx = croppedDbIdx;
                bestInsertionComplement = readAlignment[0].complement;
            }

            auto croppedLength = readCroppingSlice[1] - readCroppingSlice[0];
            auto headerBuilder = parsePacBioHeader(rawFastaEntry.header);
            headerBuilder.qualityRegionBegin = 0;
            headerBuilder.qualityRegionBegin = croppedLength;
            auto newHeader = headerBuilder.to!string;
            auto expectedFastaLength = newHeader.length + croppedLength
                + croppedLength / options.fastaLineWidth + 2;
            auto croppedFastaEntry = appender!string;
            croppedFastaEntry.reserve(expectedFastaLength);

            croppedFastaEntry ~= newHeader;
            croppedFastaEntry ~= lineSep;
            // dfmt off
            auto croppedSequence = readAlignment[0].complement
                ? rawFastaEntry[readCroppingSlice[0] .. readCroppingSlice[1]].array.reverseComplement
                : rawFastaEntry[readCroppingSlice[0] .. readCroppingSlice[1]].array;
            // dfmt on
            croppedFastaEntry ~= croppedSequence.chunks(options.fastaLineWidth).joiner(lineSep);
            croppedFastaEntry ~= lineSep;

            croppedFastaEntries ~= croppedFastaEntry.data;
        }

        string croppedPileUpDb = buildDamFile(croppedFastaEntries.data, options);

        with (CoordinateTransform) with (ReadAlignmentType)
            {
                // dfmt off
                return Hit(
                    Insertion.make(
                        Coordinate(
                            croppingRefPositions.backContigId,
                            croppingRefPositions.backIdx,
                        ),
                        Coordinate(
                            croppingRefPositions.frontContigId,
                            croppingRefPositions.frontIdx,
                        ),
                        bestInsertionLength,
                        pileUpType,
                    ),
                    bestInsertionIdx,
                    bestInsertionComplement,
                    croppedPileUpDb,
                );
                // dfmt on
            }
    }

    protected size_t getInsertionScore(in ReadAlignment readAlignment,
            in ReadAlignmentType pileUpType,
            in Flag!"preferSpanning" preferSpanning = No.preferSpanning) pure
    {
        immutable shortAnchorPenaltyMagnitude = AlignmentChain.maxScore / 512;
        immutable notSpanningPenaltyMagnitude = AlignmentChain.maxScore / 2;
        immutable improperAlignmentPenaltyMagnitude = AlignmentChain.maxScore / 8;

        long numAlignments = readAlignment.length;
        long expectedAlignmentCount = pileUpType == ReadAlignmentType.gap ? 2 : 1;
        auto alignmentAnchor = readAlignment[].map!getRegion.fold!"a | b" - repetitiveRegions;
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
        long improperAlignmentPenalty = readAlignment[].map!"a.isProper ? 0 : 1".sum * improperAlignmentPenaltyMagnitude / numAlignments;
        // dfmt off
        size_t score = max(0, (
              avgAlignmentScore
            - shortAnchorPenalty
            - notSpanningPenalty
            - improperAlignmentPenalty
        ));
        // dfmt on

        // dfmt off
        debug logJsonDebug(
            "expectedAlignmentCount", expectedAlignmentCount,
            "avgAnchorSize", avgAnchorSize,
            "avgAlignmentLength", avgAlignmentLength,
            "avgAlignmentScore", avgAlignmentScore,
            "shortAnchorPenalty", shortAnchorPenalty,
            "score", score,
        );
        // dfmt on

        return score;
    }

    protected Hit buildConsensus(Hit croppedDbResult)
    {
        static struct ConsensusAlignmentOptions
        {
            string[] daccordOptions;
            string[] dalignerOptions;
            string[] dbsplitOptions;
            string workdir;
        }

        // dfmt off
        string consensusDb = getConsensus(croppedDbResult.dbFile, croppedDbResult.readId, const(ConsensusAlignmentOptions)(
            options.daccordOptions,
            options.dalignerOptions ~ [
                format!"-l%d"(options.minAnchorLength),
                format!"-e%f"((1 - options.readsErrorRate)^^2),
            ],
            options.dbsplitOptions,
            options.workdir,
        ));
        // dfmt on

        debug logJsonDebug("consensusDb", consensusDb);

        croppedDbResult.dbFile = consensusDb;
        croppedDbResult.readId = 1;

        return croppedDbResult;
    }

    protected DJunctor fetchTracePoints(PileUp pileUp)
    {
        auto allAlignmentChains = pileUp.getAlignmentChainRefs();
        allAlignmentChains.sort!("*a < *b", SwapStrategy.stable);
        allAlignmentChains.attachTracePoints(options.refDb, options.readsDb,
                options.damapperOptions, options);

        return this;
    }

    protected DJunctor insertHits()
    {
        logJsonDiagnostic("state", "enter", "function", "djunctor.insertHits");
        foreach (hit; catHits)
        {
            insertHit(hit);
        }
        // dfmt off
        logJsonDiagnostic(
            "step", "djunctor.insertHits",
            "coordTransform", coordTransform.toJson,
            "coordTransformPython", coordTransform.toString(),
        );
        // dfmt on
        logJsonDiagnostic("state", "exit", "function", "djunctor.insertHits");

        return this;
    }

    protected DJunctor insertHit(in Hit hit)
    {
        logJsonDiagnostic("state", "enter", "function", "djunctor.insertHit");

        // dfmt off
        const(size_t)[] refContigIds = [
            hit.insertion.begin.contigId,
            hit.insertion.end.contigId,
        ];
        // dfmt on
        auto complement = hit.complement;
        size_t readId = hit.readId;
        auto read = (() => {
            auto fastaEntries = getFastaEntries(hit.dbFile, [readId], options);
            auto front = fastaEntries.front;
            auto parsedFastaRecord = front.parseFastaRecord;

            return parsedFastaRecord;
        })()();

        auto insertSequence = read[];
        CoordinateTransform.Insertion insertionInfo = hit.insertion;
        // Update insertion info with the actual consensus length.
        insertionInfo.length = read.length;
        CoordinateTransform.Insertion trInsertionInfo = insertionInfo;
        trInsertionInfo.begin = T(insertionInfo.begin);
        trInsertionInfo.end = T(insertionInfo.end);
        auto refContigFastaEntries = refContigIds.map!(
                contigId => getReferenceFastaEntry(T(contigId))).array;
        auto endContigLength = refContigFastaEntries[$ - 1].length;
        string newHeader = refContigFastaEntries[0].header;

        if (insertionInfo.isExtension)
        {
            auto effectiveExtensionLength = getEffectiveExtensionLength(insertionInfo,
                    endContigLength);

            if (effectiveExtensionLength < minExtensionLength)
            {
                // dfmt off
                logJsonDiagnostic(
                    "info", "skipping insertion of short extension",
                    "effectedContig", insertionInfo.end.contigId,
                    "extensionLength", effectiveExtensionLength,
                    "extensionType", insertionInfo.type.to!string,
                );
                // dfmt on

                return this;
            }
        }

        size_t newContigLength = trInsertionInfo.totalInsertionLength(endContigLength);
        immutable lineSep = typeof(refContigFastaEntries[0]).lineSep;
        alias ReferenceSequence = typeof(refContigFastaEntries[0][0 .. 0]);
        ReferenceSequence beforeInsertSequence;
        ReferenceSequence afterInsertSequence;
        auto fastaBuilder = appender!string;
        auto expectedFastaLength = newHeader.length + newContigLength
            + newContigLength / options.fastaLineWidth + 2;
        fastaBuilder.reserve(expectedFastaLength);

        fastaBuilder ~= newHeader;
        fastaBuilder ~= lineSep;

        with (CoordinateTransform.Insertion.Type)
        {
            final switch (insertionInfo.type)
            {
            case gap:
                beforeInsertSequence = refContigFastaEntries[0][0 .. trInsertionInfo.begin.idx];
                afterInsertSequence = refContigFastaEntries[1][trInsertionInfo.end.idx .. $];
                break;
            case front:
                beforeInsertSequence = refContigFastaEntries[0][0 .. 0];
                afterInsertSequence = refContigFastaEntries[0][trInsertionInfo.end.idx .. $];
                break;
            case back:
                beforeInsertSequence = refContigFastaEntries[0][0 .. trInsertionInfo.begin.idx],
                    afterInsertSequence = refContigFastaEntries[0][0 .. 0];
                break;
            case relabel:
                assert(0, "invalid insertion");
            }
        }

        auto joinedSequence = chain(beforeInsertSequence, insertSequence, afterInsertSequence);
        fastaBuilder ~= joinedSequence.chunks(options.fastaLineWidth).joiner(lineSep);
        fastaBuilder ~= lineSep;

        debug assert(expectedFastaLength >= fastaBuilder.data.length,
                format!"avoid reallocation: expected %d but used %d"(expectedFastaLength,
                    fastaBuilder.data.length));

        auto overlayDb = buildDamFile([fastaBuilder.data], options);

        T ~= insertionInfo;
        referenceDb.addAlias(T(refContigIds[0]), overlayDb, 1);

        immutable sequencePreviewLength = 500;
        // dfmt off
        logJsonDebug(
            "step", "insertHits",
            "type", insertionInfo.type.to!string,
            "contigIds", refContigIds.array.toJson,
            "readId", readId,
            "complement", complement.to!bool,
            "insertSequence", insertSequence.array.to!string,
            "beforeInsertSequence", beforeInsertSequence.tail(sequencePreviewLength).array.to!string,
            "afterInsertSequence", afterInsertSequence.take(sequencePreviewLength).array.to!string,
            "insertionInfo", insertionInfo.toJson(),
            "transformedInsertionInfo", trInsertionInfo.toJson(),
            "overlayDb", overlayDb,
            "expectedFastaLength", expectedFastaLength,
            "newHeaderLength", newHeader.length,
            "newContigLength", newContigLength,
            "optionsFastaLineWidth", options.fastaLineWidth,
            "fastaBuilderDataLength", fastaBuilder.data.length,
        );
        // dfmt on
        debug assert(expectedFastaLength >= fastaBuilder.data.length,
                format!"avoid reallocation: expected %d but used %d"(expectedFastaLength,
                    fastaBuilder.data.length));

        logJsonDiagnostic("state", "exit", "function", "djunctor.insertHit");

        return this;
    }

    protected auto getReferenceFastaEntry(in size_t contigId) const
    {
        auto dbRef = referenceDb[contigId];

        return getReferenceFastaEntry(dbRef);
    }

    protected auto getReferenceFastaEntry(in DBUnion.DBReference dbRef) const
    {
        return parseFastaRecord(getFastaEntries(dbRef.dbFile, [dbRef.contigId + 0], options).front);
    }

    // dfmt off
    protected size_t getEffectiveExtensionLength(
        in CoordinateTransform.Insertion insertionInfo,
        in size_t endContigLength
    ) const
    // dfmt on
    {
        // dfmt off
        return insertionInfo.isFrontExtension
            ? insertionInfo.length - insertionInfo.end.idx
            : insertionInfo.length - (endContigLength - insertionInfo.begin.idx);
        // dfmt on
    }

    protected DJunctor finish()
    {
        logJsonDiagnostic("state", "enter", "function", "djunctor.finish");

        // dfmt off
        this
            .mergeOverlayDb
            .writeCoordTransform;
        // dfmt on

        logJsonDiagnostic("state", "exit", "function", "djunctor.finish");

        return this;
    }

    protected DJunctor mergeOverlayDb()
    {
        size_t numReferenceContigs = getNumContigs(referenceDb.baseDb, options);
        // dfmt off
        auto contigSources = iota(1, numReferenceContigs + 1)
            .map!(contigId => T(contigId))
            .map!(trContigId => tuple(trContigId, referenceDb[trContigId]))
            .group!"a[0] == b[0]"
            .map!"a[0]";
        auto fastaEntries = contigSources
            .save
            .map!"a[1]"
            .map!(dbRef => getReferenceFastaEntry(dbRef));
        immutable lineSep = typeof(fastaEntries.front).lineSep;

        fastaEntries
            .map!(fastaEntry => fastaEntry.toFasta(options.fastaLineWidth))
            .joiner(lineSep)
            .write;
        // dfmt on

        // dfmt off
        logJsonDebug(
            "numReferenceContigs", numReferenceContigs,
            "contigSources", contigSources.array.toJson,
        );
        // dfmt on

        foreach (newContigId, contigSource; contigSources.enumerate(1))
        {
            size_t oldContigId = contigSource[0];

            if (oldContigId != newContigId)
            {
                T ~= CoordinateTransform.Insertion.makeRelabel(oldContigId, newContigId);
                logJsonDebug("relabel", ["oldContigId" : oldContigId,
                        "newContigId" : newContigId].toJson);
            }
        }

        return this;
    }

    protected DJunctor writeCoordTransform()
    {
        if (!(options.coordTransform is null))
        {
            auto coordTransformFile = File(options.coordTransform, "w");

            coordTransformFile.write(
                    coordTransform.toString(CoordinateTransform.ExportLanguage.python));
        }

        return this;
    }

    private void logFillingInfo(string step, string readState)(in PileUp[] pileUpsByGap)
    {
        static auto getGapInfo(in PileUp pileUp)
        {
            bool isGap = pileUp.isGap;
            auto type = isGap ? ReadAlignmentType.gap : pileUp[0].type;

            alias lengthFilter = p => isGap ? p.length == 2 : p.length == 1;

            // dfmt off
            return Json([
                "type": Json(type.to!string),
                "contigIds": pileUp
                    .filter!lengthFilter
                    .map!(readAlignment => readAlignment[]
                        .map!"a.contigA.id"
                        .array)
                    .front
                    .toJson,
                "estimateLengthMean": pileUp
                    .filter!lengthFilter
                    .map!"a.getInsertionSize"
                    .array
                    .mean
                    .toJson,
                "estimateLengthMedian": pileUp
                    .filter!lengthFilter
                    .map!"a.getInsertionSize"
                    .array
                    .median
                    .toJson,
                "numReads": Json(pileUp.length),
            ]);
            // dfmt on
        }

        // dfmt off
        logJsonDiagnostic(
            "step", step,
            "readState", readState,
            "numGaps", pileUpsByGap.length,
            "gapInfo", pileUpsByGap.map!getGapInfo.array,
        );
        // dfmt on
    }
}

interface AlignmentChainFilter
{
    AlignmentChain[] opCall(AlignmentChain[] alignmentChains);
}

abstract class ReadFilter : AlignmentChainFilter
{
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

        return alignmentChains.filter!(ac => !discardedReadIds.contains(ac.contigB.id)).array;
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

/**
    This realizes a transform which translates coordinates before insertions
    to coordinates after them. The insertions can be added sequentially. The
    transform is fully functional at any point in time and can be converted
    to a Python 2.7 script for external usage.
*/
struct CoordinateTransform
{
    static struct Coordinate
    {
        size_t contigId;
        size_t idx;
    }

    static struct Insertion
    {
        enum Type
        {
            front,
            gap,
            back,
            relabel,
        }

        Coordinate begin;
        Coordinate end;
        private size_t _length;

        @disable this(Coordinate);
        @disable this(Coordinate, Coordinate);
        private this(Coordinate begin, Coordinate end, size_t length) pure nothrow
        {
            this.begin = begin;
            this.end = end;
            this._length = length;
        }

        static Insertion makeRelabel(size_t oldContigId, size_t newContigId) pure nothrow
        {
            return Insertion(Coordinate(newContigId), Coordinate(oldContigId), 0);
        }

        static Insertion make(Coordinate begin, Coordinate end, size_t length,
                in ReadAlignmentType insertionType) pure nothrow
        {
            final switch (insertionType)
            {
            case ReadAlignmentType.front:
                begin.idx = 0;
                break;
            case ReadAlignmentType.gap:
                break;
            case ReadAlignmentType.back:
                end.idx = begin.idx + length;
                break;
            }

            return Insertion(begin, end, length);
        }

        invariant
        {
            assert(begin != end || (begin.idx == end.idx && _length > 0), "empty insertion");
            if (begin.contigId == end.contigId)
            {
                assert(end.idx >= begin.idx,
                        "insertion begin index should be before/equal to end index on same contig");
                assert(_length >= end.idx - begin.idx, "inserted sequence too small");
            }
        }

        @property length() pure const nothrow
        {
            return _length;
        }

        @property length(size_t newLength) pure nothrow
        {
            final switch (type)
            {
            case Type.back:
                end.idx = begin.idx + newLength;
                goto case;
            case Type.front:
            case Type.gap:
                _length = newLength;
                return;
            case Type.relabel:
                return;
            }
        }

        @property isExtension() pure const nothrow
        {
            return begin.contigId == end.contigId;
        }

        @property isFrontExtension() pure const nothrow
        {
            return isExtension && !isBackExtension;
        }

        @property isBackExtension() pure const nothrow
        {
            return isExtension && (begin.idx > 0 && end.idx == begin.idx + length);
        }

        @property isGap() pure const nothrow
        {
            return !isExtension && !isRelabel;
        }

        @property isRelabel() pure const nothrow
        {
            return !isExtension && begin.idx == 0 && end.idx == 0 && length == 0;
        }

        @property type() pure const nothrow
        {
            if (isExtension)
            {
                return isFrontExtension ? Type.front : Type.back;
            }
            else if (isGap)
            {
                return Type.gap;
            }
            else
            {
                return Type.relabel;
            }
        }

        size_t totalInsertionLength(in size_t endContigLength) pure const nothrow
        {
            final switch (type)
            {
            case Type.gap:
            case Type.front:
            case Type.relabel:
                return begin.idx + length + endContigLength - end.idx;
            case Type.back:
                return begin.idx + length;
            }
        }
    }

    static enum ExportLanguage
    {
        python,
    }

    const(Insertion)[] insertions;

    /**
        Transform a coordinate according to insertions. If the designated
        contig is not present in insertions then the original coordinate is
        returned. Otherwise the new coordinate is computed as follows:

        1. Find insertion where end contig matches input. Note: if begin
           contig matches input not action is needed because neither the
           `idx` nor the `contigId` is affected.
        2. Travel from that point backwards to the beginning of the chain –
           ie. until begin contig and end contig of two insertions are
           unequal – summing all the coordinate shifts. The new reference
           contig is the last one in the (reverse) chain.

        For a single insertion and input coordinate `x` the shifted
        coordinate `T(x)` is calculated as follows:

        ---
        Case 1:  begin contig == end contig
            Case 1.1:  extend begin of contig

                   begin/end contig

                  0  ie  x
                  |---+--+-------|
                 /    |
                |-----|
                0    li

                T(x) = x
                     - ie  // (1) make relative to insertion point of end contig
                     + li  // (2) make relative to begin of insertion
                // choose ib = 0
                     = x - ie + li + ib

            Case 1.2:  extend end of contig

                        begin/end contig

                0       x ib        ie
                |-------+--+---| - - +
                           |         |
                           |---------|
                           0        li

                T(x) = x
                // choose ie = li + ib
                     = x - ie + li + ib

        Case 2:  begin contig != end contig

              begin contig       end contig

            0         ib     0  ie  x
            |----------+---| |---+--+-------|
                       |         |
                       |---------|
                       0        li

                        insertion

            T(x) = x
                 - ie  // (1) make relative to insertion point of end contig
                 + li  // (2) make relative to begin of insertion
                 + ib  // (3) make relative to begin of begin contig
        ---

        Returns: transformed coordinate.
    */
    Coordinate transform(in size_t contigId, in size_t idx) pure const
    {
        return transform(Coordinate(contigId, idx));
    }

    size_t transform(in size_t contigId) pure const
    {
        return transform(Coordinate(contigId, 0)).contigId;
    }

    /// ditto
    Coordinate transform(in Coordinate input) pure const
    {
        Coordinate result = input;

        // dfmt off
        auto insertionChain = insertions[]
            .retro
            .find!"a.end.contigId == b"(0 + input.contigId);
        // dfmt on
        size_t lastContigId = input.contigId;

        foreach (insertion; insertionChain)
        {
            if (insertion.end.contigId != lastContigId)
            {
                break; // chain end reached
            }

            // Note: reverse order of calculation to prevent negative
            // intermediate results
            result.idx += insertion.begin.idx; // (3)
            result.idx += insertion.length; // (2)
            result.idx -= insertion.end.idx; // (1)
            // Begin contig is the new reference
            result.contigId = insertion.begin.contigId;
            lastContigId = insertion.begin.contigId;
        }

        return result;
    }

    /// ditto
    alias opCall = transform;

    ///
    unittest
    {
        CoordinateTransform coordTransform;
        //           0        495 500 0   5  100   600
        // reference |----------+---| |---+--+-------|
        //                      |         |
        // insert               |---------|
        //                      0       100
        // dfmt off
        coordTransform.add(Insertion.make(
            Coordinate(1, 495),
            Coordinate(2, 5),
            100,
            ReadAlignmentType.gap,
        ));
        // dfmt on

        auto inputCoord = Coordinate(2, 100);
        // dfmt off
        auto transformedCoord = Coordinate(
            1,
            (
                100
                - 5    // coord on contig 2 relative to insertion point
                + 100  // make relative to begin of insertion
                + 495  // make relative to begin of contig 1
            )
        );
        // dfmt on
        assert(coordTransform.transform(inputCoord) == transformedCoord);
        assert(coordTransform(inputCoord) == transformedCoord);
    }

    unittest
    {
        with (ReadAlignmentType)
        {
            CoordinateTransform coordTransform;
            // dfmt off
            coordTransform.add(Insertion.make(
                Coordinate(100, 0),
                Coordinate(101, 0),
                0,
                gap,
            ));
            coordTransform.add(Insertion.make(
                Coordinate(1, 495),
                Coordinate(2, 5),
                100,
                gap,
            ));
            coordTransform.add(Insertion.make(
                Coordinate(2, 490),
                Coordinate(5, 10),
                150,
                gap,
            ));
            coordTransform.add(Insertion.make(
                Coordinate(5, 480),
                Coordinate(3, 20),
                200,
                gap,
            ));
            coordTransform.add(Insertion.make(
                Coordinate(102, 0),
                Coordinate(103, 0),
                0,
                gap,
            ));
            coordTransform.add(Insertion.make(
                Coordinate(200),
                Coordinate(200, 5),
                10,
                front,
            ));
            coordTransform.add(Insertion.make(
                Coordinate(201, 90),
                Coordinate(201),
                20,
                back,
            ));
            coordTransform.add(Insertion.make(
                Coordinate(202),
                Coordinate(202, 15),
                30,
                front,
            ));
            coordTransform.add(Insertion.make(
                Coordinate(202),
                Coordinate(202, 20),
                40,
                front,
            ));
            coordTransform.add(Insertion.make(
                Coordinate(203, 65),
                Coordinate(203),
                50,
                back,
            ));
            coordTransform.add(Insertion.make(
                Coordinate(203, 60),
                Coordinate(203),
                60,
                back,
            ));
            coordTransform.add(Insertion.makeRelabel(
                1024,
                1000,
            ));
            // dfmt on

            {
                // Case: contig not present in insertions
                auto inputCoord = Coordinate(1337, 42);
                assert(coordTransform.transform(inputCoord) == inputCoord);
            }
            {
                // Case: contig is never end of an insertions
                auto inputCoord = Coordinate(1, 64);
                assert(coordTransform.transform(inputCoord) == inputCoord);
            }
            {
                // Case: contig is end of the last insertions of a reverse chain.
                auto inputCoord = Coordinate(2, 64);
                // dfmt off
                auto transformedCoord = Coordinate(
                    1,
                    (
                          64   // last and only insertion in chain
                        - 5    // (1)
                        + 100  // (2)
                        + 495  // (3)
                    )
                );
                // dfmt on
                assert(coordTransform.transform(inputCoord) == transformedCoord);
            }
            {
                // Case: contig is end of an inner insertions, ie. there are
                //       insertions before and after in the chain.
                auto inputCoord = Coordinate(5, 64);
                // dfmt off
                auto transformedCoord = Coordinate(
                    1,
                    (
                        (
                              64   // last insertion in chain
                            - 5    // (1)
                            + 100  // (2)
                            + 495  // (3)
                        )      // first insertion in chain
                        - 10    // (1)
                        + 150  // (2)
                        + 490  // (3)
                    )

                );
                // dfmt on
                assert(coordTransform.transform(inputCoord) == transformedCoord);
            }
            {
                // Case: contig is end of the first insertions of a reverse chain.
                auto inputCoord = Coordinate(3, 64);
                // dfmt off
                auto transformedCoord = Coordinate(
                    1,
                    (
                        (
                            (
                                  64   // last insertion in chain
                                - 5    // (1)
                                + 100  // (2)
                                + 495  // (3)
                            )      // second insertion in chain
                            - 10    // (1)
                            + 150  // (2)
                            + 490  // (3)
                        )      // first insertion in chain
                        - 20    // (1)
                        + 200  // (2)
                        + 480  // (3)
                    )

                );
                // dfmt on
                assert(coordTransform.transform(inputCoord) == transformedCoord);
            }
            {
                // Case: coordinate on rejected part of end contig, ie. `x < ie`.
                auto inputCoord = Coordinate(2, 3);
                // dfmt off
                auto transformedCoord = Coordinate(
                    1,
                    (
                          3    // last and only insertion in chain
                        - 5    // (1)
                        + 100  // (2)
                        + 495  // (3)
                    )
                );
                // dfmt on
                assert(coordTransform.transform(inputCoord) == transformedCoord);
            }

            {
                // Case: simple front extension
                auto inputCoord = Coordinate(200, 3);
                // dfmt off
                auto transformedCoord = Coordinate(
                    200,
                    (
                          3    // only extension (in chain)
                        - 5    // (1)
                        + 10   // (2)
                        + 0    // (3)
                    )
                );
                // dfmt on
                assert(coordTransform.transform(inputCoord) == transformedCoord);
            }
            {
                // Case: simple back extension
                auto inputCoord = Coordinate(201, 3);
                // dfmt off
                auto transformedCoord = Coordinate(
                    201,
                    (
                          3    // only extension (in chain)
                        - 110  // (1)
                        + 20   // (2)
                        + 90   // (3)
                    )
                );
                // dfmt on
                assert(coordTransform.transform(inputCoord) == transformedCoord);
            }
            {
                // Case: double front extension
                auto inputCoord = Coordinate(202, 3);
                // dfmt off
                auto transformedCoord = Coordinate(
                    202,
                    (
                        (
                              3    // last extension (in chain)
                            - 15   // (1)
                            + 30   // (2)
                            +  0   // (3)
                        )     // first extension (in chain)
                        - 20  // (1)
                        + 40  // (2)
                        +  0  // (3)
                    )
                );
                // dfmt on
                assert(coordTransform.transform(inputCoord) == transformedCoord);
            }
            {
                // Case: double back extension
                auto inputCoord = Coordinate(203, 3);
                // dfmt off
                auto transformedCoord = Coordinate(
                    203,
                    (
                        (
                              3    // last extension (in chain)
                            - 115  // (1)
                            + 50   // (2)
                            + 65   // (3)
                        )      // first extension (in chain)
                        - 120  // (1)
                        + 60   // (2)
                        + 60   // (3)
                    )
                );
                // dfmt on
                assert(coordTransform.transform(inputCoord) == transformedCoord);
            }
            {
                // Case: relabel
                auto inputCoord = Coordinate(1024, 3);
                auto transformedCoord = Coordinate(1000, 3);

                assert(coordTransform.transform(inputCoord) == transformedCoord);
            }

            // dfmt off
            coordTransform.add(Insertion.make(
                Coordinate(203, 115),
                Coordinate(202, 5),
                50,
                gap,
            ));
            // dfmt on

            {
                // Case: double front extension + double back extension + gap spanned
                auto inputCoord = Coordinate(202, 3);
                // dfmt off
                auto transformedCoord = Coordinate(
                    203,
                    (
                        (
                            (
                                  3    // last front extension (in chain)
                                - 15   // (1)
                                + 30   // (2)
                                +  0   // (3)
                            )     // first front extension (in chain)
                            - 20  // (1)
                            + 40  // (2)
                            +  0  // (3)
                        )      // only gap (in chain)
                        -   5  // (1)
                        +  50  // (2)
                        + 115  // (3)
                    )
                );
                // dfmt on
                assert(coordTransform.transform(inputCoord) == transformedCoord);
            }
        }
    }

    /**
        Add an insertion to this transform.

        Returns: this transform.
        See_Also: `CoordinateTransform.transform`.
    */
    CoordinateTransform add(in Insertion newInsertion) pure
    {
        size_t idx = getInsertionIndex(newInsertion);

        if (idx == 0)
        {
            insertions = [newInsertion] ~ insertions;
        }
        else if (idx == insertions.length)
        {
            insertions ~= newInsertion;
        }
        else
        {
            insertions = insertions[0 .. idx] ~ [newInsertion] ~ insertions[idx .. $];
        }

        return this;
    }

    /// ditto
    void opOpAssign(string op)(in Insertion newInsertion) pure if (op == "~")
    {
        this.add(newInsertion);
    }

    unittest
    {
        with (ReadAlignmentType)
        {
            Insertion getDummyInsertion(in size_t beginContigId,
                    in size_t endContigId, in ReadAlignmentType insertionType = gap)
            {
                static immutable dummyLength = 42;

                return Insertion.make(Coordinate(beginContigId, dummyLength / 2),
                        Coordinate(endContigId, dummyLength / 2), dummyLength, insertionType);
            }

            CoordinateTransform getDummyTransform()
            {
                CoordinateTransform coordTransform;

                // dfmt off
                coordTransform.insertions = [
                    getDummyInsertion(1, 2),
                    getDummyInsertion(3, 4),
                    getDummyInsertion(9, 8),
                    getDummyInsertion(6, 5),
                ];
                // dfmt on

                return coordTransform;
            }

            {
                auto coordTransform = getDummyTransform();
                // Case 1 (Gap): newInsertion fits between two existing insertions
                // dfmt off
                coordTransform.add(getDummyInsertion(2, 3));
                assert(coordTransform.insertions == [
                    getDummyInsertion(1, 2),
                    getDummyInsertion(2, 3),
                    getDummyInsertion(3, 4),
                    getDummyInsertion(9, 8),
                    getDummyInsertion(6, 5),
                ]);
                coordTransform.add(getDummyInsertion(8, 6));
                assert(coordTransform.insertions == [
                    getDummyInsertion(1, 2),
                    getDummyInsertion(2, 3),
                    getDummyInsertion(3, 4),
                    getDummyInsertion(9, 8),
                    getDummyInsertion(8, 6),
                    getDummyInsertion(6, 5),
                ]);
                // dfmt on
            }
            {
                auto coordTransform = getDummyTransform();
                // Case 2* (Gap): newInsertion extends an existing insertion chain in the front
                // dfmt off
                coordTransform.add(getDummyInsertion(42, 1));
                assert(coordTransform.insertions == [
                    getDummyInsertion(42, 1),
                    getDummyInsertion(1, 2),
                    getDummyInsertion(3, 4),
                    getDummyInsertion(9, 8),
                    getDummyInsertion(6, 5),
                ]);
                // Case 2 (Gap): newInsertion extends an existing insertion chain in the front
                coordTransform.add(getDummyInsertion(42, 9));
                assert(coordTransform.insertions == [
                    getDummyInsertion(42, 1),
                    getDummyInsertion(1, 2),
                    getDummyInsertion(3, 4),
                    getDummyInsertion(42, 9),
                    getDummyInsertion(9, 8),
                    getDummyInsertion(6, 5),
                ]);
                // dfmt on
            }
            {
                auto coordTransform = getDummyTransform();
                // Case 2* (Extension): newInsertion extends an existing insertion chain in the front
                // dfmt off
                coordTransform.add(getDummyInsertion(1, 1, front));
                assert(coordTransform.insertions == [
                    getDummyInsertion(1, 1, front),
                    getDummyInsertion(1, 2),
                    getDummyInsertion(3, 4),
                    getDummyInsertion(9, 8),
                    getDummyInsertion(6, 5),
                ]);
                coordTransform.add(getDummyInsertion(1, 1, front));
                assert(coordTransform.insertions == [
                    getDummyInsertion(1, 1, front),
                    getDummyInsertion(1, 1, front),
                    getDummyInsertion(1, 2),
                    getDummyInsertion(3, 4),
                    getDummyInsertion(9, 8),
                    getDummyInsertion(6, 5),
                ]);
                // Case 2 (Extension): newInsertion extends an existing insertion chain in the front
                coordTransform.add(getDummyInsertion(9, 9, front));
                assert(coordTransform.insertions == [
                    getDummyInsertion(1, 1, front),
                    getDummyInsertion(1, 1, front),
                    getDummyInsertion(1, 2),
                    getDummyInsertion(3, 4),
                    getDummyInsertion(9, 9, front),
                    getDummyInsertion(9, 8),
                    getDummyInsertion(6, 5),
                ]);
                coordTransform.add(getDummyInsertion(9, 9, front));
                assert(coordTransform.insertions == [
                    getDummyInsertion(1, 1, front),
                    getDummyInsertion(1, 1, front),
                    getDummyInsertion(1, 2),
                    getDummyInsertion(3, 4),
                    getDummyInsertion(9, 9, front),
                    getDummyInsertion(9, 9, front),
                    getDummyInsertion(9, 8),
                    getDummyInsertion(6, 5),
                ]);
                // dfmt on
            }
            {
                auto coordTransform = getDummyTransform();
                // Case 3 (Gap): newInsertion extends an existing insertion chain in the back
                // dfmt off
                coordTransform.add(getDummyInsertion(4, 42));
                assert(coordTransform.insertions == [
                    getDummyInsertion(1, 2),
                    getDummyInsertion(3, 4),
                    getDummyInsertion(4, 42),
                    getDummyInsertion(9, 8),
                    getDummyInsertion(6, 5),
                ]);
                coordTransform.add(getDummyInsertion(5, 42));
                assert(coordTransform.insertions == [
                    getDummyInsertion(1, 2),
                    getDummyInsertion(3, 4),
                    getDummyInsertion(4, 42),
                    getDummyInsertion(9, 8),
                    getDummyInsertion(6, 5),
                    getDummyInsertion(5, 42),
                ]);
                // dfmt on
            }
            {
                auto coordTransform = getDummyTransform();
                // Case 3 (Extension): newInsertion extends an existing insertion chain in the back
                // dfmt off
                coordTransform.add(getDummyInsertion(4, 4, back));
                assert(coordTransform.insertions == [
                    getDummyInsertion(1, 2),
                    getDummyInsertion(3, 4),
                    getDummyInsertion(4, 4, back),
                    getDummyInsertion(9, 8),
                    getDummyInsertion(6, 5),
                ]);
                coordTransform.add(getDummyInsertion(4, 4, back));
                assert(coordTransform.insertions == [
                    getDummyInsertion(1, 2),
                    getDummyInsertion(3, 4),
                    getDummyInsertion(4, 4, back),
                    getDummyInsertion(4, 4, back),
                    getDummyInsertion(9, 8),
                    getDummyInsertion(6, 5),
                ]);
                coordTransform.add(getDummyInsertion(5, 5, back));
                assert(coordTransform.insertions == [
                    getDummyInsertion(1, 2),
                    getDummyInsertion(3, 4),
                    getDummyInsertion(4, 4, back),
                    getDummyInsertion(4, 4, back),
                    getDummyInsertion(9, 8),
                    getDummyInsertion(6, 5),
                    getDummyInsertion(5, 5, back),
                ]);
                coordTransform.add(getDummyInsertion(5, 5, back));
                assert(coordTransform.insertions == [
                    getDummyInsertion(1, 2),
                    getDummyInsertion(3, 4),
                    getDummyInsertion(4, 4, back),
                    getDummyInsertion(4, 4, back),
                    getDummyInsertion(9, 8),
                    getDummyInsertion(6, 5),
                    getDummyInsertion(5, 5, back),
                    getDummyInsertion(5, 5, back),
                ]);
                // dfmt on
            }
            {
                auto coordTransform = getDummyTransform();
                // Case 4 (Gap): open new insertion chain
                // dfmt off
                coordTransform.add(getDummyInsertion(1337, 42));
                assert(coordTransform.insertions == [
                    getDummyInsertion(1, 2),
                    getDummyInsertion(3, 4),
                    getDummyInsertion(9, 8),
                    getDummyInsertion(6, 5),
                    getDummyInsertion(1337, 42),
                ]);
                coordTransform.add(getDummyInsertion(1337, 42));
                assert(coordTransform.insertions == [
                    getDummyInsertion(1, 2),
                    getDummyInsertion(3, 4),
                    getDummyInsertion(9, 8),
                    getDummyInsertion(6, 5),
                    getDummyInsertion(1337, 42),
                    getDummyInsertion(1337, 42),
                ]);
                // dfmt on
            }
            {
                auto coordTransform = getDummyTransform();
                // Case 4 (Extension): open new insertion chain
                // dfmt off
                coordTransform.add(getDummyInsertion(1337, 1337, front));
                assert(coordTransform.insertions == [
                    getDummyInsertion(1, 2),
                    getDummyInsertion(3, 4),
                    getDummyInsertion(9, 8),
                    getDummyInsertion(6, 5),
                    getDummyInsertion(1337, 1337, front),
                ]);
                coordTransform.add(getDummyInsertion(42, 42, back));
                assert(coordTransform.insertions == [
                    getDummyInsertion(1, 2),
                    getDummyInsertion(3, 4),
                    getDummyInsertion(9, 8),
                    getDummyInsertion(6, 5),
                    getDummyInsertion(1337, 1337, front),
                    getDummyInsertion(42, 42, back),
                ]);
                // dfmt on
            }
            {
                CoordinateTransform emptyCoordTransform;

                // Case 4* (Gap): open new insertion chain
                assert(emptyCoordTransform.add(getDummyInsertion(1337, 42))
                        .insertions == [getDummyInsertion(1337, 42)]);
            }
            {
                CoordinateTransform emptyCoordTransform;

                // Case 4* (Extension): open new insertion chain
                assert(emptyCoordTransform.add(getDummyInsertion(1337, 1337,
                        front)).insertions == [getDummyInsertion(1337, 1337, front)]);
            }
            {
                CoordinateTransform emptyCoordTransform;

                // Operator Style
                emptyCoordTransform ~= getDummyInsertion(1337, 42);
                assert(emptyCoordTransform.insertions == [getDummyInsertion(1337, 42)]);
            }
        }
    }

    private size_t getInsertionIndex(in Insertion newInsertion) pure const
    {
        // dfmt off
        if (
            // Case 4*: open new insertion chain
            insertions.length == 0 ||
            // Case 2*: newInsertion extends an existing insertion chain in the front
            newInsertion.end.contigId == insertions[0].begin.contigId
        )
        {
            return 0;
        }
        // dfmt on
        // Case 3*: newInsertion extends an existing insertion chain in the back
        // Case 4**: open new insertion chain
        else if (insertions.length == 1)
        {
            // Note: this case is needed in order for `slice(2)` to always yield
            //       pairs (of size two).
            return 1;
        }

        size_t case1Result = size_t.max;
        size_t case2Result = size_t.max;
        size_t case3Result = size_t.max;

        foreach (i, insAB; insertions.slide(2).enumerate)
        {
            auto insA = insAB[0];
            auto insB = insAB[1];

            // Case 1: newInsertion fits between two existing insertions; take first result
            if (case1Result == size_t.max && newInsertion.begin.contigId == insA.end.contigId
                    && newInsertion.end.contigId == insB.begin.contigId)
            {
                case1Result = i + 1;
            }

            // Case 2: newInsertion extends an existing insertion chain in the front; take first result
            if (case2Result == size_t.max && newInsertion.begin.contigId != insA.end.contigId
                    && newInsertion.end.contigId == insB.begin.contigId)
            {
                case2Result = i + 1;
            }

            // Case 3: newInsertion extends an existing insertion chain in the back; take last result
            if (newInsertion.begin.contigId == insA.end.contigId
                    && newInsertion.end.contigId != insB.begin.contigId)
            {
                case3Result = i + 1;
            }
        }

        if (case1Result != size_t.max)
            return case1Result;
        else if (case2Result != size_t.max)
            return case2Result;
        else if (case3Result != size_t.max)
            return case3Result;

        // Case 4: open new insertion chain
        return insertions.length;
    }

    unittest
    {
        Insertion getDummyInsertion(in size_t beginContigId, in size_t endContigId)
        {
            return Insertion(Coordinate(beginContigId, 0), Coordinate(endContigId, 0), 0);
        }

        CoordinateTransform coordTransform;

        // dfmt off
        coordTransform.insertions = [
            getDummyInsertion(1, 2),
            getDummyInsertion(3, 4),
            getDummyInsertion(9, 8),
            getDummyInsertion(6, 5),
        ];
        // dfmt on
        {
            // Case 1: newInsertion fits between two existing insertions
            assert(coordTransform.getInsertionIndex(getDummyInsertion(2, 3)) == 1);
            assert(coordTransform.getInsertionIndex(getDummyInsertion(8, 6)) == 3);
        }
        {
            // Case 2*: newInsertion extends an existing insertion chain in the front
            assert(coordTransform.getInsertionIndex(getDummyInsertion(42, 1)) == 0);
            // Case 2: newInsertion extends an existing insertion chain in the front
            assert(coordTransform.getInsertionIndex(getDummyInsertion(42, 9)) == 2);
        }
        {
            // Case 3: newInsertion extends an existing insertion chain in the front
            assert(coordTransform.getInsertionIndex(getDummyInsertion(4, 42)) == 2);
            assert(coordTransform.getInsertionIndex(getDummyInsertion(5, 42)) == 4);
        }
        {
            // Case 4: open new insertion chain
            assert(coordTransform.getInsertionIndex(getDummyInsertion(1337, 42)) == 4);
            assert(coordTransform.getInsertionIndex(getDummyInsertion(42, 1337)) == 4);
        }
        {
            CoordinateTransform emptyCoordTransform;

            // Case 4*: open new insertion chain
            assert(emptyCoordTransform.getInsertionIndex(getDummyInsertion(1337, 42)) == 0);
            assert(emptyCoordTransform.getInsertionIndex(getDummyInsertion(42, 1337)) == 0);
        }
        coordTransform.insertions = [getDummyInsertion(1, 2)];
        {
            // Case 3*: newInsertion extends an existing insertion chain in the front
            assert(coordTransform.getInsertionIndex(getDummyInsertion(2, 42)) == 1);
        }
        {
            // Case 4**: open new insertion chain
            assert(coordTransform.getInsertionIndex(getDummyInsertion(1337, 42)) == 1);
        }
    }

    string toString(in ExportLanguage lang = ExportLanguage.python) pure const
    {
        immutable estimatePrefaceLength = 250;
        immutable estimateLengthPerInsertion = 100;
        immutable estimateEpilogueLength = 2000;

        auto result = appender!string;
        result.reserve(estimatePrefaceLength + estimateLengthPerInsertion
                * insertions.length + estimateEpilogueLength);

        final switch (lang)
        {
        case ExportLanguage.python:
            buildPythonString(result);
            break;
        }

        return result.data;
    }

    private void buildPythonString(out Appender!string builder) pure const
    {
        immutable preface = q"EOF
            from __future__ import print_function
            from collections import namedtuple
            from sys import argv, exit

            Coordinate = namedtuple("Coordinate", ["contig_id", "idx"])
            Insertion = namedtuple("Insertion", ["begin", "end", "length"])

            INSERTIONS = [
EOF".outdent;
        immutable epilogue = q"EOF
            ]


            def transform(contig_id, idx, insertions=INSERTIONS):
                result_contig_id = contig_id
                result_idx = idx
                last_contig_id = contig_id

                for insertion in reversed(insertions):
                    if insertion.end.contig_id != last_contig_id and \
                            last_contig_id != contig_id:
                        break  # chain end reached
                    elif insertion.end.contig_id == last_contig_id or \
                            last_contig_id != contig_id:
                        # chain begin reached or inside chain

                        result_idx += insertion.begin.idx
                        result_idx += insertion.length
                        result_idx -= insertion.end.idx
                        result_contig_id = insertion.begin.contig_id
                        last_contig_id = insertion.begin.contig_id

                return Coordinate(result_contig_id, result_idx)


            def print_help():
                print("Usage: {} [-h] CONTIG_ID IDX".format(argv[0]))
                print()
                print(
                    "Transform coordinates from input coordinates to output " +
                    "coordinates."
                )
                print("Prints `NEW_CONTIG_ID NEW_IDX` to STDOUT.")
                print()
                print("Positional arguments:")
                print(" CONTIG_ID  Contig ID in .dam before transformation")
                print(" IDX        Base index (zero-based) in sequence")
                print()
                print("Optional arguments:")
                print(" -h         Prints this help.")

            if __name__ == "__main__":
                if len(argv) != 3 or "-h" in argv:
                    print_help()
                    exit(1)

                try:
                    contig_id = int(argv[1])
                    idx = int(argv[2])

                    if contig_id < 0 or idx < 0:
                        raise Exception()
                except Exception:
                    print_help()
                    exit(2)

                transformed = transform(contig_id, idx)
                print("{} {}".format(transformed.contig_id, transformed.idx))
EOF".outdent;
        immutable insertionTemplate = q"EOF
            Insertion(
                Coordinate(%d, %d),
                Coordinate(%d, %d),
                %d
            ),
EOF".outdent.indent(4);

        builder ~= preface;
        foreach (insertion; insertions)
        {
            // dfmt off
            formattedWrite!insertionTemplate(
                builder,
                insertion.begin.contigId,
                insertion.begin.idx,
                insertion.end.contigId,
                insertion.end.idx,
                insertion.length,
            );
            // dfmt on
        }
        builder ~= epilogue;
    }
}

/// Access contigs distributed over various DBs by ID.
struct DBUnion
{
    alias DBReference = Tuple!(string, "dbFile", size_t, "contigId");

    string baseDb;
    DBReference[size_t] overlayDbs;

    /// Return the `dbFile` for `contigId`.
    DBReference opIndex(in size_t contigId) pure const
    {
        return overlayDbs.get(contigId, DBReference(baseDb, contigId));
    }

    unittest
    {
        auto dbQuery = DBUnion("baseDb");
        dbQuery.overlayDbs[42] = DBReference("overlayDb42", 1);

        assert(dbQuery[0] == DBReference("baseDb", 0));
        assert(dbQuery[7] == DBReference("baseDb", 7));
        assert(dbQuery[42] == DBReference("overlayDb42", 1));
    }

    /// Set the source for `contigId`.
    void addAlias(in size_t contigId, in string dbFile, in size_t newContigId) pure
    {
        addAlias(contigId, DBReference(dbFile, newContigId));
    }

    /// ditto
    void addAlias(in size_t contigId, in DBReference dbRef) pure
    {
        overlayDbs[contigId] = dbRef;
    }

    unittest
    {
        auto dbQuery = DBUnion("baseDb");

        assert(dbQuery[42] == DBReference("baseDb", 42));

        dbQuery.addAlias(42, DBReference("overlayDb42", 1));

        assert(dbQuery[42] == DBReference("overlayDb42", 1));
    }
}

///
unittest
{
    auto dbQuery = DBUnion("baseDb");

    dbQuery.addAlias(42, "overlayDb42", 1);

    assert(dbQuery[0].dbFile == "baseDb");
    assert(dbQuery[0].contigId == 0);
    assert(dbQuery[7].dbFile == "baseDb");
    assert(dbQuery[7].contigId == 7);
    assert(dbQuery[42].dbFile == "overlayDb42");
    assert(dbQuery[42].contigId == 1);

    dbQuery.addAlias(7, "overlayDb7", 1);

    assert(dbQuery[7].dbFile == "overlayDb7");
    assert(dbQuery[7].contigId == 1);
}
