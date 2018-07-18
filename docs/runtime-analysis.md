Runtime & Memory Analysis
=========================

[TOC]

\\(
    \newcommand{\O}{\mathcal{O}}
\\)
This document describes an detailed analysis of runtime and memory usage of
this algorithm. Analysis is made separately for each step of the algorithm as
outlined below.

```d
void djunctor()
{
    init();
    assessRepeatStructure();
    filterAlignments();

    PileUp[] pileUps = buildPileUps();
    foreach (ref pileUp; pileUps)
    {
        processPileUp(pileUp);
    }

    writeNewAssembly();
    writeUnusedReadsList();
}
```

Definitions
-----------

- \\(G_r\\)    … Size of the reference genome.
- \\(G_i\\)    … Size of the reference assembly.
- \\(n_c\\)    … Number of contigs in the reference.
- \\(n_s\\)    … Number of scaffolds in the reference.
- \\(n_r\\)    … Number of reads.
- \\(C\\)      … Read coverage.
- \\(l_g\\)    … Length of longest gap that is being filled.
- \\(n_{mi}\\) … Number of mask intervals in the input mask.
- \\(n_{mr}\\) … Number of mask intervals in the repeat mask.
- \\(n_{sa}\\) … Number of local alignments in the self-alignment of the reference.
- \\(n_{ra}\\) … Number of local alignments in the ref vs. reads alignment.
- \\(n_{tp}\\) … Number of trace points in the ref vs. reads alignment.

Analysis
--------

### `init()`

* Total Runtime: \\( \O(n_c \log(n_c) + n_{mi} \log(n_{mi}) + n_{sa} \log(n_{sa}) + n_{ra} \log(n_{ra}) + n_R) \\)
* Total Memory:  \\( \O(n_c + n_{mi} + n_{sa} + n_{ra} + n_r) \\)

```d
void init()
{
    // O(1) in time and memory.
    numReferenceContigs = getNumContigs(options.refDb, options);
    // O(1) in time and memory.
    numReads = getNumContigs(options.readsDb, options);
    // O(n_c) in time and memory.
    scaffoldStructure = getScaffoldStructure(options.refDb, options).array;

    if (options.inMask !is null)
    {
        // O(n_mi log(n_mi)) in time and O(n_mi) in memory.
        this.repetitiveRegions = ReferenceMask(readMask!ReferenceInterval(options.refDb,
                options.inMask, options));
    }


    // O(n_sa log(n_sa)) in time and O(n_sa) in memory
    selfAlignment = getAlignments(options.refDb, options.selfAlignmentFile, options);
    // O(n_ra log(n_ra)) in time and O(n_ra) in memory
    readsAlignment = getAlignments(options.refDb, options.readsDb,
            options.refVsReadsAlignmentFile, options);

    // O(n_c log(n_c)) in time and O(n_c) in memory.
    initResultScaffold();
    // O(n_r) in time and O(n_r) in memory.
    initUnusedReads();
}
```

### `assessRepeatStructure()`
* Total Runtime: \\( \O(n_{sa} \log(n_c) + n_{ra} \log(n_c) + n_{mi}) \\)
* Total Memory:  \\( \O(n_{sa} + n_c + n_{ra} + n_{mi}) \\)


```d
void assessRepeatStructure()
{
    // O(1) in time and memory
    auto alphaHalf = (1 - options.confidence) / 2;
    // O(n_sa) in time and memory
    auto selfCoverage = alignmentCoverage(selfAlignment);
    // O(n_ra) in time and memory
    auto readsCoverage = alignmentCoverage(readsAlignment);
    // TODO
    auto selfCoverageConfidenceInterval = tuple(
        invPoissonCDF(alphaHalf, selfCoverage),
        invPoissonCDF(1 - alphaHalf, selfCoverage),
    );
    // TODO
    auto readsCoverageConfidenceInterval = tuple(
        invPoissonCDF(alphaHalf, readsCoverage),
        invPoissonCDF(1 - alphaHalf, readsCoverage),
    );

    // O(n_sa log(n_c)) in time and O(n_sa + n_c) in memory
    auto repetitiveRegions1 = new BadAlignmentCoverageAssessor(selfCoverageConfidenceInterval.expand)(selfAlignment);
    // O(1) in time since `repetitiveRegions1` is already sorted and no new memory is required
    this.repetitiveRegions |= repetitiveRegions1;

    // O(n_ra log(n_c)) in time and O(n_ra + n_c) in memory
    auto repetitiveRegions2 = new BadAlignmentCoverageAssessor(readsCoverageConfidenceInterval.expand)(readsAlignment);
    // O((n_sa + n_c) log(n_ra + n_c)) in time since both operands are already sorted and no new memory is required
    this.repetitiveRegions |= repetitiveRegions2;

    if (options.outMask != null)
    {
        // O(n_mi) in time and memory
        writeMask(options.refDb, options.outMask, this.repetitiveRegions.intervals, options);
    }
}
```
### `filterAlignments()`

* Total Runtime: \\( \O(n_{ra} log^2(n_{ra})) \\)
* Total Memory:  \\( \O(n_{ra}) \\)

```d
void filterAlignments()
{
    auto filters = tuple(
        // O(n_ra n_mr) in time and O(n_ra) in memory
        new WeaklyAnchoredAlignmentChainsFilter(repetitiveRegions, options.minAnchorLength),
        // O(n_ra) in time and O(n_ra) in memory
        new ImproperAlignmentChainsFilter(),
        // O(n_ra log²(n_ra)) in time and O(n_ra) in memory
        new AmbiguousAlignmentChainsFilter(&unusedReads),
        // O(n_ra log²(n_ra)) in time and O(n_ra) in memory
        new RedundantAlignmentChainsFilter(&unusedReads),
    );
    // O(1) in time and memory
    auto filterInput = readsAlignment[];

    // see above
    foreach (i, filter; filters)
    {
        auto filterOutput = filter(filterInput[]);
        filterInput = filterOutput;
    }
    // O(1) in time and memory
    auto filterPipelineOutput = filterInput;
    // O(1) in time and memory
    this.readsAlignment = filterPipelineOutput;
}
```
### `buildPileUps()`

* Total Runtime: \\( \O(n_{ra} \log(n_{ra}) + n_c \log(n_c)) \\)
* Total Memory:  \\( \O(n_{ra} + n_c) \\)

```d
PileUp[] buildPileUps()
{
    return buildPileUpsFromReadAlignments(numReferenceContigs, readsAlignment)
        .filter!(pileUp => pileUp.length >= options.minReadsPerPileUp)
        .array;
}

PileUp[] buildPileUpsFromReadAlignments(in size_t numReferenceContigs, AlignmentChain[] candidates)
{
    // O(n_ra log(n_ra)) in time and O(1) in memory
    candidates.sort!("a.contigB.id < b.contigB.id", SwapStrategy.stable);
    // O(n_ra) in time and memory
    auto readAlignmentJoins = candidates
        // O(1) in time and memory
        .chunkBy!isSameRead
        // O(1) in time and memory
        .map!collectReadAlignments
        .joiner
        // O(1) in time and memory
        .filter!"a.isValid"
        // O(1) in time and memory
        .map!"a.getInOrder()"
        .map!(to!ReadAlignmentJoin)
        .array;
    auto alignmentsScaffold = buildScaffold!(concatenatePayloads!Payload, Payload)(numReferenceContigs + 0, readAlignmentJoins)
        // O(n_c + n_ra) in time and O(1) in memory
        .discardAmbiguousJoins!Payload
        // O(n_c log(n_c)) in time and O(1) in memory
        .mergeExtensionsWithGaps!("a ~ b", Payload);
    // O(n_c) in time and memory
    auto pileUps = alignmentsScaffold
        .edges
        .filter!"a.payload.length > 0"
        .map!"a.payload"
        .filter!(pileUp => pileUp.isValid)
        .array;

    return pileUps;
}
```

### `processPileUp(pileUp)`

* Total Runtime: \\( \O(n_c (C (l_g + log(C)) + n_mi + log(n_c))) \\)
* Total Memory:  \\( \O(C l_g) \\)

```d
// This is called O(n_c) times
void processPileUp(ref PileUp pileUp)
{
    // O(C l_g + C (l_g + log(C)) + n_mi) in time and O(C l_g) in memory
    auto croppingResult = cropPileUp(pileUp, repetitiveRegions, options);
    // O(C) in time and O(1) in memory
    auto referenceReadIdx = bestReadAlignmentIndex(pileUp, Yes.preferSpanning, options);
    auto referenceRead = pileUp[referenceReadIdx];
    // TODO
    auto consensusDb = buildConsensus(croppingResult.db, referenceReadIdx + 1);

    if (pileUp.isExtension && shouldSkipShortExtension(croppingResult,
            consensusDb, referenceRead))
    {
        return;
    }

    // O(log(n_c)) in time and O(1) in memory
    addInsertionToScaffold(referenceRead, consensusDb, croppingResult.referencePositions);
    // O(log(n_c)) in time and O(1) in memory
    addFlankingContigSlicesToScaffold(croppingResult.referencePositions);
    // O(C) in time and O(C) in memory
    markReadsAsUsed(pileUp);
}
```
### `writeNewAssembly()`

* Total Runtime: \\( \O(G + n_c log(n_c))) \\)
* Total Memory:  \\( \O(G) \\)

```d
void writeNewAssembly()
{
    if (!options.shouldExtendContigs)
    {
        // O(n_c) in time and O(1) in memory
        catHits = catHits.removeExtensions!InsertionInfo();
    }
    catHits = catHits
        // O(n_c log(nc_)) in time and O(1) in memory
        .enforceJoinPolicy!InsertionInfo(options.joinPolicy)
        // O(n_c log(nc_)) in time and O(1) in memory
        .normalizeUnkownJoins!InsertionInfo()
        // O(n_c log(nc_)) in time and O(1) in memory
        .fixContigCropping();
    // O(G) in time and O(G) in memory
    foreach (startNode; contigStarts!InsertionInfo(catHits))
    {
        writeNewContig(startNode);
    }
}
```

### `writeUnusedReadsList()`

* Total Runtime: \\( \O(n_{ra})) \\)
* Total Memory:  \\( \O(n_{ra}) \\)
