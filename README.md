djunctor
=========

Close assembly gaps using long-reads with focus on correctness.


Building
--------

```sh
git clone https://github.com/a-ludi/djunctor.git
cd djunctor
dub build
```


Developer Notes
---------------

```sh
# run the program for development
dub run -- ARGS...

# run tests
dub test

# clean up
dub clean
```


Functional Description
-------------------

`djunctor` takes two inputs, a reference assembly and a set of long reads
(PacBio) and tries to close scaffold gaps or reduce their size.

### Functional Overview

```pascal
program djunctor(ReferenceAssembly, LongReads)
begin
    FindAlignments()

    repeat
        FilterUselessReads()
        BuildPileUpsFromAlignments()

        for each pileUp do
            SelectGoodReads()
            BuildConsensus()
            InsertHit()
    until (numHits == 0 or maxIterationsReached)

    OutputResult()
end
```

### Difficulties

#### Low Complexity Regions/Repeats

1. Low Complexity Regions
    ```fasta
    ccgcacctcaaatcgtcaccgttgtgtatcgaggggacttatagtgc
    tcctgtgacatgtcactgttgcggtcgaaccggtcgtgcaatccgac
    gtcccaatgcccgccgcattaacggtagccatAAAAAAAAAAAAAAA
    AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
    AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAcgc
    atcaccgatcggggtcggtaataaaaggacaaagttagtgttggcca
    cgaacttctcacgaataagttccctggttttgcgagggaatgcatct
    gctaggcgtcactggacacagtgggaaagctgccgggggcga
    ```
2. Low Complexity Regions
   ```fasta
   ccgcacctcaaatcgtcaccgttgtgtatcgaggggacttatagtgc
   tcctgtgacatgtcAGTAGTAGTAGTAGTAGTAGTAGTAGTAGTAGT
   AGTAGTAGTAGTAGTAGTAGTAGTAGTAGTAGTAGTAGTAGTAGTAG
   TAGTAGTAGTAGTAGTAGTAGTAGTAGTAGTAGTccacgagctggag
   cctaaaacaattccatgagactggtctaggttacgcagtgtagccgc
   atcaccgatcggggtcggtaataaaaggacaaagttagtgttggcca
   cgaacttctcacgaataagttccctggttttgcgagggaatgcatct
   gctaggcgtcactggacacagtgggaaagctgccgggggcga
   ```
3. Tandem Repeats
   ```fasta
   ccgcacctcaaatcgtcaccgttgtgtatcgaggggacttatagtgc
   tcctgtgacatgtcactgttgcggtcgTTGTGTATCGAGGGGACTTA
   TAGTGCTCCTGTTTGTGTATCGAGGGGACTTATAGTGCTCCTGTTTG
   TGTATCGAGGGGACTTATAGTGCTCCTGTTTGTGTATCGAGGGGACT
   TATAGTGCTCCTGTTTGTGTATCGAGGGGACTTATAGTGCTCCTGTT
   TGTGTATCGAGGGGACTTATAGTGCTCCTGTTTGTGTATCGAGGGGA
   CTTATAGTGCTCCTGTaagttccctggttttgcgagggaatgcatct
   gctaggcgtcactggacacagtgggaaagctgccgggggcga
   ```
4. Transposable elements
   ```fasta
   ccgcacctcaaatcgtcaccgttgtgtatcgaggggacttatagtgc
   tcctgtGACATGTCACTGTTGCGGTCGAACCGGTCGTGcaatccgac
   gtcccaatgcccgccgcattaacggtagccataGACATGTCACTGTT
   GCGGTCGAACCGGTCGTGtagcgcgacaaaaaccccacgagctggag
   cctaaaacaattccatgagactggtctaggGACATGTCACTGTTGCG
   GTCGAACCGGTCGTGtcgtaaaggtctgtcatagtttgtgtgtgtga
   gcggaagtataaacgaaaagaggaccagaaaaGACATGTCACTGTTG
   CGGTCGAACCGGTCGTGacagtgggaaagctgccgggggcga
   ```


License
-------

This project is licensed under [MIT License](./LICENSE).
