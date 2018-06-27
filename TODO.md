# TODO

-   [ ] Improve trimming/filtering of contigs using depth covariance information.
-   [ ] Better organization of intermediate data
    -   [ ] Drop the res/seq dichotomy
    -   [ ] Organize split versions of to-be-combined files into appropriately named
        folders.
        -   Will I know what the appropriate naming is?
-   [ ] Are the differences between the mice at UT that we KNOW have
        vA and those with vB?
-   [ ] Construct OTU-3.vB (based on the libraries chosen fro OTU-1.vB) as a
        feel for how much effect there is of poor-er sampling of vB.
-   [x] Regenerate non-reassembled, refined, depth-trimmed draft genome comparison for OTU-1
    -   The back-mapping of reads to the un-assembled genome may not have been done correctly.
-   [x] Regenerate quast reports for OTU-1-A and OTU-1-B
    -   The depth trimming parameters have been changed, and I'm hoping they result in longer
        contigs.
-   [x] Draft genome generation with reassembled scaffolds instead of contigs.
    -   I'm hoping this generates longer genome fragments.
-   [ ] Rerun refinement for OTU-2 and OTU-4
-   [ ] Generate genomes/refinement for Muribaculaceae MAGs
-   [ ] Consider writing about how most of the differences in OTU-1-A and -B
    can be described as duplication or deletion events of multi-copy genes (as
    opposed to novel genes).
    -   This would require some additional evidence that duplications aren't just
        due to bad genome building.
    -   I think I could get at this by looking at coverage.
        Bad genome building would result in lower coverage for the parts that I
        _thought_ were duplication events.
        Real duplication events would not have the same dip in coverage.
-   CheckM draft OTU-1, OTU-7 genomes.
-   [ ] Add library ID tag to the backmapping BAMs (careful with the timestamps)
-   [x] Pilon for library-set specific genome reconstruction
-   [x] Import unique sequence abundance, taxonomy, and OTU abundance from longev 16S results
-   [ ] Automate bin combination w/ 16S results
-   [x] Look for evidence of strain variation in gene content
-   [x] File rather than directory based PROKKA output
-   [x] COGs to Modules
-   [x] Check OTU-2 bin assembly against Lactobacillus johnsonii
-   [ ] "Strain resolved analysis of an abundant Muribaculaceae species across murine hosts"

