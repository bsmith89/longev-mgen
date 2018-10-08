# Requirements

-   prokka
    -   With dbCAN family level database installed see http://csbl.bmb.uga.edu/dbCAN/download.php
    -   Be sure to insteall dbCAN alphabetically after HAMAP, or else it'll supersede
        the more interprettable HAMAP annotations.
    -   TODO: Consider naming these HMM dbs as 0.HAMAP 1.dbCAN etc.
-   MinPath (v1.4, patched shebang line to 'python2' not 'python'
-   Bignorm (see https://git.informatik.uni-kiel.de/axw/Bignorm)
-   MCL (https://micans.org/mcl/)
-   TODO

Metadata:

-   `raw/longev_rrs_results.db` can be acquired from the longev
    project as `res/lite.results.db`.
    Code for producing this file is in `Makefile`.
    Current md5sum is `5e74824958c80488a523fdbd8a108e84`.
-   TODO: Archive the results of data extraction from raw/longev_rrs_results.db,
    since this data will ultimately be publicly available, while the DB itself
    will not.

    27e7cc16f81420bace9f5ef1c44a3079  res/core.r.count.tsv
    9085eff5efb579a23ba86aec6d5140d1  res/core.r.taxonomy.tsv

# Bin refinement

TODO: Encapsulate this work in a script/recipes.

Combining bins is important.
Let's ignore missassembly.
Let's assume that I've been sufficiently conservative with binning that I
don't put contigs in the same bin that are never found in the same genome.
If I just trust the clustering (based on abundance co-variance) I'll have
a few problems, _even when my CheckM scores are high_!
The key to this conclusion is that co-abundance clustering is bad at
combining sequence that differentiates strains (the parts of the pan-genome
that are not universal..."epigenome"?)
These contigs will not have as high coverage covariance as contigs in the
core will have to each other.
If they're all in a single cassette, they may still bin together, but not
with the core.
The result is that you'll have one bin with the core genome.
Intuitively, this core bin will probably have a high CheckM score.
And then a bunch of smaller bins with very low CheckM scores with epigenes.
This means that checkM's marker gene informed merge approach won't work,
since there won't be any marker genes in the other bins to inform merging.
Linkage is a better option, but bins linked to the core bin will probably
have lower coverage than the core, since they're (by definition) not
ubiquitous.

I like to think about the _real_ genome abundances as a latent variable
that we can never observe.
Then, the coverage of contigs (or bins) is a linear combination of these genomes,
since various genomes can share contigs (both core and epi).
This latent variable model might be implicit in a PCA.
(Or the compositional equivalent?
The null covariance structure is going to be important in whatever model
I use, ultimately.)
This thought process extends to the 16S sequence information, too.
Some of the genomes have the same 16S sequence, and some have different.
Even within an OTU, we have several unique sequences that probably
reflect some variation in gene content which would be lost from
normal co-abundance binning.
As a result, the relative abundance of 16S sequences is also informative
for getting at this latent genome composition.
We can combine and conquor!
The model that I just described (two multivariate datasets that
both reflect an underlying latent property of the system)
is implicit in a Cannonical Correspondance Analysis (CCA).
A related approach (with implicit shrinkage) is one version of partial
least squares (PLS).

I identified bins that covaried in whole or in part with unique V4
sequences by running a 'PLSCanonical' (correspondence analysis in sklearn)
between the V4 relative abundance and bin relative coverage (normalized
to total coverage).
This model should reflect the fact that 16S % abundance is a linear combination
of the abundance of genomes (a latent variable) with that sequence,
and the bin abundance is a the same.
Therefore, we account for redundancy, whereby multiple genomes might have the
same V4 sequence, but they also might have some subset of the same genes.
We use our 44 replicated libraries with both types of data to estimate
abundance covariances (which we would not be able to do without a
reasonable sample size).
It's also important that we use enough latent components;
these represent the actual genome abundances, so if you have too few you're
not accurately modeling the data.
Would cross-validation be a good way to pick this number?
I'm not sure I'm interpretting the results correctly, but we can
calculate the matrix product of the loadings on each term (X and Y,
`x_loadings_`, and `y_loadings_` of the fitted PLSCanonical instance),
which serves as a "importance score" and indicates the association between
each unique sequence and each bin.
This would presumably also work at the contig level (rather than bin),
but by binning first we incorporate tetranucleotide information
and reduce the dimensionality.

At first picking where to set the inclusion cutoff was difficult.
I think I settled on the following:
Include anything that is within 50% of the score of the best scoring match.

It seems to work well.  I recover the following genomes when I use 30 latent
components (this was determined by 6-fold cross validation), sqrt transform OTU
relative abundance _and_ relative coverage, and pick using the 50% of max rule.

OTU1A bin1256,bin1311,bin1379,bin503
OTU1B bin069,bin1169,bin1256,bin1311,bin1408,bin1464,bin491
OTU1C bin1169,bin1408,bin1449,bin1832,bin491
OTU2A bin306
OTU3A bin137
OTU4A bin1178,bin1257,bin1323,bin1408,bin1965,bin293,bin396,bin560,bin855,bin955
OTU4B bin1178,bin1257,bin293,bin560,bin855,bin955,bin976
OTU5A bin196
OTU6A bin1568,bin1617,bin1811,bin922
OTU6C bin1568,bin1617,bin1811,bin922
OTU7A bin1379,bin1465,bin1906,bin244,bin503,bin937
OTU8A bin1256,bin1464,bin1527,bin1552,bin1778,bin1906,bin503,bin922,bin956

These look pretty good!!!

-   OTU1
    -   For OTU1_A and B, bin1311 has 96% completeness, and _none_ of the
        matched bins (for either unique sequence) have any of the marker genes
        (at least, none of the universal marker genes).
        Our power to find matches for C is lower, because it's only found at UT
        and it's at much lower abundance.
    -   It also MAY not be at greater abundance
        in ACA mice.
    -   Can these features explain our inability to associate bin1311
        with that OTU?  It does have an association in the PLS analysis, it
        just doesn't rise to the 50% threshold.

-   OTU4
    -   For OTU4_A and B, all three of the core bins (bin560, bin293, and
        bin1257) are combined in both OTUs (something that would have otherwise
        required complementarity data, and that is well confirmed by linkage
        data.) A has some extra bins that B doesn't, while the
        inverse is only true for bin976, a single contig, 10kb bin.
    -   There is overlap with OTU1: e.g. bin1408, an 8 contig, 40kb bin
        which has linkage to bins in both OTU4 and OTU1, so who knows.

-   OTU2 / OTU3
    -   It's not apparent from these combinations, but there is a _strong_
        tendency for OTU2 and OTU3 (the spike) to get clustered together.
        No idea why this is.
    -   Tetranucleotide frequencies would definitely solve this problem.
    -   Why is this happening?  Perhaps it is because ACA reduces the abundance
        of OTU2, and also increases the density of the bacterial population.
        That means the spike is _less_ abundant.
        So the two covary because of an external common cause.
    -   When it doesn't quite rise to this threshold, OTU2 and OTU3 each
        have just a single bin assigned to them with high completeness.
    -   OTU3 has high linkage to a ton of other stuff.
        Given high abundance of the OTU, this could be because of
        random mapping, greater probability of missassembly, or greater
        probability of chimeras forming during PCR.
        In any case, it's bad.

I'm going to take these genomes as correct with _one_ adjustment;
I'll add bin1311 to OTU1C, since I'm confident it should be there.
I also won't do anything with OTUs 7 and 8, since they don't appear to be
clean (lots of overlap and they're a little on the large side).

## Commentary

OTUs are approximately the right size: all around 3 Mb

Interestingly, we end up with some inter-OTU overlap.
Some of it _could_ be real and interesting:
multiple strains of _Muribaculaceae_ with
reads mapping to shared bins...could it be parasitic sequence?
I should check bin1256.

It could also be spurious correlations.
Maybe I need some way for sharing to be _rare_, but not impossible (since we
want similar OTUs to share.)
Or maybe I should be more conservative with what shares a bin.
Combining sequence into contigs or bins that shouldn't be together
may introduce enough coverage correlation to be problematic.
This would _also_ confuse any linkage based approaches leveraging the paired ends.

Speaking of linkage, right now my best approach to combining bins is a combination
of checkM complementarity and paired-end linkage.
This can be further informed with the 16S leveraged cross-covariance structure,
although I don't have any automated way of doing this.

I wonder if my use of "2 blocks canonical PLS" is messed up.

It's also possible that site dependence (pseudoreplication) is messing
up the process.
Two OTUs that only show up at (e.g.) UT are going to be highly correlated
even though they wouldn't be conditional on site.
How do I include this correction?

How conservative should I be to answer my particular questions?
I'm increasingly interested in looking at these large genomic
differences.
Maybe I should just grab everything into one OTU1 including the whole pangenome.

## Decision

For OTU-1 and OTU-4 (here called OTU-7), I used the following to put
all relevant contigs together:

```
cat seq/core.a.bins.d/{bin00560,bin01257,bin00855,bin00955,bin00293,bin01178,bin01965,bin01811,bin01568}.fn > seq/core.a.mbins.d/Otu0007.fn
cat seq/core.a.bins.d/{bin01311,bin00503,bin01379,bin01256,bin01169,bin01408,bin00491,bin01449,bin01784,bin01832,bin01178,bin01464,bin00820}.fn > seq/core.a.mbins.d/Otu0001.fn
```

(date: 2018-05-
Binning, metabinning, and strain-variation finding all in one step?
All three of these things revolve around finding nested components.
As in: there's a main component that describes the bin, then there's
a meta-component that describes bins that go together
(just without any tetranucleotide component), and finally there's a
search bins/contigs that are less abundant than the core, but that co-vary
with the core in SOME samples.

Levels:

-   fine-grained covariance will identify contigs that ALWAYS go together
    these could be core genomes or any other fragments of sequence that
    always go together
-   Then there's nesting, where the coverage of core genomes accounts for
    the maximum abundance of fragments that belong to that species, but
    cannot predict lower-abundances, since the framgnet might be missing
    from some or all of the community.
-   I wonder if these cassettes nest at all, like is there a fractal structure?
-   This all breaks down if the same sequnce is variably present across multiple
    other genomes.  It would also break down with biased coverage estimation
    (which could be due to multiple copies).  This second option seems like
    it could be dealt with in the model.
-   I have a general sense that we get lots of extra information from the clustering
    pattern; like, if one sample has suprisingly high abundance of a particular
    sequence, we can put it in the context of its "cassette" and we don't have
    to think that that sequence is variably present, because it could just be
    measurement error.
-   Lots of opportunities for regularization.
-   The presence/absence of large sequence fragments from a core genome might
    (hypothesis!) mirror SNP variation (assuming there's no HGT)



# CAZy analysis

Plan:

-   Identify all CDS with CAZy hits
-   Cluster these bins
    -   Blast all-by-all
    -   Use bitscore to cluster proteins (k-means?)

# Proposal

The mammalian gut microbiome is a complex ecological system that influences
energy balance, pathogen resistance, and inflammation, among other processes
with importance to host health.
Understanding how the bacterial inhabitants of the gut respond to environmental
perturbations is a major step in developing a predictive framework for
microbiome-based therapies.
Acarbose (ACA) is an alpha-amylase and alpha-glucosidase inhibitor that
has been shown to increase the amount of starch entering the lower
digestive system after a meal.
In mice, rats, and humans, individuals treated with ACA have increased
concentrations of short-chain fatty acids (SCFAs) in feces, presumably as
a result of greater bacterial fermentation of this polysaccharide substrate.
Interestingly, ACA increase median lifespan by 20% in male mice, and 5% in
females.

# Paper Topics

-   Question 1: Why do OTUs 1 and 4 (7 here) respond strongly to acarbose while
    other Muribaculaceae do not? (TODO: Confirm that other Muri's do not
    respond positively to ACA.)
    -   Hypothesis 1.A: It has something to do with the carbohydrate
        utilization potential of these two strains.
        -   TODO: Follow-up on the following genes that differentate OTU-1 and
            OTU-4 from the Ormerod strains:
            -   ec3.2.1.39 Glucan endo-1,3-beta-D-glucosidase.
                (only these three taxa have it in my strains and only 10% of
                Ormerod strains)
            -   ec3.1.1.3 Triacylglycerol lipase.
                (only these strains and one more of those here have it, and 10
                of Ormerod strains)
        -   Figure 1.A.i: Venn-diagram of OTU-1 and OTU-4 COGs colored by category
        -   Figure 1.A.ii: Carb COGs that they DO have in common to the
            exclusion of most other here-strains (ordination makes them appear
            pretty different in carb space.)
    -   Hypothesis 1.B: It has something to do with the competitiveness for
        nitrogen.
        -   Figure 1.B.i: Ordination w/ Ormerod strains both panels for both
            carbohydrate and nitrogen COGs
            -   OTU-1 and OTU-4 differentiate themselves in carb-cog space.
                While they could both be alpha-glucan guild members, OTU-4
                looks _more_ like a host-glycan utilizer; it is
                borderline, though.
            -   In AA COG space, all of the guilds jumble together, although
                there does appear to be at least one principle component of
                differentiation.  It appears that OTU-1 and OTU-4 are similarly
                placed in that space.
-   Question 2: Why does OTU-1 NOT respond as strongly to ACA at UT where OTU-4
    does?
    -   Hypothesis 2.A: Competition w/ OTU-4 (for carbs or nitrogen)
    -   Hypothesis 2.B: OTU-1 is different at UT
        -   Figure 2.B.i: Evidence that OTU-1 comes in two varieties.
            -   TODO: Figure out where the 16S went.
            -   TODO: Figure out if there's nucleotide variation in the V4 for each strain
            -   TODO: Mapping of vB libraries to vA contigs and reverse to
                intended to confirm that we aren't being more strict towards
                one than the other.  (try to detect false negatives)
        -   Figure 2.B.ii: What genes differentiate the two strains?
            -   Evidence of duplication events being a major component of
                the diff?
            -   TODO: Follow-up on the genes that differ between
                OTU-1.vA and OTU-1.vB, especially ec3.2.1.31 Beta-glucuronidase,
                which is found in OTU-1.vA _and_ OTU-4, but not OTU-1.vB.
                Also see ec3.2.1.18 / COG4409 (sialidase) which is
                present in OTU-1.vA but not vB.
            -   TODO: Follow-up on the mucin associated story suggested
                by beta-glucouronidase and sialidase.
                -   Can we find evidence for differentail metabolism of
                    proline, threonine, and serine, (and glycine, which
                    can be easily interconverted with threonine/serine?)
                    TODO: Follow-up on Glycine hydroxymethyltransferase
                    (COG0112 / ec2.1.2.1) which is present in OTU-1.vA and
                    OTU-4, but not OTU-1.vB (could it be that this gene is used
                    in the conversion of serine/threonine to glycine, and isn't
                    useful if OTU-1.vB doesn't have access to mucin
                    amino-acids?)
            -   TODO: CAZy annotation of full genomes and then match them and
                see which genes don't overlap. (This is basically to catch
                what HAMAP and COGS didn't.)
                -   This seems to work well:
                -   I identify potential carb.  active enzymes by pulling
                    anything that hits dbCAN
                -   I then do an all-by-all blastp (diamond) to calculate
                    bitscores on alignments.
                -   Then I cluster everything using AverageNeighbor (sklearn)
                -   I can show that even clusters that include hits to the same
                    CAZy domains are quite distinct (but internally consistent)
                -   See Opu0085, Opu0095, Opu0219, Opu0239 (which all contain
                    GH31) but where alignment shows that they are clearly
                    different proteins.

        -   TODO: What differentiates mice with one strain versus the other?
            -   I think there's some evidence that sex and treatment
                differentiate mice with the two strains.
            -   This would be a particularly interesting result, since OTU-4
                seems to favor females and OTU-1 males.
        -   TODO: Any evidence for HGT between OTU-4 and OTU-1?
            -   See `OTU-7_202_pilon_0_3743`, `OTU-1_273_pilon_0_4562`, which
                have much overlapping length.
-   Question 3: Can we add anything to the discussion from Ormerod?
    -   Hypothesis 3.A: dN/dS will identify genes under strong selection in the
        Muribaculaceae.
    -   Hypothesis 3.B: Any evidence for HGT?
    -   Hypothesis 3.C: Can we get more detailed in our analysis of carb active
        enzymes?

# Paper Outline

## Introduction
The mammalian gut microbiome is a complex ecological system that influences
energy balance [@TODO], pathogen resistance [@TODO], and inflammation [@TODO],
among other processes with importance to host health.
Understanding how the bacterial inhabitants of the gut respond to environmental
perturbations is a major step in developing a predictive framework for
microbiome-based therapies.
Acarbose (ACA) is an alpha-glucosidase inhibitor prescribed for the
treatment of type 2 diabetes mellitus because it reduces the absorption of
glucose from starch in the small intestine [@TODO].
In rodents, ACA has been shown to increase the amount of starch entering the
lower digestive system after a meal [@TODO], resulting
in changes to the composition of the gut microbiota and its fermentation
products [@TODO].
Interestingly, long-term treatment with ACA has been shown to
substantially increase longevity in male mice and to a lesser extent in females
[@TODO].

<!--
TODO: Something big picture can go here.
-->

In a recent experimental study we showed that the relative abundance of a
number of bacterial taxa as well as the concentrations of propionate and
butyrate respond to long term treatment with ACA [@TODO].
This study was notable in being replicated across three sites: The University
of Michigan (UM) in Ann Arbor, The University of Texas Health Science
Center at San Antonio (UT), and The Jackson Laboratory (TJL) in Bar Harbor,
Maine.
At UM and TJL one highly abundant bacterial taxon,
classified as a member of the Bacteroides family Muribaculaceae and here
designated as B1, was found to be enriched nearly 4-fold in ACA treated mice.
B1 was also present at UT but was not found to be significantly more abundant
in ACA treated mice relative to controls.
Instead, a different member of the Muribaculaceae, designated B2, was found to
be highly abundant and 4-fold enriched in ACA-treated mice, but was nearly
absent at UM and TJL.
Five other Muribaculaceae, designated B3 through B7, were also identified
as among the most abundant members of the mouse gut microbiota across the three
sites, although none of these were found to be enriched in ACA treatment.

<!--
This site-dependence illustrates the importance of the local
bacterial metacommunity in determining the response of the gut bacterial
community to ACA, and the dangers of generalizing microbiome studies in one
population of mice to others.
-->

Family Muribaculaceae&mdash;formerly the S24-7 and sometimes referred to as
Homeothermaceae&mdash;has only one published cultivar [@TODO] despite being a
common and abundant inhabitant of the mammalian gut, especially in mice.
Previous studies have suggested that the Muribaculaceae specialize on
the fermentation of complex polysaccharides, much like members of the genus
Bacteroides also in order Bacteroidales.
Based on genomes reconstructed from metagenomes, Ormerod et al. [@Ormerod2016]
proposed three distinct carbohydrate utilization guilds, which they describe
as specialists on alpha-glucans, plant glycans, and host glycans, respectively.
While it is reasonable to expect that alpha-glucan specialists would be most
benefited by the large influx of starch to the gut resulting from ACA
treatment, this prediction has not been tested,
and physiological inferences based on the genome content of members of
this clade have been largely divorced from biological observations.

<!--
TODO: Introduce MAGs?
-->

Here we analyze genomes assembled from fecal metagenomes for mice at UT and UM
and explore two closely related questions about the niche of B1 and B2 in the
lower digestive system.
First, why do B1 and B2 each increase with ACA treatment, while other
Muribaculaceae do not?
And second, why is the response of B1 site specific;
it is _not_ found to be enriched in treated mice at UT, but is at high relative
abundance in the local community?
Despite similar patterns of abundance at their respective sites, these two taxa
seem to be only distantly related, sharing just 90% of positions in their 16S
rRNA gene V4 hypervariable region.
We nonetheless find evidence that B1 and B2 occupy overlapping niches,
specializing in the degradation of alpha-glucans, a role not held by the other
Muribaculaceae described in our study.
In addition, we identify two distinct strains of B1, with functionally relevant
differences in gene content, one of which is specific to the samples collected
at UT.


<!--
increasing
4-fold from 10% to approximately 40% relative abundance.
Estimates of the absolute density of 16S suggest that
the change in proportion were due to an increase in M1 abundance
rather than decreases in others.
M1 was similarly increased in ACA treated mice at TJL.
At UT, however, a different member of the Muribaculaceae, designated M2,
was 4-fold more abundant in ACA treated mice, going from a median
relative abundance of 6% to 26%.
M2 was not found or at very low levels in the vast majority of mice at UM
and TJL.

-->



<!--
Alternatively, the excess of carbohydrate may lead to selection for other
traits that provide access to key resources (e.g. nitrogen).

Given the increased abundance observed for both B1 and B2, but not other
Muribaculacae in acarbose treated mice, we identify features that may underlie
a differential response

Two alternative may account for the differential pattern of abundance
across sites.
Competition with M2 may inhibit the dramatic expansion of M1 at UT.
B2 is not found at UM or TJL where B1 responds to ACA.
Alternatively, despite being the same OTU across sites, hidden strain variation
not paralleled by variation in the V4 region may differentiate
B1 at UT from the other two sites.

Interestingly, these two OTUs are not closely related, sharing only about
90% of their 16S rRNA gene V4 hypervariable region sequence.
Identifying shared functional potential between these two responding taxa may
therefore present an opportunity to better understand the traits that poise
some bacterial populations to expand due to ACA treatment, while others
are unaffected or reduced.
-->


In this study we apply strain-resolved genome reconstruction from metagenomes
to compare and contrast the functional potential among Muribaculaceae at UM and
UT.
We demonstrate the utility of culture-free genomics to understand
the ecological role of key members of the mouse gut microbial community,
and explore several hypotheses that may explain differences in the distribution
and response of bacteria to perturbations.
Hypotheses derived from this analysis provide a foundation for future
physiological studies in recently available cultivars.
While the vast majority of host-associated bacterial species have never been
isolated, let alone characterized, combining experimental data from complex
communities with the analysis of reconstructed genomes provides a powerful tool
for expanding understanding to these potentially important taxa.

## Recovery of Genome

Metagenomes were sequenced for a set of samples containing an average of X
percent of N OTUs classified (based on the V4 region) as Muribaculaceae.
We roughly calculate a TODO table of expected coverages for these genomes,
suggesting that we might be able to pick up all of them.
Genomes were recovered through a process of
assembly, binning, bin refinement, reassembly, reassembly refinement, and trimming.
Raw assembly statistics TODO.
Mapping statistics TODO.
Genomes for all 7 OTUs were recovered with less than 5% contamination and
in 5 cases >95% completeness.
Two versions of OTU-1 were obtained; this is further discussed in a later
section.
The remainder of this section will focus on the genome recovery method,
and the evidence that these genomes are an accurate reflection of the true
genomes.

A fraction of the samples for which we obtained sequence data were spiked
with S. alaskensis, a marine bacterium that is highly unlikely to be found
in the mouse gut microbiome.
This allowed us to asses the quality of genome recovery and reconstruction
in an internal control.
However, this is still an idealized scenario, since recovery was not
hindered by strain variation or the presence of closely related species.
While binning contigs based on a mixture model of high-dimensional gaussians
(as in CONCOCT [@TODO]) produced feasibly complete and low contamination
bins in some cases, it was not very successful in others, including
taxa that we had the most interest in.



## Recovery of Genomes (2)

TODO: It seems like we're having special trouble recovering OTU-1 compared
to other genomes.
This could possibly because of strain heterogeneity, both that differentating
variant A from variant B, but also heterogeneity represented in the V4 region
variation (which seems to be somewhat independent of the rest of the genome).

-   S. alaskensis recovery was very good.
-   CheckM results suggest success.
-   TABLE of genome results
    -   Length
    -   Num contigs (N50?)
    -   GC content
    -   checkM results
        -   Completeness contamination
    -   Gene count and annotation rate
    -   What is the sequencing depth? Variance?
    -   MetaCyc pathways as a more robust measure of contamination?? (no, probably not)
-   TODO: Standardize the choice of contigs for each genome.
    -   Consider not reassembling
-   SUPPLEMENTAL FIGURE of genome recovery workflow
-   FIGURE phylogenetics with Ormerod
-   TODO: Add back 16S genes and identify variation
-   TODO: Do I discuss recovery of strain variants here? (probably)
-   TODO: Unique genes that differentiate OTU-1 from the Ormerod strains?
    (answer the question: why is it special that I recovered them?)
-   TODO: Figure mapping OTU-1 strains against each other
    -   This will also include depth lines for alignments from
        libraries positive for each.
-   TODO: Mapping from OTU-7 to M6 (not exactly pretty, but I could quote the percent of genomes that align)

## Shared niche of OTU-1 and OTU-7

A number of genes are shared by both OTUs to the exclusion of most
other Muribaculaceae in my study.
Both OTU-1 strains fall cleanly into the alpha-glucan guild.
TODO: Is this still as clean, with the most recent annotations?
TODO: Which annotation do we use to assert this?
OTU-7 is found much closer to the host-glycan folks.
OTU-5 is the only other OTU to fall out with the alpha-glucans.

Looking for genes (especially those that contain a GH13) that
are specific to OTU-1 _and_ OTU-7 provides one way to identify functions
that provide an advantage in acarbose treated mice.
It makes sense that that advantage would also be experienced in
the absence of ACA, since starch is probably still very common
in the gut.
I believe that both OTU-1 and OTU-7 are particularly unlikely to have
amino-acid synthesis machinery (see the minpath results, these
should be MUCH more robust when predicting the absence of pathways
compared to the presence...maybe?  TODO: check this intuition)

TODO: Why does OTU-5 _not_ behave like OTU-1 and OTU-7?

This section is written without consideration for Otu0001\_vB, since
we're not talking about why that strain fails to bloom at UT.

TODO: A grand survey of genes that are shared by OTU-1 and OTU-7, but that
are rare in other strains (especially if they're common in either
the starch or the host-glycan ormerod guilds)

TODO: Some in-depth amino-acid auxotrophy analysis.
See <https://doi.org/10.1371/journal.pgen.1007147>.

## Strain differences between OTU-1 at UM and UT

There are two strains of OTU-1:
one that's found at both sites (vA) and one found only at UT (vB).
These two strains share approximately 80% of their length.
(CheckM says that each is 90% complete, so it's possible that they are
actually identical, but we just failed to sequence/assemble/confirm/bin a different 10%
of the each genome, but this seems unlikely for a number of reasons:
-   the two assemblies were not independent
-   the coverage is precipitously lower in libraries lacking that strain for
    the strain specific parts of each genome.
    It's enough lower that it's not just random noise, and it includes large
    contigs which we would expect to get _some_ coverage on, regardless of
    stochasticity.
-   the unshared parts have some overlap in gene content with different context
    (TODO: confirm this)
The fact that they differ could explain the lack of response of OTU-1 at UT.
It is unkown which strain is found at TJL (TODO: Is it worth considering
getting some sequence?)
It is interesting to ask what genes differentiate the two strains.
It is interesting that no other Muri's have multiple strains. (perhaps
they just aren't the same OTU, but some strains are specific to sites?)
There is also good evidence of nucleotide (synonymous and non-synonymous)
variation between strains.  TODO: How do I demonstrate this?
TODO: What is the average amino-acid identity in shared genes?

The 16S gene recovered with each genome, does it share a V4 region?
I believe the answer is yes, but that there are variants at other positions.
Interestingly, this would suggest that the V4 region has shared
variation at several sites, despite differing outside of the V4
(TODO: can we see this in the samtools pileup results?)
Does it make sense for the V4 regions to not differ in the recovered sequence
while they clearly differ in the 16S amplicon results at each site?

Can we assert that any of the strain differences in gene content are real
and not just the result of incomplete recovery?
Does the sequencing depth lend credence?

TODO: Should we talk about the apparent presence of a phage in the OTU-1 genome
at UT?

# TODO

-   Consider the newly identified secreted GH13 containing enzyme
    -   What is its genomic context?  Part of a PUL?
    -   Is there any evidence for a similar enzyme in OTU-7?
        -   If not, how is OTU-7 blooming with starch??
    -   Who else has the same enzyme?
        -   Nothing in blast has higher than ~50% id.  Is that weird??
        -   Nothing with the same OPF, but is that because they're too granular?
        -   Any other taxa not in this study?  Are they starch degraders?
-   Figure out what's going on with Opf01942:
    -   This GH13 containing enzyme is found in OTU-1 (both sites),
        OTU-9,
        several Ormerod strains (M11, H7, M10, all annotated as starch degraders),
        as well as B. ovatus (but not B. theta),
    -   It's annotated as 'Isomalto-dextranase' in all cases.
    -   In most cases it includes a cysteine just after a signal-peptide.
        -   One notable exception (the only exception?) is in B. ovatus, which apparently lacks the signal peptide (SignalP score approx. 0.36)
    -   Interestingly, OTU-1 actually has THREE copies of this OPF (no other mag has copies), but they're substantially different.
        -   Two of the copies are in successive positions in the genome
            -   Of these, one copy is annotated as a GH97 domain, not GH13 (despite being clustered as the same OPF)
            -   and is fused with a ["SusF_SusE" domain](https://www.ebi.ac.uk/interpro/entry/IPR032187); no other MAG has this fusion.
            -   When I mask using Gblocks (so we're just looking at the GH13) and
                estimate a phylogeny, the two copies in OTU-1 are more closely
                related to each other than to any of the other sequences, suggesting
                duplication since divergence from the other lineages.
                -   TODO: dN/dS analysis to show positive selection on one or both copies?
        -   The third copy is more similar to the other MAGs sequences than the first two.
-    OTU-7, does it have starch activity on the outer membrane??
     -  In OTU-7, I get just one hit, even using a liberal filter, for proteins
        containing GH13, GH97, or GH31.  (Are there other GH domains I should
        include in my search? perhaps GH66)  Otu0007_vA_01521, is super short (75 AAs),
        contains a partial GH97 domain.  I think it's either a pseudogene or an
        assembly error, because the following ORF, which also contains a GH97,
        is missing a ~80 residues from the N-terminus (based on an alignment to
        other proteins in the same cluster).  I'm assuming the actual protein
        is their concatenation.
    -   So there is enzymatic breakdown of starch in OTU-7!  But it's a GH97,
        not GH13.  It's therefore not at all clear if OTU-7 would be
        particularly good at growing on extracellular starch.  That protein
        family (Opf01276) is found across Muribaculaceae (including OTU-1) and
        in both B. theta and B. ovatus,
        (TODO: check the following) but OTU-7's is the only sequence with
        evidence for secretion.
        -   OTU-1 and OTU-7 are the only Muri's in the ITP mice that have a secreted GH97.
            Opf01297 in OTU-7 and one of the versions of the duplicated
            Opf01942 in OTU-1.
    -   Widening the search to CBMs 20, 25, 26, and 69 I get three lipoproteins
        (Aside: should I include other CBMs as well?  I'm having a hard time
        finding a list of those that might bind starch.):
        -   Otu0007_vA_00858 containing one (or maybe two) CBM26s, preceded by
            a 'P_gingi_FimA' domain.  I'm imagining the structure as a long
            protein "stick" with a starch "hook" on the end.  It's genomically
            adjacent to another fimbrial-looking protein with a CBM16 on the
            end and an 'Mfa2' containing protein.  None of the other species
            have homologous proteins, although I do find it in the Ormerod
            reconstruction of the same species.
        -   Otu0007_vA_01525 with repeated SusE and SusF domains.
        -   And Otu0007_vA_01826 containing a large fraction, but not all of a
            SusE domain.  From my limited reading of the literature, this would
            be a putative inner-membrane lipoprotein (since there is an Asp, 2
            residues from the cleavage site).
    -   I get one more hit when I widen the search still further to CBM48.
        -   Otu0007_vA_00489: besides the CBM48, it also contains a 'PCMD'
            which is a Pfam domain hypothesized to have GH activity...
            Homologous sequences are found in B. ovatus, OTU-5, OTU-9, and a
            whole bunch of the Ormerod MAGs.
-   I've filtered my OPF list by those in which most of the member features
    are annotated with CAZY domains (to remove everything but the polysaccharide active proteins).
    -   In this subset I find an OPF (Opf01517) that has 18 members, mostly
        OTU-1 duplicates (3 copies in vB and 4 copies in vC !!)
    -   OTU-7 also has one representative, as does OTU-9.
    -   This OPF is interesting because
        -   It's GH31 containing
        -   It has multiple copies in OTU-1 and is also found in OTU-7
        -   It's found in B. ovatus, but not B. theta
        -   In both OTU-1 _and_ M11 Opf01517 is found in the same PUL as
            Opf01942!  (and Opf01276 is in that same PUL in OTU-1)
        -   _However_: Opf01517 is apparently cytoplasmic in OTU-7 soooo...
-   OTU-7 really does have just one PUL with extracellular alpha-1-4 activity.
    -   And it doesn't even have any GH13 domains.
    -   Is it still able to consume starch, or is it relegated to scavenging malto-oligos?
    -   It _does_ have periplasmic GH13 activity (shared functionality with
        OTU-1), so these could be acting on longer poly-alpha-1,4-glucans
    -   TODO: Consider including OTU-2 in this analysis, since it is less
        abundant with ACA.  Does it only have GH13/31 rather than GH97?

-   TODO: SecretomeP v2.0 ??
-   TODO: Correct Otu0007_vA_01521 sequence (it get's cleaved from the rest of the
    protein)  How should I do this?
        -   I could manually correct the "rfn" version of the genome.
        -   I could try to automate gene combining, look for cases of split domains
            or something?
        -   I could _choose_ a different refinement; give up on some accuracy
            to get that one gene.
        -   I could try to improve the assembly in general by using MEGAHIT
            instead of SPADES (maybe this would be better because of
            highly variable coverage...?)
        -   I could split out the ORF generation step from prokka; what's it
            doing at this point anyway??
        -   The full protein sequence (which certainly COULD be a pseudogene),
            is in core-k161_1293133_pilon_0_116895 and runs from 61922 to 63915
            with a putative deletion just after 62149.
            -   The sequence has two codons on either side of that nonsense
                mutation/error is 'CGA|AAG|-TG|AGC|GGT', with a '-' indicating
                the putative deletion.
            -   I could write a correction script that replaces that (plus
                extra sequence length) with a sequences that has an N in the
                deletion.  This would shift everything back into frame.
            -   That's one way to correct things, but we could also add the N
                anywhere upstream of that point (until it introduces a stop
                codon) That point happens at position approx. 135 in the
                protein (or 62057), where ACG|TAT|GAA|GTG| becomes
                AC|GTA|TGA|AGT|G where TGA is a stop codon.
            -   Basically any erroneous deletion between 62057 and 62149 would
                produce a putative truncation that shouldn't actually be there;
                so we have some options.
            -   Conversely, any REAL deletion in that range would produce a
                REAL trunctation that looks exactly the same.
            -   I'm really counting on the alignments for evidence of this.
            -   If an erroneous deletion shows up somewhere before position
                62149 (but after 62057) then that frameshift will greatly
                change the amino-acids present between those two coords.
                It might therefore form an especially bad match with the HMM.
                We know that the GH97 domain bridges this frameshift.
                The fact that the rest of the protein is "pristine" suggests
                that this isn't a real nonsense mutation and is instead a
                sequencing/assembly error.
            -   When I run a 3-frame tblast on the nucleotide sequence, all
                of the alignments are split at the truncation
                -   Matches to Bacteroides acidifaciens and xylanisolvans (sp?).
                -   Interestingly, the front piece (first 75 AAs) includes
                    a match at the very last position.  This suggests
                    that the frameshift (mutation or error) is actually
                    right ON the stop codon itself.
                -   Conversely, when I run a 3-frame blastx, I get more AAs
                    matching when the new frame starts at position +210 (62267)
                -   This would be odd under either model.
            -   As it currently stands (2018-09-14) I settled on adding an
                'N' before the stop codon in a "final" genome file.
            -   Correction (2018-09-17): This doesn't really help a ton because
                prodigal now fails to include the signal peptide in the protein
                when I make this change.  I think I'll take the BLASTx
                chimeric translation as the correct one, but do all
                my analyses with the original sequence data, just manually
                adding the corrected protein when necessary.
            -   Assuming this protein exists is risky:
                I see no evidence from the sequence data
                that this erroneous deletion happened during assembly.
                It's unclear to me if an error here would be more surprising
                than a pseudogene without any OTHER mutations (since I expect
                mutation to accumulate quickly in pseudogenes).

Corrected sequence region (this may not be completely accurate; for
illustration only; look for 'N'):

```
>core-k161_188734_pilon_0_17445:61921-63916
AAAAAAAGAACCTTAGCTTTGGCCATGATACCTGCCATGACGCTGTGCGGCATCTCGGCA
AAGGAATACCGGGTGACCTCGCCCGACGGGCAACTTAGCGCAATCGTCGAGACCGGGGGG
AAACTGACGTATGAAGTGCAGCTCGAAGGGCACACCGTAATTTCACCTTCCGCCATAGGG
ATGGAGCTTTCGAACGGTGTGAAACTCGGGAAAAACCGAAAGNTGAGCGGTGTGAAACGG
AACAGTGTGGATGAGATGGTGCCATCGCCGTTCTATCGTGCCGAGAGCATACGCGACAAT
TACAACGAGTTGGCGATGAATATAGGGAAGGGGTGGACTTTGGTGTTCCGCGCCTTCAAT
GATGCGGTGGCCTATCGTTTCGTGTGTAAGGATAAGAAACCGTTCGAAATTATAAATGAG
ACTGTTGAATACCGGTTCCCCGCCGATTTCACAGCCACTGCACCTTATGTGCGCTCCGGC
ACCTCCGGCGATTTCGAATCGCAGTTTATGAATTCGTTTGAAAACACCTACACTGTGGCG
CCGGTTTCGCAACTCGACAACGGGCGGCTGGTGTTTCTGCCTATGGCAGTGGAGACTCCC
GAAGGCGTTACGGTGGCTATGACAGAGGCCGCGCTCGAAAATTACCCCGGTCTTTACCTC
AACAATACCGGCATGGGCACAGGCATGAAGGGGGTGTTTGCCAAGCGTCCGAAAACCATG
GAACAGGGTGGCCATAACCGCCTGCAGATGGTTGTGAAAGACCGTGAAAACTTTATTGCG
AAGGTTGAGGGGCCGAGGAGTTTCCCGTGGCGTGTCGCCGTGGTCACCAAAAACGATAAA
GACCTTGCGGCAAGCAATATCAGCTATCTGCTTTCCGACCCTTCACGTGTTGCCGATACC
TCCTGGATCAAGCCCGGAAAGGTGGCATGGGAATGGTGGAACGACTGGAACCTGGAAGGA
GTTGATTTCAAAACCGGTGTCAATAACGAGACCTATAAGGCTTATATCGACTTTGCTGCT
GAGAAAGGTATCGAATATGTGATTCTTGACGAAGGATGGGCTGTGAACCTCCAGGCCGAC
CTGATGCAGGTGGTTGACAACATCAACCTGCCGGAACTCGTGGATTATGCTGCCAAAAAG
GGGGTGGGCTTGATTTTGTGGGCGGGATATTATGCATTCGACCGTGATATGGAGAATGTG
TGCCGACATTATGCCGATATGGGTATAAAGGGGTTCAAGGTCGACTTCATGGACCGCGAT
GACCAGATTATGACCGATTTCAATTATAGGGCAGCCGAAACAGCCGCGCGTCATCACCTT
GTGCTCGACCTTCACGGGACCCATAAACCTGCAGGGATGAACCGCACTTGGCCGAATGTC
CTGAATTTCGAGGGTGTCAACGGCCTTGAGCAACTGAAATGGGCCCAGGAAAAGGATGTG
GACCAGATTAAAAACGATGTTACCATCCCGTTCCTGCGTCAGCTTGCCGGGCCGATGGAT
TATACCCAGGGGGCTATGCTCAACGGGGTGCCGGGCAGCTACCGTCCCAGCAATTCCGAG
CCGATGAGCCAGGGTACACGTTGCCGCCAGCTGGCCTTGTATATGATTTTCGATTCCCCG
CTAAACATGCTTTGCGACAGCCCTAGCAATTACCGTCGCGAGCCGGAATGCACGGGGTTT
ATTGCCGAAGTACCTACGACATGGGATGAGACCCGTGTGCTCAATGCCTCGATGGGGGAA
TATATAGTGACGGCACGCCGCAAGGGCAATTCATGGTATGTGGGCGGTATCACCGACCGC
ACACCACGCGACATCAAAGTTGACCTTTCCTTCCTGGCTCCGGGCAGCCATAAGGCGGTT
ATATTCCGCGACGGGGTTAACGCCCACCGCAAGGGGAGCGATTACAAGCGTGTTGCGAAG
GATGTCAATGCCGATATGGTTCTTGACATGCACCTTGCCCCAGGCGGAGGCTTCGCAATC
AGGATTGACTAA
```
