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
### OTU-1 and OTU-4 (summary of longev paper)
The mammalian gut microbiome is a complex ecological system that influences
energy balance [@TODO], pathogen resistance [@TODO], and inflammation [@TODO],
among other processes with importance to host health.
Understanding how the bacterial inhabitants of the gut respond to environmental
perturbations is a major step in developing a predictive framework for
microbiome-based therapies.
Acarbose (ACA) is an alpha-glucosidase inhibitor that
has been shown to increase the amount of starch entering the lower
digestive system in (TODO: mice/rats?) after a meal [@TODO].

In a recent multi-site experimental study on the effect of ACA on the gut
microbiota and fermentation products, the relative abundance of a number of
bacterial taxa were found to respond to treatment with the drug, as well as the
concentrations of propionate and butyrate.
At two of the three study sites, UM and TJL, one OTU responded dramatically
<!--
TODO: Introduce UT, TJL, UM, OTU abbrevs.
-->
while a different OTU responded at UT.
Both OTUs were classified as members of the recently cultured family
Muribaculaceae&mdash;formerly the S24-7 and in some literature referred to as
Homeothermaceae&mdash;which has been shown to be a common and abundant
inhabitant of the mammalian gut, especially in mice.
Interestingly, these two OTUs are not closely related, sharing only about
90% of their 16S rRNA gene V4 hypervariable region.
While several other Muribaculaceae were also among the most abundant OTUs,
there is not compelling evidence that these taxa also increased in abundance
with treatment.
Identifying shared functional potential between these two responding taxa may
therefore present an opportunity to better understand the traits that poise
some bacterial populations to expand due to ACA treatment, while others
are unaffected or reduced.

Previous studies have suggested that the Muribaculaceae specialize on
the fermentation of complex polysaccharides, much like members of the genus
Bacteroides also in order Bacteroidales.
Based on genomes reconstructed from metagenomes, Ormerod et al. [@Ormerod2016]
identified three distinct carbohydrate utilization guilds, which they describe
as specialists on alpha-glucans, plant glycans, and host glycans, respectively.
While it is reasonable to expect that alpha-glucan specialists would be
most benefitted by the large influx of starch to the gut caused by ACA
treatment, this prediction has not been tested.
Alternatively, the excess of carbohydrate may lead to selection for other
traits that provide access to key resources (e.g. nitrogen).

Interestingly, the OTU that responds at UM and TJL is also common
at UT.
Why it does not show a significant increase in population size at UT,
where the second OTU responds instead, is not clear.
Two alternative may account for the differential pattern of abundance
across sites.
Competition with OTU-4 at UT may inhibit the dramatic expansion
of OTU-1 at that site.
OTU-4 is not found at UM or TJL where OTU-1 responds to ACA.
Alternatively, despite being the same OTU across sites, hidden strain variation
not paralleled by variation in the V4 region may differentiate
OTU-1 at UT from the other two sites.

In this study we apply strain-resolved genome reconstruction from metagenomes
to compare and contrast the functional potential among Muribaculaceae at UM and
UT.
We demonstrate the utility of culture-free genomics to understand
the ecological role of key members of the mouse gut microbial community,
and explore several hypotheses that may explain differences in the distribution
and response of bacteria to perturbations.
<!--
TODO: Figure out better ways to refer to the two OTUs.
-->

## Recovery of Genomes

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
