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
