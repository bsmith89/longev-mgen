---
title: '_Muribaculaceae_ genomes assembled from metagenomes suggest
    genetic drivers of differential response to acarbose treatment in mice'

---

# Abstract

The drug acarbose (ACA) is used to treat
diabetes, and by inhibiting alpha-amylase in the small intestine increases the amount of
starch entering the lower digestive tract.
This results in
changes to the composition of the microbiota and its fermentation products.
Acarbose also increases longevity in mice, an effect that could be related to
increased production of the
short-chain fatty acids propionate and butyrate.
In experiments replicated across three study sites, two distantly related
species in the bacterial family _Muribaculaceae_ were dramatically more
abundant in ACA treated mice, distinguishing these responders from other
members of the family.
Members of the _Muribaculaceae_ likely produce propionate and are abundant and
diverse in the guts of mice, although few isolates are
available.
We reconstructed genomes from metagenomes (MAGs) for eight populations of
_Muribaculaceae_ to examine what distinguishes species that respond positively
to acarbose.
We found two closely related MAGs (B1-A and B1-B) from one responsive species
that both contain a polysaccharide utilization locus with a predicted
extracellular alpha-amylase.
<!--
TODO: Check the following and add if true: "...and these were the only two to
posses a predicted OM-localized GH13 domain containing gene product.
-->
These also shared a periplasmic neopullulanase with another, distantly related
MAG (B2) representative of the only other responsive species.
This gene differentiated these three MAGs from MAGs representative of
non-responding species.
Differential gene content in B1-A and B1-B may be associated with the
inconsistent response of this species to acarbose across study sites.
This work demonstrates the utility of culture-free genomics for infering
the ecological roles of gut bacteria including their response to pharmeceutical
perturbations.

<!-- TODO: Ensure abstract not over word limit. -->


# Importance


The drug acarbose is used to treat diabetes by preventing the breakdown of
starch in the small intestine, resulting in dramatic changes in the abundance
of some members of the gut microbiome and its fermentation products.
In mice, several of the bacteria that respond most positively are classified in
the family _Muribaculaceae_, members of which produce propionate as a primary
fermentation product.
Propionate has been associated with gut health and increased longevity in mice.
We found that genomes of the most responsive _Muribaculaceae_ showed signs of
specialization for starch fermentation, presumably providing them a competitive
advantage in the large intestine of animals consuming acarbose.
Comparisons among genomes support existing models for the ecological niches
occupied by members of this family.
In addition, genes encoding one type of enzyme known to participate in starch
breakdown were found in all three genomes from responding species, but none of
the others.

# Background

The mammalian gut microbiome is a complex ecological system that
influences energy balance [@Turnbaugh2006], pathogen
resistance [@Britton2012], and inflammation [@Syal2018], among other
processes with importance to host health.
Understanding how the
bacterial inhabitants of the gut respond to pharmaceutical and dietary
perturbations is a major step in developing a predictive framework for
microbiome-based therapies.
Acarbose (ACA) is an alpha-glucosidase inhibitor
prescribed for the treatment of type 2 diabetes mellitus because it
reduces the absorption of glucose from starch in the small
intestine [@Hiele1992].
In rats, ACA has been shown to increase the
amount of starch entering the lower digestive system after a
meal [@Dehghan-Kooshkghazi2004].
ACA treatment also changes the composition
of the gut microbiota and its fermentation products in many rodents.
[@Zhao2018; @Holt1996; @Wolever2000; @Dehghan-Kooshkghazi2004; @Weaver1997; @Weaver2000; @Wolin1999; @Zhang2017; @Baxter2019; @Smith2019].
Interestingly, long-term treatment with ACA has been shown to
substantially increase longevity in male mice and to a lesser extent in
females [@Harrison2014; @Strong2016; @Harrison2019].

Previously we have shown that the relative abundance of a
number of bacterial species as well as the concentrations of propionate and
butyrate respond to long term treatment with ACA [@Smith2019].
This study was
notable in being replicated across three sites: The University of Michigan (UM)
in Ann Arbor, The University of Texas Health Science Center at San
Antonio (UT), and The Jackson Laboratory (TJL) in Bar Harbor, Maine.
At
UM and TJL one highly abundant bacterial species was enriched nearly
4-fold in ACA treated mice.
This species, defined at a 97% identity threshold of the 16S rRNA gene V4
region and designated as OTU-1, was classified as a member of the family
_Muribaculaceae_ in order _Bacteroidales_.
OTU-1 was also present and abundant at UT but was not
significantly more abundant in ACA treated mice relative to
controls.
Instead, a different _Muribaculaceae_ species,
designated OTU-4, was found to be highly abundant and 4-fold enriched in
ACA-treated mice, but was nearly absent at UM and TJL.
Other _Muribaculaceae_ were also identified as among the most abundant members
of the mouse gut microbiota across the three sites, although none of
these were found to be enriched in ACA treatment.

The family _Muribaculaceae_---previously referred to as the S24-7 after an
early clone [@s247clone; @Salzman2002],
<!--
TODO: Convert s247clone reference to an inlined URL.
-->
or sometimes as _Candidatus Homeothermaceae_ [@Ormerod2016]---has only
a handful of published cultivars
[@Lagkouvardos2016; @Lagkouvardos2019; @Miyake2020]
despite being a common and abundant
inhabitant of the mammalian gut, especially in mice [@Ormerod2016].
Previous studies have suggested that the _Muribaculaceae_ specialize on
the fermentation of complex polysaccharides [@Ormerod2016], much like
members of the genus _Bacteroides_, which is also a member of the order
_Bacteroidales_.
Genomic analysis has also suggested that the capacity for propionate production
is widespread in the family [@Ormerod2016].
<!--
TODO: Add the Baxter 2019 finding of reduced muri relative abundance
with ACA.
TODO: Add functional findings from other Muri MAG/cultivation studies here?
TODO: Other literature on propionate production?
-->

Recently, techniques have been developed to reconstruct
genomes of uncultivated members of bacterial
communities [@Parks2017; @Lee2017].
<!-- TODO: Add additional Muri mags to this paragraph.
[@Lesker2020; @Lagkouvardos2019; @Miyake2019] -->
Based on 30 such metagenome
assembled genomes (MAGs) they reconstructed using this approach,
Ormerod and colleagues [@Ormerod2016] proposed that the _Muribaculaceae_ fall
into three distinct carbohydrate utilization guilds, which they describe
as specialists on alpha-glucans, plant glycans, and host glycans,
respectively.
While it is reasonable to expect that alpha-glucan specialists
would benefit the most from the large influx of starch to the gut
resulting from ACA treatment, this prediction has not been tested, and
physiological inferences based on the genome content of members of this
clade have been largely divorced from biological observations.

Experimental perturbations of complex microbial communities present an
opportunity to observe ecological features of many bacterial taxa
without cultivated members and generate hypotheses about their
physiology.
Given the observed, dramatically increased relative
abundance of OTU-1 and OTU-4 (here referred to as "responders") in mice
treated with ACA, we hypothesize that these species are capable of robust
growth on starch, while the other _Muribaculaceae_ found in the study
("non-responders"), lack the genomic features necessary for the
utilization of polysaccharides that reach the colon in
greater quantities following ACA treatment.
Alternatively, responders may be
resistant to the inhibitory effects of ACA, or benefit from elevated
levels of intermediate starch degradation products.
Since isolates of
the _Muribaculaceae_ strains in these mice are not available for
characterization, a comparative genomic approach is taken to explore
their functional potential.

Most of the research on the genomic components of polysaccharide
degradation in gram negative bacteria has been carried out in the genus
_Bacteroides_, and in particular _B. thetaiotaomicron_ [@Martens2009].
Starch utilization in _B. thetaiotaomicron_ is dependent on an ensemble
of eight proteins, SusRABCDEFG that enable recognition, binding,
hydrolysis, and import of starch and related
polysaccharides [@Foley2016].
<!--
TODO: Do I need to be more explicit about what polysaccharides make up the
category "starch"?
-->
Homologs of SusC and SusD characterize all
known polysaccharide utilization systems in this clade [@Grondin2017],
are encoded in Sus-like genomic regions known as polysaccharide
utilization loci (PULs), and are widespread in the phylum
_Bacteroidetes_ [@Fernandez-Gomez2013].
The molecular range of these
systems is determined by the carbohydrate-active enzymes and structural
proteins they encode, based on the specificity of glycoside hydrolase
(GH) and carbohydrate binding module (CBM) domains, which have been
extensively cataloged in the dbCAN database [@Yin2012; @Zhang2018].

Here MAGs from the feces of mice at UT and UM are analyzed to explore
two closely related questions about the niche of OTU-1 and OTU-4 in the lower
digestive system.
First, why do these species each increase in relative abundance with ACA
treatment, while other species of _Muribaculaceae_ do not?
And second, why is the
response of OTU-1 site specific?
Despite similar patterns of abundance at
their respective sites, the two responding species seem to be only distantly
related, sharing just 90% of nucleotides in their 16S rRNA gene V4
hypervariable region [@Smith2019].
We nonetheless find genomic evidence that OTU-1
and OTU-4 occupy overlapping niches, specializing in the degradation of
alpha-glucans, a role not held by the other _Muribaculaceae_ described in
this study.
In addition, we identify two distinct genomic variants of OTU-1,
referred to as B1-A and B1-B, which are differentially distributed
between UM and UT and have functionally relevant differences in gene
content.

Reconstructing genomes from metagenomes allows for the comparison of the
functional potential of _Muribaculaceae_ at UM and UT.
This work
demonstrates the utility of culture-free genomics to understand the
potential ecological roles of these key members of the mouse gut microbial
community and explore several hypotheses that may explain differences in
the distribution and response of these bacteria to ACA treatment.
Hypotheses
derived from this analysis provide a foundation for future physiological
studies in recently obtained cultivars.
While a large fraction of
host-associated bacterial species are without isolated representatives
[@Lagkouvardos2017],
let alone characterized [@Stewart2012],
combining experimental data from complex
communities with the analysis of reconstructed genomes provides a
powerful tool for expanding understanding to these understudied taxa.

# Results

## Recovered population genomes are of high quality and resemble other _Muribaculaceae_ genomes

MAGs were constructed for 8 populations classified as members of the family
_Muribaculaceae_,
including for two species, OTU-1 and OTU-4, previously shown to
respond positively to ACA.
For OTU-1, two closely related genomic variants were recovered,
here designated B1-A and B1-B,
possessing 0.63 and 0.36 Mbp of unshared sequence, respectively
(Table 2).
We designate the MAG constructed for OTU-4 as B2.
MAGs obtained from non-responding species are designated B3 through B7.
All 8 novel MAGs are estimated to be
mostly complete and all had less than 1% estimated contamination
based on the recovery of ubiquitous, single-copy genes (Table 1).
The median N50
statistic was approximately 71 kbp, suggesting that assembly was
suitable for inferring the genomic context of functional genes.
Estimated genome sizes, GC%, and number of predicted genes are
all similar to previously published MAGs belonging to the
family _Muribaculaceae_,
as well as the finished genome for _Muribaculum intestinale_ strain YL27.

+----------+-----------------+-------------+-----------+----------+---------+---------+-------------------------+
| Genome   | Completeness^1^ | Scaffolds   | Length^2^ | N50      | GC      | in      | Top nr BLAST            |
|  / MAG   |                 |             |           |          |         | Smith   | Hit                     |
|          |                 |             |           |          |         | et al., | (Identity)              |
|          |                 |             |           |          |         | 2019    |                         |
+:=========+================:+============:+==========:+=========:+========:+:========+=========================+
| YL-27^3^ | 99%             | 1           | 3.3       | 33,070   | 50.1%   | _n/a_   | _n/a_                   |
+----------+-----------------+-------------+-----------+----------+---------+---------+-------------------------+
| B1-A     | 97%             | 228         | 3.2       | 41,412   | 46.6%   | OTU-1   | WP_123406077.1 (99.92%) |
+----------+-----------------+-------------+-----------+----------+---------+---------+-------------------------+
| B1-B     | 97%             | 152         | 3.0       | 59,916   | 46.9%   | OTU-1   | WP_123406077.1 (100%)   |
+----------+-----------------+-------------+-----------+----------+---------+---------+-------------------------+
| B2       | 98%             | 65          | 2.6       | 79,454   | 50.5%   | OTU-4   | WP_128713622.1 (92.52%) |
+----------+-----------------+-------------+-----------+----------+---------+---------+-------------------------+
| B3       | 86%             | 98          | 2.2       | 63,818   | 54.0%   | OTU-6   | OKY86749.1 (90.84%)     |
+----------+-----------------+-------------+-----------+----------+---------+---------+-------------------------+
| B4       | 98%             | 31          | 2.7       | 148,039  | 55.2%   | OTU-5   | WP_123486179.1 (91.80%) |
+----------+-----------------+-------------+-----------+----------+---------+---------+-------------------------+
| B5       | 86%             | 50          | 2.5       | 78,179   | 55.7%   | OTU-8   | WP_123486179.1 (92.43%) |
+----------+-----------------+-------------+-----------+----------+---------+---------+-------------------------+
| B6       | 99%             | 110         | 3.2       | 87,115   | 48.3%   | OTU-30  | WP_123613567.1 (100%)   |
+----------+-----------------+-------------+-----------+----------+---------+---------+-------------------------+
| B7       | 98%             | 97          | 2.5       | 59,037   | 53.9%   | OTU-39  | WP_123541885.1 (100%)   |
+----------+-----------------+-------------+-----------+----------+---------+---------+-------------------------+

: Table 1: Summary of novel MAGs compared to the genome of _Muribaculum intestinale_ YL27

^1^ Estimated by CheckM [@Parks2015] \
^2^ Total length in Mbp \
^3^ _Muribaculum intestinale_ YL-27 reference genome

To confirm the assertion that each of the reconstructed genomes
is representative of a previously described _Muribaculaceae_ species identified
in these mice [@Smith2019],
the per-library mapping rates of each genome were
compared to the relative abundance of the associated 16S rRNA gene in
amplicon libraries.
Pearson correlation coefficients between the fraction of reads
mapped and species relative abundance were above 0.86 for all MAGs,

![
Figure 1: Comparison of novel and previously described
_Muribaculaceae_ genomes.  Novel MAGs ("B1-A", "B1-B", and "B2" through "B7")
are combined with the finished genome for _M. intestinale_ strain YL27, as well
as 30 previously described MAGs that are
hypothesized to reflect three polysaccharide utilization guilds: specializing on
alpha-glucans (points and labels colored blue), host glycans (violet), and plant glycans
(green) [@Ormerod2016].
(**A**) MAGs described in this study were placed in a phylogenetic context using a
maximum-likelihood concatenated gene tree based on an amino acid alignment of 9 shared,
single-copy genes, and four other _Bacteroidales_ species as an outgroup (not shown).
Nodes with
less than 70% confidence are collapsed into polytomies and topological support greater
than 95% is indicated (black dots).
Branch lengths indicate an estimate of expected
substitutions per site.
(**B**, **C**) Functional comparisons were visualized by plotting
the first two principal components of an ordination on annotation counts of either
(**B**) eight COGs previously identified as maximally discriminatory between
hypothesized guilds [@Ormerod2016], or (**C**) _de novo_ clusters based on sequence similarity of GH
domain containing proteins.
PCA was performed with the 30 previously constructed MAGs from family _Muribaculaceae_
and the percent of variation described by the first two
components is included in the axis labels.
All genomes were then projected onto that space.
Novel MAGs (black triangles) are labeled, as are the previously described MAGs M1,
M6, and the proposed H. arabinoxylanisolvens (Ha),
<!-- TODO: Discuss w/ Tom H. arabinoxylanisolvens is a species name without
standing given by Ormerod et al. -->
and the finished genome of M.
intestinale (Mi, gray circle).
](fig/muri_comparison.pdf)

### Phylogenetics

To better understand the evolutionary relationships between these
organisms, a concatenated gene tree was constructed to assess relationships
among the 8 new MAGs and those previously described [@Ormerod2016], and
_M. intestinale_ YL27.
The tree was rooted by four other _Bacteroidales_
species: _Bacteroides ovatus_ (ATCC-8483), _Bacteroides
thetaiotaomicron_ VPI-5482, _Porphyromonas gingivalis_ (ATCC-33277), and
_Barnesiella viscericola_ (DSM-18177).
Most internal nodes were supported with
high topological confidence (>95% bootstrap support), and the placement of the
MAGs derived by Ormerod and colleagues was highly consistent with their
published tree.
To further check the robustness of our phylogeny, a second maximum
likelihood tree was constructed based on the _rpoB_ gene, which is generally
not thought to be transmitted horizontally (despite exceptions [@Kim2013]),
This approach also recapitulated the published topology (See Supplementary
Results).
<!--
TODO: Supplementary results
-->
The estimated phylogeny shows
that the newly generated MAGs encompass most of the
documented diversity of _Muribaculaceae_.
Two of the new MAGs, B2 and B6, appear to be closely related to previous
MAGs M6, and M1, respectively [@Ormerod2016].
Nonetheless, this phylogenetic analysis suggests that the majority of the
MAGs derived here had not been described at the time of this study.
To further link our MAGs with recently generated MAGs and isolate
genomes, we BLASTed translated _rpoB_ sequences against
the GenBank non-redundant protein database (updated 2020-06-22).
Accession numbers for the closest hits as well as the amino acid
identity are reported in Table 1.
Despite a growing number of _Muribaculaceae_ genomes
deposited in public repositories, four of the eight MAGs reconstructed
here have rpoB genes with less than 93% amino acid identity to their highest
scoring match.
<!--
TODO: Add information about
 <https://www.ncbi.nlm.nih.gov/assembly/GCF_004102775.1/> which has perfect
 matches to (at least) B1A_279 - 284 (one of the starch PULs).
 This means that I have to repeat that part of the analysis, because
 I'll definitely have to include this genome (and other isolates/MAGs)
 in all of my comparisons.
TODO: Do all of the above statements still hold?
-->

### Novel protein families

Annotations based on alignment to a database of previously characterized
sequences may provide only limited insight into the function of gene products,
in particular for genes
from largely unstudied families of bacteria.
In order to identify
previously uncharacterized orthologous groups of genes, _de novo_
clustering [@Schloss2008] was carried out based on amino acid similarity
of all putative genes found in the 8 novel MAGs, 30 previously
generated MAGs, _M. intestinale_, four publicly available draft
genomes from the family, and four reference _Bacteroidales_
(see Fig. 1).
<!-- TODO: Long term: Add all available other genomes/MAGs to this analysis. -->
The
resulting clusters are referred to as operational protein families
(OPFs).
While a fraction of the 12,648 resulting OPFs may be due to
spurious sequence similarity and without biological relevance, 5,767 had
representatives in at least three genomes, increasing the likelihood
that these reflect evolutionarily conserved protein sequences.
Of these,
only 2,404 had members annotated with any COG, KO, or putative function.
The remaining 3,363 unannotated OPFs encompass 17,831 predicted protein
sequences across the 47 genomes.
Annotations of predicted genes in MAG and reference genomes with OPFs, COGs,
KOs, and Pfam and CAZy domains are all available in the Supplementary Results.
<!--
TODO: Supplementary Results
-->>

### Ordination of gene content

To compare the novel MAGs to other available MAGs and reference genomes,
a previous published
analysis was recreated, harnessing a set of 8 COGs found by Ormerod and
colleagues to
maximally differentiate the three hypothesized guilds [@Ormerod2016].
We replicated the original ordination analysis using our annotation of the
publicly available MAGs, and then projected all other genomes onto this same
space (see Fig. 1).
Newly available genomes were compared
to the three clusters hypothesized to represent specialization on
alpha-glucans, plant glycans, and host glycans.
While the 8 novel MAGs
inhabit approximately the same volume as those in the original analysis,
and some could be plausibly classified based on these criteria, the
ambiguous placement of B4 and _M. intestinale_ suggests that new genomes
will present additional exceptions to the three-guild model.

It is notable that MAGs from responding species cluster with the proposed
alpha-glucan guild, consistent with a functional potential for starch
utilization absent in the non-responders.
To expand on this descriptive analysis
and to leverage the more comprehensive view provided by _de novo_
clustering to explore differences and similarities in carbohydrate
utilization potential, a second ordination of genomes was performed,
this time based on OPF labels of predicted genes found to contain GH
domains (Fig. 1).
Similar to the previous ordination based
on COGs, three groups of genomes approximately reflecting those proposed
by Ormerod and colleagues are apparent.
However, the placement of B2 (as well
as the closely related M6) relative to the proposed guilds are
substantially different.

<!--
TODO: Consider adding some results on the ability of these organisms to produce
propionate.
-->

## Analysis of MAGs from species responsive to ACA treatment suggests genes involved in starch utilization.

Based on the characterization of genes and genomic regions with a role
in starch utilization in the closely related genus _Bacteroides_, it is
plausible that alpha-amylase localized to the outer membrane may be common to
starch utilizing bacteria in the order _Bacteroidales_ [@Shipman1999].
Indeed, B1-A and B1-B both have three OM-localized genes predicted to code for
GH13 containing lipoproteins (B1A280, B1A301, B1A333 in B1-A and B1B...TODO...in B1-B),
<!--
TODO: Add gene numbers for B1B GH13 containing lipoproteins.
-->
each in a separate PUL (see Fig. 2).
<!-- TODO: Update figure references and be consistent. -->
While it also includes members without
this activity, GH13 is the main family of alpha-amylases [@Janecek2014].
These genomic regions also possess additional genes with
carbohydrate-active domains that are expected to interact with alpha-glucans.

![
Figure 2: Polysaccharide utilization loci in _Bacteroidales_.
Diagrams of the Sus
operon (**A**) and the dextran associated PUL (**B**) of _B. thetaiotaomicron_
VPI-5482 along with
five putative starch-associated PULs identified in _Muribaculaceae_ MAGs
B1-A and B4
(**C**-**G**).
B1-A PULs shown here are syntenic in B1-B (not shown).
Predicted protein coding sequences are shown as boxes pointed in the direction
of transcription.
Homology to SusC, SusD, and SusEF is indicated.
Protein regions
with homology to starch-associated GHs, as well as GH66, and CBMs are shown
as shallow rectangles, and are colored as indicated in the legend.
Several OPFs
are noted with members in multiple genomes, including clusters that contain SusR
(Opf01144), SusA (Opf01391), and SusB (Opf00018).
The inferred localization of
each protein product is also indicated: cytoplasmic (genes labeled C), periplasmic
(**P**), outer membrane (**O**), or inner membrane (**I**).
](fig/pul_diagrams.pdf)


Besides B1-A and B1-B, B5 is the only other MAG to possess a putative PUL
coding for a full complement of predicted starch-active proteins.
Several of the OPF annotations found in these presumptive starch PULs are
shared by _B. thetaiotaomicron_, suggesting shared function.
This set including SusC-homologs
Opf01277, Opf02066, which includes relatives of SusD, and Opf02791 whose
members possess CBM20 starch-binding domains.
However, while B5 also has
a GH13 containing lipoprotein (B51713), its predicted localization is on
the inner membrane.
It is unclear whether this explains B5's
non-response in ACA-treated mice.
Plausible OM-localized, GH13
containing proteins are not found in any non-responders.
While this
characteristic does not seem to perfectly discriminate responder from
non-responders---B2 also lacks such a gene---it nonetheless
demonstrates concordance between inferred genomic features and observed
population dynamics of the corresponding species.

Despite the absence of a GH13 domain on the outer-membrane, it is
plausible that B2 is capable of degrading starch using other enzymatic
machinery.
We speculate about one putative locus
(see Fig. 2F),
which has a similar gene content to characterized
[@Ravcheev2013; @Rogers2013; @VanBueren2015]
dextran PULs in _B. thetaiotaomicron_ and _B. ovatus_.

To expand the search for relevant genetic features, _de novo_ protein
clusters were filtered to those with members in B1-A, B1-B, and
B2.
Of these OPFs, several stood out as particularly relevant.
Opf01144
includes SusR, the regulator of transcription of the starch utilization
system in _B. thetaiotaomicron_, as well as its homolog in _B. ovatus_.
It is an apparent subcluster of the larger family defined by K21557, and
in many cases is encoded directly upstream of _susC_ in putative PULs
that are considered likely to have affinity for alpha-glucans.
In B1-A and B1-B, two of the
three putative starch PULs encode a member of Opf01144, and it is
similarly located in PULs with starch-active CBM and GH domains in B2
and B5.
In addition, of the seven MAGs constructed by Ormerod _et al._
that encode a member of this cluster, five of them are classified to the
alpha-glucan guild.
It is plausible that members of Opf01144 share a
functional role regulating transcriptional responses to alpha-glucans.

Opf01391, which recapitulates K21575, includes SusA: the periplasmic
neopullulanase of _B. thetaiotaomicron_ and an important component of
starch utilization in that organism [@DElia1996].
This family is found
in the MAGs associated with responders, B1-A, B1-B, and B2, and none of the
other MAGs generated in this study.
What's more, it's found in twelve of the thirteen
alpha-glucan and a minority of the plant glycan guild members.
Interestingly,
although it is encoded by the Sus operon in _B. thetaiotaomicron_ and
its homologous locus in _B. ovatus_, in the _Muribaculaceae_ members of
Opf01391 do not in general appear to be encoded in PULs.

## Unshared gene content in B1-A and B1-B

Two distinct genomic variants were associated with OTU-1
with one found in a majority
of the UT mouse metagenomes, and the other ubiquitous at UM.
Using the
nucmer tool for genome alignment [@Delcher2002], 19.6% of the B1-A MAG
sequence and 12.2% of B1-B were found to not align to the other.
While
these hundreds of kbp may in part reflect errors in genome recovery,
much of the unaligned length suggests differences in gene content
between these populations of OTU-1.
This observation was confirmed
by assessing the mapping of metagenomic reads against predicted protein
coding genes in each variant.
For each pairing of metagenomic read
library to genomic variant, gene coverage was normalized by the median
gene coverage in order to identify genes with conspicuously fewer reads
in particular subsets of the mice.
Libraries have low coverage of large
portions of either the B1-A or B1-B MAG (see
Fig. 3), suggesting that mice are primarily inhabited
by one of the two variants, and that a portion of genes are variant
specific.

![
Figure 3: Visualization of differential gene content in two OTU-1 populations.
Heatmaps
depict mapping coverage of metagenomes against putative protein coding genes in
MAGs B1-A or B1-B normalized to the median coverage.
Rows represent one or
more pooled libraries for each mouse included in the study and columns
represent individual genes.
The site at which each mouse was housed is
indicated by triangles in the far left column: UT (green, left pointing) or UM
(blue, right).
Filled triangles correspond to those mice
representative of just B1-A or just B1-B (not a mixture) flagged for downstream
analysis.
Genes are shown
only where the median normalized coverage ratio between these B1-A and B1-B
specific metagenomes is greater than 1.5.
Rows and columns are arbitrarily
ordered to maximize visual distinction between variants.
](fig/b1_vars.pdf)

Metagenomic libraries manually chosen as unambiguous representatives of
either B1-A or B1-B were used to systematically identify genes
differentiating the two.
The median normalized mapping depths in each
set of libraries against predicted genes in each MAG were compared,
providing a measure of the relative enrichment or depletion of genomic
sequences between the two populations of OTU-1 (see Supplementary Results).
<!--
TODO: Supplementary Results
-->
This analysis found 12.8%
of predicted genes in B1-A were depleted at least 5-fold in B1-B
populations, and 12.4% the reverse.
While this observed depletion could
indicate variation in copy number, differential gene content between
variants is a more parsimonious explanation for most loci.
These
predicted genes reflect 2.7% of unique KOs in B1-A and 1.9% in B1-B.
Interestingly, the fraction of variant specific OPFs is greater, 7.5%
and 7.1% respectively, suggesting that _de novo_ clustering could be
more sensitive to potential differences in physiology.

+:-----------------------------------------------+:------+:---------+:------+:---------+
|                                                |B1-A              |B1-B              |
+:-----------------------------------------------+:------+:---------+:------+:---------+
|                                                | Total | Specific | Total | Specific |
+------------------------------------------------+-------+----------+-------+----------+
| Nucleotide length^1^                           | 3.23  | 0.63     | 2.96  | 0.36     |
+------------------------------------------------+-------+----------+-------+----------+
| Genes                                          | 2,710 | 348      | 2,496 | 309      |
+------------------------------------------------+-------+----------+-------+----------+
| OPFs^2^                                        | 2,308 | 173      | 2,202 | 157      |
+------------------------------------------------+-------+----------+-------+----------+
| KOs^2^                                         | 1,056 | 29       | 1,033 | 20       |
+------------------------------------------------+-------+----------+-------+----------+
| COGs^2^                                        | 716   | 8        | 709   | 3        |
+------------------------------------------------+-------+----------+-------+----------+
|                                                |       |          |       |          |
+------------------------------------------------+-------+----------+-------+----------+

: Table 2: Summary of variant specific features in two highly similar MAGs

^1^ in Mbp \
^2^ unique

Given the observation that the relative abundance of OTU-1 was dramatically
increased with ACA treatment at UM, while not being significantly
affected at UT, and that B1-B was not found in metagenomes at UM, we
hypothesized that differences in the genomic potential
of B1-A and B1-B could explain the different response to ACA
at the two sites.

Genomic regions apparently specific to B1-A---defined as an at least
5-fold enrichment in B1-A specific libraries relative to
B1-B specific libraries---include just one PUL
(SusC-homolog encoded by B1A00048).
This locus includes a predicted outer membrane localized GH30-containing
protein.
Proteins that contain a GH30 domain have
beta-glucosylceramidase, beta-1,6-glucanase, or beta-xylosidase
activity [@StJohn2010].
Given that this PUL also encodes a periplasmic,
GH3 containing protein, it appears to be unlikely that it has
specificity for starch.
The B1-A also possesses numerous phage
insertions not seen in B1-B.
Conversely, a CRISPR
operon including 25 repeat units (Cas9 encoded by B1B01367) appears to
be specific to B1-B.

Most strikingly, a 16 kbp region (from B1A01498 to B1A01514) specific to
B1-A was found to contain several genes with homology to cell capsule
and exopolysaccharide synthesizing enzymes.
Based on annotations with
KEGG orthologous groups, these include homologs of _tuaG_ (K16698),
_tagE_ (K00712), _gmhB_ (K03273), _gmhA_/_lpcA_ (K03271),
_hddA_ (K07031), _exoO_ (K16555), _waaH_ (K19354), and _tagF_ (K09809).
Interestingly, the B1-B MAG contains a different such region of about
6.5 kbp (B1B00851 to B1B00856) with _wfeD_ (K21364), _pglJ_ (K17248),
and _epsH_ (K19425).
For each, several of the OPFs in the respective
regions were not found anywhere in the opposing genome, suggesting that
the makeup of each variant's exterior surface might be distinctly
different.

# Discussion

Mice are a key model system for study of the mammalian gut microbiome,
with an outsized importance in testing mechanistic hypotheses for the
roles of this community in host health [@Nguyen2015].
The
generalizability of observations made in mice is a constant
concern [@Nguyen2015], in part due to extensive difference in taxonomic
composition compared to humans [@Lagkouvardos2016].
Bacteria classified in the family
_Muribaculaceae_ are abundant in the murine gut
microbiome [@Ormerod2016].
While these bacteria are also found in humans
(although at lower abundance), only a few members of this clade
have been cultivated and described
[@Lagkouvardos2016; @Lagkouvardos2019; @Miyake2020].
As a result, the ecological
roles of these bacteria have not yet been characterized,
limiting the value of the mouse as a model system.
Better understanding the ecology of _Muribaculaceae_ in the murine gut
will increase the transferability of microbiome findings from mice to humans.
Attempts to study these
organisms make use of genomes reconstructed from metagenomic reads, and
have suggested---in the absence of experimental data---that members of
the family consume a diversity of polysaccharides in the lower gut.

Here we have extended that approach to eight new genomes, and associated
those with species for which changes in relative abundance in response to
ACA treatment have been experimentally assessed.
This enabled us to
explore why two responding species, represented by MAGs B1-A, B1-B, and B2,
increase with ACA treatment, while
the other species of _Muribaculaceae_ do not.
Annotations of reconstructed genomes
suggest that the responders may possess starch degradation
capabilities absent in the non-responders.

We examine the three-guild model proposed by Ormerod and
colleagues [@Ormerod2016] by reproducing their dimensional reduction approach
with the addition of these new genomes.
In this analysis, annotations of B1-A, B1-B, and B2
are consistent with a hypothesized guild of alpha-glucan
degrading species, supporting their interpretation.
A more nuanced
approach to annotation was also applied by constructing _de novo_
clusters of proteins based on homology.
Interestingly, this analysis
indicates that B2, and the closely related M6, share physiological
potential with MAGs in the host-glycan guild, suggesting that a more
detailed examination can identify specific functions that discriminate
responders from non-responders.
This approach is bolstered by the
phylogenetic and genomic distinction between B2 and both B1-A and B1-B,
reducing the confounding effects of shared evolutionary history.


By including otherwise unannotated genes, genomic comparisons based on OPFs may
reflect shared functional potential better than applying previously defined
orthologies.
Besides the identification of potentially
novel gene families, _de novo_ homology clustering [@Schloss2008] also
enables differentiation of sub-groups not captured by standard
annotations.
For instance, hypothetical genes annotated as homologs of
SusC, SusD, and SusEF, were clustered into 119, 162, and 33 different
OPFs respectively.
It is plausible that this sub-clustering captures
differences in protein structure with importance in oligo- and
polysaccharide recognition, import, and binding.
Combined with annotation of characterized functional domains, these clusters
more narrowly predict the polysaccharide utilization ranges of uncultured
organisms.
Testing these predictions will require characterization of the metabolic
potential of these genes after obtaining cultivars or through heterologous
expression in appropriate hosts.
<!-- TODO: Eventually, I'd like to add references here to the B1-like isolate
that I just found: <https://www.ncbi.nlm.nih.gov/assembly/GCF_004102775.1/> -->

A detailed analysis of PULs identified multiple loci shared in both B1-A and
B1-B that appear
to be adapted to the degradation of starch or related carbohydrates, due
to the presence of an OM localized GH13 containing
protein [@Koropatkin2010].
Counterintuitively, B2 had no such PUL,
suggesting that its response to ACA may result from other enzymatic
capabilities.
Of particular interest is a PUL encoding proteins with
GH97, CBM20, and CBM69 domains, all of which have documented activity on
starch [@Naumoff2005; @Boraston2004].
While the only outer-membrane
localized hydrolase in this PUL is a GH66, and members of this family
have characterized activity on the alpha-1,6 linkages between glucose
monomers in dextran [@Kim2012],
it is plausible that this PUL can be
repurposed and confers some ability to grow on starch.
<!-- TODO: Figure out if [@Chaudet2016] show that the PUL in B.t. that is
closely related to B2's putative starch PUL is active on starch.
I don't THINK so.  I think it's a different PUL.
How to reference this study?
-->

In addition, a gene encoding a SusA homolog was identified in B1-A, B1-B,
and B2 but in none of the non-responders.
While it is unclear how
expression of this important component of starch utilization might be
regulated, given that it is not located in a PUL in any of the
responding populations, SusA is important for growth on amylopectin in
_B.  thetaiotaomicron_ [@DElia1996].
Since inhibition by ACA is variable
across enzymes [@Kim1999], it is possible that ACA treatment
results in elevated production of dextrin and maltooligosaccharides in the
lower guts of mice due to residual alpha-amylase activity, even at levels
sufficient to prohibit host digestion.
Periplasmic hydrolysis of these
starch breakdown products may be sufficient for increased abundance of
these species in ACA treated mice.

It is notable that of the closely related variants, B1-A and B1-B
associated with OTU-1,
B1-B is found at UT and not UM.
We previously observed site-specificity of the ACA response of
this species, in which OTU-1 did not have a significantly
increased abundance in treated mice at UT, while it was the most
dramatic change at UM.
Differences in the functional potential due to differences in gene content of
populations found at each of the sites is one possible explanation for this
pattern.
Intriguingly, while we do not conjecture a mechanistic link,
an ACA-by-site interaction effect on longevity has been
previously observed in the mouse colonies sampled here,
with male mice at UT showing a larger increase in median longevity with ACA
treatment than those at UM [@Harrison2014; @Harrison2019].

Despite evidence that
large differences in gene content can be found between
even closely related populations [@Rasko2008; @Medini2005],
studies reconstructing genomes from
metagenomes have just started to consider these pangenome dynamics
[@Scholz2016; @Truong2017; @Delmont2018; @Segata2018a].
The discovery of two
populations of OTU-1 therefore demonstrates the value of considering
pangenome dynamics, and presents a potential explanation for the
observed site-specific response of that taxon.
The finding that both
variants have the same complement of three PULs apparently specializing
in starch utilization and the same SusA homolog does not support the
hypothesis that differences in starch utilization potential account for
these abundance patterns.
We did, however, identify numerous differences
in the gene content of B1-A and B1-B, including variant specific loci
that may influence the structure and function of the outer surface of
the cell.
Capsule variation is known to greatly affect both ecological
and host interactions [@Merino2015].

While these results do not establish a mechanistic explanation for
differences in the response of B1-A, B1-B at UM and UT, nor conclusively
identify the pathways that enable starch utilization in B2,
they do suggest a number of
genomic features that likely contribute to previously observed patterns
in taxon abundance.
Future studies utilizing metatranscriptomic analysis
might demonstrate active expression of these genes, or differential
expression in mice treated with acarbose compared to controls.
Likewise,
even in the absence of a B2 cultivar, the potential role of its dextran PUL
in enrichment under acarbose treatment could be tested using
available cultivars, like _B. thetaiotaomicron_, that posses a homologous gene
cluster.

# Conclusions

In this study we have reconstructed and described genomes representing 7
species in the family _Muribaculaceae_ from the mouse fecal microbiome, and
have found features that differentiate those bacterial species
that respond positively to ACA treatment from those that do not.
This analysis suggests that
utilization of starch and related polysaccharides enables increased
population size in mice treated with ACA---an alpha-amylase inhibitor.
In
addition, two distinct genomic variants of one species were identified
that differ in functional gene content, potentially explaining
site-specific differences in response.
By combining observed changes in
relative abundance during experimental manipulation with inferred
functional gene content, we are able to study mammalian symbionts in the
absence of cultured representatives.
This sequence-based approach is
broadly applicable in microbial ecology and enables improved
understanding of _in situ_ dynamics within complex microbial
communities.

# Methods

## Mouse treatment, sample collection, extraction and sequencing

Mice were bred, housed, and treated as described in [@Harrison2014].
Briefly, genetically heterogeneous UM-HET3 mice at each study site were
produced by the four-way cross between (BALB/cByJ x C57BL/6J) F1 mothers and
(C3H/HeJ x DBA.2J) F1 fathers, as detailed in [@Miller2011].
Mice were fed LabDiet 5LG6 (TestDiet Inc.) from weaning onwards.
Starting at 8 months
of age, mice randomly assigned to treatment were fed chow with 1,000 ppm
ACA (Spectrum Chemical Manufacturing Corporation).
Mice were housed 3
males or 4 females to a cage.
Colonies were assessed for infectious
agents every 3 months, and all tests were negative.

Individual fecal pellets were collected from a single mouse per cage.
16S rRNA gene libraries and metabolite analyses of these samples are
as described previously [@Smith2019].
From this collection, a subset of samples were
non-randomly selected for metagenomic sequencing in order to test
several unrelated hypotheses about SCFA production.
<!-- TODO: Is this still too ambiguous? Too specific? -->
Samples were from 54 mice, with at least six treated and
control representatives of both males and females at each site.

Fecal samples were slurried with nuclease free water at a 1:10 (w/v)
ratio, and most samples were spiked with _Sphingopyxis alaskensis_
RB2256 prepared as described previously [@Smith2019] before DNA extraction and
sequencing.
Based on alignment to the reference genome, sequenced reads from
_S. alaskensis_ can be distinguished from all endogenous bacteria in mouse
feces.
This spike was added as an internal standard to quantify total 16S rRNA
gene abundance, and also provides a benchmark for the reconstruction
of bacterial genomes.
<!-- TODO: Also discuss S.a. assembly quality? Perhaps skim over this by
explaining that spiking was done for experimental consideration not dealt
with in this manuscript. -->
A small number of these were split for both spiked and unspiked samples,
which we used to validate this procedure.
For each, 150 uL of this
sample was transferred for extraction using the MoBio PowerMag
Microbiome kit.
Metagenomic libraries were prepared using standard
procedures sequenced on the Illumina HiSeq 400 platform using the v4
paired-end 2x150 bp.

## Assembly, binning, and MAG refinement

Raw metagenomic reads were deduplicated using FastUniq [@Xu2012],
adapters trimmed using Scythe [@Buffalo2018], and quality trimmed using
Sickle [@Joshi2011] to produce processed reads for all downstream
analyses.
The resulting paired-end reads were assembled into primary
contigs using MEGAHIT [@Li2014].
Reads were then mapped back to these
contigs with Bowtie2 [@Langmead2012], and per-library coverage was
estimated for each contig.

For all contigs \>1000 bp in length, dimensional reductions built into
CONCOCT [@Alneberg2014] were applied to produce input data for a
Gaussian mixture model (GMM) similar to the procedure used by that
program for binning.
However, unlike CONCOCT---due to computational
limitations---the model was trained on just 10% of the input data,
sampled randomly, before assigning bins to all contig.
While this may
have reduced the accuracy of the binning procedure, we believe that
subsequent refinement steps mitigated the impact of this decision.

OTUs were classified taxonomically and relative abundance was calculated
for matched libraries as described in [@Smith2019].
Bins were then recruited to one or more OTUs by
calculating a Canonical partial least squares between OTU abundance and
bin coverage as implemented in the scikit-learn machine learning library
for Python [@Pedregosa2012].
For bins recruited to OTUs classified as
_Muribaculaceae_, contigs were re-clustered based on coverage across
samples.
First "trusted contigs" were manually selected which correlated
closely with OTU abundance.
The mean coverage of these was used to
normalize the per-library coverage of all other contigs.
Then, using a
GMM, groups of contigs were clustered such that the normalized coverage
across samples was consistent.
These groups were used to inform the
manual assignment of contigs to MAGs.
Libraries in which MAGs had
non-negligible coverage were identified and used in subsequent
refinements.
While clustering contigs associated with OTU-1 a number of
groups containing on the order of $10^5$ bp were found with a bimodal
coverage distribution, low normalized coverage in a subset of libraries,
and normal coverage in others.
By this criterion, contigs in these
"variable" groups were partitioned into two MAG variants, A and B, with
the non-variable contig groups shared by both.
To avoid challenges associated with admixture, only libraries that appeared on
further inspection to have just one of the two variants were considered
in downstream refinement steps.
The mice matching these libraries are
highlighted in Fig. 3.
Genomic variants were not found associated with any of the other
_Muribaculaceae_ OTUs described in this study.

For each MAG, several alternative refinement procedures were performed
from which the best quality result was selected.
Reads mapping to the
curated contigs were digitally
normalized [@Wedemeyer2017; @Brown2012a; @Zhang2014a] and reassembled
with SPAdes [@Bankevich2012].
This reassembly as well as the original
contigs were cleaned using a single pass of the Pilon assembly
refinement tool [@Walker2014a].
Finally, the per-library mapping depths
of each position in these assemblies were compared to the average
mapping depth of the "trusted contigs" selected earlier, and regions
with low cosine similarity were excised from the final assemblies.

Genome completeness and contamination estimates were calculated based on
ubiquitous single-copy genes using the program CheckM [@Parks2015].
Based on these results, the final assembly with the highest completeness
and with contamination \< 1% was selected from the various refinements.

## Reference genomes

The _Muribaculum intestinale_ genome sequence was obtained from GenBank
(accession GCA002201515.1), as well as four additional draft genomes
(GCA003024805.1, GCA003024815.1, GCA002633305.1, GCA002633115.1)
enabling comparison of MAGs to cultured isolates.
<!-- TODO: Does the analysis need to be updated due to the addition of other
isolate genomes? -->
While other genomes labeled as _Muribaculaceae_ had also been deposited at the
time of this work, they were excluded from this analysis due to redundancy or
apparent misidentification to the family.
More recently, several other isolate genomes have become available but have
also not been included here.
Thirty previously constructed MAGs [@Ormerod2016] were obtained from the SRA.
For
comparison, nucleotide sequences for _B. thetaiotaomicron_ VPI-5482
(AE015928.1), _B. ovatus_ (CP012938.1), _Barnesiella viscericola_
(GCA000512915.1), and _Porphyromonas gingivalis_ (GCA000010505.1), were
also downloaded from GenBank.

## Genome annotation

All genomes were initially annotated with Prokka [@Seemann2014] version
1.13, which uses Prodigal [@Hyatt2010] for gene finding.
Putative
protein sequences were additionally annotated with domains from both the
dbCAN database [@Yin2012] release 6 of carbohydrate-active domains and
Pfam [@Punta2012] release 31.0, using HMMER3 [@Eddy2011; @Eddy2009]
version 3.1b2.
Protein sequences were also annotated with KO numbers by
BLAST using the KEGG database as of March 2018 as the reference and
taking the best hit with a maximum E-value of 1e-10.

Lipoproteins were predicted using LipoP [@Juncker2003] (version 1.0a)
and a score cutoff of 5 and a margin cutoff of 2.
Lipoproteins with an
arginine at position +2 relative to the cleavage site were labeled as
localized to the inner membrane.
Periplasmic proteins were identified
with SignalP [@Petersen2011] (version 4.1).
Predicted protein sequences
from all annotated genomes were locally all-by-all aligned using the
DIAMOND implementation of the BLAST algorithm [@Buchfink2014].
Each pair
was then assigned a similarity value as the bitscore of their best local
alignment normalized by the greater of the two self-alignments.
This
results in a matrix of pairwise scores reflecting the proximity to
perfect homology.
Scores less than 0.2 were replaced with 0.
Clusters
were formed using the MCL algorithm [@Enright2002] with an inflation
parameter of 5.

SusCDEF homologs were identified based on relatively relaxed criteria,
harnessing OPF assignments, domain predictions, and Prokka annotations
to avoid false negatives while maintaining specificity.
For each OPF,
all KOs assigned to members were collected as plausible KOs for the
cluster.
Protein sequences in OPF clusters which included K21572 were
flagged as putative SusC-homologs, as were sequences directly annotated
as such by Prokka.
Using a similar approach, proteins in clusters tagged
with K21571 or with any of domains PF12771, PF14322, PF12741, PF07980
were identified as putative SusD.
Proteins in clusters tagged with
K21571, or with either PF14292 or PF16411, were considered SusEF
homologs.
PULs were identified by a SusC-homolog with its start codon
within 5 kbp of a SusD-homolog's start on the same strand.
Manual
inspection supported the vast majority of these identifications.

## Phylogenetics

Predicted amino acid sequences for ORFs from all MAGs and all reference genomes
were search for homology to
TIGRFAM protein clusters [@Selengut2007] using `hmmsearch` in the HMMER3
software package version 3.2.1 [@Eddy2011].
Hits were filtered at the "trusted-cutoff" score threshold
defined separately for each protein model.
Sequences found in just one copy in every genome were used as taxonomic
marker genes for phylogenetic analysis.
Marker gene sequences were aligned to their respective models using `hmmalign`
dropping unaligned residues.
Aligned markers were then concatenated for each MAG and reference genome,
and an approximate maximum likelihood phylogeny was estimated using the
FastTree software version 2.1.10 [@Price2010] with the default parameters.

## Data Availability

Metagenomic libraries will be uploaded to the short read archive,
and the SRA Accession will be added in an upcoming version of this manuscript.
All code and additional metadata needed to reproduce this analysis are
available at <https://github.com/bsmith89/longev-mgen>.

# References

<!-- TODO: Move websites to links in the text rather than in the refs. -->
<!-- TODO: Fix capitalization of Sp. Nov. in references. -->
<!-- TODO: Fix capitalization of species names (due to italization) in references -->
<!-- TODO: Figure out appropriate use of journal title capitalization. -->
