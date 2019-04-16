---
title: '_Muribaculaceae_ genomes assembled from metagenomes suggest
    genetic drivers of differential response to acarbose treatment in mice'

---

# Background

The mammalian gut microbiome is a complex ecological system that
influences energy balance [@Turnbaugh2006], pathogen
resistance [@Britton2012], and inflammation [@Syal2018], among other
processes with importance to host health. Understanding how the
bacterial inhabitants of the gut respond to pharmaceutical and dietary
perturbations is a major step in developing a predictive framework for
microbiome-based therapies. Acarbose (ACA) is an -glucosidase inhibitor
prescribed for the treatment of type 2 diabetes mellitus because it
reduces the absorption of glucose from starch in the small
intestine [@Hiele1992]. In rodents, ACA has been shown to increase the
amount of starch entering the lower digestive system after a
meal [@Dehghan-Kooshkghazi2004], resulting in changes to the composition
of the gut microbiota and its fermentation
products [@Zhao2018; @Holt1996; @Wolever2000; @Dehghan-Kooshkghazi2004; @Weaver1997; @Weaver2000; @Wolin1999; @Zhang2017].
Interestingly, long-term treatment with ACA has been shown to
substantially increase longevity in male mice and to a lesser extent in
females [@Harrison2014; @Strong2016].

In Chapter [\[ch:longev\]](#ch:longev){reference-type="ref"
reference="ch:longev"} it was shown that the relative abundance of a
number of bacterial taxa as well as the concentrations of propionate and
butyrate respond to long term treatment with ACA. This study was notable
in being replicated across three sites: The University of Michigan (UM)
in Ann Arbor, The University of Texas Health Science Center at San
Antonio (UT), and The Jackson Laboratory (TJL) in Bar Harbor, Maine. At
UM and TJL one highly abundant operational taxonomic unit (OTU),
classified as a member of the *Bacteroidales* family *Muribaculaceae*
and here designated as B1, was found to be enriched nearly 4-fold in ACA
treated mice. B1 was also present and abundant at UT but was not found
to be significantly more abundant in ACA treated mice relative to
controls. Instead, a different member of the *Muribaculaceae*,
designated B2, was found to be highly abundant and 4-fold enriched in
ACA-treated mice, but was nearly absent at UM and TJL. Other
*Muribaculaceae* were also identified as among the most abundant members
of the mouse gut microbiota across the three sites, although none of
these were found to be enriched in ACA treatment.

Family *Muribaculaceae*---formerly the S24-7 and sometimes referred to
as *Candidatus Homeothermaceae*---has only one published
cultivar [@Lagkouvardos2016] despite being a common and abundant
inhabitant of the mammalian gut, especially in mice [@Ormerod2016].
Previous studies have suggested that the *Muribaculaceae* specialize on
the fermentation of complex polysaccharides [@Ormerod2016], much like
members of the genus *Bacteroides* also in order *Bacteroidales*.

Recently, techniques have been developed for the reconstruction of
genomes of uncultivated members of bacterial
communities [@Parks2017; @Lee2017]. Based on 30 such metagenome
assembled genomes (MAGs) they reconstructed using this approach,
Ormerod *et al.* [@Ormerod2016] proposed that the *Muribaculaceae* fall
into three distinct carbohydrate utilization guilds, which they describe
as specialists on -glucans, plant glycans, and host glycans,
respectively. While it is reasonable to expect that -glucan specialists
would be most benefited by the large influx of starch to the gut
resulting from ACA treatment, this prediction has not been tested, and
physiological inferences based on the genome content of members of this
clade have been largely divorced from biological observations.

Experimental perturbations of complex microbial communities present an
opportunity to observe ecological features of many bacterial taxa
without cultivated members and generate hypotheses about their
physiology. Given the observed, dramatically increased relative
abundance of B1 and B2 (here referred to as "responders") in mice
treated with ACA, we hypothesize that these OTUs are capable of robust
growth on starch, while the other *Muribaculaceae* found in the study
("non-responders"), lack the genomic features necessary for the
utilization of the polysaccharide. Alternatively, responders may be
resistant to the inhibitory effects of ACA, or benefit from elevated
levels of intermediate starch degradation products. Since isolates of
the *Muribaculaceae* species in these mice are not available for
characterization, a comparative genomic approach is taken to explore
their functional potential.

Most of the research on the genomic components of polysaccharide
degradation in gram negative bacteria has been carried out in the genus
*Bacteroides*, and in particular *B. thetaiotaomicron* [@Martens2009].
Starch utilization in *B. thetaiotaomicron* is dependent on an ensemble
of eight proteins, SusRABCDEFG that enable recognition, binding,
hydrolysis, and import of starch and related
polysaccharides [@Foley2016]. Homologs of SusC and SusD characterize all
known polysaccharide utilization systems in this clade [@Grondin2017],
are encoded in Sus-like genomic regions known as polysaccharide
utilization loci (PULs), and are widespread in the
*Bacteroidetes* [@Fernandez-Gomez2013]. The molecular range of these
systems is determined by the carbohydrate-active enzymes and structural
proteins they encode, based on the specificity of glycoside hydrolase
(GH) and carbohydrate binding module (CBM) domains, which have been
extensively cataloged in the dbCAN database [@Yin2012; @Zhang2018].

Here MAGs from the feces of mice at UT and UM are analyzed to explore
two closely related questions about the niche of B1 and B2 in the lower
digestive system. First, why do B1 and B2 each increase with ACA
treatment, while other *Muribaculaceae* do not? And second, why is the
response of B1 site specific? Despite similar patterns of abundance at
their respective sites, these two OTUs seem to be only distantly
related, sharing just 90% of nucleotides in their 16S rRNA gene V4
hypervariable region (see
Appendix [\[app:otu\_tax\]](#app:otu_tax){reference-type="ref"
reference="app:otu_tax"}). We nonetheless find genomic evidence that B1
and B2 occupy overlapping niches, specializing in the degradation of
-glucans, a role not held by the other *Muribaculaceae* described in
this study. In addition, we identify two distinct variants of B1,
referred to as B1-A and B1-B, which are differentially distributed
between UM and UT and have functionally relevant differences in gene
content.

Reconstructing genomes from metagenomes allow for the comparison of the
functional potential of *Muribaculaceae* at UM and UT. This work
demonstrates the utility of culture-free genomics to understand the
ecological role of these key members of the mouse gut microbial
community and explore several hypotheses that may explain differences in
the distribution and response of bacteria to perturbations. Hypotheses
derived from this analysis provide a foundation for future physiological
studies in recently obtained cultivars. While a preponderance of
host-associated bacterial species have never been isolated, let alone
characterized [@Stewart2012], combining experimental data from complex
communities with the analysis of reconstructed genomes provides a
powerful tool for expanding understanding to these understudied taxa.

# Results

## Recovered population genomes are of high quality and resemble other *Muribaculaceae* genomes

MAGs were constructed for 7 populations in the family *Muribaculaceae*,
including ACA responders B1 and B2, and non-responders B3 through B7.
For B1, two genomic variants were recovered, B1-A and B1-B, MAGs that
possess 0.63 and 0.36 Mbp of unshared sequence, respectively (additional
details about these variants are in
Section [2.3](#subsec:mags-b1-var){reference-type="ref"
reference="subsec:mags-b1-var"}). All 8 novel MAGs are estimated to be
of high completeness and all had less than 1% estimated contamination
based on the recovery of ubiquitous, single-copy genes. The median N50
statistic was approximately 71 kbp, indicating successful assembly, and
suggesting that inferences based on genomic context are generally
possible. Estimated genome sizes, GC%, and number of predicted genes are
all similar to previously published MAGs as well as the finished
*Muribaculum intestinale* YL27 genome.

+---------+---------+---------+---------+---------+---------+---------+
| Taxon   | Complet | Scaffol | Length^ | N50     | GC      | in      |
|         | eness^1 | ds      | 2^      |         |         | Chapter |
|         | ^       |         |         |         |         |  [\[ch: |
|         |         |         |         |         |         | longev\ |
|         |         |         |         |         |         | ]](#ch: |
|         |         |         |         |         |         | longev) |
|         |         |         |         |         |         | {refere |
|         |         |         |         |         |         | nce-typ |
|         |         |         |         |         |         | e="ref" |
|         |         |         |         |         |         | referen |
|         |         |         |         |         |         | ce="ch: |
|         |         |         |         |         |         | longev" |
|         |         |         |         |         |         | }       |
+:========+:========+========:+:========+========:+:========+:========+
| YL-27^3^ | 99%     | 1       | 3.3     | 3,307,0 | 50.1%   |         |
+---------+---------+---------+---------+---------+---------+---------+
| B1-A    | 97%     | 228     | 3.2     | 41,412  | 46.6%   | OTU-1   |
+---------+---------+---------+---------+---------+---------+---------+
| B1-B    | 97%     | 152     | 3.0     | 59,916  | 46.9%   | OTU-1   |
+---------+---------+---------+---------+---------+---------+---------+
| B2      | 98%     | 65      | 2.6     | 79,454  | 50.5%   | OTU-4   |
+---------+---------+---------+---------+---------+---------+---------+
| B3      | 86%     | 98      | 2.2     | 63,818  | 54.0%   | OTU-6   |
+---------+---------+---------+---------+---------+---------+---------+
| B4      | 98%     | 31      | 2.7     | 148,039 | 55.2%   | OTU-5   |
+---------+---------+---------+---------+---------+---------+---------+
| B5      | 86%     | 50      | 2.5     | 78,179  | 55.7%   | OTU-8   |
+---------+---------+---------+---------+---------+---------+---------+
| B6      | 99%     | 110     | 3.2     | 87,115  | 48.3%   | OTU-30  |
+---------+---------+---------+---------+---------+---------+---------+
| B7      | 98%     | 97      | 2.5     | 59,037  | 53.9%   | OTU-39  |
+---------+---------+---------+---------+---------+---------+---------+

^7^ _Muribaculum intestinale_ YL-27 reference genome 

: [\[tbl:mag\_qual\_summary\]]{#tbl:mag_qual_summary
label="tbl:mag_qual_summary"} Summary of novel MAGs compared to the
genome of *Muribaculum intestinale* YL27

\vspace{2mm}
In order to confirm the assertion that each of the reconstructed genomes
is representative of *Muribaculaceae* OTU described in
Chapter [\[ch:longev\]](#ch:longev){reference-type="ref"
reference="ch:longev"}, per library mapping rates of each genome were
compared to the relative abundance of the associated 16S rRNA gene in
amplicon libraries. Despite the biases and technical variability
inherent to both sequencing methods, and the limitations of mapping
software, Pearson correlation coefficients between the fraction of reads
mapped and OTU relative abundance were above 0.86 for all MAGs,

![image](fig/muri_comparison.pdf){width="\textwidth"}

### Phylogenetics

To better understand the evolutionary relationships between these
organisms, a concatenated gene tree was constructed for all 8 novel
MAGs, as well as 30 publicly available MAG sequences [@Ormerod2016], and
*M. intestinale* YL27. The tree was rooted by four other *Bacteroidales*
species: *Bacteroides ovatus* (ATCC-8483), *Bacteroides
thetaiotaomicron* VPI-5482, *Porphyromonas gingivalis* (ATCC-33277), and
*Barnesiella viscericola* (DSM-18177). Most internal nodes were found to
have high topological confidence, and the placement of the MAGs
reconstructed by Ormerod *et al.* was highly consistent with their
published tree. To check that this concatenated approach is reflective
of the organismal evolutionary history, a second maximum likelihood tree
was constructed based on the *rpoB* gene, which is generally not thought
to be transmitted horizontally, (despite exceptions [@Kim2013]), also
recapitulating the published topology. The estimated phylogeny shows
that the 8 OTUs with newly reconstructed MAGs encompass most of the
documented diversity of *Muribaculaceae*. Two of our taxa, B2 and B6,
appear to be closely related to taxa with genomes reconstructed by
Ormerod *et al.*: M6, and M1, respectively. Nonetheless, this
phylogenetic analysis suggests that many of the genomes reconstructed
here have not been described previously.

### Novel protein families

Annotations based on alignment to a database of previously characterized
sequences may provide only limited insight, in particular for genomes
from largely unstudied families of bacteria. In order to identify
previously uncharacterized orthologous groups, *de novo*
clustering [@Schloss2008] was carried out based on amino acid similarity
of all putative genes found in the 8 novel MAGs, 30 previously
reconstructed MAGs, *M. intestinale*, four publicly available draft
genomes from the family, and the four reference *Bacteroidales*. The
resulting clusters are referred to as operational protein families
(OPFs). While a fraction of the 12,648 resulting OPFs may be due to
spurious sequence similarity and without biological relevance, 5,767 had
representatives in at least three genomes, increasing the likelihood
that these reflect evolutionarily conserved protein sequences. Of these,
only 2,404 had members annotated with any COG, KO, or putative function.
The remaining 3,363 OPFs include 17,831 predicted proteins across the 47
genomes

### Annotation ordination

To compare novel MAGs to other available genomes, a previous published
analysis was recreated, harnessing a set of 8 COGs found by Ormerod *et
al.* to maximally differentiate the three hypothesized guilds. By
projecting genome annotations onto a reproduction of this previously
defined space (see
Figure [\[fig:muri\_comparison\]](#fig:muri_comparison){reference-type="ref"
reference="fig:muri_comparison"}), newly available genomes were compared
to the three clusters hypothesized to represent specialization on
-glucans, plant glycans, and host glycans. While the 8 novel MAGs
inhabit approximately the same volume as those previously reconstructed,
and some could be plausibly classified based on these criteria, the
ambiguous placement of B4 and *M. intestinale* suggests that new genomes
will present additional exceptions to the three-guild model.

It is notably that both responders cluster with the proposed -glucan
guild, consistent with a functional potential for starch utilization not
present in the non-responders. To expand on this descriptive analysis
and to leverage the more comprehensive view provided by *de novo*
clustering to explore differences and similarities in carbohydrate
utilization potential, a second ordination of genomes was performed,
this time based on OPF labels of predicted genes found to contain GH
domains
(Figure [\[fig:pul\_diagrams\]](#fig:pul_diagrams){reference-type="ref"
reference="fig:pul_diagrams"}). Similar to the previous ordination based
on COGs, three groups of genomes approximately reflecting those proposed
by Ormerod *et al.* are apparent. However, the placement of B2 (as well
as the closely related M6) relative to the proposed guilds are
substantially different.

## Comparison of responder and non-responder MAGs suggest genomic features with role in starch utilization

Based on the characterization of genes and genomic regions with a role
in starch utilization in the closely related genus *Bacteroides*, it is
plausible that -amylase localized to the outer membrane may be common to
starch utilizing bacteria in the order *Bacteroidales* [@Shipman1999].
Indeed, B1 has three OM-localized genes predicted to code for GH13
containing lipoproteins (B1A280, B1A301, B1A333), each in a separate PUL
(see
Figure [\[fig:pul\_diagrams\]](#fig:pul_diagrams){reference-type="ref"
reference="fig:pul_diagrams"}). While it also includes members without
this activity, GH13 is the main family of -amylases [@Janecek2014].
These genomic regions also possess additional genes with
carbohydrate-active domains that are expected to interact with -glucans.

![image](fig/pul_diagrams.pdf){width="\textwidth"}

Besides B1, B5 is the only other OTU to possess a putative PUL coding
for a full complement of predicted starch-active proteins. Several OPFs
have members in both this region and either B1 or *B. thetaiotaomicron*
PULs, suggesting shared function. This set including SusC-homologs
Opf01277, Opf02066, which includes relatives of SusD, and Opf02791 whose
members possess CBM20 starch-binding domains. However, while B5 also has
a GH13 containing lipoprotein (B51713), its predicted localization is on
the inner membrane. It is unclear whether this explains B5's
non-response in ACA-treated mice. Plausible OM-localized, GH13
containing proteins are not found in any non-responders. While this
characteristic does not seem to perfectly discriminate responder from
non-responder OTUs---B2 also lacks such a gene---it nonetheless
demonstrates concordance between inferred genomic features and observed
population dynamics.

Despite the absence of a GH13 domain on the outer-membrane, it is
plausible that B2 is capable of degrading starch using other enzymatic
machinery. We speculate about one putative locus (see
Figure [\[fig:pul\_diagrams\]](#fig:pul_diagrams){reference-type="ref"
reference="fig:pul_diagrams"} panel F), which has a similar gene content
to characterized [@Ravcheev2013; @Rogers2013; @VanBueren2015] dextran
PULs in *B. thetaiotaomicron* and *B. ovatus*.

To expand the search for relevant genetic features, *de novo* protein
clusters were filtered to those with members in the MAGs for both B1 and
B2. Of these OPFs, several stood out as particularly relevant. Opf01144
includes SusR, the regulator of transcription of the starch utilization
system in *B. thetaiotaomicron*, as well as its homolog in *B. ovatus*.
It is an apparent subcluster of the larger family defined by K21557, and
in many cases is encoded directly upstream of *susC* in putative PULs
which consider likely to have affinity for -glucans. In B1, two of the
three putative starch PULs encode a member of Opf01144, and it is
similarly located in PULs with starch-active CBM and GH domains in B2
and B5. In addition, of the seven MAGs reconstructed by Ormerod *et al.*
that encode a member of this cluster, five of them are classified to the
-glucan guild. It is plausible that members of Opf01144 share a
functional role regulating transcriptional responses to -glucans.

Opf01391, which recapitulates K21575, includes SusA: the periplasmic
neopullulanase of *B. thetaiotaomicron* and an important component of
starch utilization in that organism [@DElia1996]. This family is found
in the MAGs of both responders, B1 and B2, and none of the
non-responders. What's more, it's found in twelve of the thirteen
-glucan and a minority of the plant glycan guild members. Interestingly,
although it is encoded by the Sus operon in *B. thetaiotaomicron* and
its homologous locus in *B. ovatus*, in the *Muribaculaceae* members of
Opf01391 do not in general appear to be encoded in PULs.

## Genomic variation in B1 {#subsec:mags-b1-var}

Two distinct variants of B1 were identified with one found in a majority
of the UT mouse metagenomes, and the other ubiquitous at UM. Using the
nucmer tool for genome alignment [@Delcher2002], 19.6% of B1-A MAG
sequence and 12.2% of B1-B were found to not align to the other. While
these hundreds of kbp may in part reflect errors in genome recovery,
much of the unaligned length suggests differences in gene content
between distinct sub-populations of B1. This observation was confirmed
by assessing the mapping of metagenomic reads against predicted protein
coding genes in each variant. For each pairing of metagenomic read
library to genomic variant, gene coverage was normalized by the median
gene coverage in order to identify genes with conspicuously fewer reads
in particular subsets of the mice. Libraries have low coverage of large
portions of either the B1-A or B1-B MAG (see
Figure [\[fig:b1\_vars\]](#fig:b1_vars){reference-type="ref"
reference="fig:b1_vars"}), suggesting that mice are primarily inhabited
by one of the two variants, and that a portion of genes are variant
specific.

![[\[fig:b1\_vars\]]{#fig:b1_vars label="fig:b1_vars"} Visualization of
differential gene content in two B1 populations. Heatmaps depict mapping
coverage of metagenomes against putative protein coding genes in the
B1-A or B1-B MAG normalized to the median coverage. Rows represent one
or more pooled libraries for each mouse included in the study and
columns represent individual genes. The site at which each mouse was
housed is indicated by triangles in the far left column: UT (green, left
pointing) or UM (blue, right). Filled triangles correspond to those mice
flagged as representative of a single B1 variant for downstream
analysis. Genes are shown only where the median normalized coverage
ratio between these B1-A and B1-B specific metagenomes is greater than
1.5. Rows and columns are arbitrarily ordered to maximize visual
distinction between variants. ](fig/b1_vars.pdf){width="\textwidth"}

Metagenomic libraries manually chosen as unambiguous representatives of
a single B1 MAG were used to systematically identify genes
differentiating the two. The median normalized mapping depths in each
set of libraries against predicted genes in each MAG were compared,
providing a measure of the relative enrichment or depletion of genomic
sequences between the two populations of B1. This analysis found 12.8%
of predicted genes in B1-A were depleted at least 5-fold in B1-B
populations, and 12.4% the reverse. While this observed depletion could
indicate variation in copy number, differential gene content between
variants is a more parsimonious explanation for most loci. These
predicted genes reflect 2.7% of unique KOs in B1-A and 1.9% in B1-B.
Interestingly, the fraction of variant specific OPFs is greater, 7.5%
and 7.1% respectively, suggesting that *de novo* clustering could be
more sensitive to potential differences in physiology.

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

^2^ unique

: [\[tbl:b1\_vars\]]{#tbl:b1_vars label="tbl:b1_vars"}Summary of variant
specific features in two B1 MAGs

\vspace{2mm}
Given the observation that the relative abundance of B1 was dramatically
increased with ACA treatment at UM, while not being significantly
affected at UT, and that B1-B was not found in metagenomes at UM, we
searched for differences in functional potential between the two
variants that could explain this pattern.

Genomic regions apparently specific to B1-A---defined as an at least
5-fold enrichment---include just one PUL (SusC-homolog encoded by
B1A00048). This locus includes a predicted outer membrane localized GH30
containing protein. Characterized GH30 containing proteins have
-glucosylceramidase, -1,6-glucanase, or -xylosidase
activity [@StJohn2010]. Given that this PUL also encodes a periplasmic,
GH3 containing protein, it appears to be unlikely that it has
specificity for starch. The B1-A MAG also possesses numerous phage
insertions not seen in the B1-B reconstruction. Conversely, a CRISPR
operon including 25 repeat units (Cas9 encoded by B1B01367) appears to
be specific to B1-B.

Most strikingly, a 16 kbp region (from B1A01498 to B1A01514) specific to
B1-A was found to contain several genes with homology to cell capsule
and exopolysaccharide synthesizing enzymes. Based on annotations with
KEGG orthologous groups, these include homologs of *tuaG* (K16698),
*tagE* (K00712), *gmhB* (K03273), *gmhA*/*lpcA* (K03271),
*hddA* (K07031), *exoO* (K16555), *waaH* (K19354), and *tagF* (K09809).
Interestingly, the B1-B MAG contains a different such region of about
6.5 kbp (B1B00851 to B1B00856) with *wfeD* (K21364), *pglJ* (K17248),
and *epsH* (K19425). For each, several of the OPFs in the respective
regions were not found anywhere in the opposing genome, suggesting that
the makeup of each variant's exterior surface might be distinctly
different.

# Discussion

Mice are a key model system for study of the mammalian gut microbiome,
with an outsized importance in testing mechanistic hypotheses for the
role of this community on host health [@Nguyen2015]. The
generalizability of observations made in mice is a constant
concern [@Nguyen2015], in part due to extensive difference in taxonomic
composition compared to humans [@Lagkouvardos2016]. The members of the
*Muribaculaceae* are abundant in the murine gut
microbiome [@Ormerod2016]. While these bacteria are also found in humans
(although at lower abundance), only one cultivated member of this clade
has been described [@Lagkouvardos2016]. As a result, the ecological
roles of these taxa have not been characterized, and observations in
mouse model systems are therefore less valuable for understanding
related processes in the human gut microbiome. Attempts to study these
organisms leverage genomes reconstructed from metagenomic reads, and
have proposed---in the absence of experimental data---that members of
the family consume a diversity of polysaccharides in the lower gut.

Here we have extended that approach to eight new genomes, and associated
those with taxa for which changes in relative abundance in response to
ACA treatment have been experimentally assessed. This enabled us to
explore why responders B1 and B2 each increase with ACA treatment, while
the other *Muribaculaceae* do not. Annotations of reconstructed genomes
suggest that these may possess starch degradation capabilities absent in
the non-responders.

We examine the three-guild model proposed by Ormerod *et
al.* [@Ormerod2016] by reproducing their dimensional reduction approach
with the addition of these new genomes. In this analysis, B1 and B2
annotations appear to be consistent with a hypothesized -glucan
degradation guild, supporting their interpretation. A more nuanced
approach to annotation was also applied by constructing *de novo*
clusters of proteins based on homology. Interestingly, this analysis
indicates that B2, and the closely related M6, share physiological
potential with taxa in the host-glycan guild, suggesting that a more
detailed examination can identify specific functions that discriminate
responders from non-responders. This approach is bolstered by the
phylogenetic and genomic distinction between B1 and B2, reducing the
confounding effects of shared evolutionary history.

By including otherwise unannotated genes, genomic comparisons based on
OPFs instead of previously defined gene orthologies may better reflect
shared functional potential. Besides the identification of potentially
novel gene families, *de novo* homology clustering [@Schloss2008] also
enables differentiation of sub-groups not captured by standard
annotations. For instance, hypothetical genes annotated as homologs of
SusC, SusD, and SusEF, were clustered into 119, 162, and 33 different
OPFs respectively. It is plausible that this sub-clustering captures
differences in protein structure with importance in oligo- and
polysaccharide recognition, import, and binding. Combined with
annotation of characterized functional domains, these clusters may
better predict the polysaccharide utilization ranges of uncultured
organisms.

A detailed analysis of PULs identified multiple loci in B1 that appear
to be adapted to the degradation of starch or related carbohydrates, due
to the presence of an OM localized GH13 containing
protein [@Koropatkin2010]. Counterintuitively, B2 had no such PUL,
suggesting that its response to ACA may result from other enzymatic
capabilities. Of particular interest is a PUL encoding proteins with
GH97, CBM20, and CBM69 domains, all of which have documented activity on
starch [@Naumoff2005; @Boraston2004]. While the only outer-membrane
localized hydrolase in this PUL is a GH66, and members of this family
have characterized activity on the -1,6 linkages between glucose
monomers in dextran [@Kim2012]. It is plausible that this PUL can be
repurposed and confers some ability to grow on starch.

In addition, a gene encoding a SusA homolog was identified in both B1
and B2 but in none of the non-responders. While it is unclear how
expression of this important component of starch utilization might be
regulated, given that it is not located in a PUL in either of the
responders, SusA is important for growth on amylopectin in *B.
thetaiotaomicron* [@DElia1996]. Since inhibition by acarbose is variable
across enzymes [@Kim1999], it is possible that acarbose treatment
results in elevated levels of dextrin and maltooligosaccharides in the
lower guts of mice due to residual -amylase activity, even at levels
sufficient to prohibit host digestion. Periplasmic hydrolysis of these
starch breakdown products may be sufficient for increased abundance of
these taxa in acarbose treated mice.

It is notable that two distinct variants of B1 were identifiable in
these metagenomes, and that the distribution of B1-A and B1-B are
reminiscent of the previously observed site-specificity of ACA response.
Despite evidence that genomic variation is common in the bacterial
world [@Rasko2008; @Medini2005], studies reconstructing genomes from
metagenomes often ignore this possibility (with a few notably
exceptions [@Truong2017a; @Delmont2018]). The discovery of two
subpopulations of B1 therefore demonstrates the value of considering
pangenome dynamics, and presents a potential explanation for the
observed site-specific response of that taxon. The finding that both
variants have the same complement of three PULs apparently specializing
in starch utilization and the same SusA homolog does not support the
hypothesis that differences in starch utilization potential account for
these abundance patterns. We did, however, identify numerous differences
in the gene content of B1-A and B1-B, including variant specific loci
that may influence the structure and function of the outer surface of
the cell. Capsule variation is known to greatly affect both ecological
and host interactions [@Merino2015].

While these results do not establish a mechanistic explanation for
differences in the response of B1 at UM and UT, nor conclusively
identify starch utilization pathways in B2, they do suggest a number of
genomic features that likely contribute to previously observed patterns
in taxon abundance. Future studies utilizing metatranscriptomic analysis
might demonstrate active expression of these genes, or differential
expression in mice treated with acarbose compared to controls. Likewise,
even in the absence of a B2 cultivar, the sufficiency of the dextran PUL
for increased growth with acarbose treatment could be tested using
available cultivars, including *B. thetaiotaomicron*.

# Conclusions

In this study we have reconstructed and described genomes representing 7
OTUs in the family *Muribaculaceae* from the mouse fecal microbiome, and
have found features that differentiate those that respond positively to
ACA treatment from those that do not. This analysis suggests that
utilization of starch and related polysaccharides enables increased
population size in mice treated with the -amylase inhibitor. In
addition, two distinct genomic variants of one taxon were identified
that differ in functional gene content, potentially explaining
site-specific differences in response. By combining observed changes in
relative abundance during experimental manipulation with inferred
functional gene content, we are able to study mammalian symbionts in the
absence of cultured representatives. This sequence-based approach is
broadly applicable in microbial ecology and enables improved
understanding of *in situ* dynamics within complex microbial
communities.

# Methods

## Mouse treatment, sample collection, extraction and sequencing

Mice were bred, housed, and treated as described in [@Harrison2014].
Briefly, genetically heterogeneous UM-HET3 mice at each study site were
produced by the four-way cross detailed in [@Miller2011]. Mice were fed
LabDiet (TestDiet Inc.) 5LG6 from weaning onwards. Starting at 8 months
of age, mice randomly assigned to treatment were fed chow with 1,000 ppm
ACA (Spectrum Chemical Manufacturing Corporation). Mice were housed 4
males or 5 females to a cage. Colonies were assessed for infectious
agents every 3 months, and all tests were negative.

Individual fecal pellets were collected from a single mouse per cage.
16S rRNA gene libraries and metabolite analyses of these samples are
described in Chapter [\[ch:longev\]](#ch:longev){reference-type="ref"
reference="ch:longev"}. From this collection, a subset of samples were
non-randomly selected for metagenomic sequencing based on various
criteria. Samples were from 54 mice, with at least six treated and
control representatives of both males and females at each site.

Fecal samples were slurried with nuclease free water at a 1:10 (w/v)
ratio, and most samples were spiked with *Sphingopyxis alaskensis*
RB2256 prepared as described in
Chapter [\[ch:longev\]](#ch:longev){reference-type="ref"
reference="ch:longev"} before DNA extraction and sequencing. Based on
alignment to the reference genome, sequenced reads from *S. alaskensis*
can be distinguished from all endogenous bacteria in mouse feces. A
small number of these were split for both spiked and unspiked samples,
which we used to validate this procedure. For each, 150 $\mu$Lof this
sample was transferred for extraction using the MoBio PowerMag
Microbiome kit. Metagenomic libraries were prepared using standard
procedures sequenced on the Illumina HiSeq 400 platform using the v4
paired-end 2x150 bp.

## Assembly, binning, and MAG refinement

Raw metagenomic reads were deduplicated using FastUniq [@Xu2012],
adapters trimmed using Scythe [@Buffalo2018], and quality trimmed using
Sickle [@Joshi2011] to produce processed reads for all downstream
analyses. The resulting paired-end reads were assembled into primary
contigs using MEGAHIT [@Li2014]. Reads were then mapped back to these
contigs with Bowtie2 [@Langmead2012], and per-library coverage was
estimated for each contig.

For all contigs \>1000 bp in length, dimensional reductions built into
CONCOCT [@Alneberg2014] were applied to produce input data for a
Gaussian mixture model (GMM) similar to the procedure used by that
program for binning. However, unlike CONCOCT---due to computational
limitations---the model was trained on just 10% of the input data,
sampled randomly, before assigning bins to all contig. While this may
have reduced the accuracy of the binning procedure, we believe that
subsequent refinement steps mitigated the impact of this decision.

OTUs were classified taxonomically and relative abundance was calculated
for matched libraries as described in
Chapter [\[ch:longev\]](#ch:longev){reference-type="ref"
reference="ch:longev"}. Bins were then recruited to one or more OTUs by
calculating a Canonical partial least squares between OTU abundance and
bin coverage as implemented in the scikit-learn machine learning library
for Python [@Pedregosa2012]. For bins recruited to OTUs classified as
*Muribaculaceae*, contigs were re-clustered based on coverage across
samples. First "trusted contigs" were manually selected which correlated
closely with OTU abundance. The mean coverage of these was used to
normalize the per-library coverage of all other contigs. Then, using a
GMM, groups of contigs were clustered such that the normalized coverage
across samples was consistent. These groups were used to inform the
manual assignment of contigs to MAGs. Libraries in which MAGs had
non-negligible coverage were identified and used in subsequent
refinements. For the B1 reconstruction, but no other MAGs, a number of
groups containing on the order of $10^5$ bp were found with low coverage
in just a subset of libraries. By this criterion, contigs in these
"variable" groups were partitioned into two MAG variants, A and B, with
non-variable groups shared by both. Only libraries that appeared on
further inspection to have just one of the two variants were considered
in downstream refinement steps. The mice matching these libraries are
highlighted in
Figure [\[fig:b1\_vars\]](#fig:b1_vars){reference-type="ref"
reference="fig:b1_vars"}.

For each MAG, several alternative refinement procedures were performed
from which the best quality result was selected. Reads mapping to the
curated contigs were digitally
normalized [@Wedemeyer2017; @Brown2012a; @Zhang2014a] and reassembled
with SPADES [@Bankevich2012]. This reassembly as well as the original
contigs were cleaned using a single pass of the Pilon assembly
refinement tool [@Walker2014a]. Finally, the per-library mapping depths
of each position in these assemblies were compared to the average
mapping depth of the "trusted contigs" selected earlier, and regions
with low cosine similarity were excised from the final assemblies.

Genome completeness and contamination estimates were calculated based on
ubiquitous single-copy genes using the program CheckM [@Parks2015].
Based on these results, the final assembly with the highest completeness
and with contamination \< 1% was selected from the various refinements.

## Reference genomes

The *Muribaculum intestinale* genome sequence was obtained from GenBank
(accession GCA002201515.1), as well as four additional draft genomes
(GCA003024805.1, GCA003024815.1, GCA002633305.1, GCA002633115.1). While
other genomes labeled as *Muribaculaceae* have also been deposited, they
were excluded from this analysis due to redundancy or apparent
misidentification to the family. The 30 MAGs reconstructed by
Ormerod *et al.* [@Ormerod2016] were obtained from the SRA. For
comparison, nucleotide sequences for *B. thetaiotaomicron* VPI-5482
(AE015928.1), *B. ovatus* (CP012938.1), *Barnesiella viscericola*
(GCA000512915.1), and *Porphyromonas gingivalis* (GCA000010505.1), were
also downloaded from GenBank.

## Genome annotation

All genomes were initially annotated with Prokka [@Seemann2014] version
1.13, which uses Prodigal [@Hyatt2010] for gene finding. Putative
protein sequences were additionally annotated with domains from both the
dbCAN database [@Yin2012] release 6 of carbohydrate-active domains and
Pfam [@Punta2012] release 31.0, using HMMER3 [@Eddy2011; @Eddy2009]
version 3.1b2. Protein sequences were also annotated with KO numbers by
BLAST using the KEGG database as of March 2018 as the reference and
taking the best hit with a maximum E-value of 1e-10.

Lipoproteins were predicted using LipoP [@Juncker2003] (version 1.0a)
and a score cutoff of 5 and a margin cutoff of 2. Lipoproteins with an
arginine at position +2 relative to the cleavage site were labeled as
localized to the inner membrane. Periplasmic proteins were identified
with SignalP [@Petersen2011] (version 4.1). Predicted protein sequences
from all annotated genomes were locally all-by-all aligned using the
DIAMOND implementation of the BLAST algorithm [@Buchfink2014]. Each pair
was then assigned a similarity value as the bitscore of their best local
alignment normalized by the greater of the two self-alignments. This
results in a matrix of pairwise scores reflecting the proximity to
perfect homology. Scores less than 0.2 were replaced with 0. Clusters
were formed using the MCL algorithm [@Enright2002] with an inflation
parameter of 5.

SusCDEF homologs were identified based on relatively relaxed criteria,
harnessing OPF assignments, domain predictions, and Prokka annotations
to avoid false negatives while maintaining specificity. For each OPF,
all KOs assigned to members were collected as plausible KOs for the
cluster. Protein sequences in OPF clusters which included K21572 were
flagged as putative SusC-homologs, as were sequences directly annotated
as such by Prokka. Using a similar approach, proteins in clusters tagged
with K21571 or with any of domains PF12771, PF14322, PF12741, PF07980
were identified as putative SusD. Proteins in clusters tagged with
K21571, or with either PF14292 or PF16411, were considered SusEF
homologs. PULs were identified by a SusC-homolog with its start codon
within 5 kbp of a SusD-homolog's start on the same strand. Manual
inspection supported the vast majority of these identifications.

# References