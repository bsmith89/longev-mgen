## Table Formats ##

```
== res/*.cvrg.tsv ==
library_id	contig_id	coverage	hit_fraction	left	right	contiguous

== res/*.depth.tsv ==
contig_id	position	depth

== res/*.contigs.tsv ==
contig_id	flag	multi	length
```

```
core-k161_147711	S. alaskensis
core-k161_23672	OTU-1
```

## Thoughts on what to present at lab meeting ##

Assembly statistics and quality check using meta-quast

Also include the Salask genome as a reference here

Binning discussiona and the approaches I ultimately tried

Single-seed, euclidean distance, manually select threshold
Seed set by choosing the most-mapped contigs and blasting to NR.
The very top one (a BIG contig) blasted to _Muribaculum intestinale_.

I chose a very complete genome without too much contamination (based on checkm using Bacteroidales db).
This was done at a cutoff of 5 (distance)

Prokka to annotate this bin.

Quast to look at SAlask and OTU1 bins, but then also MetaQuast to look
at my full assembly and compare it to Niel's.
