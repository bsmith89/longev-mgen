.bail ON
.read schema.sql
PRAGMA cache_size = 1000000;
PRAGMA foreign_keys = TRUE;
.separator \t
.import meta/library.noheader.tsv library
.import res/core.a.proc.contigs.nlength.noheader.tsv contig
.import res/core.a.proc.contigs.bins.noheader.tsv contig_bin
.import res/core.a.proc.contigs.cvrg.noheader.tsv contig_coverage
.import res/core.a.proc.contigs.bins.checkm_details.noheader.tsv bin_checkm
ANALYZE;
