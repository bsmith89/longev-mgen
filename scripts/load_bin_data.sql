.read scripts/schema.sql
PRAGMA cache_size = 1000000;
.separator \t
.import res/core.a.proc.scontigs.cvrg.noheader.tsv coverage
.import res/core.a.proc.scontigs.bins.noheader.tsv bin
.import res/core.a.proc.scontigs.nlength.noheader.tsv contig
ANALYZE;
