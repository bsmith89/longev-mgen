-- {{{1 Denormalize Views
-- TODO: Do I need to denormalize more tables?  Fewer?
-- TODO: Do I need to add any indices?

CREATE TABLE __bin_coverage AS SELECT * FROM bin_coverage;
DROP VIEW bin_coverage;
ALTER TABLE __bin_coverage RENAME TO bin_coverage;

CREATE TABLE __contig_linkage AS SELECT * FROM contig_linkage;
DROP VIEW contig_linkage;
DROP TABLE _contig_linkage;
ALTER TABLE __contig_linkage RENAME TO contig_linkage;

CREATE TABLE __bin_linkage AS SELECT * FROM bin_linkage;
DROP VIEW bin_linkage;
ALTER TABLE __bin_linkage RENAME TO bin_linkage;
