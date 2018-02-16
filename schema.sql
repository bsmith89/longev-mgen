-- {{{1 Tables

CREATE TABLE library
  ( library_id TEXT PRIMARY KEY
  , extraction_id TEXT
  , run TEXT
  , file_r1 TEXT
  , file_r2 TEXT
  );

CREATE TABLE contig
  ( contig_id TEXT PRIMARY KEY
  , length INT
  );

CREATE TABLE contig_coverage
  ( library_id TEXT REFERENCES library(library_id)
  , contig_id TEXT REFERENCES contig(contig_id)
  , coverage FLOAT

  , PRIMARY KEY (contig_id, library_id)
  );
CREATE INDEX idx_contig_coverage__library_id ON contig_coverage(library_id);

CREATE TABLE bin_checkm
  ( bin_id TEXT PRIMARY KEY
  , markers_used INT
  , completeness FLOAT
  , contamination FLOAT
  , heterogeneity FLOAT
  );

CREATE TABLE contig_bin
  ( contig_id TEXT PRIMARY KEY
                   REFERENCES contig(contig_id)
  , bin_id TEXT
  );
CREATE INDEX idx_contig_bin__bin_id ON contig_bin(bin_id);

CREATE TABLE rrs_taxon_rabund
  ( extraction_id TEXT
  , taxon_id TEXT
  , relative_abundance FLOAT

  , PRIMARY KEY (extraction_id, taxon_id)
  );
CREATE INDEX idx_rrs_taxon_rabund__taxon_id ON rrs_taxon_rabund(taxon_id);

CREATE TABLE contig_linkage
  ( contig_id_1 TEXT
  , contig_id_2 TEXT
  , tally INT

  , PRIMARY KEY (contig_id_1, contig_id_2)
  );
CREATE TRIGGER trigger_contig_linkage_add_flipped_records
  AFTER INSERT ON contig_linkage FOR EACH ROW
  BEGIN
    INSERT INTO contig_linkage
      VALUES (NEW.contig_id_2, NEW.contig_id_1, NEW.tally)
    ;
  END
;
CREATE INDEX idx_contig_linkage__contig_id_2 ON contig_linkage(contig_id_2);


-- {{{1 Views

CREATE VIEW bin_length AS
  SELECT bin_id, SUM(length) AS length, COUNT(contig_id) AS n_contigs
  FROM contig
  JOIN contig_bin USING (contig_id)
  GROUP BY bin_id
;

CREATE VIEW total_nucleotides_mapping AS
  SELECT contig_id, library_id, coverage * length AS mapping_count
  FROM contig_coverage
  JOIN contig USING (contig_id)
;

CREATE VIEW bin_coverage AS
  SELECT
      bin_id
    , library_id
    , SUM(mapping_count) / length AS coverage
  FROM total_nucleotides_mapping
  JOIN contig_bin USING (contig_id)
  JOIN bin_length USING (bin_id)
  GROUP BY bin_id, library_id
;
