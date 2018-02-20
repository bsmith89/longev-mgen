-- {{{1 Tables

CREATE TABLE library
  ( library_id TEXT PRIMARY KEY
  , extraction_id TEXT
  , run TEXT
  , file_r1 TEXT
  , file_r2 TEXT
  );

CREATE TABLE library_asmbl_group
  ( library_id TEXT REFERENCES library
  , asmbl_group TEXT
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

CREATE TABLE _contig_linkage
  ( contig_id_1        TEXT REFERENCES contig(contig_id)
  , contig_id_2        TEXT REFERENCES contig(contig_id)
  , library_id         TEXT REFERENCES library(library_id)
  , tally              INT

  , PRIMARY KEY (contig_id_1, contig_id_2, library_id)
  );
CREATE INDEX idx_contig_linkage__contig_id_2 ON _contig_linkage(contig_id_2);
CREATE INDEX idx_contig_linkage__library_id ON _contig_linkage(library_id);

CREATE TABLE _bin_complementarity
  ( bin_id_1 TEXT
  , bin_id_2 TEXT
  , score FLOAT

  , PRIMARY KEY (bin_id_1, bin_id_2)
  );
CREATE INDEX idx__bin_complementarity__bin_id_2 ON _bin_complementarity(bin_id_2);

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

CREATE VIEW contig_linkage AS
  SELECT contig_id_1 AS contig_id, contig_id_2 AS contig_id_linked, library_id, tally
  FROM (SELECT * FROM _contig_linkage
        UNION
        SELECT contig_id_2 AS contig_id_1, contig_id_1 AS contig_id_2, library_id, tally FROM _contig_linkage)
;

CREATE VIEW bin_complementarity AS
  SELECT
      bin_id_1 AS bin_id
    , bin_id_2 AS bin_id_combined
    , score
  FROM (SELECT * FROM _bin_complementarity
        UNION
        SELECT bin_id_2 AS bin_id_1, bin_id_1 AS bin_id_2, score FROM _bin_complementarity)
;

CREATE VIEW bin_linkage AS
SELECT
    b1.bin_id AS bin_id_1
  , b2.bin_id AS bin_id_2
  , COUNT(DISTINCT b1.contig_id) AS contig_count_1
  , COUNT(DISTINCT b2.contig_id) AS contig_count_2
  , SUM(tally) AS tally
FROM (SELECT contig_id, contig_id_linked, SUM(tally) AS tally
      FROM contig_linkage
      GROUP BY contig_id, contig_id_linked
     )
JOIN contig_bin AS b1 USING (contig_id)
JOIN contig_bin AS b2 ON contig_id_linked = b2.contig_id
GROUP BY bin_id_1, bin_id_2
;

CREATE VIEW contig_total_coverage AS
  SELECT contig_id, SUM(coverage) AS coverage
  FROM contig_coverage
  GROUP BY contig_id
;

CREATE VIEW bin_total_coverage AS
  SELECT bin_id, SUM(coverage) AS coverage
  FROM bin_coverage
  GROUP BY bin_id
;

CREATE VIEW contig_total_linkage AS
  SELECT
      contig_id
    , contig_id_linked
    , COUNT(DISTINCT library_id) AS library_count
    , SUM(tally) AS tally
  FROM contig_linkage
  GROUP BY contig_id, contig_id_linked
;
