-- {{{1 Tables

-- general metadata {{{2
CREATE TABLE mouse
  ( mouse_id            TEXT PRIMARY KEY
  , cohort              TEXT
  , sex                 TEXT
  , treatment           TEXT       -- {acarbose, inulin, control, young (control),
                                   -- + cr (calorie-restricted), (17-alpha-)estradiol
                                   -- + rapamycin}
  , site                TEXT       -- {JL (Jackson Labs), UM (University of Michigan),
                                   -- + UT (University of Texas, San Antonio)}
  , cage_id                TEXT
  , date_of_birth       DATETIME
  , age_at_death_or_censor  INTEGER    -- in days
  , censored                BOOLEAN
  );

CREATE TABLE sample
  ( sample_id                   TEXT PRIMARY KEY
  , mouse_id             TEXT REFERENCES mouse(mouse_id)
  , collection_date      DATETIME
  , collection_time      DATETIME    -- Hour out of 24 local time
  , empty_tube_weight    FLOAT       -- in g
  , full_tube_weight     FLOAT       -- in g
  , hydrated_tube_weight FLOAT       -- in g of the sample after
                                     -- + homogenizing in water, spinning, and
                                     -- + pulling off the accessible supernatant
  , comments             TEXT
  );

CREATE TABLE extraction
  ( extraction_id                      TEXT PRIMARY KEY
  , sample_id               TEXT REFERENCES sample(sample_id)
  , weight                  FLOAT    -- of the sample extracted from in g
  , hydrated_weight         FLOAT    -- in g of the sample after
                                     -- + homogenizing in water, spinning, and
                                     -- + pulling off the accessible supernatant
  , volume                  FLOAT    -- in mL of water used to extract the sample
  , spike_id        TEXT             -- FIXME: This needs to be constrained to
                                     -- + values in spike(spike_id)
  , spike_volume    FLOAT            -- in uL
  , dna_concentration  FLOAT   -- ng/ul
  , comments                TEXT
  );

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

CREATE TABLE rrs_taxon_count
  ( extraction_id TEXT
  , sequence_id TEXT
  , otu_id TEXT
  , tally INT

  , PRIMARY KEY (extraction_id, sequence_id)
  );
CREATE INDEX idx_rrs_taxon_count__sequence_id ON rrs_taxon_count(sequence_id);

CREATE TABLE taxonomy
  ( sequence_id TEXT PRIMARY KEY
  , domain_ TEXT
  , phylum_ TEXT
  , class_ TEXT
  , order_ TEXT
  , family_ TEXT
  , genus_ TEXT
  );

CREATE TABLE _contig_linkage
  ( contig_id_1        TEXT REFERENCES contig(contig_id)
  , contig_id_2        TEXT REFERENCES contig(contig_id)
  , library_id         TEXT REFERENCES library(library_id)
  , read_count              INT

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
  SELECT contig_id_1 AS contig_id, contig_id_2 AS contig_id_linked, library_id, read_count
  FROM (SELECT * FROM _contig_linkage
        UNION
        SELECT contig_id_2 AS contig_id_1, contig_id_1 AS contig_id_2, library_id, read_count FROM _contig_linkage)
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
  , SUM(read_count) AS read_count
FROM (SELECT contig_id, contig_id_linked, SUM(read_count) AS read_count
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
    , COUNT(DISTINCT extraction_id) AS extraction_count
    , COUNT(DISTINCT library_id) AS library_count
    , SUM(read_count) AS read_count
  FROM contig_linkage
  JOIN library USING (library_id)
  GROUP BY contig_id, contig_id_linked
;
