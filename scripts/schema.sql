CREATE TABLE coverage
  ( library_id TEXT
  , contig_id TEXT
  , coverage FLOAT

  , PRIMARY KEY (contig_id, library_id)
  );
CREATE INDEX idx_coverage__library_id ON coverage(library_id);

CREATE TABLE contig
  ( contig_id TEXT PRIMARY KEY
  , length INT
  );

CREATE TABLE bin
  ( contig_id TEXT PRIMARY KEY REFERENCES contig(contig_id)
  , bin_id TEXT
  );
CREATE INDEX idx_bin__bin_id ON bin(bin_id);

CREATE VIEW bin_length AS
  SELECT bin_id, SUM(length) AS length, COUNT(contig_id) AS n_contigs
  FROM contig
  JOIN bin USING (contig_id)
  GROUP BY bin_id
;

CREATE VIEW total_mapping AS
  SELECT contig_id, library_id, coverage * length AS mapping_count
  FROM coverage
  JOIN contig USING (contig_id)
;

CREATE VIEW bin_coverage AS
  SELECT
      bin_id
    , library_id
    , SUM(mapping_count) / length AS coverage
  FROM total_mapping
  JOIN bin USING (contig_id)
  JOIN bin_length USING (bin_id)
  GROUP BY bin_id, library_id
;
