-- {{{1 Tables

-- CREATE TABLE mag
-- ( mag_id PRIMARY KEY
-- , description
-- );

CREATE TABLE checkm
( mag_id PRIMARY KEY --REFERENCES mag(mag_id)
, n_markers INT
, completeness FLOAT
, contamination FLOAT
, strain_heterogeneity FLOAT
);

CREATE TABLE quast
( mag_id PRIMARY KEY REFERENCES checkm(mag_id)
, n_contigs INT
, n_contigs_gt_1000 INT
, n_contigs_gt_5000 INT
, n_contigs_gt_10000 INT
, n_contigs_gt_25000 INT
, n_contigs_gt_50000 INT
, total_length INT
, total_length_gt_1000 INT
, total_length_gt_5000 INT
, total_length_gt_10000 INT
, total_length_gt_25000 INT
, total_length_gt_50000 INT
, n_contigs_ INT
, length_largest_contig INT
, total_length_ INT
, gc_percent FLOAT
, n50 INT
, n75 INT
, l50 INT
, l75 INT
, ambiguous_bases_per_100kbp FLOAT
);

CREATE TABLE _sequence
( sequence_id PRIMARY KEY
, mag_id REFERENCES checkm(mag_id)
);

CREATE TABLE _sequence_length
( sequence_id PRIMARY KEY REFERENCES _sequence(sequence_id)
, nlength INT
);

CREATE TABLE ko
( ko_id PRIMARY KEY
, description
);

CREATE TABLE cog
( cog_id PRIMARY KEY
, function_category
, description
);

CREATE TABLE domain
( domain_id PRIMARY KEY
, pfam_description
, cazy_description
);

CREATE TABLE feature
( feature_id PRIMARY KEY
, sequence_id REFERENCES _sequence(sequence_id)
, left INT
, right INT
);
CREATE INDEX idx_feature__sequence_id ON feature(sequence_id);

CREATE TABLE _feature_details
( feature_id PRIMARY KEY REFERENCES feature(feature_id)
, ftype
, nlength INT
, product_description
);

CREATE TABLE feature_to_cog
( feature_id PRIMARY KEY REFERENCES feature(feature_id)
, cog_id REFERENCES cog(cog_id)
);

CREATE TABLE feature_to_ko
( feature_id REFERENCES feature(feature_id)
, ko_id REFERENCES ko(ko_id)
);

CREATE TABLE feature_to_opf
( feature_id PRIMARY KEY REFERENCES feature(feature_id)
, opf_id
);
CREATE INDEX idx_feature_to_opf__opf_id ON feature_to_opf(opf_id);

CREATE TABLE feature_domain
( feature_id REFERENCES feature(feature_id)
, domain_id
, score FLOAT
, left INT
, right INT
);
CREATE INDEX idx_feature_domain__domain_id ON feature_domain(domain_id);

CREATE TABLE feature_to_architecture
( feature_id PRIMARY KEY REFERENCES feature(feature_id)
, architecture
);

CREATE TABLE feature_signal_peptide
( feature_id PRIMARY KEY REFERENCES feature(feature_id)
, cleavage_position INT
, score FLOAT
, closest_cysteine INT
);

CREATE TABLE feature_tmh
( feature_id PRIMARY KEY REFERENCES feature(feature_id)
, tmhelix_count INT
);

CREATE TABLE feature_lipop
( feature_id PRIMARY KEY REFERENCES feature(feature_id)
, lipop_type
, score FLOAT
, margin FLOAT
, cleavage_position INT
, aa_at_position_plus_2
);

-- {{{1 Views

CREATE VIEW sequence AS
SELECT * FROM _sequence JOIN _sequence_length USING (sequence_id)
;

CREATE VIEW feature_neighborhood AS
SELECT
    a.feature_id AS seed_id
  , b.feature_id AS feature_id
  , ABS(((a.left + a.right) / 2) - ((b.left + b.right) / 2)) AS distance
 FROM feature AS a
 JOIN feature AS b
   ON a.sequence_id = b.sequence_id
;

CREATE VIEW opf_architecture AS
SELECT
    opf_id
  , architecture
  , COUNT(feature_id) AS tally
FROM feature_to_opf
LEFT JOIN feature_to_architecture USING (feature_id)
GROUP BY opf_id, architecture
ORDER BY opf_id, tally DESC
;

-- Match each OPF to it's most common architecture, with metadata.
CREATE VIEW opf_to_architecture AS
SELECT
    a.opf_id
  , a.architecture
  , a.tally * 1.0 / m.total_tally AS fraction
  , m.total_tally AS out_of
FROM opf_architecture AS a
JOIN (SELECT
          opf_id
        , architecture
        , MAX(tally) AS max_tally
        , SUM(tally) AS total_tally
      FROM opf_architecture
      GROUP BY opf_id
     ) AS m
  ON a.opf_id = m.opf_id AND a.tally = max_tally
GROUP BY a.opf_id
;

CREATE VIEW feature_localization AS
SELECT
    feature_id
  , CASE
    WHEN lp.score > 8 AND lp.aa_at_position_plus_2 != 'D' AND lp.margin > 4 THEN 'OM'
    WHEN lp.score > 5 AND lp.aa_at_position_plus_2 != 'D' AND lp.margin > 2 THEN 'OM?'
    WHEN lp.score > 5 AND lp.aa_at_position_plus_2 != 'D' THEN 'OM??'
    WHEN lp.score > 8 AND lp.aa_at_position_plus_2 = 'D' AND lp.margin > 4 THEN 'IM'
    WHEN lp.score > 5 AND lp.aa_at_position_plus_2 = 'D' AND lp.margin > 2 THEN 'IM?'
    WHEN lp.score > 5 AND lp.aa_at_position_plus_2 = 'D' THEN 'IM??'
    WHEN sp.score > 0.5 AND sp.closest_cysteine < 4 THEN 'PP?/OM?'
    WHEN sp.score > 0.7 THEN 'PP'
    WHEN sp.score > 0.5 THEN 'PP?'
    WHEN sp.score > 0.4 THEN 'PP??'
    ELSE ''
    END AS localization
  , lp.score AS lp_score
  , lp.margin AS lp_margin
  , sp.score AS sp_score
  , sp.closest_cysteine AS closest_cysteine
  , tmhelix_count
FROM feature
LEFT JOIN feature_lipop AS lp USING (feature_id)
LEFT JOIN feature_signal_peptide AS sp USING (feature_id)
LEFT JOIN feature_tmh USING (feature_id)
;

CREATE VIEW feature_details AS
SELECT *
FROM feature
LEFT JOIN _feature_details USING (feature_id)
LEFT JOIN feature_localization USING (feature_id)
LEFT JOIN feature_to_opf USING (feature_id)
LEFT JOIN feature_to_architecture USING (feature_id)
LEFT JOIN feature_to_cog USING (feature_id)
LEFT JOIN feature_to_ko USING (feature_id)
;
