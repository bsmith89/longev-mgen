-- {{{1 Tables

CREATE TABLE genome
( genome_id PRIMARY KEY
, is_mag BOOL
, genome_type
, taxon_family
, ormerod_guild
);

CREATE TABLE checkm
( genome_id PRIMARY KEY REFERENCES genome(genome_id)
, n_markers INT
, completeness FLOAT
, contamination FLOAT
, strain_heterogeneity FLOAT
);

CREATE TABLE quast
( genome_id PRIMARY KEY REFERENCES genome(genome_id)
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
, genome_id REFERENCES genome(genome_id)
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

CREATE TABLE pfam_domain
( domain_id PRIMARY KEY
, pfam_description
);

CREATE TABLE cazy_domain
( domain_id PRIMARY KEY
, cazy_description
);

CREATE TABLE tigr_domain
( domain_id PRIMARY KEY
, tigr_description
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

CREATE TABLE feature_x_ko
( feature_id REFERENCES feature(feature_id)
, ko_id REFERENCES ko(ko_id)
);

CREATE TABLE feature_to_opf
( feature_id PRIMARY KEY REFERENCES feature(feature_id)
, opf_id
);
CREATE INDEX idx_feature_to_opf__opf_id ON feature_to_opf(opf_id);

CREATE TABLE feature_x_pfam_domain
( feature_id REFERENCES feature(feature_id)
, domain_id REFERENCES pfam_domain(domain_id)
, score FLOAT
, left INT
, right INT
);
CREATE INDEX idx_feature_x_pfam_domain__domain_id ON feature_x_pfam_domain(domain_id);

CREATE TABLE feature_x_cazy_domain
( feature_id REFERENCES feature(feature_id)
, domain_id REFERENCES cazy_domain(domain_id)
, score FLOAT
, left INT
, right INT
);
CREATE INDEX idx_feature_x_cazy_domain__domain_id ON feature_x_cazy_domain(domain_id);

-- Like the above, but excluding less-ideal hits.
CREATE TABLE feature_x_cazy_minimal_domain
( feature_id REFERENCES feature(feature_id)
, domain_id REFERENCES cazy_domain(domain_id)
, score FLOAT
-- TODO: is_minimal BOOL
, left INT
, right INT
);
CREATE INDEX idx_feature_x_cazy_minimal_domain__domain_id ON feature_x_cazy_minimal_domain(domain_id);

CREATE TABLE feature_x_tigr_domain
( feature_id REFERENCES feature(feature_id)
, domain_id REFERENCES tigr_domain(domain_id)
, score FLOAT
, left INT
, right INT
);
CREATE INDEX idx_feature_x_tigr_domain__domain_id ON feature_x_tigr_domain(domain_id);

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

CREATE TABLE variant_cross_coverage
( feature_id
, genome_id
, coverage_ratio FLOAT
, PRIMARY KEY (feature_id, genome_id)
);

CREATE TABLE library_size
( library_id
, tally INT
);

CREATE TABLE feature_library_coverage
( feature_id REFERENCES feature(feature_id)
, library_id
, coverage INT
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

CREATE VIEW opf_ko AS
SELECT
    opf_id
  , ko_id
  , COUNT(feature_id) AS tally
FROM feature_to_opf
LEFT JOIN feature_x_ko USING (feature_id)
GROUP BY opf_id, ko_id
ORDER BY opf_id, tally DESC
;

-- Match each OPF to it's most common KO, with metadata.
CREATE VIEW opf_to_ko AS
SELECT
    a.opf_id
  , a.ko_id
  , a.tally * 1.0 / m.total_tally AS fraction
  , m.total_tally AS out_of
FROM opf_ko AS a
JOIN (SELECT
          opf_id
        , ko_id
        , MAX(tally) AS max_tally
        , SUM(tally) AS total_tally
      FROM opf_ko
      GROUP BY opf_id
     ) AS m
  ON a.opf_id = m.opf_id AND a.tally = max_tally
GROUP BY a.opf_id
;

CREATE VIEW opf_cog AS
SELECT
    opf_id
  , cog_id
  , COUNT(feature_id) AS tally
FROM feature_to_opf
LEFT JOIN feature_to_cog USING (feature_id)
GROUP BY opf_id, cog_id
ORDER BY opf_id, tally DESC
;

-- Match each OPF to it's most common COG, with metadata.
CREATE VIEW opf_to_cog AS
SELECT
    a.opf_id
  , a.cog_id
  , a.tally * 1.0 / m.total_tally AS fraction
  , m.total_tally AS out_of
FROM opf_cog AS a
JOIN (SELECT
          opf_id
        , cog_id
        , MAX(tally) AS max_tally
        , SUM(tally) AS total_tally
      FROM opf_cog
      GROUP BY opf_id
     ) AS m
  ON a.opf_id = m.opf_id AND a.tally = max_tally
GROUP BY a.opf_id
;

CREATE VIEW feature_localization AS
SELECT
    feature_id
  , CASE
    WHEN lp.score > 5 AND lp.aa_at_position_plus_2 != 'D' AND lp.margin > 3 THEN 'OM'
    WHEN lp.score > 5 AND lp.aa_at_position_plus_2 = 'D' AND lp.margin > 3 THEN 'IM'
    WHEN sp.score > 0.4 THEN 'PP'
    ELSE 'CY'
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
-- DO NOT join KOs, because feature_x_ko is many-to-many:
;

--  Many-to-many relationship between features and plausible KO membership.
CREATE VIEW feature_possible_ko AS
SELECT DISTINCT *
FROM (
  SELECT *
  --  All KOs assigned to the feature
  FROM feature_x_ko AS k1

  UNION

  SELECT
  --  All KOs assigned to any feature in the same OPF
      f1.feature_id AS feature_id
    , k2.ko_id AS ko_id
  FROM feature AS f1
  JOIN feature_to_opf AS o1 USING (feature_id)
  JOIN feature_to_opf AS f2 USING (opf_id)
  JOIN feature_x_ko AS k2
    ON k2.feature_id = f2.feature_id
     )
WHERE feature_id != ''
  AND ko_id != ''
;

CREATE VIEW starch_active_gh_hit AS
SELECT
    feature_id
  , domain_id
  , score
FROM feature_x_cazy_domain
WHERE
    ( domain_id LIKE 'GH13|_%' ESCAPE '|' OR domain_id = 'GH13'
   OR domain_id LIKE 'GH97|_%' ESCAPE '|' OR domain_id = 'GH97'
   OR domain_id LIKE 'GH31|_%' ESCAPE '|' OR domain_id = 'GH31'
   OR domain_id LIKE 'GH57|_%' ESCAPE '|' OR domain_id = 'GH57'
   OR domain_id LIKE 'GH70|_%' ESCAPE '|' OR domain_id = 'GH70'
   OR domain_id LIKE 'GH77|_%' ESCAPE '|' OR domain_id = 'GH77'
    )
;

CREATE VIEW starch_active_cbm_hit AS
SELECT
    feature_id
  , domain_id
  , score
FROM feature_x_cazy_domain
WHERE domain_id IN ('CBM20', 'CBM21', 'CBM25',
                    'CBM26', 'CBM41', 'CBM48',
                    'CBM53', 'CBM58', 'CBM74',
                    'CBM82', 'CBM83', 'CBM69')
;

CREATE VIEW starch_active_domain_hit AS
SELECT feature_id, domain_id, score
FROM starch_active_gh_hit
UNION
SELECT feature_id, domain_id, score
FROM starch_active_cbm_hit
;

CREATE VIEW starch_active_domain_best_hit AS
SELECT
    feature_id
  , domain_id
  , score
FROM starch_active_domain_hit
JOIN (SELECT feature_id, MAX(score) AS max_score
      FROM starch_active_domain_hit
      GROUP BY feature_id
     ) USING (feature_id)
WHERE score = max_score
;

CREATE VIEW putative_susC AS
SELECT DISTINCT feature_id
FROM feature_details
LEFT JOIN feature_possible_ko USING (feature_id)
WHERE ko_id = 'K21573'
   OR product_description = 'TonB-dependent receptor SusC'
;

CREATE VIEW putative_susD AS
SELECT DISTINCT feature_id
FROM feature
LEFT JOIN feature_x_pfam_domain USING (feature_id)
LEFT JOIN feature_possible_ko USING (feature_id)
WHERE ko_id = 'K21572'
   OR domain_id LIKE 'SusD-like%'
;

CREATE VIEW putative_susEF AS
SELECT DISTINCT feature_id
FROM feature
LEFT JOIN feature_possible_ko USING (feature_id)
JOIN feature_x_pfam_domain USING (feature_id)
WHERE ko_id = 'K21571'
   OR domain_id LIKE '%SusE%'
   OR domain_id LIKE '%SusF%'
;

CREATE VIEW putative_susG AS
SELECT DISTINCT feature_id
FROM feature_details
JOIN feature_x_cazy_domain USING (feature_id)
WHERE lp_score > 5
  AND domain_id LIKE 'GH%'
;

CREATE VIEW putative_PUL_susC AS
SELECT DISTINCT
    seed_id AS feature_id
  , SUM(susC) AS tally_susC
  , SUM(susD) AS tally_susD
  , SUM(susEF) AS tally_susEF
  , SUM(susG) AS tally_susG
  , SUM(starch_active_gh) AS tally_starch_active_gh
  , SUM(starch_active_cbm) AS tally_starch_active_cbm
FROM feature_neighborhood
JOIN feature_details USING (feature_id)
LEFT JOIN (SELECT feature_id, 1 AS susC FROM putative_susC) AS c USING (feature_id)
LEFT JOIN (SELECT feature_id, 1 AS susD FROM putative_susD) AS d USING (feature_id)
LEFT JOIN (SELECT feature_id, 1 AS susEF FROM putative_susEF) AS e USING (feature_id)
LEFT JOIN (SELECT feature_id, 1 AS susG FROM putative_susG) AS g USING (feature_id)
LEFT JOIN (SELECT DISTINCT feature_id, 1 AS starch_active_gh
           FROM starch_active_gh_hit
          ) AS gh USING (feature_id)
LEFT JOIN (SELECT DISTINCT feature_id, 1 AS starch_active_cbm
           FROM starch_active_cbm_hit
          ) AS cbm USING (feature_id)
WHERE DISTANCE < 15000
  AND seed_id IN (SELECT * FROM putative_susC)
GROUP BY seed_id
HAVING tally_susC > 0 AND tally_susD > 0
;
