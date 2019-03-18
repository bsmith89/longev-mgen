CREATE TEMP VIEW pul_x_domain AS
SELECT
    seed_id
  , feature_id
  , domain_id
  , domain_start
  , domain_stop
FROM
  ( SELECT feature_id AS seed_id, strand
    FROM pul_susC
    JOIN feature_to_strand USING (feature_id)
  )
JOIN feature_distance USING (seed_id)
JOIN feature_to_strand USING (feature_id, strand)
JOIN feature_x_cazy_minimal_domain USING (feature_id)
WHERE distance < 15000
;

CREATE TEMP VIEW starch_active_pul AS
SELECT DISTINCT
    feature_id AS seed_id
FROM pul_x_domain
LEFT JOIN starch_active_cbm_domain AS cbm USING (domain_id)
LEFT JOIN starch_active_gh_domain AS gh USING (domain_id)
GROUP BY seed_id
HAVING COUNT(cbm.domain_id) > 0 AND COUNT(gh.domain_id) > 0
;

CREATE TEMP VIEW opf_starch_pul_count AS
SELECT opf_id, COUNT(feature_id) AS tally
FROM ( SELECT DISTINCT feature_id
       FROM starch_active_pul
       JOIN feature_distance USING (seed_id)
       WHERE distance < 10000
     )
JOIN feature_to_opf USING (feature_id)
WHERE opf_id NOT NULL
GROUP BY opf_id
;

SELECT opf_id, s.tally, a.tally
FROM opf_starch_pul_count AS s
JOIN ( SELECT
           opf_id
         , COUNT(feature_id) AS tally
       FROM feature_to_opf
       GROUP BY opf_id
     ) AS a USING (opf_id)
ORDER BY s.tally DESC
;


-- CREATE TEMP VIEW susEF AS
-- SELECT DISTINCT feature_id
-- FROM feature
-- LEFT JOIN feature_possible_ko USING (feature_id)
-- JOIN feature_x_pfam_domain USING (feature_id)
-- WHERE ko_id = 'K21571'
--    OR domain_id LIKE '%SusE%'
--    OR domain_id LIKE '%SusF%'
-- ;
--
-- CREATE TEMP VIEW susG AS
-- SELECT DISTINCT feature_id
-- FROM feature_details
-- JOIN feature_x_cazy_domain USING (feature_id)
-- WHERE lp_score > 5
--   AND domain_id LIKE 'GH%'
-- ;

-- CREATE TEMP VIEW pul AS
-- SELECT DISTINCT
--     seed_id AS seed_id
-- FROM feature_neighborhood
-- LEFT JOIN (SELECT feature_id, 1 AS susC, strand AS susC_strand FROM susC JOIN feature_to_strand USING (feature_id)) AS c USING (feature_id)
-- LEFT JOIN (SELECT feature_id, 1 AS susD, strand AS susD_strand FROM susD JOIN feature_to_strand USING (feature_id)) AS d USING (feature_id)
-- WHERE DISTANCE < 10000
--   AND seed_id IN (SELECT * FROM susC)
-- GROUP BY seed_id
-- HAVING SUM(susC) > 0 AND SUM(susD) > 0
-- ;
--
-- CREATE TEMP VIEW starch_pul AS
-- SELECT
--     seed_id
-- FROM pul
-- JOIN feature_neighborhood USING (seed_id)
-- JOIN (SELECT feature_id, 1 AS starch_active_gh FROM starch_active_gh_hit) USING (feature_id)
-- JOIN (SELECT feature_id, 1 AS starch_active_cbm FROM starch_active_cbm_hit) USING (feature_id)
-- WHERE DISTANCE < 10000
-- GROUP BY seed_id
-- HAVING SUM(starch_active_gh) > 0 AND SUM(starch_active_cbm) > 0
-- ;
--
-- -- CREATE TEMP VIEW putative_PUL AS
-- SELECT DISTINCT
--     genome_id
--   , putative_PUL_susC.*
-- FROM putative_PUL_susC
-- JOIN feature USING (feature_id)
-- JOIN sequence USING (sequence_id)
-- ;
--
-- CREATE TEMP VIEW putative_PUL_count AS
-- SELECT genome_id, SUM(MIN(tally_susC, tally_susD)) AS tally
-- FROM putative_PUL
-- JOIN sequence USING (sequence_id)
-- GROUP BY genome_id
-- ;
--
-- SELECT
--     genome_id
--   , tally
--   , completeness
--   , contamination
--   , (tally * (1 - (contamination / 100))) / (completeness / 100) AS adjusted_tally
-- FROM putative_PUL_count
-- JOIN checkm USING (genome_id)
-- ORDER BY adjusted_tally DESC
-- ;
