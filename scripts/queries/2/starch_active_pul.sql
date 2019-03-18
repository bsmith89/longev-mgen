CREATE TEMP VIEW pul_x_domain AS
SELECT
    seed_id
  , feature_id
  , ss.strand
  , fs.strand
  , domain_id
  , domain_start
  , domain_stop
FROM
  ( SELECT feature_id AS seed_id, strand
    FROM pul_susC
    JOIN feature_to_strand USING (feature_id)
  ) AS ss
JOIN feature_distance USING (seed_id)
JOIN feature_to_strand AS fs USING (feature_id)
JOIN feature_x_cazy_minimal_domain USING (feature_id)
WHERE distance < 15000
AND ss.strand == fs.strand
;

SELECT DISTINCT
    seed_id
FROM pul_x_domain
LEFT JOIN starch_active_cbm_domain AS cbm USING (domain_id)
LEFT JOIN starch_active_gh_domain AS gh USING (domain_id)
GROUP BY seed_id
HAVING COUNT(cbm.domain_id) > 0 AND COUNT(gh.domain_id) > 0
;
