SELECT DISTINCT
    genome_id
  , sequence_id
  , feature_id
  , feature_start
  , feature_stop
  , (feature_start < feature_stop) AS positive_strand
  , susC
  , susD
  , susEF
  , starch_active_domain
  , localization
  , opf_id
  , architecture
  , product_description
  , cazy_domain_list
FROM (SELECT feature_id AS seed_id FROM pul_susC)
JOIN feature_distance USING (seed_id)
JOIN feature_details USING (feature_id)
LEFT JOIN (SELECT feature_id, 1 AS susC FROM susC) USING (feature_id)
LEFT JOIN (SELECT feature_id, 1 AS susD FROM susD) USING (feature_id)
LEFT JOIN (SELECT feature_id, 1 AS susEF FROM susEF) USING (feature_id)
LEFT JOIN (SELECT DISTINCT feature_id, domain_id AS starch_active_domain
           FROM starch_active_domain_best_hit
          ) USING (feature_id)
LEFT JOIN (SELECT
               feature_id
             , GROUP_CONCAT(domain_id, ',') AS cazy_domain_list
           FROM feature_x_cazy_minimal_domain
           GROUP BY feature_id) USING (feature_id)
JOIN sequence USING (sequence_id)
WHERE DISTANCE < 25000
;
