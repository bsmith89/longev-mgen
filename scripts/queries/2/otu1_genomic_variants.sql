-- Find PULs specific to one strain of OTU-1 or the other.
-- Find all of the putative SusD features while being robust to annotation errors.

SELECT DISTINCT
    feature_id
  , left
  , right
  , susC
  , susD
  , susEF
  , susG
  , starch_active_domain
  , localization
  , coverage_ratio AS cross_coverage_ratio
  , opf_id
  , architecture
  , product_description
  , cazy_domain_list
FROM (SELECT feature_id, sequence_id FROM feature)  -- So I don't lose non-coding features
JOIN sequence USING (sequence_id)
LEFT JOIN feature_details USING (feature_id)
LEFT JOIN (SELECT feature_id, 1 AS susC FROM putative_susC) AS c USING (feature_id)
LEFT JOIN (SELECT feature_id, 1 AS susD FROM putative_susD) AS d USING (feature_id)
LEFT JOIN (SELECT feature_id, 1 AS susEF FROM putative_susEF) AS e USING (feature_id)
LEFT JOIN (SELECT feature_id, 1 AS susG FROM putative_susG) AS g USING (feature_id)
LEFT JOIN (SELECT DISTINCT feature_id, domain_id AS starch_active_domain
           FROM starch_active_domain_best_hit
          ) AS s USING (feature_id)
LEFT JOIN (SELECT
               feature_id
             , GROUP_CONCAT(domain_id, ',') AS cazy_domain_list
           FROM feature_x_cazy_minimal_domain
           GROUP BY feature_id) AS ca USING (feature_id)
LEFT JOIN (SELECT feature_id, coverage_ratio FROM variant_cross_coverage) USING (feature_id)
WHERE mag_id IN ('Otu0001_vB', 'Otu0001_vC', 'Otu0007_vA', 'Bacteroides_ovatus_ATCC_8483', 'Bacteroides_thetaiotaomicron_VPI5482')
;
