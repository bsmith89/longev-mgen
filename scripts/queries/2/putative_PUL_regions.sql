-- Find all of the putative SusD features while being robust to annotation errors.
--
-- This is a pretty liberal filter.

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
  , opf_id
  , architecture
  , product_description
  , cazy_domain_list
FROM feature_neighborhood
JOIN feature_details USING (feature_id)
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
WHERE DISTANCE < 25000
  AND seed_id IN (SELECT feature_id FROM putative_PUL_susC)
;
