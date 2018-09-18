-- Find all of the potentially secreted proteins that are potentially active on starch.

CREATE TEMP VIEW domain_best_hits AS
SELECT starch_active_domain_hits.*
FROM starch_active_domain_hits
JOIN ( SELECT feature_id, MAX(matched_domain_score) AS max_score
       FROM starch_active_domain_hits
       GROUP BY feature_id
     ) AS m USING (feature_id)
WHERE matched_domain_score = max_score
;

SELECT
    feature_id
  , sequence_id
  , left
  , right
  , opf_id
  , localization
  , architecture
  , matched_domain
FROM domain_best_hits
LEFT JOIN feature_to_opf USING (feature_id)
LEFT JOIN feature_localization USING (feature_id)
LEFT JOIN feature_to_architecture USING (feature_id)
JOIN feature USING (feature_id)
WHERE localization IN ('OM', 'OM?')
ORDER BY feature_id
;
