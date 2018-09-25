-- Find all of the potentially secreted proteins that are potentially active on starch.

SELECT
    feature_id
  , sequence_id
  , left
  , right
  , opf_id
  , feature_lipop.score AS lipop_score
  , feature_lipop
  , architecture
  , domain
FROM starch_active_domain_best_hit
LEFT JOIN feature_to_opf USING (feature_id)
LEFT JOIN feature_lipop USING (feature_id)
LEFT JOIN feature_to_architecture USING (feature_id)
JOIN feature USING (feature_id)
WHERE feature_lipop.score > 5
ORDER BY feature_id
;
