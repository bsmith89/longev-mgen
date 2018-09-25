-- Find all of the potentially secreted proteins that are potentially active on starch.

CREATE TEMP VIEW proximate_susC AS
SELECT DISTINCT n.feature_id
FROM feature_neighborhood AS n
JOIN putative_susC AS s ON s.feature_id = n.seed_id
WHERE distance < 15000
;

SELECT
    feature_id
  , sequence_id
  , left
  , right
  , opf_id
  , localization
  , proximate_susC
  , architecture
  , domain_id
FROM starch_active_domain_best_hit
JOIN (SELECT DISTINCT feature_id FROM starch_active_gh_hit) USING (feature_id) -- Only select those with hits to GHs.
LEFT JOIN feature_to_opf USING (feature_id)
LEFT JOIN feature_localization USING (feature_id)
LEFT JOIN feature_to_architecture USING (feature_id)
LEFT JOIN (SELECT feature_id, 1 AS proximate_susC FROM proximate_susC) USING (feature_id)
JOIN feature USING (feature_id)
WHERE lp_score > 5
ORDER BY feature_id
;
