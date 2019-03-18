SELECT DISTINCT
    feature_id
  , left
  , right
  -- , distance
  , opf_id
  , localization
  , tmhelix_count
  , architecture
  , product_description
FROM feature_neighborhood
LEFT JOIN (SELECT feature_id AS seed_id, opf_id AS seed_opf_id FROM feature_to_opf) USING (seed_id)
LEFT JOIN feature_localization USING (feature_id)
LEFT JOIN feature_to_opf USING (feature_id)
LEFT JOIN feature_to_architecture USING (feature_id)
JOIN feature_details USING (feature_id)
WHERE distance < 15000
  AND seed_opf_id IN ('Opf01517', 'Opf01942', 'Opf01276', 'Opf01351')
ORDER BY feature_id
;
