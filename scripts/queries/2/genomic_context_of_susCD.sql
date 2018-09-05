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
LEFT JOIN (SELECT
               feature_id AS seed_id
             , product_description AS seed_product
           FROM feature_details
          ) USING (seed_id)
LEFT JOIN feature_localization USING (feature_id)
LEFT JOIN feature_to_opf USING (feature_id)
LEFT JOIN feature_to_architecture USING (feature_id)
LEFT JOIN feature_details USING (feature_id)
WHERE distance < 15000
  AND seed_product = 'TonB-dependent receptor SusC'
ORDER BY feature_id
;
