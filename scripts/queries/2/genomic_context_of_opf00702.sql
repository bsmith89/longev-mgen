SELECT
    feature_id
  , left
  , right
  , opf_id
  , localization
  , architecture
  , product_description
FROM feature_neighborhood
LEFT JOIN (SELECT
               feature_id AS seed_id
             , architecture AS seed
           FROM feature_to_architecture
          ) USING (seed_id)
LEFT JOIN feature_details USING (feature_id)
WHERE distance < 10000
  AND seed LIKE '%Porph_ging%'
;
