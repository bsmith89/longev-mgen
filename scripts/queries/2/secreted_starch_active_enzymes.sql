-- Find all of the potentially secreted proteins that are potentially active on starch.

CREATE TEMP VIEW domain_hits AS
SELECT
    feature_id
  , domain_id AS matched_domain
  , feature_domain.score AS matched_domain_score
FROM feature_domain
WHERE
   -- domain_id LIKE 'CBM%'
    ( domain_id LIKE 'GH13|_%' ESCAPE '|' OR domain_id = 'GH13'
   OR domain_id LIKE 'GH97%'
   OR domain_id LIKE 'GH31%'
   OR domain_id LIKE 'GH57%'
   OR domain_id LIKE 'GH70%'
   OR domain_id LIKE 'GH77%'
   OR domain_id IN ('CBM20', 'CBM21', 'CBM25',
                    'CBM26', 'CBM41', 'CBM48',
                    'CBM53', 'CBM58', 'CBM74',
                    'CBM82', 'CBM83')
    )
;

CREATE TEMP VIEW domain_best_hits AS
SELECT domain_hits.*
FROM domain_hits
JOIN ( SELECT feature_id, MAX(matched_domain_score) AS max_score
       FROM domain_hits
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
WHERE localization IN ('OM', 'OM?', 'OM??', 'PP?/OM?')
ORDER BY feature_id
;
