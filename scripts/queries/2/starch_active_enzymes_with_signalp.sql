-- Find all of the potentially secreted proteins that are potentially active on starch.

CREATE TEMP VIEW domain_hits AS
SELECT
    feature_id
  , domain_id AS matched_domain
  , feature_x_cazy_domain.score AS matched_domain_score
FROM feature_x_cazy_domain
WHERE
    ( domain_id LIKE 'GH13|_%' ESCAPE '|' OR domain_id = 'GH13'
   OR domain_id LIKE 'GH97%'
   OR domain_id LIKE 'GH31%'
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
  , feature_signal_peptide.score AS score
  , closest_cysteine
  , cleavage_position
  , opf_id
  , architecture
  , matched_domain
  , matched_domain_score
FROM feature_signal_peptide
JOIN feature USING (feature_id)
JOIN sequence USING (sequence_id)
LEFT JOIN feature_to_opf USING (feature_id)
LEFT JOIN feature_to_architecture USING (feature_id)
JOIN domain_best_hits USING (feature_id)
WHERE ABS(closest_cysteine) <= 4 AND feature_signal_peptide.score > 0.5
ORDER BY genome_id, feature_id, matched_domain_score DESC
;
