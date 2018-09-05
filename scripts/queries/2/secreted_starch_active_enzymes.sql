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
   OR domain_id IN ('CBM26', 'CBM25', 'CBM20', 'CBM69', 'CBM48', 'CBM58')
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
  , feature_signal_peptide.score AS signalp_score
  , feature_lipop.score AS lipop_score
  , opf_id
  , architecture
  , closest_cysteine AS signalp_closest_cysteine
  , feature_signal_peptide.cleavage_position AS signalp_cleavage_position
  , feature_lipop.margin AS lipop_score_margin
  , feature_lipop.cleavage_position AS lipop_cleavage_position
  , matched_domain
  , matched_domain_score
FROM feature_signal_peptide
JOIN feature USING (feature_id)
JOIN sequence USING (sequence_id)
LEFT JOIN feature_to_opf USING (feature_id)
LEFT JOIN feature_to_architecture USING (feature_id)
JOIN domain_best_hits USING (feature_id)
LEFT JOIN feature_lipop USING (feature_id)
WHERE ( (feature_signal_peptide.score > 0.5 AND ABS(signalp_closest_cysteine) < 4)
     OR feature_lipop.score > 1
      )
  -- AND mag_id IN ('Otu0007_vA', 'Otu0001_vB', 'Otu0001_vC')
ORDER BY mag_id, feature_id, matched_domain_score DESC
;
