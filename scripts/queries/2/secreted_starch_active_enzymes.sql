-- Find all of the potentially secreted proteins that are potentially active on starch.
--
-- Potentially starch-active is defined as containing a GH13, GH31, or GH97 domain.
-- Potentially secreted is defined as having a signal peptide region which is
-- followed by a cysteine directly after the cleavage site, suggesting the
-- possibility of export to fully extracellular space.
--
-- I include signal peptide scores that are close to 0.5, despite this increasing
-- false positives, in part because Muribaculaceae wouldn't be used in training
-- SignalP, so they might be missed.
--
-- A larger radius around the putative cleavage site (+/- 4) is used, because
-- there can be ambiguity in cleavage position.
-- In cases where the cysteine is not at position 0, the full SignalP
-- output should be followed-up on.

CREATE TEMP VIEW domain_hits AS
SELECT
    feature_id
  , feature_signal_peptide.score AS score
  , closest_cysteine
  , cleavage_position
  , opf_id
  , architecture
  , domain_id AS matched_domain
  , feature_domain.score AS matched_domain_score
FROM feature_signal_peptide
JOIN feature USING (feature_id)
JOIN sequence USING (sequence_id)
LEFT JOIN feature_to_opf USING (feature_id)
LEFT JOIN feature_domain USING (feature_id)
LEFT JOIN feature_to_architecture USING (feature_id)
WHERE ABS(closest_cysteine) <= 4 AND feature_signal_peptide.score > 0.5
  AND ( domain_id LIKE 'GH13|_%' ESCAPE '|' OR domain_id = 'GH13'
     OR domain_id LIKE 'GH97%'
     OR domain_id LIKE 'GH31%'
      )
ORDER BY mag_id, feature_id, feature_domain.score DESC
;

SELECT domain_hits.*
FROM domain_hits
JOIN ( SELECT feature_id, MAX(matched_domain_score) AS max_score
       FROM domain_hits
       GROUP BY feature_id
     ) AS m USING (feature_id)
WHERE matched_domain_score = max_score
;
