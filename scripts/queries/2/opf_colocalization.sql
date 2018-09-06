-- Identify likely associations between OPFs based on finding them as frequent
-- neighbors.
--
-- Counts may be off in cases where there are multiple of the same OPF within
-- the search radius.

SELECT
  o1.opf_id AS opf_id_1
, o2.opf_id AS opf_id_2
, COUNT(*) / 2 AS tally
FROM feature_neighborhood AS n
JOIN feature_to_opf AS o1 USING (feature_id)
JOIN feature_to_opf AS o2 ON n.seed_id = o2.feature_id
WHERE distance < 10000
  AND o1.feature_id != o2.feature_id
GROUP BY opf_id_1, opf_id_2
ORDER BY tally DESC
;
