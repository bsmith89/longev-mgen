-- Identify likely associations between OPFs based on finding them as frequent
-- neighbors.
--
-- Counts may be off in cases where there are multiple of the same OPF within
-- the search radius.

SELECT
  a.opf_id AS feature_annot
, b.opf_id AS neighbor_annot
, COUNT(*) / 2 AS count_neighboring_pairs
FROM feature_neighborhood AS n
JOIN feature_to_opf AS a
  ON n.feature_id = a.feature_id
JOIN feature_to_opf AS b
  ON n.neighbor_id = b.feature_id
WHERE distance < 10000
GROUP BY feature_annot, neighbor_annot
ORDER BY count_neighboring_pairs DESC
;
