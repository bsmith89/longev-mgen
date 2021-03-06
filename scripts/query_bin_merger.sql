SELECT
    b.*
  , 1.0 * b.read_count / (n1.length + n2.length) AS linkage_score_1
  , (b.read_count * b.read_count) / (c1.coverage * c2.coverage) / (n1.length + n2.length) AS linkage_score_2
  , (b.read_count * b.read_count) / (c1.coverage * c2.coverage) AS linkage_score_3
  , n1.n_contigs AS ncontigs_1
  , n2.n_contigs AS ncontigs_2
  , n1.length AS length_1
  , n2.length AS length_2
  , d1.completeness AS completeness_1
  , d2.completeness AS completeness_2
  , c1.coverage AS coverage_1
  , c2.coverage AS coverage_2
  , m.score AS complementarity_score
FROM bin_linkage AS b
JOIN bin_length AS n1 ON n1.bin_id = bin_id_1
JOIN bin_length AS n2 ON n2.bin_id = bin_id_2
JOIN bin_checkm AS d1 ON d1.bin_id = bin_id_1
JOIN bin_checkm AS d2 ON d2.bin_id = bin_id_2
JOIN (SELECT bin_id, SUM(coverage) AS coverage FROM bin_coverage GROUP BY bin_id) AS c1 ON c1.bin_id = bin_id_1
JOIN (SELECT bin_id, SUM(coverage) AS coverage FROM bin_coverage GROUP BY bin_id) AS c2 ON c2.bin_id = bin_id_2
LEFT JOIN bin_complementarity AS m ON m.bin_id = bin_id_1 AND m.bin_id_combined = bin_id_2
ORDER BY linkage_score_1 DESC, complementarity_score DESC
;
