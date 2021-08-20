CREATE TEMP VIEW proximate_cazy_domain_tally AS
SELECT seed_id, COUNT(*) AS tally
FROM (SELECT feature_id AS seed_id FROM pul_susC)
JOIN feature_distance USING (seed_id)
JOIN feature_x_cazy_domain USING (feature_id)
WHERE distance < 10000
GROUP BY seed_id
;

SELECT genome_id, seed_id, tally
FROM (SELECT feature_id AS seed_id FROM pul_susC)
LEFT JOIN proximate_cazy_domain_tally USING (seed_id)
JOIN feature ON (seed_id = feature_id)
JOIN sequence USING (sequence_id)
;
