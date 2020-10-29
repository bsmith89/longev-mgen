CREATE TEMP VIEW feature_with_function_annotation AS
SELECT DISTINCT feature_id
FROM feature_details
LEFT JOIN feature_x_ko USING (feature_id)
WHERE (product_description NOT LIKE '%hypothetical protein%')
   OR (ko_id NOT NULL)
   OR (cog_id NOT NULL)
;

CREATE TEMP VIEW opf_details AS
SELECT
  opf_id,
  COUNT(DISTINCT feature_id) AS feature_tally,
  COUNT(DISTINCT genome_id) AS genome_tally,
  feature_id IN (SELECT feature_id FROM feature_with_function_annotation) AS has_function_annotation
FROM feature_to_opf
JOIN feature USING (feature_id)
JOIN sequence USING (sequence_id)
GROUP BY opf_id
HAVING opf_id NOT NULL
;

SELECT COUNT(DISTINCT feature_id) AS tally_all_features
FROM feature_details
;

SELECT COUNT(DISTINCT opf_id) AS tally_all_opfs
FROM opf_details
WHERE genome_tally > 0
;

SELECT SUM(feature_tally) AS tally_all_features_in_opfs
FROM opf_details
;

SELECT COUNT(DISTINCT opf_id) AS tally_opf_in_gt2_genomes
FROM opf_details
WHERE genome_tally > 2
;

SELECT COUNT(DISTINCT opf_id) AS tally_opfs_in_gt2_genomes_and_with_annotated_members
FROM opf_details
WHERE genome_tally > 2
AND has_function_annotation = 1
;

SELECT COUNT(DISTINCT opf_id) AS tally_opfs_in_gt2_genomes_and_without_annotated_members
FROM opf_details
WHERE genome_tally > 2
AND has_function_annotation = 0
;

SELECT SUM(feature_tally) AS tally_features_in_opfs_in_gt2_genomes_and_without_annotated_members
FROM opf_details
WHERE genome_tally > 2
AND has_function_annotation = 0
;

SELECT SUM(feature_tally) AS tally_features_in_opfs_in_gt5_genomes_and_without_annotated_members
FROM opf_details
WHERE genome_tally > 5
AND has_function_annotation = 0
;
