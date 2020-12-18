CREATE TEMPORARY VIEW feature_to_ko_list AS
SELECT feature_id, group_concat(ko_id, '|') AS ko_list
FROM feature_x_ko
GROUP BY feature_id
;

CREATE TEMPORARY VIEW feature_to_cazy_domain_list AS
SELECT feature_id, group_concat(domain_id, '|') AS cazy_domain_list
FROM feature_x_cazy_domain
GROUP BY feature_id
;

CREATE TEMPORARY VIEW details AS
SELECT *
FROM feature_details
JOIN sequence USING (sequence_id)
LEFT JOIN (SELECT feature_id, 1 AS susC FROM susC) USING (feature_id)
LEFT JOIN (SELECT feature_id, 1 AS susD FROM susD) USING (feature_id)
LEFT JOIN (SELECT feature_id, 1 AS susEF FROM susEF) USING (feature_id)
LEFT JOIN feature_to_ko_list USING (feature_id)
LEFT JOIN feature_to_cazy_domain_list USING (feature_id)
;

SELECT
  genome_id,
  sequence_id,
  feature_id,
  feature_start,
  feature_stop,
  nlength,
  localization,
  susC,
  susD,
  susEF,
  opf_id,
  architecture,
  cog_id,
  ko_list,
  cazy_domain_list,
  product_description
FROM details
ORDER BY genome_id, sequence_id, feature_start
;
