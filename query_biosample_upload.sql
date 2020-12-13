CREATE TEMPORARY TABLE _extraction_to_biosample ( extraction_id PRIMARY KEY REFERENCES extraction(extraction_id), biosample_id );
.import extraction_to_biosample.tsv _extraction_to_biosample

CREATE TEMPORARY VIEW sample_metadata AS
  SELECT
    sample_id
  , mouse_id
  , collection_date
  , collection_time
  , full_tube_weight - empty_tube_weight AS weight
  , hydrated_tube_weight - empty_tube_weight AS hydrated_weight
  , (full_tube_weight - empty_tube_weight) /
      (hydrated_tube_weight - empty_tube_weight)
      AS hydration_factor
  , sample.comments AS comments
  , cohort
  , sex
  , treatment
  , site
  , cage_id
  , date_of_birth
  , JULIANDAY(collection_date) - JULIANDAY(date_of_birth) AS age_at_collect
  , age_at_death_or_censor
  , censored
  FROM sample
  LEFT JOIN mouse
    USING (mouse_id)
;

CREATE TEMPORARY VIEW extraction_weight AS
  SELECT
    extraction_id
  , CASE WHEN extraction.weight IS NOT NULL
      THEN extraction.weight
      ELSE sample.weight
    END AS weight
  , CASE WHEN extraction.hydrated_weight IS NOT NULL
      THEN extraction.hydrated_weight
      ELSE sample.hydrated_weight
    END AS hydrated_weight
  FROM extraction
  LEFT JOIN sample_metadata AS sample
    USING (sample_id)
;

CREATE TEMPORARY VIEW extraction_metadata AS
  SELECT
    extraction_id
  , sample_id
  , weight.weight AS weight
  , weight.hydrated_weight AS hydrated_weight
  , sample.weight AS sample_weight
  , weight.weight / weight.hydrated_weight AS hydration_factor
  , volume
  , spike_id
  , spike_volume
  , dna_concentration
  , extraction.comments
  , mouse_id
  , collection_date
  , collection_time
  , sample.comments AS sample_comments
  , cohort
  , sex
  , treatment
  , site
  , cage_id
  , date_of_birth
  , age_at_collect
  , age_at_death_or_censor
  , censored
  FROM extraction
  LEFT JOIN sample_metadata AS sample
    USING (sample_id)
  LEFT JOIN extraction_weight AS weight
    USING (extraction_id)
;

-- CREATE TEMPORARY TABLE _mgen_library_to_extraction ( mgen_library_id PRIMARY KEY, extraction_id REFERENCES extraction(extraction_id));
-- .import library_to_extraction.tsv _mgen_library_to_extraction

-- SELECT * FROM extraction_metadata;

SELECT DISTINCT
  -- library_id,
  -- _extraction_to_biosample.biosample_id,
  extraction_id AS sample_name,
  'PRJNA448009' AS bioproject_accession,
  'mouse gut metagenome' AS organism,
  collection_date,
  collection_time || ':00' AS collection_time,
  'not applicable' AS env_broad_scale,
  'not applicable' AS env_local_scale,
  'feces' AS env_medium,
  CASE
    WHEN site = 'UM' THEN 'USA: Ann Arbor, MI'
    WHEN site = 'UT' THEN 'USA: San Antonio, TX'
    WHEN site = 'TJL' THEN 'USA: Bar Harbor, ME'
  END AS geo_loc_name,
  'Mus musculus' AS host,
  CASE
    WHEN site = 'UM' THEN '42.28 N 83.74 W'
    WHEN site = 'UT' THEN '29.51 N 98.58 W'
    WHEN site = 'TJL' THEN '44.37 N 68.20 W'
  END AS lat_lon,
  CAST(age_at_collect AS INT) || ' days' AS host_age,
  age_at_death_or_censor AS host_age_at_death,
  'stool' AS host_body_product,
  'LabDiet 5LG6' AS host_diet,
  'UM-HET3' AS host_genotype,
  sex AS host_sex,
  sample_id AS host_subject_id,
  CASE WHEN treatment IS 'acarbose' THEN 'acarbose' ELSE '' END AS host_treatment,
  CASE WHEN treatment IS 'acarbose' THEN '1000 ppm' ELSE '' END AS host_treatment_dosage,
  CASE WHEN treatment IS 'acarbose' THEN '8 months' ELSE '' END AS host_treatment_start,
  '10090' AS host_taxid,
  ROUND(sample_weight, 3) || ' g' AS samp_size,
  'less than 6 months' AS samp_store_dur,
  '-80 C' AS samp_store_temp,
  sample_id AS source_material_id,
  mouse_id AS host_id,
  CASE WHEN spike_volume NOT NULL AND spike_volume > 0 AND spike_volume IS NOT ''
    THEN 'Homogenizing in water; centrifuged to separate solid fraction; spiked with external standard (S. alaskensis)'
    ELSE 'Homogenizing in water; centrifuged to separate solid fraction'
  END AS samp_mat_process,
  CASE WHEN spike_volume NOT NULL AND spike_volume > 0 AND spike_volume IS NOT ''
    THEN 'Sphingopyxis alaskensis stationary phase culture'
    ELSE ''
  END AS external_standard,
  CASE
    WHEN spike_volume IS NOT '' AND weight IS NOT '' THEN weight || ' g'
    ELSE ''
  END AS samp_weight_external_standard_added,
  spike_volume AS external_standard_volume
  -- , *
FROM extraction_metadata
-- JOIN library USING (extraction_id)
LEFT JOIN _extraction_to_biosample USING (extraction_id)
JOIN library USING (extraction_id)
  -- WHERE biosample_id IS NULL
ORDER BY biosample_id
;
