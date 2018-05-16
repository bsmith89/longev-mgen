-- Query longev results DB to get long-form table of 16S based relative
-- abundance estimates.
SELECT
    extraction_id
  , taxon_id AS sequence_id
  , taxon_id_b AS otu_id
  , SUM(tally)
FROM rrs_library
JOIN _rrs_library_taxon_count USING (rrs_library_id)
JOIN taxonomy USING (taxon_id, taxon_level)
JOIN extraction USING (extraction_id)
WHERE taxon_level = 'unique'
GROUP BY extraction_id, sequence_id
;
