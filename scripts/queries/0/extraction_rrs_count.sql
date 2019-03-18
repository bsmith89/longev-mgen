SELECT
    extraction_id
  , otu_id
  , SUM(tally) AS tally
FROM rrs_taxon_count
GROUP BY extraction_id, otu_id
;
