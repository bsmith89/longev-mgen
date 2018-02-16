-- Query longev results DB to get long-form table of 16S based relative
-- abundance estimates.
SELECT extraction_id, taxon_id, relative_abundance AS relative_abundance
FROM rrs_library
JOIN mgen_library USING (extraction_id)
JOIN rrs_library_taxon_relative_abundance USING (rrs_library_id)
WHERE taxon_level = 'otu-0.03'
GROUP BY extraction_id
;
