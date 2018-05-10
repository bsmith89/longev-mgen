-- Query longev results DB to get long-form table of 16S based relative
-- abundance estimates.
SELECT
    sequence_id
  , taxon_id AS otu_id
  , domain_.taxon_id_b AS domain_
  , phylum_.taxon_id_b AS phylum_
  , class_.taxon_id_b AS class_
  , order_.taxon_id_b AS order_
  , family_.taxon_id_b AS family_
  , genus_.taxon_id_b AS genus_
  FROM (SELECT taxon_id AS sequence_id, taxon_id_b AS taxon_id FROM taxonomy WHERE taxon_level = 'unique' AND taxon_level_b = 'otu-0.03') AS s
  JOIN (SELECT taxon_id, taxon_id_b FROM taxonomy WHERE taxon_level_b = 'domain') AS domain_ USING (taxon_id)
  JOIN (SELECT taxon_id, taxon_id_b FROM taxonomy WHERE taxon_level_b = 'phylum') AS phylum_ USING (taxon_id)
  JOIN (SELECT taxon_id, taxon_id_b FROM taxonomy WHERE taxon_level_b = 'class') AS class_ USING (taxon_id)
  JOIN (SELECT taxon_id, taxon_id_b FROM taxonomy WHERE taxon_level_b = 'order') AS order_ USING (taxon_id)
  JOIN (SELECT taxon_id, taxon_id_b FROM taxonomy WHERE taxon_level_b = 'family') AS family_ USING (taxon_id)
  JOIN (SELECT taxon_id, taxon_id_b FROM taxonomy WHERE taxon_level_b = 'genus') AS genus_ USING (taxon_id)
;
