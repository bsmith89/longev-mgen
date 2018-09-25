CREATE TEMP VIEW feature_has_domain AS
SELECT DISTINCT feature_id, domain_id
FROM feature_cazy_domain
;

-- How many features are clustered into this OPF?
CREATE TEMP VIEW opf_count AS
SELECT opf_id, COUNT(feature_id) AS tally
FROM feature_to_opf
GROUP BY opf_id
;

-- How many of the features clustered into this OPF have the domain?
CREATE TEMP VIEW opf_has_domain_count AS
SELECT opf_id, domain_id, COUNT(feature_id) AS tally
FROM feature_to_opf
JOIN feature_has_domain USING (feature_id)
GROUP BY opf_id, domain_id
;

CREATE TEMP VIEW opf_domain AS
SELECT opf_id, domain_id, c1.tally * 1.0 / c2.tally AS frequency
FROM opf_has_domain_count AS c1
JOIN opf_count AS c2 USING (opf_id)
;

-- Which OPFs are associated with dbCAN domains?
SELECT opf_id, domain_id
FROM cazy_domain
JOIN opf_domain USING (domain_id)
WHERE domain.cazy_description != ''
  AND frequency > 0.5
;
