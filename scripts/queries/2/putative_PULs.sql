-- Find all the likely PULs by finding proximate susC, susD, and susEF loci.
--
-- TODO: Fix this.  There should be far more found in B. theta.
-- Probably means that the filters for putative susCDEF are too tight/wrong?
-- Try using a different search string (KOs not working for some reason?)

CREATE TEMP VIEW putative_susC AS
SELECT DISTINCT feature_id
FROM feature_details
WHERE opf_id IN (SELECT opf_id FROM feature_details WHERE architecture = 'CarbopepD_reg_2:Plug:TonB_dep_Rec')
;

CREATE TEMP VIEW putative_susD AS
SELECT DISTINCT feature_id
FROM feature_details
WHERE opf_id IN (SELECT opf_id FROM feature_details WHERE ko_id = 'K21572')
;

-- CREATE TEMP VIEW putative_susCD AS
-- SELECT * FROM putative_susC
-- UNION
-- SELECT * FROM putative_susD
-- ;

CREATE TEMP VIEW putative_susEF AS
SELECT DISTINCT feature_id
FROM feature_details
WHERE opf_id IN (SELECT opf_id FROM feature_details WHERE ko_id = 'K21571')
;

SELECT DISTINCT
    seed_id AS susC
  , SUM(susC) AS tally_susC
  , SUM(susD) AS tally_susD
  , SUM(susEF) AS tally_susEF
FROM feature_neighborhood
JOIN feature_details USING (feature_id)
LEFT JOIN (SELECT feature_id, 1 AS susC FROM putative_susC) AS c USING (feature_id)
LEFT JOIN (SELECT feature_id, 1 AS susD FROM putative_susD) AS d USING (feature_id)
LEFT JOIN (SELECT feature_id, 1 AS susEF FROM putative_susEF) AS e USING (feature_id)
WHERE DISTANCE < 20000
  AND seed_id IN (SELECT * FROM putative_susC)
GROUP BY seed_id
HAVING tally_susC > 0 AND tally_susD > 0 AND tally_susEF > 0
;
