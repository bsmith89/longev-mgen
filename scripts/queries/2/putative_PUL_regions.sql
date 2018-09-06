-- Find all of the putative SusD features while being robust to annotation errors.
--
-- This is a pretty liberal filter.

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
    feature_id
  , left
  , right
  , susC
  , susD
  , susEF
  , localization
  , opf_id
  , architecture
  , product_description
FROM feature_neighborhood
JOIN feature_details USING (feature_id)
LEFT JOIN (SELECT feature_id, 1 AS susC FROM putative_susC) AS c USING (feature_id)
LEFT JOIN (SELECT feature_id, 1 AS susD FROM putative_susD) AS d USING (feature_id)
LEFT JOIN (SELECT feature_id, 1 AS susEF FROM putative_susEF) AS e USING (feature_id)
WHERE DISTANCE < 10000
  AND seed_id IN (SELECT * FROM putative_susC)
;
