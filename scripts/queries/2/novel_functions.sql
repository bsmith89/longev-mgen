
-- Interface
-- CREATE TEMP VIEW feature_x_function AS
-- SELECT feature_id, cog_id AS function_id
-- FROM feature_to_cog
-- ;
--
-- CREATE TEMP VIEW function_ AS
-- SELECT cog_id AS function_id, description
-- FROM cog
-- ;

CREATE TEMP VIEW feature_x_function AS
SELECT feature_id, opf_id AS function_id FROM feature_to_opf
;

CREATE TEMP VIEW function_ AS
SELECT opf_id AS function_id, architecture AS description
FROM opf_to_architecture
;

-- CREATE TEMP VIEW feature_x_function AS
-- SELECT feature_id, ko_id AS function_id FROM feature_x_ko
-- ;
--
-- CREATE TEMP VIEW function_ AS
-- SELECT ko_id AS function_id, description
-- FROM ko
-- ;

-- Subqueries
CREATE TEMP VIEW mag_having_function AS
SELECT mag_id, function_id, COUNT(feature_id) > 0 AS present
from feature
JOIN sequence USING (sequence_id)
JOIN feature_x_function USING (feature_id)
GROUP BY mag_id, function_id
;

CREATE TEMP VIEW mag_class_tally AS
SELECT mag_class, COUNT(mag_id) AS tally
FROM mag
GROUP BY mag_class
;

CREATE TEMP VIEW mag_class_function_frequency AS
SELECT mag_class, function_id, SUM(present) * 1.0 / tally AS frequency
FROM mag_having_function
JOIN mag USING (mag_id)
JOIN mag_class_tally USING (mag_class)
GROUP BY mag_class, function_id
;

-- Query
SELECT function_id
, CASE WHEN frequency_here IS NULL THEN 0 ELSE frequency_here END AS frequency_here
, CASE WHEN frequency_ormerod IS NULL THEN 0 ELSE frequency_ormerod END AS frequency_ormerod
, description
FROM function_
LEFT JOIN (SELECT function_id, frequency AS frequency_here FROM mag_class_function_frequency WHERE mag_class = 'here') USING (function_id)
LEFT JOIN (SELECT function_id, frequency AS frequency_ormerod FROM mag_class_function_frequency WHERE mag_class = 'ormerod') USING (function_id)
WHERE frequency_here > 0
ORDER BY frequency_ormerod, frequency_here DESC
;
