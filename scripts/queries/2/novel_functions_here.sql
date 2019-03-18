CREATE TEMP VIEW mag_having_ko
SELECT mag_id, ko_id, COUNT(feature_id) > 0 AS present
from feature
JOIN sequence USING (feature_id)
JOIN feature_to_ko USING (feature_id)
GROUP BY mag_id, ko_id
;

CREATE TEMP VIEW mag_class_tally
SELECT mag_class, COUNT(mag_id) AS tally
FROM mag
GROUP BY mag_class
;

CREATE TEMP VIEW ko_frequency
SELECT mag_class, ko_id, SUM(present) * 1.0 / mag_class_tally AS frequency
FROM mag_having_ko
JOIN mag USING (mag_id)
JOIN mag_class_tally USING (mag_class)
GROUP BY (mag_class, ko_id)
;

SELECT * FROM ko_frequency LIMIT 5;
