-- Find all of the putative SusC features while being robust to annotation errors.
--
-- This is a pretty liberal filter.

SELECT opf_id, COUNT(feature_id) AS tally
FROM feature_details
WHERE opf_id IN (SELECT opf_id FROM feature_details WHERE architecture = 'CarbopepD_reg_2:Plug:TonB_dep_Rec')
GROUP BY opf_id
ORDER BY tally DESC
;
