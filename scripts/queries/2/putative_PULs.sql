-- Find all the likely PULs by finding proximate susC, susD, and susEF loci.
--
-- TODO: Consider using a susEF filter as well.
-- But I think that these are hard to identify consistently.

-- CREATE TEMP VIEW putative_PUL AS
SELECT DISTINCT
    mag_id
  , putative_PUL_susC.*
FROM putative_PUL_susC
JOIN feature USING (feature_id)
JOIN sequence USING (sequence_id)
;

-- CREATE TEMP VIEW putative_PUL_count AS
-- SELECT mag_id, SUM(MIN(tally_susC, tally_susD)) AS tally
-- FROM putative_PUL
-- JOIN sequence USING (sequence_id)
-- GROUP BY mag_id
-- ;
--
-- SELECT
--     mag_id
--   , tally
--   , completeness
--   , contamination
--   , (tally * (1 - (contamination / 100))) / (completeness / 100) AS adjusted_tally
-- FROM putative_PUL_count
-- JOIN checkm USING (mag_id)
-- ORDER BY adjusted_tally DESC
-- ;
