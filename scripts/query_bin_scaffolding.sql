SELECT
    s.contig_id AS contig_id_1
  , s.contig_id_linked AS contig_id_2
  , b1.bin_id AS bin_id_1
  , b2.bin_id AS bin_id_2
  , n1.length AS contig_length_1
  , n2.length AS contig_length_2
  , tally
FROM contig_linkage AS s
JOIN contig_bin AS b1 USING (contig_id)
JOIN contig_bin AS b2 ON b2.contig_id = s.contig_id_linked
JOIN contig AS n1 USING (contig_id)
JOIN contig AS n2 ON n2.contig_id = s.contig_id_linked
WHERE b1.bin_id = 'bin353'
  AND (b2.bin_id != 'bin353' OR b1.contig_id < b2.contig_id)
  AND tally > 50
