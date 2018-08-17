-- {{{1 Tables

-- CREATE TABLE mag
-- ( mag_id PRIMARY KEY
-- , description
-- );

CREATE TABLE checkm
( mag_id PRIMARY KEY --REFERENCES mag(mag_id)
, n_markers INT
, completeness FLOAT
, contamination FLOAT
, strain_heterogeneity FLOAT
);

CREATE TABLE quast
( mag_id PRIMARY KEY REFERENCES checkm(mag_id)
, n_contigs INT
, n_contigs_gt_1000 INT
, n_contigs_gt_5000 INT
, n_contigs_gt_10000 INT
, n_contigs_gt_25000 INT
, n_contigs_gt_50000 INT
, total_length INT
, total_length_gt_1000 INT
, total_length_gt_5000 INT
, total_length_gt_10000 INT
, total_length_gt_25000 INT
, total_length_gt_50000 INT
, n_contigs_ INT
, length_largest_contig INT
, total_length_ INT
, gc_percent FLOAT
, n50 INT
, n75 INT
, l50 INT
, l75 INT
, ambiguous_bases_per_100kbp FLOAT
);

CREATE TABLE sequence
( sequence_id PRIMARY KEY
, mag_id REFERENCES checkm(mag_id)
);

CREATE TABLE feature
( feature_id PRIMARY KEY
, sequence_id REFERENCES sequence(sequence_id)
, left INT
, right INT
);
CREATE INDEX idx_feature__sequence_id ON feature(sequence_id);

CREATE TABLE feature_details
( feature_id PRIMARY KEY REFERENCES feature(feature_id)
, ftype
, nlength INT
, product_description
);

CREATE TABLE ko
( ko_id PRIMARY KEY
, description
);

CREATE TABLE feature_to_ko
( feature_id REFERENCES feature(feature_id)
, ko_id REFERENCES ko(ko_id)
);

CREATE TABLE cog
( cog_id PRIMARY KEY
, function_category
, description
);

CREATE TABLE feature_to_cog
( feature_id REFERENCES feature(feature_id)
, cog_id REFERENCES cog(cog_id)
);


CREATE TABLE feature_to_opf
( feature_id REFERENCES feature(feature_id)
, opf_id
);

-- {{{1 Views

CREATE VIEW opf AS
SELECT
    opf_id
  , COUNT(opf_id) AS feature_count
FROM feature_to_opf
GROUP BY opf_id
;
