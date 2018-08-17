-- {{{1 Tables

-- general metadata {{{2
CREATE TABLE mouse
  ( mouse_id            TEXT PRIMARY KEY
  , cohort              TEXT
  , sex                 TEXT
  , treatment           TEXT       -- {acarbose, inulin, control, young (control),
                                   -- + cr (calorie-restricted), (17-alpha-)estradiol
                                   -- + rapamycin}
  , site                TEXT       -- {JL (Jackson Labs), UM (University of Michigan),
                                   -- + UT (University of Texas, San Antonio)}
  , cage_id                TEXT
  , date_of_birth       DATETIME
  , age_at_death_or_censor  INTEGER    -- in days
  , censored                BOOLEAN
  );

CREATE TABLE sample
  ( sample_id                   TEXT PRIMARY KEY
  , mouse_id             TEXT REFERENCES mouse(mouse_id)
  , collection_date      DATETIME
  , collection_time      DATETIME    -- Hour out of 24 local time
  , empty_tube_weight    FLOAT       -- in g
  , full_tube_weight     FLOAT       -- in g
  , hydrated_tube_weight FLOAT       -- in g of the sample after
                                     -- + homogenizing in water, spinning, and
                                     -- + pulling off the accessible supernatant
  , comments             TEXT
  );

CREATE TABLE extraction
  ( extraction_id                      TEXT PRIMARY KEY
  , sample_id               TEXT REFERENCES sample(sample_id)
  , weight                  FLOAT    -- of the sample extracted from in g
  , hydrated_weight         FLOAT    -- in g of the sample after
                                     -- + homogenizing in water, spinning, and
                                     -- + pulling off the accessible supernatant
  , volume                  FLOAT    -- in mL of water used to extract the sample
  , spike_id        TEXT             -- FIXME: This needs to be constrained to
                                     -- + values in spike(spike_id)
  , spike_volume    FLOAT            -- in uL
  , dna_concentration  FLOAT   -- ng/ul
  , comments                TEXT
  );

CREATE TABLE library
  ( library_id TEXT PRIMARY KEY
  , extraction_id TEXT
  , run TEXT
  , file_r1 TEXT
  , file_r2 TEXT
  );

CREATE TABLE library_asmbl_group
  ( library_id TEXT REFERENCES library
  , asmbl_group TEXT
  );

CREATE TABLE rrs_taxon_count
  ( extraction_id TEXT
  , sequence_id TEXT
  , otu_id TEXT
  , tally INT

  , PRIMARY KEY (extraction_id, sequence_id)
  );
CREATE INDEX idx_rrs_taxon_count__sequence_id ON rrs_taxon_count(sequence_id);

CREATE TABLE rrs_taxonomy
  ( sequence_id TEXT PRIMARY KEY
  , otu_id TEXT
  , domain_ TEXT
  , phylum_ TEXT
  , class_ TEXT
  , order_ TEXT
  , family_ TEXT
  , genus_ TEXT
  );
