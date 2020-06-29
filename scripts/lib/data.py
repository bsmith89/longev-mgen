import sqlite3
import pandas as pd

def load_data(db_path):
    con = sqlite3.connect(db_path)
    conc = (pd.read_sql(
        """
        SELECT mouse_id, molecule_id, concentration
        FROM mouse_concentration
        WHERE sample_count > 0
        AND (    (channel = 'RI' AND molecule_id = 'butyrate')
                OR (channel = 'RI' AND molecule_id = 'acetate')
                OR (channel = 'RI' AND molecule_id = 'propionate')
                OR (channel = 'RI' AND molecule_id = 'glucose')
                OR (channel = 'RI' AND molecule_id = 'lactate')
                OR (channel = 'RI' AND molecule_id = 'succinate')
            )
        """,
        con=con, index_col=['mouse_id', 'molecule_id'])
            ['concentration'].unstack('molecule_id', fill_value=0)
        )

    mols = ['acetate', 'butyrate', 'succinate', 'lactate', 'propionate', 'glucose']
    mol_c_count = {'acetate': 2, 'butyrate': 4, 'succinate': 4,
                'lactate': 3, 'propionate': 3, 'glucose': 6}

    conc[[m + '_c' for m in mols]] = conc[mols].apply(lambda x: x * mol_c_count[x.name])
    conc['fermented_c'] = conc[[m + '_c' for m in mols]].drop('glucose_c', axis='columns').sum(axis='columns')
    conc[[m + '_frac' for m in mols]] = conc[mols].apply(lambda x: x / conc.fermented_c)
    conc['total_scfa'] = conc.acetate + conc.butyrate + conc.propionate
    conc['propionate_scfa_frac'] = conc.propionate / conc.total_scfa
    conc['butyrate_scfa_frac'] = conc.butyrate / conc.total_scfa
    conc['acetate_scfa_frac'] = conc.acetate / conc.total_scfa

    taxonomy = (pd.read_sql(
        """
        SELECT taxon_id AS otu_id, taxon_level_b AS level, taxon_id_b
        FROM taxonomy
        WHERE taxon_level = 'otu-0.03'
        AND taxon_level_b IN ('phylum', 'class', 'order', 'family', 'genus')
        AND confidence > 0.7
        AND taxon_id_b NOT LIKE '%_unclassified'
        """,
        con=con, index_col=['otu_id', 'level'])
                        .taxon_id_b.unstack('level')
                        [['phylum', 'class', 'order', 'family', 'genus']]
                        .apply(lambda x: x.str.replace('[-.]', '_')))

    abund = (pd.read_sql("""
        SELECT mouse_id, rrs_library_id, taxon_id, absolute_abundance
        FROM rrs_library_taxon_absolute_abundance
        JOIN rrs_library_metadata USING (rrs_library_id)
        WHERE taxon_level = 'otu-0.03'
                        """, con=con,
                        index_col=['mouse_id', 'rrs_library_id', 'taxon_id'])
            # Reshape into wide-format
            ['absolute_abundance'].unstack().fillna(0)
            # Drop libraries without an associated mouse_id
            .reset_index().dropna(subset=['mouse_id']).set_index('mouse_id')
            .drop('rrs_library_id', axis='columns')
            )
    assert abund.index.is_unique

    rabund = abund.apply(lambda x: x / x.sum(), axis='columns')

    abund_family = abund.groupby(taxonomy.family, axis=1).sum()
    abund_family['unclassified'] = (abund.sum(axis=1) - abund_family.sum(axis=1))
    abund_family.unclassified[abund_family.unclassified < 0] = 0
    rabund_family = rabund.groupby(taxonomy.family, axis=1).sum()
    rabund_family['unclassified'] = 1 - rabund_family.sum(axis=1)
    rabund_family.unclassified[rabund_family.unclassified < 0] = 0

    families = list(abund_family.columns)

    mouse = pd.read_sql(
        """
        SELECT
            mouse_id
        , cohort
        , sex
        , treatment
        , site
        , cage_id
        , age_at_death_or_censor
        , censored
        FROM mouse
        WHERE mouse_id NOT NULL
        """, con=con, index_col='mouse_id')
    mouse.site = mouse.site.apply(lambda x: 'TJL' if x == 'JL' else x)

    sample = pd.read_sql(
        """
        SELECT
            sample_id
        , mouse_id
        , age_at_collect
        , weight AS sample_weight
        , hydrated_weight AS sample_hydrated_weight
        , weight / hydrated_weight AS hydration_factor
        FROM sample_metadata
        JOIN (SELECT mouse_id, MAX(age_at_collect) AS max_age_at_collect
            FROM sample_metadata
            GROUP BY mouse_id
            ) USING (mouse_id)
        WHERE age_at_collect = max_age_at_collect
        """, con=con, index_col='mouse_id')

    mouse = mouse.join(sample, how='left')
    assert mouse.index.is_unique
    assert len(mouse.censored.unique()) <= 3  # Only 1.0, 0.0, and NaN
    mouse['dead'] = mouse['censored'].map({1.0: False, 0.0: True})
    mouse.drop('censored', axis='columns', inplace=True)
    mouse['dens'] = abund.sum(axis='columns')
    count = (pd.read_sql("""
        SELECT mouse_id, rrs_library_id, taxon_id, tally
        FROM rrs_library_taxon_count
        JOIN rrs_library_metadata USING (rrs_library_id)
        WHERE taxon_level = 'otu-0.03'
        AND spike_id NOT NULL
                        """, con=con,
                        index_col=['mouse_id', 'rrs_library_id', 'taxon_id'])
            # Reshape into wide-format
            ['tally'].unstack().fillna(0).astype(int)
            # Drop libraries without an associated mouse_id
            .reset_index().dropna(subset=['mouse_id']).set_index('mouse_id')
            .drop('rrs_library_id', axis='columns')
            )
    assert count.index.is_unique

    count_family = count.groupby(taxonomy.family, axis=1).sum()
    count_family['unclassified'] = (count.sum(axis=1) - count_family.sum(axis=1))

    return { 'con': con
           , 'conc': conc
           , 'mols': mols
           , 'mol_c_count': mol_c_count
           , 'taxonomy': taxonomy
           , 'abund': abund
           , 'rabund': rabund
           , 'abund_family': abund_family
           , 'rabund_family': rabund_family
           , 'families': families
           , 'mouse': mouse
           , 'sample': sample
           , 'count': count
           , 'count_family': count_family
           }
