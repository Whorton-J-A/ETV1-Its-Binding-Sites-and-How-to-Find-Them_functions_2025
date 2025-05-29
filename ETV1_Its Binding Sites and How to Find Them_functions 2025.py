import pandas as pd
from pathlib import Path



def tf_chrom_locations(tfs=None, score_threshold=400, chr_num='chrall', save_dir='.'):# Filters TFBS based of motif similarity to JASPAR 2020; deterimed by UCSC GENENOME BOWSER.
    if tfs is None:
        raise ValueError("You must provide a list of TF DataFrames.")
        
    save_path = Path(save_dir)
    save_path.mkdir(parents=True, exist_ok=True)
    
    tfs_combined = pd.concat(tfs)
    tfs_combined_threshold = tfs_combined[tfs_combined['score'] >= score_threshold]
    tfs_chrom_groupings = dict(tuple(tfs_combined_threshold.groupby('chrom')))
    
    saved_files = []

    if chr_num == 'chrall':
        print(f'The chromosomes with a TF score threshold of at least {score_threshold} are {list(tfs_chrom_groupings.keys())}')
        for chr_key, tfs_chromdf in tfs_chrom_groupings.items():
            tfs_chromdf = tfs_chromdf.copy()
            tfs_chromdf['threshold'] = score_threshold
            file_path = save_path / f'locations_{chr_key}.csv'
            print(f'Saving TFs chromosome locations to machine as {file_path}')
            tfs_chromdf.to_csv(file_path, index=False)
            saved_files.append(str(file_path))
        return saved_files
    else:
        if chr_num not in tfs_chrom_groupings:
            raise ValueError(f"Chromosome '{chr_num}' not found with score threshold {score_threshold}")
        tfs_chromdf = tfs_chrom_groupings[chr_num].copy()
        tfs_chromdf['threshold'] = score_threshold
        file_path = save_path / f'locations_{chr_num}.csv'
        print(f'Saving TFs chromosome locations to machine as {file_path}')
        tfs_chromdf.to_csv(file_path, index=False)
        return str(file_path)



def tf_binding_by_gene(chip_df, tf_dict, save_dir='tf_gene_outputs'):# Finds TFBS in genes acrossed CHIP-SEQ dataset
    save_path = Path(save_dir)
    save_path.mkdir(parents=True, exist_ok=True)

    all_hits = []
    chromosomes = chip_df['chrom'].dropna().unique()

    for chrom in chromosomes:
        chip_chr = chip_df[chip_df['chrom'] == chrom]
        genes = chip_chr['Gene_Symbol'].dropna().unique()

        if chrom not in tf_dict:
            continue

        tf_chr_df = tf_dict[chrom]

        for gene in genes:
            gene_df = chip_chr[chip_chr['Gene_Symbol'] == gene]
            if gene_df.empty:
                continue

            gene_start = gene_df['txStart'].min()
            gene_end = gene_df['txEnd'].max()

            tf_hits = tf_chr_df[
                (tf_chr_df['txStart'] >= gene_start) &
                (tf_chr_df['txEnd'] <= gene_end)
            ].copy()

            if not tf_hits.empty:
                tf_hits['Gene_Symbol'] = gene
                tf_hits['chrom'] = chrom
                all_hits.append(tf_hits)

    if all_hits:
        combined = pd.concat(all_hits, ignore_index=True)
        tf_counts = combined.groupby(['chrom', 'Gene_Symbol', 'TF']).size().reset_index(name='Count')

        combined_file = save_path / 'tf_gene_raw_overlaps.csv'
        counts_file = save_path / 'tf_gene_tf_counts.csv'

        combined.to_csv(combined_file, index=False)
        tf_counts.to_csv(counts_file, index=False)

        return str(combined_file), str(counts_file), combined, tf_counts
    else:
        return None, None, pd.DataFrame(), pd.DataFrame()


def regulatory_units_by_chrom(df, save_dir=None, return_concat=False):# Finds candidate cis-regulatory elements (cCREs) from the ENCODE database acrossed chromosomes based on UCSC browser track.
    reg_units_by_chrom = dict(tuple(df.groupby('chrom')))

    if save_dir:
        save_path = Path(save_dir)
        save_path.mkdir(parents=True, exist_ok=True)
        for chrom, sub_df in reg_units_by_chrom.items():
            file_path = save_path / f"{chrom}_regulatory_units.csv"
            sub_df.to_csv(file_path, index=False)

    if return_concat:
        return reg_units_by_chrom, pd.concat(reg_units_by_chrom.values(), ignore_index=True)
    else:
        return reg_units_by_chrom
