#!/usr/bin/env python

# mamba install -c bioconda autometa -y

from autometa.common.kmers import embed,normalize
import pandas as pd
import argparse
import logging
import os

logging.basicConfig(level=logging.DEBUG)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", help="Path to kofamscan_results_matrix.tsv (assumes matrix samples have .tsv in their column names)", required=True)
    parser.add_argument("--output", help="Output directory path to write embeddings", required=True)
    parser.add_argument("--secondary-kos", help="Path to secondary_kos.txt (generated with get_secondary_kos.py)", required=False)
    args = parser.parse_args()
    df = pd.read_table(args.input, index_col='ko_num')
    sample_cols = [col for col in df.columns if ".tsv" in col]
    if args.secondary_kos:
        sdf = pd.read_table(args.secondary_kos, header=None, names=['ko_num'])
        df = df.loc[~df.index.isin(sdf.ko_num)]
    dff = df[sample_cols].T.copy()
    methods = ["bhsne", "densmap", "umap"]
    norm_method = "ilr"
    for method in methods:
        norm_df = normalize(dff, method=norm_method)
        out = os.path.join(os.path.realpath(args.output), f'kofamscan_results_{method}.tsv')
        embed_df = embed(dff, method=method)
        embed_df.index.name = 'filename'
        embed_df.to_csv(out, sep='\t', index=True, header=True)
        # Create embedding with normalized profiles...
        normed_out = os.path.join(os.path.realpath(args.output), f'kofamscan_results_{norm_method}_{method}.tsv')
        nembed_df = embed(norm_df, method=method)
        nembed_df.index.name = 'filename'
        nembed_df.to_csv(normed_out, sep='\t', index=True, header=True)


if __name__ == '__main__':
    main()