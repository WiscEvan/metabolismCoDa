#!/usr/bin/env python
# Add Organism Name column for easy conversion when visualizing...

import os
import pandas as pd
import glob
import argparse


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", help="Path to kegg module totals table or kegg ortholog totals table", required=True)
    parser.add_argument("--out", help="Path to write table with database metadata added", required=True)
    parser.add_argument("--dbdir", help="Path to directory containing ncbi_dataset.tsv files", required=True)
    args = parser.parse_args()

    search_str = os.path.join(args.dbdir, "*.tsv")
    df = pd.read_table(args.input)
    db_df = pd.concat(pd.read_table(fp) for fp in glob.glob(search_str)).drop_duplicates(subset=['Assembly Accession'])
    db_df['filename'] = db_df['Assembly Accession'].map(lambda x: f"{x}.kofamscan.tsv")

    main_df = pd.merge(df, db_df, how='left', left_on='filename', right_on='filename')
    main_df.to_csv(args.out, sep='\t', index=False, header=True)
    print(f"Wrote table with metadata added to {args.out}")

if __name__ == '__main__':
    main()