#!/usr/bin/env python

import glob
import pandas as pd


def main():
    secondary_metabolites_fpaths = glob.glob("processed/Metabolism/Biosynthesis_of_other_secondary_metabolites/*.tsv")
    not_included_secondary_fpath = glob.glob("processed/Not_Included_in_Pathway_or_Brite/Unclassified_metabolism/Secondary_metabolism.tsv")

    secondary_fpaths = secondary_metabolites_fpaths + not_included_secondary_fpath

    df = pd.concat(pd.read_table(fp) for fp in secondary_fpaths)


    secondary_kos = df.ko_num.unique()

    print(f"recovered {df.ko_num.nunique():,} unique KOs to {df.cdesc.nunique():,} kegg categories ({len(secondary_fpaths):,} total files)")

    kos_lines = "\n".join(secondary_kos) + "\n"

    with open("secondary_kos.txt", "w") as fh:
        fh.write(kos_lines)


if __name__ == '__main__':
    main()