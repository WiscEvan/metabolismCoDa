#!/usr/bin/env python

import argparse
from typing import List
import pandas as pd
import glob
import os


def read_kofamscan_detail_tsv(fpath: str) -> pd.DataFrame:
    filename = os.path.basename(fpath)
    kos = []
    with open(fpath) as fh:
        lines = fh.readlines()
    for line in lines:
        # header
        # gene name	KO	threshold	score	E-value	"KO definition"
        # any other genes are annotated as not-significant w/o the asterisk
        if not line.startswith("*"):
            continue
        (
            __,
            gene_name,
            ko_num,
            threshold,
            score,
            evalue,
            ko_def,
        ) = line.strip().split("\t")
        ko_def = ko_def.replace("'", "").replace('"', "")
        threshold, score, evalue = map(float, [threshold, score, evalue])
        kos.append(
            {
                "filename": filename,
                "gene_name": gene_name,
                "ko_num": ko_num,
                "threshold": threshold,
                "score": score,
                "evalue": evalue,
                "ko_def": ko_def,
            }
        )
    return pd.DataFrame(kos)


def add_ko_metadata(df: pd.DataFrame, db_fpaths: List[str]) -> pd.DataFrame:
    # br:ko2ec        EC number
    # br:ko2rn        Reaction
    # br:ko2cog       COG
    # br:ko2go        GO
    # br:ko2tc        TC number
    # br:ko2cazy      CAZy
    for db_fpath in db_fpaths:
        if not db_fpath:
            continue
        db_df = pd.read_table(db_fpath, index_col="#KO")
        df = pd.merge(df, db_df, left_on="ko_num", right_index=True, how="left")
    return df


def get_kegg_join_relations_format(
    df: pd.DataFrame,
    out: str,
    filenames: List[str],
    col_names: List[str] = ['gene_name','ko_num'],
) -> str:
    dff = df.loc[df.filename.isin(filenames)]
    dff = dff.set_index('filename').loc[filenames]
    lines = ""
    for filename, filename_df in dff.groupby(["filename"], sort=False)[col_names]:
        genome = filename.replace(".kofamscan.tsv", "")
        lines += f"# {genome}\n"
        lines += filename_df[col_names].to_csv(sep='\t', header=None, index=False) + "\n"
    with open(out, "w") as outfh:
        outfh.write(lines)
    return out


def kegg_decoder_format(df: pd.DataFrame)->pd.DataFrame:
    dff = df[['filename', 'gene_name', 'ko_num']].copy()
    # Example format:
    # {genome}_{gene_name}\t{ko}
    # Comamonas-testosteroni_1        K02313

    dff["genome"] = dff.filename.map(lambda fn: fn.replace(".kofamscan.tsv", "").replace("_", "-"))
    genome_gene_name_df = dff[['genome','gene_name']].astype(str).apply('_'.join, axis=1)
    genome_gene_name_df.name = 'genome_gene_name'
    return pd.concat([genome_gene_name_df, dff[['ko_num']]], axis=1)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--dirpath",
        help="Path to output directory of kofamscan with tables containing '*.kofamscan.tsv' in their filename",
        required=True,
    )
    parser.add_argument(
        "--ko2cazy",
        help="Path to CAZy (ko2cazy) database file",
        required=False,
    )
    parser.add_argument(
        "--ko2cog",
        help="Path to COG (ko2cog) database file",
        required=False,
    )
    parser.add_argument(
        "--ko2ec",
        help="Path to EC number (ko2ec) database file",
        required=False,
    )
    parser.add_argument(
        "--ko2go",
        help="Path to GO (ko2go) database file",
        required=False,
    )
    parser.add_argument(
        "--ko2rn",
        help="Path to reaction (ko2rn) database file",
        required=False,
    )
    parser.add_argument(
        "--ko2tc",
        help="Path to transporter classification (ko2tc) database file",
        required=False,
    )
    parser.add_argument(
        "--join-orgs",
        help="filenames corresponding to genomes to retrieve for kegg join relationship file",
        nargs="+",
        required=False,
    )
    parser.add_argument(
        "--keggdecoder-orgs",
        help="filenames corresponding to genomes to retrieve for kegg-decoder input",
        nargs="+",
        required=False,
        )
    parser.add_argument(
        "--out-table",
        help="Path to write merged counts table of kofamscan results with metadata",
        required=False,
        default="ko_counts.tsv",
    )
    parser.add_argument(
        "--out-matrix",
        help="Path to write merged counts matrix of kofamscan results",
        required=False,
        default="ko_counts_matrix.tsv",
    )
    parser.add_argument(
        "--out-kegg-join",
        help="Path to write merged kegg join_relations",
        required=False,
        default="kegg_join_relations.tsv",
    )
    parser.add_argument(
        "--out-keggdecoder",
        help="File path to write table formatted for input to KEGG-decoder",
        required=False,
        default="keggdecoder.tsv",
    )
    args = parser.parse_args()
    kofam_search_str = os.path.join(args.dirpath, "*kofamscan.tsv")
    df = pd.concat(
        [read_kofamscan_detail_tsv(fp) for fp in glob.glob(kofam_search_str)]
    )

    if args.join_orgs:
        join_relations_fpath = get_kegg_join_relations_format(df, out=args.out_kegg_join, filenames=args.join_orgs)
        print(f"Wrote kegg join relations file to {join_relations_fpath}")
    
    # This will create the full table of annotations
    db_fpaths = [
        args.ko2ec,
        args.ko2cog,
        args.ko2cazy,
        args.ko2go,
        args.ko2tc,
        args.ko2rn,
    ]
    df = add_ko_metadata(df, db_fpaths)

    df.to_csv(args.out_table, sep="\t", index=False, header=True)
    print(
        f"Wrote {df.shape[0]:,} knums of {df.filename.nunique():,} kofamscan results to {args.out_table}"
    )

    # Format for KEGG-decoder
    # get_kegg_colormap(cdf, outdir=args.out_colormap)
    if args.keggdecoder_orgs:
        cdf = df.loc[df.filename.isin(args.keggdecoder_orgs)].copy()
        keggdecoder_df = kegg_decoder_format(cdf)
        keggdecoder_df.to_csv(args.out_keggdecoder, index=False, header=False, sep='\t')
        print(f"Wrote KEGG-decoder formatted table to {args.out_keggdecoder}")


    # Now we need to create a separate matrix of KO x genome where each cell are the KO counts
    ko_counts = df[["filename", "ko_num"]].groupby("filename")["ko_num"].value_counts()
    ko_counts.name = "n"
    ko_matrix = (
        ko_counts.to_frame()
        .reset_index()
        .pivot(index="ko_num", columns="filename", values="n")
    )

    ko_matrix.to_csv(args.out_matrix, sep="\t", index=True, header=True)
    print(f"Wrote matrix ({ko_matrix.shape[0]:,} unique KO ids) to {args.out_matrix}")


if __name__ == "__main__":
    main()
