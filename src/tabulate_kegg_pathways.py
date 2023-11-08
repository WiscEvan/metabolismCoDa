#!/usr/bin/env python

import os
import re
import argparse
import glob
from typing import Dict, List, Set
import pandas as pd

def get_pathway_maps(line):
    pattern = re.compile(r"\[PATH:(.*)\]")
    maps = pattern.findall(line)
    if maps:
        maps = maps[0].split()
    else:
        maps = []
    return maps

def read_brite_ko(fpath: str)->pd.DataFrame:
    data = []
    with open(fpath) as f:
        for line in f:
            if line.startswith("!"):
                continue
            elif line.startswith("#"):
                continue
            elif line.startswith("A"):
                acode = line.strip().split()[0].replace("A", "")
                adesc = " ".join(line.strip().split(" ")[1:])
            elif line.startswith("B"):
                # 'B  09194 Poorly characterized\n'
                bcode = line.strip().split()[1]
                bdesc = " ".join(line.strip().split()[2:])
            elif line.startswith("C"):
                # 'C    99992 Structural proteins\n'
                ccode = line.strip().split()[1]
                cdesc = " ".join(line.strip().split()[2:])
            elif line.startswith("D"):
                __,ko_num,ko_acronym_defition = line.strip().split(maxsplit=2)
                # ['D', 'K01810', 'GPI, pgi; glucose-6-phosphate isomerase [EC:5.3.1.9]']
                data.append({
                    "acode": acode,
                    "adesc": adesc,
                    "bcode": bcode,
                    "bdesc": bdesc,
                    "ccode": ccode,
                    "cdesc": cdesc,
                    "KO": ko_num,
                    "definition": ko_acronym_defition,
                })
    df = pd.DataFrame(data).set_index("KO")
    return df


def read_brite_komodule(fpath: str)->pd.DataFrame:
    # line = "D      M00309  Non-phosphorylative Entner-Doudoroff pathway, gluconate/galactonate => glycerate [PATH:map00030 map00052 map01200 map01100 map01120]"
    data = []
    with open(fpath) as f:
        for line in f:
            if line.startswith("!"):
                continue
            elif line.startswith("#"):
                continue
            elif line.startswith("A"):
                adesc = line.strip().replace("A", "", 1)
            elif line.startswith("B"):
                # 'B  09194 Poorly characterized\n'
                bdesc = " ".join(line.strip().split()[1:])
            elif line.startswith("C"):
                # 'C    99992 Structural proteins\n'
                cdesc = " ".join(line.strip().split()[1:])
            elif line.startswith("D"):
                __,module,module_definition = line.strip().split(maxsplit=2)
                pathway_maps = get_pathway_maps(line)
                # D      M00309  Non-phosphorylative Entner-Doudoroff pathway, gluconate/galactonate => glycerate [PATH:map00030 map00052 map01200 map01100 map01120]
                # [PATH:map00030 map00052 map01200 map01100 map01120]
                # Need to write regex to get pathways
                data.append({
                    "A": adesc,
                    "B": bdesc,
                    "C": cdesc,
                    "module": module,
                    "definition": module_definition,
                    "pathway_maps":pathway_maps,
                })
    df = pd.DataFrame(data).set_index("module")
    return df


def clean_string(input_string):
    match = re.search(r'\[(BR|PATH):([^\]]+)\]', input_string)
    if match:
        id_type, id_value = match.group(1), match.group(2)
        description = re.sub(r'[\/,:*?"()<>|\[\]]', '', input_string)
        description = description.replace(f"{id_type}{id_value}", '').strip()
        description = description.replace(" ", "_").replace("_-_", "-")
        return f"{id_value}_{description}"
    else:
        return re.sub(r'[\/,:*?"<>|()\[\]]', '', input_string).replace(" ", "_").replace("_-_", "-")


def get_ko_category_totals(dff: pd.DataFrame, brite_df: pd.DataFrame, outdir: str)->pd.DataFrame:
    brite_ko_indices = brite_df.groupby(["adesc", "bdesc", "cdesc"]).indices
    dfs = []
    for a_b_c,c_brite_index in brite_ko_indices.items():
        A_category,B_category,C_category = a_b_c
        brite_c_df = brite_df.iloc[c_brite_index]
        brite_c_ko_index = brite_c_df.index
        c_df = dff.loc[dff.index.isin(brite_c_ko_index)]
        if c_df.empty:
            continue
        clean_A = clean_string(A_category)
        clean_B = clean_string(B_category)
        out_dirpath = os.path.join(outdir, clean_A, clean_B)
        clean_C = clean_string(C_category)
        out_filename = f"{clean_C}.tsv"
        out_filepath = os.path.join(out_dirpath, out_filename)
        os.makedirs(out_dirpath, exist_ok=True)
        out_df = pd.merge(c_df, brite_df, how='left', left_index=True, right_index=True)
        if not os.path.exists(out_filepath):
            out_df.index.name = "ko_num"
            out_df.to_csv(out_filepath, sep='\t', index=True, header=True)
        total_unique = brite_c_df.index.nunique()
        c_nunique_counts = c_df.ge(1).sum()
        c_total_df = c_nunique_counts.to_frame(name='n')
        c_total_df["A"] = A_category
        c_total_df["B"] = B_category
        c_total_df["C"] = C_category
        c_total_df["KO unique count"] = total_unique
        dfs.append(c_total_df)
    ko_totals = pd.concat(dfs).sort_values(["KO unique count"], ascending=False)
    ko_totals = ko_totals.assign(percent_complete=lambda row: row["n"] / row["KO unique count"] * 100)
    ko_totals.index.name = 'filename'
    ko_totals = ko_totals.reset_index().set_index(["A", "B", "C"])
    return ko_totals

def get_ko_set_from_pathway_maps(fpaths: List[str])->Set[str]:
    ko_set = set()
    for fpath in fpaths:
        with open(fpath) as fh:
            for line in fh:
                ko = line.strip().split('\t')[-1].replace("ko:", "")
                ko_set.add(ko)
    return ko_set

def get_kmod_totals(dff: pd.DataFrame, kmod_df: pd.DataFrame, pathway_maps: Dict[str,str])->pd.DataFrame:
    module_dfs = []
    for module, map_ids in kmod_df.groupby("module"):
        map_fpaths = [pathway_maps.get(map_id) for map_id in map_ids.pathway_maps.explode().tolist() if pathway_maps.get(map_id)]
        if not map_fpaths:
            # print(f"{module} does not contain pathway maps, skipping...")
            continue
        try:
            ko_set = get_ko_set_from_pathway_maps(map_fpaths)
        except TypeError:
            print(f"module:{module}, {map_fpaths}")
            continue
        # Could go one by one and compute totals by pathway map?
        module_unique_count = len(ko_set)
        pathway_df = dff.loc[dff.index.isin(ko_set)]
        pathway_nunique_counts = pathway_df.ge(1).sum()
        pathway_total_df = pathway_nunique_counts.to_frame(name='n')
        pathway_total_df['module unique KO count'] = module_unique_count
        pathway_total_df = pathway_total_df.assign(percent_complete=lambda row: row["n"] / module_unique_count * 100)
        pathway_total_df['module'] = module
        module_dfs.append(pathway_total_df)
    total_df = pd.concat(module_dfs)
    total_df = pd.merge(total_df, kmod_df, how='left', left_on='module', right_index=True)
    total_df.index.name = 'filename'
    return total_df


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--kofamscan-matrix", help="Path to kofamscan results matrix with ko_num as index and samples as columns, values are counts", required=True)
    parser.add_argument("--outdir", help="Output directory path to write tables", required=True)
    parser.add_argument("--brite-ko1", help="Path to brite ko00001.txt file", default="/media/BRIANDATA2/databases/kegg/brite/ko00001.txt")
    parser.add_argument("--brite-ko2", help="Path to brite ko00001.txt file", default="/media/BRIANDATA2/databases/kegg/brite/ko00001.txt")
    parser.add_argument("--pathways", help="Directory path to pathway maps", default="/media/BRIANDATA2/databases/kegg/pathways")
    args = parser.parse_args()
    # ko_file = "/media/BRIANDATA2/databases/kegg/brite/ko00001.txt"
    # brite_df = read_brite_ko(ko_file)
    brite_ko_df = read_brite_ko(args.brite_ko1)
    # brite_ko_df['pathway_map'] = brite_ko_df.cdesc.map(lambda x: get_pathway_maps(x)).explode()
    brite_kmod_df = read_brite_komodule(args.brite_ko2)
    # To unroll pathway_maps list to each have its own row w.r.t. module
    # brite_kmodule_df.explode("pathway_maps")
    pathway_maps = {
        os.path.basename(fp).replace(".txt", ""):fp
        for fp in glob.glob(os.path.join(args.pathways, "map*.txt"))
    }

    # for map_id, fpath in pathway_maps.items():
    #     # if map_id in maps
    #     ko_set = get_ko_set_from_pathway_map(fpath)
    

    # Read in counts matrix and remove metadata cols
    df = pd.read_table(args.kofamscan_matrix, index_col='ko_num')
    # df = pd.read_table("kofamscan_results/kofamscan_results_matrix.tsv", index_col='ko_num')
    sample_cols = [col for col in df.columns if ".tsv" in col]
    dff = df[sample_cols].copy()
    
    kmod_totals = get_kmod_totals(dff, brite_kmod_df, pathway_maps)
    # median completeness by module
    orgs = [
        "Gammaproteobacterium_IMCC1989.kofamscan.tsv",
        "AB1_endobugula_sertula.kofamscan.tsv",
        "UBA1124_sp002313265.kofamscan.tsv",
        "pacifica_merged_3.kofamscan.tsv",
        "simplex_merged_1.kofamscan.tsv",
        "Vibrio_parahaemolyticus.kofamscan.tsv",
    ]
    kmod_totals_out = os.path.join(args.outdir, "kegg_module_totals.tsv")
    kmod_totals.to_csv(kmod_totals_out, sep='\t', index=True, header=True)
    print(f"wrote kegg module counts to {kmod_totals_out}")
    # NOTE: orgs other than Ca. E. sp. provided above 
    # were identified from 03-kegg-exploratory.Rmd
    kmod_medians = kmod_totals.loc[orgs].groupby("module")['percent_complete'].median()
    # Filter any modules that are 0% complete
    kmod_medians = kmod_medians.loc[kmod_medians.gt(1)].sort_values(ascending=False)
    # Add metadata annotations back in...
    kmod_medians = pd.merge(kmod_medians, brite_kmod_df, how='left', left_index=True, right_index=True)
    kmod_medians_out = os.path.join(args.outdir, "kegg_module_median_percent_complete.tsv")
    kmod_medians.to_csv(kmod_medians_out, sep='\t', index=True, header=True)
    print(f"wrote kegg module median completeness to {kmod_medians_out}")

    # Get ko category totals
    ko_totals = get_ko_category_totals(dff, brite_ko_df, args.outdir)
    medians = ko_totals.groupby(["A","B","C"])["percent_complete"].median().to_frame(name='median_percent_completeness').copy()
    ko_medians_out = os.path.join(args.outdir, "kegg_orthology_median_percent_complete.tsv")
    medians.to_csv(ko_medians_out, sep='\t', index=True, header=True)
    print(f"wrote kegg orthology median completeness to {ko_medians_out}")
    ko_totals_out = os.path.join(args.outdir, "kegg_orthology_totals.tsv")
    ko_totals.to_csv(ko_totals_out, sep='\t', index=True, header=True)
    print(f"wrote kegg orthology counts to {ko_totals_out}")

if __name__ == "__main__":
    main()
