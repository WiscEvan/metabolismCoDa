#!/usr/bin/env bash
#SBATCH --partition=queue
#SBATCH --time=UNLIMITED
#SBATCH -N 1 # Nodes
#SBATCH --ntasks-per-node=1
#SBATCH --error=%J.tabulate_kofamscan.err
#SBATCH --output=%J.tabulate_kofamscan.out

# ENV SETUP
# mamba create -c bioconda -n autometa autometa -y
# mamba activate autometa
##
# USER input
outdir="."
kofamscan_results_dir="${outdir}/kofamscan_results"
# KEGG database files
KO2CAZY="./databases/kegg/brite/ko2cazy.txt"
KO2COG="./databases/kegg/brite/ko2cog.txt"
KO2EC="./databases/kegg/brite/ko2ec.txt"
KO2GO="./databases/kegg/brite/ko2go.txt"
KO2RN="./databases/kegg/brite/ko2rn.txt"
KO2TC="./databases/kegg/brite/ko2tc.txt"
PATHWAYS_DIR="./databases/kegg/pathways"
BRITE_KO_DB="./databases/kegg/brite/ko00001.txt"
BRITE_KMODULE_DB="./databases/kegg/brite/ko00002.txt"
# END User input
# ---------------

# BEGIN
tabulate_outdir="${outdir}/processed"
OUT_TABLE="${tabulate_outdir}/kofamscan_results_table.tsv"
OUT_MATRIX="${tabulate_outdir}/kofamscan_results_matrix.tsv"

if [ ! -d $tabulate_outdir ];then
    mkdir -p $tabulate_outdir
fi

# Alternatively:
# kegg_decoder_orgs=($(cut -d" " -f1 kegg_decoder_orgs.txt))
# where kegg_decoder_orgs has filenames in first column
# Merge kofamscan results and add metadata from kegg dbs
python ./src/merge_kegg_tables.py \
    --dirpath $kofamscan_results_dir \
    --ko2cazy $KO2CAZY \
    --ko2cog $KO2COG \
    --ko2ec $KO2EC \
    --ko2go $KO2GO \
    --ko2rn $KO2RN \
    --ko2tc $KO2TC \
    --out-table $OUT_TABLE \
    --out-matrix $OUT_MATRIX

# Now tabulate completeness of pathways according to KOs
python ./src/tabulate_kegg_pathways.py \
    --kofamscan-matrix $OUT_MATRIX \
    --outdir $tabulate_outdir \
    --brite-ko1 $BRITE_KO_DB \
    --brite-ko2 $BRITE_KMODULE_DB \
    --pathways $PATHWAYS_DIR

python ./src/embed_results.py \
    --input $OUT_MATRIX \
    --output $tabulate_outdir

