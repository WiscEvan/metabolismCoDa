#!/usr/bin/env bash
#SBATCH --partition=queue
#SBATCH --time=UNLIMITED
#SBATCH -N 1 # Nodes
#SBATCH --ntasks-per-node=8
#SBATCH --error=%J.kofamscan.err
#SBATCH --output=%J.kofamscan.out

# SETUP
# compute environment...
# mamba create -n kofamscan -c bioconda kofamscan pandas 
# mamba activate kofamscan
#  ... and databases...
# 00-download-kegg-dbs.sh
# 00-download-kofamscan-dbs.sh
# ... (optional) download reference genomes from NCBI ...
# 00-download-ncbi-reference-genomes.sh
# 

input_dir="./ncbi" # should contain MAG amino-acid fasta files
outdir="."
cpu=4

kofamscan_results_dir="${outdir}/kofamscan_results"
# These files are downloaded using 00-download-kofamscan-dbs.sh
profiles="metabolismCoDa/databases/kofamscan/profiles"
ko_list="metabolismCoDa/databases/kofamscan/ko_list"
# create outdirs
tmp="${kofamscan_results_dir}/tmp"
mkdir -p $tmp

# Retrieve all genomes' proteins
all_faas=(`find $input_dir -name "protein.faa"`)

for query in ${all_faas[@]};do
    if [ $(basename $query) == "protein.faa" ];then
        # When using 00-download-ncbi-reference-genomes.sh
        # Filepaths resemble: 
        # /home/user/metabolismCoDa/ncbi/test/data/GCF_001922985.1/protein.faa
        accession=$(basename $(dirname $query))
    else
        # e.g. <your/input/dir/MAG_name.faa>
        accession=$(basename ${query/.faa/})
    fi
    outfname="${accession}.kofamscan.tsv"
    tmp_dir="${tmp}/${accession}_tmp"
    out="${kofamscan_results_dir}/${outfname}"
    if [ -f $out ];then
        echo "$out already exists, skipping..."
        continue
    fi
    echo "annotating ${accession}"
    # exec_annotation \
    srun --mincpus $cpu --job-name="kofamscan ${accession}" exec_annotation \
        $query \
        -o $out \
        --profile $profiles \
        --ko-list $ko_list \
        --tmp-dir $tmp_dir \
        --cpu $cpu \
        --format detail-tsv &

done

wait
