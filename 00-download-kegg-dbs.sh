#!/bin/bash

# Retrieve kegg pathway maps with their respective KO identifiers.
# BEGIN USER INPUT
OUTDIR="databases"
# END USER INPUT
# ----------------------


dbdir="${OUTDIR}/kegg"
pathway_list="${dbdir}/pathway_list.txt"
pathway_dir="${dbdir}/pathways"
mkdir -p $pathway_dir
brite_list="${dbdir}/brite_list.txt"
brite_xl_list="${dbdir}/brite_xl_list.txt"
brite_dir="${dbdir}/brite"
mkdir -p $brite_dir

if [ ! -f $brite_list ];then
    wget --quiet https://rest.kegg.jp/list/brite -O $brite_list
fi
if [ ! -f $brite_xl_list ];then
    wget --quiet https://rest.kegg.jp/list/brite/xl -O $brite_xl_list
fi
if [ ! -f $pathway_list ];then
    echo "Downloading pathway list -> ${pathway_list}"
    wget --quiet http://rest.kegg.jp/list/pathway -O $pathway_list
    # easy read in pandas
    # df = pd.read_table("databases/kofam/kegg_pathways/kegg_pathway_list.txt", header=None)
fi
if [ ! -f "${dbdir}/kegg_db_info.txt" ];then
    wget --quiet https://rest.kegg.jp/info/kegg -O "${dbdir}/kegg_db_info.txt"
fi
if [ ! -f "${dbdir}/kegg_pathway_db_info.txt" ];then
    wget --quiet https://rest.kegg.jp/info/pathway -O "${dbdir}/kegg_pathway_db_info.txt"
fi
if [ ! -f "${dbdir}/kegg_organism_map.txt" ];then
    # df = pd.read_table("databases/kofam/kegg_pathways/kegg_organism_map.tsv", header=None)
    wget --quiet https://rest.kegg.jp/list/organism -O "${dbdir}/kegg_organism_map.tsv"
fi

# NOTE the organism map may be used to retrieve the gene entries for a specific reference organism
# i.e. homosapien entry id is T01001 and the abbreviated id is hsa:
# so the list API urls would be https://rest.kegg.jp/list/T01001 or https://rest.kegg.jp/list/hsa
# Can even get a subset pathway map https://www.kegg.jp/pathway/hsa00010
# You may also retrieve a list of pathways from the organism:
# i.e. https://rest.kegg.jp/list/pathway/hsa
for pathway in $(cut -f1 $pathway_list); do
    pathway_map_id="${pathway/path:/}"
    pathway_map_fpath="${pathway_dir}/${pathway_map_id}.txt"
    pathway_map_img_fpath="${pathway_dir}/${pathway_map_id}.png"
    if [ ! -f ${pathway_map_fpath} ];then
        echo "getting ${pathway_map_id} KOs"
        wget --quiet http://rest.kegg.jp/link/ko/$pathway -O $pathway_map_fpath
    fi
    # Will retrieve the associated pathway map image (2x size)
    if [ ! -f ${pathway_map_img_fpath} ];then
        echo "getting ${pathway_map_id} image"
        wget --quiet https://rest.kegg.jp/get/${pathway_map_id}/image2x -O $pathway_map_img_fpath
    fi
done

# NOTE: Can linke to associated EC info with url: https://www.kegg.jp/entry/ec:{ec_number}
# i.e. https://www.kegg.jp/entry/ec:2.7.10.1

for brite in $(cut -f1 $brite_list);do
    brite_id=${brite/br:/}
    brite_fpath="${brite_dir}/${brite_id}.txt"
    if [ ! -f $brite_fpath ];then
        echo "getting ${brite_id} file"
        wget --quiet https://rest.kegg.jp/get/br:${brite_id} -O $brite_fpath
    fi
done

for brite in $(cut -f1 $brite_xl_list);do
    brite_id=${brite/br:/}
    brite_fpath="${brite_dir}/${brite_id}.txt"
    if [ ! -f $brite_fpath ];then
        echo "getting ${brite_id} file"
        wget --quiet https://rest.kegg.jp/get/br:${brite_id} -O $brite_fpath
    fi
done
