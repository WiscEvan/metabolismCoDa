#!/usr/bin/env bash
## Downloads databases to use kofamscan to annotate sequences with KO identifiers
# After mapping, the pathways and brite hierarchies may be inferred from the other KEGG databases either online
# or downloaded using retrieve_kegg_dbs.sh

OUTDIR="databases"


dbdir="${OUTDIR}/kofamscan"
if [ ! -d $dbdir ];then
	mkdir -p $dbdir
fi

cd $dbdir

if [ ! -f ko_list ];then
	wget ftp://ftp.genome.jp/pub/db/kofam/ko_list.gz
	gzip -d ko_list.gz
fi

if [ ! -f profiles.tar.gz ];then
	wget ftp://ftp.genome.jp/pub/db/kofam/profiles.tar.gz
	tar xf profiles.tar.gz
fi

if [ ! -f README ];then
	wget ftp://ftp.genome.jp/pub/db/kofam/README
fi

