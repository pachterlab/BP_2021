#!/bin/bash

CONFIG='config.json'
if [[ $# -ge 1 ]]; then
    CONFIG=$1
fi

REFDIR=$(jq -r '.ref_dir' ${CONFIG})
SCRIPTDIR=$(jq -r '.script_dir' ${CONFIG})

for species in 'arabidopsis-tair10' 'fly-dm6' 'human-grch38' 'human_mouse-hg19_mm10' 'mouse-mm10' 'rat-rnor6' 'worm-ws260' 'zebrafish-dr82'; do
    outdir=$REFDIR/$species
    mkdir -p $outdir/kallisto
    cmd="bash $SCRIPTDIR/mkref_kallisto.sh -o $outdir/kallisto -t $outdir/ref/transcriptome.fa"
    echo $cmd
    eval $cmd
done
