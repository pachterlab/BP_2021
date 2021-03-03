#!/bin/bash

CONFG="config_all.json"
if [[ $# -ge 1 ]]; then
    CONFIG=$1
fi

READSDIR=$(jq -r '.reads_dir' ${CONFIG})
OUTDIR=$(jq -r '.out_dir' ${CONFIG})
REFDIR=$(jq -r '.ref_dir' ${CONFIG})
CONFDIR=$(jq -r '.config_dir' ${CONFIG})
SCRIPTDIR=$(jq -r '.script_dir' ${CONFIG})
V2WL=$(jq -r '.v2_whitelist' ${CONFIG})
V3WL=$(jq -r '.v3_whitelist' ${CONFIG})

mkdir -p $OUTDIR/alevin_out
mkdir -p $OUTDIR/alevin_out_transcriptome

for config in $CONFDIR/*; do
    ref=$(jq -r '.ref' ${config})
    sample=$(jq -r '.sample' ${config})
    tech=$(jq -r '.technology' ${config})
    species=$(jq -r '.species' ${config})
    
    fastqs="$READSDIR/$species-$sample/"

    TECH="chromium"
    WL=$V2WL
    if [[ $tech = *3 ]]
    then
        TECH="chromiumV3"
        WL=$V3WL
    fi

    index="$REFDIR/$species-$ref/salmon_transcriptome/salmon_index"
    t2g="$REFDIR/$species-$ref/ref/salmon_t2g.txt"
    mkdir -p $OUTDIR/alevin_out_transcriptome/$species-$sample

    cmd="bash $SCRIPTDIR/salmon_alevin-fry.sh -o $OUTDIR/alevin_out_transcriptome/$species-$sample/ -i $index -w $WL -g $t2g -x $TECH -f $fastqs"

    echo $cmd
    eval $cmd
    
    index="$REFDIR/$species-$ref/salmon/salmon_index"
    t2g="$REFDIR/$species-$ref/salmon/decoy_t2g.txt"
    mkdir -p $OUTDIR/alevin_out/$species-$sample

    cmd="bash $SCRIPTDIR/salmon_alevin-fry.sh -o $OUTDIR/alevin_out/$species-$sample/ -i $index -w $WL -g $t2g -x $TECH -f $fastqs"

    echo $cmd
    eval $cmd
    
done
