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

mkdir -p $OUTDIR/kallisto_out

for config in $CONFDIR/*; do
    ref=$(jq -r '.ref' ${config})
    sample=$(jq -r '.sample' ${config})
    tech=$(jq -r '.technology' ${config})
    species=$(jq -r '.species' ${config})
    
    kb_index="$REFDIR/$species-$ref/kallisto/index.idx"
    t2g="$REFDIR/$species-$ref/ref/t2g.txt"
    fastqs="$READSDIR/$species-$sample"

    TECH="10xv2"
    WL=$V2WL
    if [[ $tech = *3 ]]
    then
        TECH="10xv3"
        WL=$V3WL
    fi

    mkdir -p $OUTDIR/kallisto_out/$species-$sample

    cmd="bash $SCRIPTDIR/kallisto_bus.sh -o $OUTDIR/kallisto_out/$species-$sample/ -i $kb_index -w $WL -g $t2g -x $TECH -f $fastqs"

    echo $cmd
    eval $cmd
done
