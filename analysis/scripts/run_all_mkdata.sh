#!/bin/bash

CONFIG='config.json'
if [[ $# -ge 1 ]]; then
    CONFIG=$1
fi

REFDIR=$(jq -r '.ref_dir' ${CONFIG})
SCRIPTDIR=$(jq -r '.script_dir' ${CONFIG})

CONFDIR=$(jq -r '.config_dir' ${CONFIG})
OUTDIR=$(jq -r '.out_dir' ${CONFIG})

for config in $CONFDIR/*.json; do
    sample=$(jq -r '.sample' ${config})
    outdir=$OUTDIR/plotting/$sample
    mkdir -p $outdir
    cmd="python $SCRIPTDIR/mkdata.py -o $outdir -d $sample"
    echo $cmd
done
