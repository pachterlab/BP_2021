#!/bin/bash

CONFIG='config.json'
if [[ $# -ge 1 ]]; then
    CONFIG=$1
fi

REFDIR=$(jq -r '.ref_dir' ${CONFIG})
SCRIPTDIR=$(jq -r '.script_dir' ${CONFIG})

for species in 'arabidopsis-tair10' 'fly-dm6' 'human-grch38' 'human_mouse-hg19_mm10' 'mouse-mm10' 'rat-rnor6' 'worm-ws260' 'zebrafish-dr82'; do
    outdir=$REFDIR/$species
    genome=$REFDIR/$species/ref/genome.fa
    transcriptome=$REFDIR/$species/ref/transcriptome.fa
    gtf=$REFDIR/$species/ref/genes.gtf
    t2g=$REFDIR/$species/ref/t2g.txt
    
    mkdir -p $outdir/salmon
    cmd="bash $SCRIPTDIR/mkref_salmon.sh -o $outdir/salmon -g $genome -a $gtf -t $transcriptome -m $t2g"
    echo $cmd
    eval $cmd

    mkdir -p $outdir/salmon_transcriptome
    cmd="bash $SCRIPTDIR/mkref_salmon_transcriptome.sh -o $outdir/salmon_transcriptome -t $transcriptome"
    echo $cmd
    eval $cmd

    cut -d$'\t' -f1,2 $t2g > $REFDIR/$species/ref/salmon_t2g.txt; done
done
