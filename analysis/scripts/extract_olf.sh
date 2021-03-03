#!/bin/bash

OUTDIR="../../data/pbmc_nose/extract/"
T2G="../../reference/human-grch38/ref/t2g.txt"
EC="../../data/pbmc_nose/matrix.ec"
TXN="../../data/pbmc_nose/transcripts.txt"
BUS="../../data/pbmc_nose/sc.bus"
CAP="../../reference/human-grch38/olf_OR51L1.txt"
FASTQDIR="../../data/fastqs/human-pbmc10k_v3"
mkdir -p $OUTDIR

# bustools capture -o $OUTDIR/olf.bus -c $CAP -e $EC -t $TXN -s $BUS
bustools sort -o $OUTDIR/solf.bus -t 8 --flags $OUTDIR/olf.bus
bustools extract -o $OUTDIR/ -f $(ls $FASTQDIR/*R* | tr '\n' ',') -N 2 $OUTDIR/solf.bus
bustools text -fp $OUTDIR/solf.bus | cut -d$'\t' -f3 > $OUTDIR/olf_ec_from_bus.txt
cat $OUTDIR/2.fastq.gz $OUTDIR/4.fastq.gz > $OUTDIR/2_4.fastq.gz
