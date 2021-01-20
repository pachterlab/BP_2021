#!/bin/bash

usage () {
    echo "Usage: $0 [options]
    
    Options:
    -o, --output            output folder
    -g, --genome            genome (fasta)
    -a, --gtf               gtf file 
    -t, --transcriptome     transcriptome (fasta)
    -m, --t2g               transcripts to genes map
    "
    exit 1
}

while getopts ":o:g:a:t:m:" opt; do
    case $opt in
        o|--output)
            OUTDIR=$OPTARG
            ;;
        g|--genome)
            GENOME=$OPTARG
            ;;
        a|--gtf)
            GTF=$OPTARG
            ;;
        t|--transcriptome)
            TRANSCRIPTOME=$OPTARG
            ;;
        m|--t2g)
            T2G=$OPTARG
            ;;
        h)
            usage
            ;;
        \?)
            echo "Invalid argument"
            usage
            ;;
        :)
            echo "Add arguments"
            usage
            ;;
    esac
done

# check options        
if [ -z "$OUTDIR" -o -z "$GENOME" -o -z "$GTF" -o -z "$TRANSCRIPTOME" -o -z "$T2G" ]
then
    echo "Error"
    usage
fi

mkdir -p $OUTDIR

echo "[bash] making gentrome and decoys"

grep "^>" <(cat $GENOME) | cut -d " " -f 1 > $OUTDIR/decoys.txt
sed -i.bak -e 's/>//g' $OUTDIR/decoys.txt
cat $TRANSCRIPTOME $GENOME > $OUTDIR/gentrome.fa

echo "[salmon] building decoy index"
/usr/bin/time --output $OUTDIR/index.log -v \
salmon index -i $OUTDIR/salmon_index \
-d $OUTDIR/decoys.txt \
-t $OUTDIR/gentrome.fa \
-p 10

echo "[bash] making decoy t2g"
cat $T2G <(paste -d$'\t' $OUTDIR/decoys.txt $OUTDIR/decoys.txt) | cut -d$'\t' -f1,2 > $OUTDIR/decoy_t2g.txt

echo "Done."