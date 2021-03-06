#!/bin/bash

usage () {
    echo "Usage: $0 [options]
    
    Options:
    -o, --output            output folder
    -t, --transcriptome     transcriptome (fasta)
    "
    exit 1
}

while getopts ":o:g:a:t:m:" opt; do
    case $opt in
        o|--output)
            OUTDIR=$OPTARG
            ;;
        t|--transcriptome)
            TRANSCRIPTOME=$OPTARG
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
if [ -z "$OUTDIR" -o -z "$TRANSCRIPTOME" ]
then
    echo "Error"
    usage
fi

mkdir -p $OUTDIR

echo "[salmon] building transcriptome index"
/usr/bin/time --output $OUTDIR/index.log -v \
salmon index -i $OUTDIR/salmon_index \
-t $TRANSCRIPTOME \
-p 10

echo "Done."
