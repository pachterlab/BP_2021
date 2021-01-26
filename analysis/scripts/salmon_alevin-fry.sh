#!/bin/bash

usage () {
    echo "Usage: $0 [options]
    
    Options:
    -o, --output            output folder
    -i, --index             pseudoalignment index directory
    -w, --whitelist         10x barcode whitelist
    -g, --genemap           transcripts to genes map
    -x, --technology        single-cell tech 'chromium' or 'chromiumV3'
    -f, --fastqdir          folder containing fastqs
    "
    exit 1
}

while getopts ":o:i:w:g:x:f:" opt; do
    case $opt in
        o|--output)
            OUTDIR=$OPTARG
            ;;
        i|--index)
            INDEX=$OPTARG
            ;;
        w|--whitelist)
            WHITELIST=$OPTARG
            ;;
        g|--genemap)
            T2G=$OPTARG
            ;;
        x|--technology)
            TECH=$OPTARG
            ;;
        f|--fastqdir)
            FASTQDIR=$OPTARG
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
if [ -z "$OUTDIR" -o -z "$INDEX" -o -z "$WHITELIST" -o -z "$T2G" -o -z "$TECH" -o -z "$FASTQDIR" ]
then
    echo "Error"
    usage
fi

# begin workflow
mkdir -p $OUTDIR

echo ""
echo '[salmon] pseudoaligning reads..'

/usr/bin/time --output $OUTDIR/pseudoalignment.log -v \
salmon alevin \
-l ISR \
-i $INDEX \
-1 $(ls $FASTQDIR | awk -v p=$FASTQDIR '{print p$0}' | grep R1) \
-2 $(ls $FASTQDIR | awk -v p=$FASTQDIR '{print p$0}' | grep R2) \
--tgMap $T2G \
-o $OUTDIR \
--rad \
-p 10 \
--$TECH # Tech goes here
  

echo ""
echo '[alevin-fry] correcting barcodes..'

GEN="$OUTDIR/generate"

/usr/bin/time --output $OUTDIR/correct.log -v \
alevin-fry generate-permit-list \
-d either \
-i $OUTDIR \
-o $GEN \
-b $WHITELIST

echo ""
echo '[alevin-fry] generating whitelist..'

WL="$OUTDIR/whitelist"

/usr/bin/time --output $OUTDIR/whitelist.log -v \
alevin-fry generate-permit-list \
-d either \
-i $OUTDIR \
-o $WL \
--knee-distance
   
echo ""
echo '[alevin-fry] sorting rad file..'


/usr/bin/time --output $OUTDIR/sort.log -v \
alevin-fry collate \
-i $GEN \
-r $OUTDIR \
-t 10
   
echo ""
echo '[alevin-fry] counting umis..'


QNT="$OUTDIR/quant"

/usr/bin/time --output $OUTDIR/count.log -v \
alevin-fry quant \
--use-mtx \
-i $GEN \
-o $QNT \
-m $T2G \
-t 10

echo ""
echo "[alevin-fry] converting to text"
/usr/bin/time --output $OUTDIR/text.log -v \
alevin-fry view \
--rad $OUTDIR/map.rad \
--header \
> /dev/null

echo "Done."