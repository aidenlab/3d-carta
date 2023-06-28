#!/bin/bash

USAGE="
*****************************************************
Usage: ./short-assemble-region-from-sam.sh [-h] [-c min_coverage] [-q min_mapq] [-r chr:start-end] path_to_sam_file(s)

ARGUMENTS:
path_to_sam_file        Path to sam file describing read alignment to sequence used as assisting reference

OPTIONS:
-h                      Shows this help
-c min_coverage         Minimum read coverage required for the reference fragment to participate in genome reconstruction [default is 3]
-q min_mapq             Minimum mapping quality of the read to be considered in the recontruction process [default is 1]
-r chr:start-end        Reconstruct fasta from a limited genomic region of a reference. The option can be an explicit listing of the chromosome with genomic start and end location in 1-base coordinate system or just the name of a chromosome for whole-chromosome reconstruction. This option is used for filtering out the sam file [default is empty for no filtering].

Uses map-sam-to-pos-array.awk and reduce-map-into-fasta.awk that should be in the same folder as this wrapper script.
*****************************************************
"

############### preliminary work ################
#reference="/Users/olga/WORK/CAMBRIDGE/Genetics/ASSEMBLY/Assisted_assembly/anonamos/panTro4.fa"
#panTro4=${reference}
# for QC
#hg38="/Users/olga/WORK/CAMBRIDGE/Genetics/ASSEMBLY/Assisted_assembly/anonamos/hg38.fa"
############### preliminary work ################


## SET DEFAULTS
min_coverage=1
min_mapq=0

## HANDLE OPTIONS
while getopts "c:q:r:h" opt; do
case $opt in
    h) echo "$USAGE" 
        exit 0
    ;;
    c) re='^[-0-9]+$'
        if ! [[ $OPTARG =~ $re ]] ; then
            echo "): Error: Wrong syntax for coverage threshold. Using default settings min_coverage=$min_coverage."
        else
            echo "(: -c flag was triggered. Coverage threshold is set to $OPTARG."
            min_coverage=$OPTARG
        fi
    ;;
    q)  re='^[0-9]+$'
        if ! [[ $OPTARG =~ $re ]] ; then
            echo "): Error: Wrong syntax for mapping quality. Using default settings min_mapq=$min_mapq."
        else
            echo "(: -q flag was triggered. Read mapping quality threshold is set to $OPTARG."
            min_mapq=$OPTARG
        fi
    ;;
    r)  re1='^[A-Za-z0-9\.]+$'; re2='^[A-Za-z0-9\.]+\:[0-9]+\-[0-9]+$'
        if ! [[ $OPTARG =~ $re1 ]] && ! [[ $OPTARG =~ $re2 ]] ; then
            echo "): Error: Wrong syntax for region option."
        else
            IFS=":-" read chr region_start region_end <<< "$OPTARG"
            if [ "$region_start" == "" ] || [ "$region_end" == "" ] || ! [ $region_end -gt $region_start ]; then
                echo "(: -r flag was triggered, but options for region start and end are not set or are ambiguous. Doing reconstruction of whole chr $chr."
                unset region_start
                unset region_end
            else
                echo "(: -r flag was triggered. Sam file will be prefiltered for $chr:$region_start-$region_end region."
            fi
        fi
    ;;
    *) echo "$USAGE"
        exit 1
    ;;
esac
done

cov_threshold=$((min_coverage-1)) # ge rather than gt criteria
mapq_threshold=$((min_mapq-1)) # ge rather than gt criteria

shift $(( OPTIND-1 ))

## HANDLE ARGUMENTS, TODO: check that sam file

if [ $# -lt 1 ]; then
    echo "): The set of arguments provided is incomplete. Please double-check your input, exiting!"
    echo "$USAGE" 
    exit 1
fi

## CHECK FOLDER FOR RELEVANT SCRIPTS
basepath=`cd "$( dirname $0)" && pwd`
echo "(: Checking that relevant scripts are present in $basepath folder."

map_script="$basepath""/map-sam-to-pos-array.awk"
reduce_script="$basepath""/reduce-map-into-fasta-stringent.awk"
flag_script="$basepath""/list-active-regions.awk"
extract_script="$basepath""/cut-reads-from-sam.awk"
reassemble_script="$basepath""/de-bruijn-assemble-DFS-both-strands.awk"
break_script="$basepath""/break-scaffold-into-contigs.awk"
local_fasta="$basepath""/cut-seq-from-fasta.awk"
map_gaps="$basepath""/make-gap-bed-for-local-ref.awk"

print_contig_N50="$basepath""/print-contig-n50.awk"
txt_to_bed_script="$basepath""/indel-txt-to-bed.awk"
prep_for_igv="$basepath""/prep-for-igv.sh"

awk -v region="$chr $region_start $region_end" -f ${extract_script} $@ | tee test | awk -v region="$chr $region_start $region_end" -v anchored="1" -f ${reassemble_script}

exit

#cmd="awk -v region=\"$chr $region_start $region_end\" -v anchored=\"1\" -f ${reassemble_script} test"
#echo $cmd
#eval $cmd
#exit

echo ">assembly" > reassembled_fragment.fasta
output=`awk -v region="$chr $region_start $region_end" -f ${extract_script} $@ | awk -v region="$chr $region_start $region_end" -v anchored="1" -f ${reassemble_script}`
read chr region_start length cov fasta message <<< "$output"
echo $fasta >> reassembled_fragment.fasta
echo $fasta
echo $message

################ supp for visualization ################
#bash ${prep_for_igv} reassembled_fragment.fasta ${hg38}
#bash ${prep_for_igv} reassembled_fragment.fasta ${panTro4}

################ supp for visualization ################
