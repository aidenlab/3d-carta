#!/bin/bash -x

USAGE="
*****************************************************
This is a wrapper for a the Anonamos comparative genome assembler which analyzes the alignment data of the dna sequences derived from species A to a genome of closely related species B to create ref-assisted assembly for species A. (This is a short workflow version of Anonamos with no provisional genome creation and a single map-reduce procedure.) The input to the assembler are single-end SAM alignment files, typically ones produced as part of the Juicer workflow. (Note that performance might increase if the files are preprocessed: split at chimeric junctions, dedupped, filtered for overlaps etc.). Pipeline includes the following steps:

	(0) In the preliminary step the work space is organized and the B reference is preprocessed to create the map of nucleotides.
	(1) During the next step (-S sort) the sam file(s) are merge-sorted into a master sorted SAM file. The master SAM is then split into manageable chunks. The split size is adjusted dynamically to result in a predefined number of files determined by jobcount parameter (unprompted).
	(3) In the mapping step (-S map) a per base nucleotide map is constructed from split sam data. Indels and read breaks positions with respect to the referene are also recorded.
	(4) The results from individual chunks are merged during the merge step (-S merge). Reference map data is also added during this step. [[Right now also the clip data is filtered here (homozygous restriction-site associated only), but this might be reconsidered.]]
	(5) The next step (-S reduce) is preliminary consolidation of map data.
	(6) Indel and clip features that were not incorporated during preliminary data consolidation are examined and candidate regions for reassembly are listed during reassembly stage (-S reassemble). Read subsets pertaining to candidate regions are extracted from fasta and local de Bruijn graph tracing is performed. This step is skipped in the current version (180326 from March 26, 2018).
	(7) In the last step (-S finalize) the final fasta is constructed. By default the output will contain as many sequences as the reference fasta, i.e. there will be one line of output sequence per every line of input fasta sequence. If however a -a|--assembly flag is triggered, the actions encoded by the assembly file for the B reference passed by the user will be applied to produce the final fasta.

Version date: March 26, 2018.

Usage: ./run-anonamos-assembler.sh [options] <path_to_reference_fasta> <path_to_sam_file(s)>

ARGUMENTS:

path_to_reference				Path to an assisting species reference fasta.
path_to_sam_file(s)        		Path to sam file(s) containing alignments of data generated using a targe species sample to an assisting species reference.

MAIN OPTIONS:

-h|--help                      	
								Shows this help.
-a|--assembly		path_to_assembly_file	
								Path to assembly file (typically from 3D-DNA and JBAT), describing a set of large-scale assembly actions (cut, anchor, order and orient) that are required to create a genome for species of interest from the assisting reference.
-c|--consensus		consensus_threshold 
								The number of bases in the alignment data accepted to call a basepair in the final sequence, i.e. how many times the base needs to be seen across the read set (and reference if c>1) to constitute a valid call. Default: 1.
-q|--mapq			mapq_threshold  
								Minimum mapping quality of the read to be considered in the recontruction process. Default: 0.
-r|--restriction	restriciton_sequence
								Restriction sequence for the enzyme. Expected if using Hi-C read alignments as input. Default: GATC.
								

SUPPLEMENTARY OPTIONS:

** worflow **

-t|--threads		number_of_threads
								Specify the number of threads to run.
-s|--stage			stage
								Start from a particular stage. Can be sort, map, merge, reduce, reassemble and finalize.
--single-stage
								Run only one stage, whichever is passed on with -s|--stage (if -s|--stage is not triggered only preliminary stage is run).

** tuning **

--trust-reference				
								Consider reference as contributing to consensus [applicable only for c>1, if c=1 trust_reference is ignored].

*****************************************************
"

pipeline=`cd "$( dirname $0)" && pwd`


###################### SET DEFAULTS #############################

shopt -s extglob # enable pathname expansions

# problem-specific

consensus_threshold=1 			# trust single reads
ignore_reference=1				# do not use reference in consensus calling
mapq_threshold=0 				# include mapq0 reads
restriction_sequence="GATC"		# assume MboI Hi-C library

skip_prep=false
skip_sort=false
skip_map=false
skip_merge=false
skip_reduce=false
skip_reassembly=true	# this step is disabled
skip_finalize=false
single_stage=false

jobcount=`grep -c ^processor /proc/cpuinfo`
jobcount=$((jobcount*80/100))	# use 80% of available threads

# organizational

topDir=$(pwd)

debugDir=${topDir}"/debug"
tmpDir=${topDir}"/tmp"
samDir=${topDir}"/sam"
splitDir=${topDir}"/split"
mapDir=${topDir}"/map"
reduceDir=${topDir}"/reduce"

## TODO: I should probably get rid of log file.
logfile=$debugDir/log.txt

safegetfullpath(){
echo `cd $(dirname $1) && pwd -P`"/"`basename $1`
}
export -f safegetfullpath


################### WRITE COMMAND #################

optandarg="$*"

cmd=$(safegetfullpath $0) 
printf "\n-------------------------\n%s %s\n" "$cmd" "$optandarg" | tee /dev/stderr


###################### HANDLE OPTIONS ###############################

while :; do
	case $1 in
	
############ MAIN ############

		-h|--help)
			echo "$USAGE" >&1
			exit 0
			shift
        ;;
        -a|--assembly) OPTARG=$2	# TODO: check extension to match .assembly
			[ -s $OPTARG ] && echo >&1 ":) -a|--assembly flag was triggered. Will apply actions encoded in ${OPTARG} to fasta." || (echo >&2 ":( File $OPTARG not found. Exiting!" && exit)
			assembly=$OPTARG
			shift
		;;
        -c|--consensus) OPTARG=$2
        	re='^[-0-9]+$'
			if ! [[ $OPTARG =~ $re ]] || [ "$OPTARG" -eq 0 ]; then
				echo ":( Error: Wrong syntax for coverage threshold. Using default settings consensus_threshold=${consensus_threshold}." >&2
			else
				echo ":) -c flag was triggered. Coverage threshold is set to $OPTARG." >&1
				consensus_threshold=$OPTARG
			fi
			shift
		;;
		-q|--mapq) OPTARG=$2
			re='^[0-9]+$'
			if ! [[ $OPTARG =~ $re ]] ; then
				echo ":( Error: Wrong syntax for mapping quality. Using default settings mapq_threshold=${mapq_threshold}." >&2
			else
				echo ":) -q|--mapq flag was triggered. Read mapping quality threshold is set to $OPTARG." >&1
				mapq_threshold=$OPTARG
			fi
			shift
		;;
		-r|--restriction) OPTARG=$2
			re='^[ACGTNX]+$'
			if ! [[ $OPTARG =~ $re ]] ; then
				echo ":( Error: Wrong syntax for restriction sequence. Exiting." >&2
				exit 1
			else
				echo ":) -r|--restriction flag was triggered. Restriction sequence is set to $OPTARG." >&1
				restriction_sequence=$OPTARG
			fi
			shift
		;;
		
############ STAGING ############

		-s|--stage) OPTARG=$2
			if [ "$OPTARG" == "sort" ] || [ "$OPTARG" == "map" ] || [ "$OPTARG" == "merge" ] || [ "$OPTARG" == "reduce" ] || [ "$OPTARG" == "reassemble" ] || [ "$OPTARG" == "finalize" ]; then
				echo >&1 ":) -s|--stage flag was triggered. Fast-forwarding to $OPTARG stage."
				stage=$OPTARG
			else
				echo >&2 ":( Unnown option for -s|--stage parameter. Can be only sort, map, merge, reduce, and reassemble. Exiting!"
				exit 1
			fi
			shift
		;;
		--single-stage)
			echo ": --single-stage flag was triggered. Will execute only a single pipeline step prompted by --stage or the preparatory stage if no explicit --stage option was listed." >&1
			single_stage=true
		;;

############ MISC ############

		-t|--threads) OPTARG=$2
			jobcount=$OPTARG
			shift
		;;

		--trust-reference)
			(echo ":) --use-reference flag was triggered. Reference will be taken into account when making consensus calls.") >&1
			ignore_reference=0
		;;

		--) # End of all options
			shift
			break
		;;
		-?*)
			echo ":| WARNING: Unknown option. Ignoring: ${1}" >&2
		;;
		*) # Default case: If no more options then break out of the loop.
			break
	esac
	shift
done

############### DO THINGS BASED ON OPTIONS ###############

if [ ! -z "$stage" ]; then
	case $stage in
		finalize)	skipt_reassemble=true
		;&
		reassemble)	skip_reduce=true
		;&
		reduce)		skip_merge=true
		;&
		merge)		skip_map=true
		;&
		map)		skip_sort=true
		;&
		sort)		skip_prep=true
		;;
		*)			echo "$USAGE"
					exit 1
	esac
fi

# check parameter compatibility: force to ignore reference if consensus is set to 1
([ "$consensus_threshold" -eq 1 ]  && [ "$ignore_reference" -eq 0 ]) && (ignore_reference=1; echo ":| WARNING: Incompatible options. Will ignore --trust-reference, which is available only with coverage threshold > 1." | tee -a /dev/stderr)

############### HANDLE EXTERNAL DEPENDENCIES ###############

##	GNU Parallel Dependency
parallel="false"
if hash parallel 2>/dev/null; then
ver=`parallel --version | awk 'NR==1{print \$3}'`
[ $ver -ge 20150322 ] && parallel="true"
fi
[ $parallel == "false" ] && echo ":| WARNING: GNU Parallel version 20150322 or later not installed. We highly recommend to install it to increase performance. Starting pipeline without parallelization!" | tee -a /dev/stderr


###################   HANDLE ARGUMENTS  ####################

([ -z $1 ] || (! [[ $1 =~ \.fasta$ ]] && ! [[ $1 =~ \.fa$ ]] && ! [[ $1 =~ \.fna$ ]])) && echo >&2 "Not sure how to parse your input for reference: file not found at expected locations or has an unrecognized extension. Exiting!" && echo >&2 "$USAGE" && exit 1

if [ `dirname $1` != `pwd` ] ; then
	cmp --silent ${1} `basename $1` || ln -sf $1
fi
reference=`basename $1`

num_sam_files=`echo $@ | xargs -n 1 echo | grep -c ".sam$"`

[ ${num_sam_files} -lt 1 ] && echo ":( Alignment file extentions not recognized. Exiting!" >&2 && exit 1

[ $((${num_sam_files}+1)) -ne $# ] && echo ":| WARNING: some input files have unexpected extensions. These will be ignored!" | tee -a /dev/stderr


			################# START OF MAIN WORKFLOW ###############
			

################# STEP 0: ORGANIZE FOLDER, PREPARE REFERENCE  #################

if [ "$skip_prep" = false ]; then

	printf "\n:) Starting preparatory work: organizing workspace, analyzing and prepping the reference: %s.\n" "`date`" | tee -a /dev/stderr

    [ ! -d $debugDir ] && mkdir -m 777 ${debugDir} || ([[ $(ls -A $debugDir) ]] && rm ${debugDir}/*)
    [ ! -d $samDir ] && mkdir -m 777 ${samDir} || ([[ $(ls -A $samDir) ]] && rm ${samDir}/*)
    [ ! -d $tmpDir ] &&  mkdir -m 777 ${tmpDir} || ([[ $(ls -A $tmpDir) ]] && rm ${tmpDir}/*)
	[ ! -d $splitDir ] && mkdir -m 777 ${splitDir} || ([[ $(ls -A $splitDir) ]] && rm ${splitDir}/*)
	[ ! -d $mapDir ] && mkdir -m 777 ${mapDir} || ([[ $(ls -A $mapDir) ]] && rm ${mapDir}/*)
	[ ! -d $reduceDir ] && mkdir -m 777 ${reduceDir} || ([[ $(ls -A $reduceDir) ]] && rm ${reduceDir}/*)
	
    rm -f $samDir/* && echo $@ | xargs -n 1 | grep ".sam$" | parallel --will-cite safegetfullpath | xargs -I % ln -sf % sam
	[ -f ${logfile} ] && rm -f ${logfile}
	
	awk -f ${pipeline}/utils/generate-cprops-file.awk ${reference} | LC_ALL=C sort -k 1,1 | awk '{$2=NR}1' > ${reference}.sorted.cprops
	
	#bash ${pipeline}/finalize/construct-fasta-from-asm.sh ${reference}.cprops <(awk '{print $2}' ${reference}.sorted.cprops) ${reference} > ${reference}.sorted

	bash ${pipeline}/finalize/construct-fasta-from-asm.sh ${reference}.sorted.cprops <(awk '{print $2}' ${reference}.sorted.cprops) ${reference} | parallel --will-cite --pipe -k awk -f ${pipeline}/map/map-fasta-to-pos-array.awk > ${mapDir}/ref.map

	[ "$single_step" = "true" ] && (printf ":) Done and exiting: %s!\n" "`date`" | tee -a /dev/stderr) && exit
	
	printf ":) Finished with preparatory steps: %s.\n" "`date`" | tee -a /dev/stderr
fi

##################### SORT AND SPLIT SAM ####################

cd ${topDir}

if [ "$skip_sort" = false ]; then

	printf "\n:) Starting to sort sam file(s): %s.\n" "`date`" | tee -a /dev/stderr

	# check that folder structure and expected files are in place:
    if [ ! -d "$samDir" ] || [ ! -d "$debugDir" ] || [ ! -d "$tmpDir" ] || [ ! -d "$splitDir" ]; then
		printf >&2 "\t:( Can't find necessary files associated with previous steps! Please make sure you have preprocessed data or start from scratch. Exiting!\n" && exit 1
        exit 1
    fi
    	
	cmd="LC_ALL=C sort -T ${tmpDir} -k3,3 -k 4,4n -S8G --parallel=48 -s ${samDir}/* > ${samDir}/master.sam &&
	unaligned=\$(awk '{if(\$3==\"*\"){counter++}else{exit}}END{print counter+0}' ${samDir}/master.sam) &&
	(head -\${unaligned} > $splitDir/job_000.sam; split -a 3 --number=l/${jobcount} --numeric-suffixes=1 --additional-suffix=.sam - ${splitDir}/job_) < ${samDir}/master.sam"

	printf "\t%s\n" "$cmd" 
	eval ${cmd}
	
	[ "$single_step" = "true" ] && (printf ":) Done and exiting: %s!\n" "`date`" | tee -a /dev/stderr) && exit
	
	printf ":) Finished sorting sam file(s): %s.\n" "`date`" | tee -a /dev/stderr

fi

# just in case get jobcount if this is not an end-to-end run on a computer with different resources. Parameter jobcount is the actual job count - 1. TODO: adjust the naming.
jobcount=`find $splitDir -name "job_*.sam" | wc -l`
jobcount=$((jobcount-1))
if [ $jobcount -eq -1 ]; then
	printf >&2 "\t:( Can't find necessary files associated with previous steps! Please make sure you have preprocessed data or start from scratch. Exiting!\n" && exit 1
fi

#########################  MAPPING STEP  ################################
cd ${topDir}

if [ "$skip_map" = false ]; then

	printf "\n:) Starting to map sam file(s): %s.\n" "`date`" | tee -a /dev/stderr

	# check that stuff exists from prev steps, mostly for stage relaunch
	if [ ! -d "$samDir" ] || [ ! -d "$debugDir" ] || [ ! -d "$tmpDir" ] || [ ! -d "$splitDir" ]
	then
		printf "\t:( Can't find necessary files associated with previous steps! Please make sure you have preprocessed data or start from scratch. Exiting!" >&2
		exit 1
	fi
		
	for i in `seq -f "%03g" 0 $jobcount`
	do
		if [ ! -f ${splitDir}/"job_"$i".sam" ]
		then 
		printf >&2 "\t:( Can't find necessary files associated with previous steps! Please make sure you have preprocessed data or start from scratch. Exiting!\n" && exit 1
		fi
	done
	
	# handle mapping with parallel to have proper logging
	cmd="seq -f \"%03g\" 1 $jobcount | parallel --will-cite -j ${jobcount} \"printf \"\t:) Starting to map $splitDir/job_{}.sam: \" && date && awk -v mapq_threshold=${mapq_threshold} -v filename=$splitDir/job_{} -f ${pipeline}/map/map-sam-to-pos-array.awk ${reference}.sorted.cprops $splitDir/job_{}.sam && printf \"\t:) Finished mapping $splitDir/job_{}.sam: \" && date\""
	printf "\t%s\n" "$cmd"
	
	seq -f "%03g" 1 $jobcount | parallel --will-cite -j ${jobcount} "printf \"\t:) Starting to map $splitDir/job_{}.sam: \" && date && awk -v mapq_threshold=${mapq_threshold} -v filename=$splitDir/job_{} -f ${pipeline}/map/map-sam-to-pos-array.awk ${reference}.sorted.cprops $splitDir/job_{}.sam && printf \"\t:) Finished mapping $splitDir/job_{}.sam: \" && date"
		
	[ "$single_step" = "true" ] && (printf ":) Done and exiting: %s!\n" "`date`" | tee -a /dev/stderr) && exit

	printf ":) Finished mapping sam file(s): %s.\n" "`date`" | tee -a /dev/stderr

fi

######################### MERGE MAPPING RESULTS ################################

cd ${topDir}

if [ "$skip_merge" = false ]
then
	
	printf "\n:) Starting to merge map file(s): %s.\n" "`date`" | tee -a /dev/stderr
	
	if [ ! -f $mapDir/ref.map ] || [ ! -d $splitDir ] || [ ! -d $mapDir ] || [ ! -d $tmpDir ]
	then
		printf >&2 "\t:( Can't find necessary files associated with previous steps! Please make sure you have preprocessed data or start from scratch. Exiting!\n" && exit 1
	fi 	

	[[ $(ls -A $tmpDir) ]] && rm ${tmpDir}/*

	# First accumulate mapping data from all chunks.
	printf "\tMerging maps: %s.\n" "`date`" | tee -a /dev/stderr

	for i in `seq -f "%03g" 1 $jobcount`
	do
		[ -f ${splitDir}/"job_"$i".map.txt" ] || continue
		mkfifo $tmpDir/"job_"$i
		awk 'BEGIN{counter=0}$0~/@/{while(counter!=substr($1,2)){print ""; counter++}; next}1' ${splitDir}/"job_"$i".map.txt" > $tmpDir/"job_"$i &
	done
	
	## some thread racing at the beginning is going to happen... 
	
	# log
	cmd="paste $mapDir/ref.map $tmpDir/* | parallel --will-cite --pipe -k -j $((jobcount/2 + 1)) \"awk '{for(i=1;i<=NF;i++){entry[(i+4) % 5]+=\$i}; print entry[0], entry[1], entry[2], entry[3], entry[4]; entry[0]=0; entry[1]=0; entry[2]=0; entry[3]=0; entry[4]=0}'\" > ${mapDir}/master.map"
	printf "\t%s" "$cmd"
	
	#eval $cmd # check and substitute
	
	paste $mapDir/ref.map $tmpDir/* | parallel --will-cite --pipe -k -j $((jobcount/2 + 1)) "awk '{for(i=1;i<=NF;i++){entry[(i+4) % 5]+=\$i}; print entry[0], entry[1], entry[2], entry[3], entry[4]; entry[0]=0; entry[1]=0; entry[2]=0; entry[3]=0; entry[4]=0}'" > ${mapDir}/master.map
		
	rm 	$tmpDir/*
			
	# Once main mapping is finished merge indels.
	printf "\tMerging indels: %s.\n" "`date`" | tee -a /dev/stderr
	
	# log
	cmd="seq -f \"%03g\" 1 $jobcount | xargs -i bash -c \"test -f ${splitDir}/job_{}.indel.txt && cat ${splitDir}/job_{}.indel.txt\" | sort --parallel=16 -k 1,1n -k 2,2 | awk 'NR==1{prev1=$1; prev2=$2; prev3=$3; next}$1==prev1&&$2==prev2{prev3+=$3;next}{print prev1, prev2, prev3; prev1=$1; prev2=$2; prev3=$3}END{print prev1, prev2, prev3}' > ${mapDir}/master.indel.txt"
	printf "\t%s" "$cmd"
	
	seq -f "%03g" 1 $jobcount | xargs -i bash -c "test -f ${splitDir}/job_{}.indel.txt && cat ${splitDir}/job_{}.indel.txt" | sort --parallel=16 -k 1,1n -k 2,2 | awk 'NR==1{prev1=$1; prev2=$2; prev3=$3; next}$1==prev1&&$2==prev2{prev3+=$3;next}{print prev1, prev2, prev3; prev1=$1; prev2=$2; prev3=$3}END{print prev1, prev2, prev3}' > ${mapDir}/master.indel.txt
		
	# Once main merging of indels is finished merge clips.
	printf "\tMerging clips: %s.\n" "`date`" | tee -a /dev/stderr
	
	# log
	cmd="seq -f \"%03g\" 1 $jobcount | xargs -i bash -c \"test -f ${splitDir}/job_{}.clip.txt && cat ${splitDir}/job_{}.clip.txt\" | sort --parallel=16 -k 1,1n | awk 'NR==1{prev1=$1; prev2=$2; next}$1==prev1{prev2+=$2;next}{print prev1, prev2; prev1=$1; prev2=$2}END{print prev1, prev2}' > ${mapDir}/master.clip.txt"
	printf "\t%s" "$cmd"
	
	seq -f "%03g" 1 $jobcount | xargs -i bash -c "test -f ${splitDir}/job_{}.clip.txt && cat ${splitDir}/job_{}.clip.txt" | sort --parallel=16 -k 1,1n | awk 'NR==1{prev1=$1; prev2=$2; next}$1==prev1{prev2+=$2;next}{print prev1, prev2; prev1=$1; prev2=$2}END{print prev1, prev2}' > ${mapDir}/master.clip.txt
		
	
	# Filter obvious clips leaving only ones that cannot be explained by restriction site [Note that heterozygosity will cause problems here. Perhaps will be done during mapping data in some future iterations. TODO: also filter obvious clips due to gaps in reference.]
	printf "\tFiltering clips: %s.\n" "`date`" | tee -a /dev/stderr
	
	cmd="(paste ${mapDir}/master.map ${mapDir}/ref.map | parallel --will-cite --pipe -j ${jobcount} -k -N 100000 --block 1G \"var=`echo {#}` && var=\$(((var-1)*100000)) && awk -f ${pipeline}/reduce/basic-reduce-map.awk - \" | awk -v seq=${restriction_sequence} -v t=1 -f ${pipeline}/reduce/basic-filter-clips.awk ${mapDir}/master.clip.txt - | sort -k 1,1n > ${mapDir}/master.clip.filtered.txt) 2> >(tee -a /dev/stderr 2>&1)"
	printf "\t%s" "$cmd"
	
	(paste ${mapDir}/master.map ${mapDir}/ref.map | parallel --will-cite --pipe -j ${jobcount} -k -N 100000 --block 1G "var=`echo {#}` && var=\$(((var-1)*100000)) && awk -f ${pipeline}/reduce/basic-reduce-map.awk - " | awk -v seq=${restriction_sequence} -v t=1 -f ${pipeline}/reduce/basic-filter-clips.awk ${mapDir}/master.clip.txt - | sort -k 1,1n > ${mapDir}/master.clip.filtered.txt) 2> >(tee -a /dev/stderr 2>&1)
	
	[ "$single_step" = "true" ] && (printf ":) Done and exiting: %s!\n" "`date`" | tee -a /dev/stderr) && exit
	
	printf ":) Finished merging map file(s): %s.\n" "`date`" | tee -a /dev/stderr
	
fi

######################### BASIC REDUCE STEP ###############################

cd ${topDir}

if [ "$skip_reduce" = false ]
then
	
	printf "\n:) Starting basic reduce map into fasta: %s.\n" "`date`" | tee -a /dev/stderr

	# check that stuff is in from previous steps
	if [ ! -f ${mapDir}/ref.map ] || [ ! -f ${mapDir}/master.map ] || [ ! -f ${mapDir}/master.indel.txt ] || [ ! -f ${mapDir}/master.clip.filtered.txt ] || [ ! -d ${splitDir} ] || [ ! -d ${mapDir} ] || [ ! -d ${tmpDir} ] || [ ! -d ${reduceDir} ]
	then
		printf >&2 "\t:( Can't find necessary files associated with previous steps! Please make sure you have preprocessed data or start from scratch. Exiting!\n" && exit 1
	fi

#  	(paste ${mapDir}/master.map ${mapDir}/ref.map | parallel --tmpdir=${tmpDir} --will-cite --pipe -j ${jobcount} -k -N 100000 --block 1G "var=`echo {#}` && var=\$(((var-1)*100000)) && awk -v consensus_threshold=${consensus_threshold} -v shift=\${var} -v ignore_reference=${ignore_reference} -f ${pipeline}/reduce/indel-reduce-map.awk ${mapDir}/master.indel.txt ${mapDir}/master.clip.filtered.txt - " > ${reduceDir}/master.reduce.txt) 2> "${reduceDir}/unaccounted.txt"
  	#>(tee -a ${logfile} 2>&1)
  	
  	(paste ${mapDir}/master.map ${mapDir}/ref.map | awk -v consensus_threshold=${consensus_threshold} -v ignore_reference=${ignore_reference} -f ${pipeline}/reduce/indel-reduce-map.awk ${mapDir}/master.indel.txt ${mapDir}/master.clip.filtered.txt - > ${reduceDir}/master.reduce.txt) 2> "${reduceDir}/unaccounted.txt"
  	
  	sort -k 1,1n ${reduceDir}/unaccounted.txt > ${reduceDir}/unaccounted.sorted.txt && rm ${reduceDir}/unaccounted.txt
  	
  	[ "$single_step" = "true" ] && (printf ":) Done and exiting: %s!\n" "`date`" | tee -a /dev/stderr) && exit

	printf ":) Finished basic reduce map into fasta: %s.\n" "`date`" | tee -a /dev/stderr

fi

########################## PATCH WITH REASSEMBLY ########################

cd ${topDir}

touch ${reduceDir}/reassembly_output.txt

if [ "$skip_reassembly" = false ]
then

	printf "\n:) Starting local reassembly to improve fasta: %s.\n" "`date`" | tee -a /dev/stderr

	# check that stuff is in from previous steps
	if [ ! -f $mapDir/ref.map ] || [ ! -f $mapDir/master.map ] || [ ! -f ${reduceDir}/unaccounted.sorted.txt ] || [ ! -d $splitDir ] || [ ! -d $mapDir ] || [ ! -d $tmpDir ] || [ ! -d $reduceDir ]
	then
		printf >&2 "\t:( Can't find necessary files associated with previous steps! Please make sure you have preprocessed data or start from scratch. Exiting!\n" && exit 1
	fi
	
	# list problematic regions as active_regions.bed file and extract associated sequences from sam file
	awk -f ${pipeline}/reassemble/list-active-regions-weighted-by-cov.awk ${reduceDir}/unaccounted.sorted.txt ${mapDir}/master.map | \
	awk 'FNR==1{start=$1; end=$2;next}$1<=end{end=$2;next}{print start, end; start=$1; end=$2}' | tee >(awk 'BEGIN{print "header"; OFS="\t"}{$9=$1; $10=$2}1' - | \
	awk -f ${pipeline}/utils/lift-asm-annotations-to-input-annotations.awk ${reference}.sorted.cprops <(awk '{print $2}' ${reference}.sorted.cprops) - | \
	awk 'BEGIN{FS="\t";OFS="\t"}NR>1{print $1, $2, $3}' > active_regions.bed) | \
	awk -f ${pipeline}/reassemble/cut-multiple-anchored-regions-from-sam.awk ${reference}.sorted.cprops - ${samDir}/master.sam > anchored_cut_output.txt
	#| parallel --will-cite --pipe --jobs 80% --rrs -N1 --regex --recend '---' "read h; read h && printf \"\$h \" && awk -v anchored=1 -f ${pipeline}/reassemble/de-bruijn-assemble-DFS-both-strands.awk -" > ${reduceDir}/reassembly_output.txt

  	[ "$single_step" = "true" ] && (printf ":) Done and exiting: %s!\n" "`date`" | tee -a /dev/stderr) && exit

	printf ":) Finished local reassembly to improve fasta: %s.\n" "`date`" | tee -a /dev/stderr

fi

########################### FINALIZE TO FASTA ########################### 	

cd ${topDir}

if [ "$skip_finalize" = false ]
then

	printf "\n:) Starting to finalize output: %s.\n" "`date`" | tee -a /dev/stderr

	# check that stuff is in from previous steps
	if [ ! -f ${reference}.sorted.cprops ] || [ ! -f ${reduceDir}/master.reduce.txt ] || [ ! -f ${reduceDir}/reassembly_output.txt ] || [ ! -d $splitDir ] || [ ! -d $mapDir ] || [ ! -d $tmpDir ] || [ ! -d $reduceDir ]
	then
		printf >&2 "\t:( Can't find necessary files associated with previous steps! Please make sure you have preprocessed data or start from scratch. Exiting!\n" && exit 1
	fi

	if [ -z ${assembly} ]; then
		
		# if no assembly file is defined dump fasta without modification, i.e. one ref-assisted fasta sequence for every fasta sequence in the reference (except for those that end up empty)
		ln -sf ${reference}.sorted.cprops ${reference}.final.cprops
		awk '{print $2}' ${reference}.final.cprops > ${reference}.final.asm
	else		
		if ! [ -f `basename $assembly` ]; then	
			ln -sf $assembly .
			assembly=`basename $assembly`
		fi
		awk -f ${pipeline}/utils/convert-assembly-to-cprops-and-asm.awk ${assembly}		
		bash ${pipeline}/utils/reconstruct-edits-from-cprops.sh `basename ${assembly} .assembly`.cprops | awk -f ${pipeline}/utils/edit-cprops-according-to-annotations.awk - ${reference}.sorted.cprops > ${reference}.final.cprops		
		awk -f ${pipeline}/utils/from-cnames-to-onames.awk `basename ${assembly} .assembly`.cprops `basename ${assembly} .assembly`.asm | awk -f ${pipeline}/utils/from-onames-to-cnames.awk ${reference}.final.cprops - > ${reference}.final.asm
		rm `basename ${assembly} .assembly`.cprops `basename ${assembly} .assembly`.asm
	fi

  	awk -f ${pipeline}/reassemble/patch.awk ${reduceDir}/reassembly_output.txt ${reduceDir}/master.reduce.txt | awk -f ${pipeline}/finalize/split-lines-according-to-cprops.awk ${reference}.final.cprops - | awk -f ${pipeline}/utils/wrap-fasta-sequence.awk - > ${reference}.final.fasta
  	
  	awk -f ${pipeline}/utils/generate-cprops-file.awk ${reference}.final.fasta > ${reference}.final.cprops
  	bash ${pipeline}/finalize/finalize-output.sh -l "assisted" ${reference}.final.cprops ${reference}.final.asm ${reference}.final.fasta final
		
fi	

	[ "$single_step" = "true" ] && (printf ":) Done and exiting: %s!\n" "`date`" | tee -a /dev/stderr) && exit

	printf ":) Finished finalizing output: %s.\n" "`date`" | tee -a /dev/stderr
	
########################### THE END ########################### 	

	printf "\n:) Finished assisted assembly pipeline: %s!\n" "`date`" | tee -a /dev/stderr
	printf "\n-------------------------\n" | tee -a /dev/stderr



##
##
##
##
##
##
##
