#!/bin/bash -x

USAGE="
*****************************************************
This is a wrapper for a the Anonamos comparative genome assembler which analyzes the alignment data of the dna sequences derived from species A to a genome of closely related species B to create ref-assisted contigs for species A. (This is a short workflow version of Anonamos with no provisional genome creation and a single map-reduce procedure.) The input to the assembler are single-end SAM alignment files (note that performance might increase if the files are preprocessed: split at chimeric junctions, dedupped, filtered for overlaps). At this point we just use raw sam output.

	The first step (-S sort) the sam file(s) are merge-sorted into a master sorted SAM file. The master SAM is then split into manageable chunks. The split size is adjusted dynamically to result in a predefined number of files determined by jobcount parameter (unprompted).
	In the second step (-S map) a per base nucleotide map is constructed for each scaffold in the reference assembly. Indels and read breaks positions with respect to the referene are also recorded. The step includes merging of the chunk map results. The result is a master map file.
	The third (optional?) step (-S addref) is adding reference nucleotides to the map and reference gap data.
	The forth step (-S reduce) is consolidation of preliminary fasta. Should probably include some reassembly but right now is just heuristics.

Version date: October 12, 2017.

Usage: ./run-anonamos-assembler.sh [-h] [-c|--consensus consensus_threshold] [-q|--mapq mapq_threshold] [-s|--stage pipeline_stage] reference_fasta path_to_sam_file(s)

ARGUMENTS:
path_to_sam_file        Path to sam file describing read alignment to sequence used as assisting reference

MAIN OPTIONS:
-h                      			Shows this help
-c|--consensus	consensus_threshold The level of required consensus required for a base call, i.e. how many times the base needs to be seen across reads (and reference if c>1) to constitute a valid call [default is 1].
-q|--mapq		mapq_threshold      Minimum mapping quality of the read to be considered in the recontruction process [default is 0, i.e. use all reads]
-s|--stage		stage				Start from a particular stage. Can be split, map, addref, consol etc.

SUPPLEMENTARY OPTIONS:
--ignore-reference					Do not consider reference as contributing to consensus [applicable only for c>1, if c=1 ignore_reference is forced].

*****************************************************
"

pipeline=`cd "$( dirname $0)" && pwd`


###################### SET DEFAULTS #############################

shopt -s extglob # enable pathname expansions

# problem-specific

ignore_reference=0
consensus_threshold=1 # trust single reads
mapq_threshold=0 # include mapq0
restriction_sequence="GATC"

skip_prep=false
skip_sort=false
skip_map=false
skip_merge=false
skip_reduce=false

jobcount=`grep -c ^processor /proc/cpuinfo`
jobcount=$((jobcount*80/100))

# organizational

topDir=$(pwd)

debugDir=${topDir}"/debug"
tmpDir=${topDir}"/tmp"
samDir=${topDir}"/sam"
splitDir=${topDir}"/split"
mapDir=${topDir}"/map"

logfile=$debugDir/log.txt


###################### HANDLE OPTIONS ###############################

while :; do
	case $1 in
		-h)
			echo "$USAGE" >&1
			exit 0
			shift
        ;;
        -c|--coverage) OPTARG=$2
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
				echo ":) -q flag was triggered. Read mapping quality threshold is set to $OPTARG." >&1
				mapq_threshold=$OPTARG
			fi
			shift
		;;
		-s|--stage) OPTARG=$2
			stage=$OPTARG
			shift
		;;
		-j|--jobs) OPTARG=$2
			jobcount=$OPTARG
			shift
		;;
		--ignore-reference)
			(echo ":) --ignore-reference flag was triggered. Reference will not be considered when making consensus calls.") >&1
			ignore_reference=1
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

# force to ignore reference is consensus is set to 1
[ "$consensus_threshold" -eq 1 ] && ignore_reference=1

############### HANDLE EXTERNAL DEPENDENCIES ###############

##	GNU Parallel Dependency
parallel="false"
if hash parallel 2>/dev/null; then
ver=`parallel --version | awk 'NR==1{print \$3}'`
[ $ver -ge 20150322 ] && parallel="true"
fi
[ $parallel == "false" ] && echo ":| WARNING: GNU Parallel version 20150322 or later not installed. We highly recommend to install it to increase performance. Starting pipeline without parallelization!" >&2


###################   HANDLE ARGUMENTS  ####################

safegetfullpath(){
echo `cd $(dirname $1) && pwd -P`"/"`basename $1`
}
export -f safegetfullpath

([ -z $1 ] || (! [[ $1 =~ \.fasta$ ]] && ! [[ $1 =~ \.fa$ ]] && ! [[ $1 =~ \.fna$ ]])) && echo >&2 "Not sure how to parse your input for reference: file not found at expected locations or has an unrecognized extension. Exiting!" && echo >&2 "$USAGE" && exit 1

if [ `dirname $1` != `pwd` ] ; then
	cmp --silent ${1} `basename $1` || ln -sf $1
fi
reference=`basename $1`

num_sam_files=`echo $@ | xargs -n 1 echo | grep -c ".sam$"`

[ ${num_sam_files} -lt 1 ] && echo ":( Alignment file extentions not recognized. Exiting!" >&2 && exit 1

[ $((${num_sam_files}+1)) -ne $# ] && echo ":| WARNING: some input files have unexpected extensions. These will be ignored!" >&2


################# STEP 0: ORGANIZE FOLDER, PREPARE REFERENCE  #################

if [ "$skip_prep" = false ]; then

	(printf "\n:)Organizing workspace.\n:)Analyzing and prepping the reference.\n") >&1

    [ ! -d $debugDir ] && mkdir -m 777 ${debugDir} || ([[ $(ls -A $debugDir) ]] && rm ${debugDir}/*)
    [ ! -d $samDir ] && mkdir -m 777 ${samDir} || ([[ $(ls -A $samDir) ]] && rm ${samDir}/*)
    [ ! -d $tmpDir ] &&  mkdir -m 777 ${tmpDir} || ([[ $(ls -A $tmpDir) ]] && rm ${tmpDir}/*)
	[ ! -d $splitDir ] && mkdir -m 777 ${splitDir} || ([[ $(ls -A $splitDir) ]] && rm ${splitDir}/*)
	[ ! -d $mapDir ] && mkdir -m 777 ${mapDir} || ([[ $(ls -A $mapDir) ]] && rm ${mapDir}/*)
	
    rm -f $samDir/* && echo $@ | xargs -n 1 | grep ".sam$" | parallel --will-cite safegetfullpath | xargs -I % ln -sf % sam
	[ -f ${logfile} ] && rm -f ${logfile}
	
	awk -f ${pipeline}/utils/generate-cprops-file.awk ${reference} > ${reference}.cprops
	
	LC_ALL=C sort -k 1,1 ${reference}.cprops > ${reference}.sorted.cprops
	
	bash ${pipeline}/finalize/construct-fasta-from-asm.sh ${reference}.cprops <(awk '{print $2}' ${reference}.sorted.cprops) ${reference} > ${reference}.sorted

	bash ${pipeline}/finalize/construct-fasta-from-asm.sh ${reference}.cprops <(awk '{print $2}' ${reference}.sorted.cprops) ${reference} | parallel --will-cite --pipe -k awk -f ${pipeline}/utils/parse-reference.awk > ${mapDir}/ref.map


fi

################### WRITE COMMAND TO LOG  #################

echo $* >> ${logfile}
echo >> ${logfile}

##################### SORT AND SPLIT SAM ####################

cd ${topDir}

if [ "$skip_sort" = false ]; then

	# check that folder structure and expected files are in place:
    if [ ! -d "$samDir" ] || [ ! -d "$debugDir" ] || [ ! -d "$tmpDir" ] || [ ! -d "$splitDir" ]; then
        echo ":( Can't find expected folders in working folder! Please make sure you have preprocessed data or start from scratch." >&2
        exit 1
    fi

	(printf "\n:) Starting to sort sam file(s): " && date) >&1 | tee -a ${logfile}
	
	cmd="LC_ALL=C sort -T ${tmpDir} -k3,3 -k 4,4n -S8G --parallel=48 -s ${samDir}/* > ${samDir}/master.sam &&
	unaligned=\$(awk '{if(\$3==\"*\"){counter++}else{exit}}END{print counter+0}' ${samDir}/master.sam) &&
	(head -\${unaligned} > $splitDir/job_000.sam; split -a 3 --number=l/${jobcount} --numeric-suffixes=1 --additional-suffix=.sam - ${splitDir}/job_) < ${samDir}/master.sam"

	echo ${cmd} >> ${logfile}
	eval ${cmd} | tee -a ${logfile}
	echo "" >> ${logfile}
	
	(printf ":) Finished sorting sam file(s): " && date) >&1 | tee -a ${logfile}
fi

#########################  MAPPING STEP  ################################

cd ${topDir}

if [ "$skip_map" = false ]; then
	printf "\n:) Starting to map sam files: " && date >&1 | tee -a ${logfile}

	# check that stuff exists from prev steps, mostly for stage relaunch
	if [ ! -d "$samDir" ] || [ ! -d "$debugDir" ] || [ ! -d "$tmpDir" ] || [ ! -d "$splitDir" ]
	then
		echo ":( Can't find necessary files associated with previous steps! Please make sure you have preprocessed data or start from scratch. Exiting!" >&2
		exit 1
	fi
		
	for i in `seq -f "%03g" 0 $jobcount`
	do
		if [ ! -f ${splitDir}/"job_"$i".sam" ]
		then 
			echo >&2 ":( Can't find necessary files associated with previous steps! Please make sure you have preprocessed data or start from scratch. Exiting!" && exit 1
		fi
	done
	
	# handle mapping with parallel to have proper logging

	seq -f "%03g" 1 $jobcount | parallel --will-cite -j ${jobcount} "printf \"\t:) Starting to map $splitDir/job_{}.sam: \" && date && awk -v mapq_threshold=${mapq_threshold} -v filename=$splitDir/job_{} -f ${pipeline}/map-sam-to-pos-array.awk ${reference}.sorted.cprops $splitDir/job_{}.sam && printf \"\t:) Finished mapping $splitDir/job_{}.sam: \" && date" >> ${logfile}
	
	(printf ":) Finished maping sam files: " && date) >&1 | tee -a ${logfile}

fi
######################### MERGE MAPPING RESULTS ################################

if [ "$skip_merge" = false ]
then
	
	if [ ! -f $mapDir/ref.map ] || [ ! -d $splitDir ] || [ ! -d $mapDir ] || [ ! -d $tmpDir ]
	then
		echo >&2 ":( Can't find necessary files associated with previous steps! Please make sure you have preprocessed data or start from scratch. Exiting!" && exit 1
	fi 
	
	(printf "\n:) Starting to merge map files: " && date) >&1 | tee -a ${logfile}

	[[ $(ls -A $tmpDir) ]] && rm ${tmpDir}/*

		
	for i in `seq -f "%03g" 1 $jobcount`
	do
		[ -f ${splitDir}/"job_"$i".map.txt" ] || continue
		mkfifo $tmpDir/"job_"$i
		awk 'BEGIN{counter=0}$0~/@/{while(counter!=substr($1,2)){print ""; counter++}; next}1' ${splitDir}/"job_"$i".map.txt" > $tmpDir/"job_"$i &
		
# 		might change to +1 global positioning, not sure
# 		start=`awk 'NR==1{print substr($1,2); exit}' ${splitDir}/"job_"$i".map.txt"`
# 		end=`wc -l < ${splitDir}/"job_"$i".map.txt"`
# 		mkfifo $tmpDir/head $tmpDir/body $tmpDir/tail
# 		( head -$start > $tmpDir/head ; head -$end | paste - <(tail -n +2 $splitDir/"job_"$i".map.txt") > $tmpDir/body; cat > $tmpDir/tail ; ) < $mapDir/master.map
# 		
# 		cat $tmpDir/head <(parallel --pipepart -a $tmpDir/body --will-cite -k "awk 'NF>5{\$1+=\$6; \$2+=\$7; \$3+=\$8; \$4+=\$9; \$5+=\$10}{print \$1, \$2, \$3, \$4, \$5}'") $tmpDir/tail > ${splitDir}/"job_"$i".map.txt"
# 				
# 		mv ${splitDir}/"job_"$i".map.txt" $mapDir/master.map
# 		rm 	$tmpDir/head $tmpDir/body $tmpDir/tail
	
	done
	
	# some thread racing at the beginning is going to happen... 
	
	paste $mapDir/ref.map $tmpDir/* | parallel --will-cite --pipe -k -j $((jobcount/2)) "awk '{for(i=1;i<=NF;i++){entry[(i+4) % 5]+=\$i}; print entry[0], entry[1], entry[2], entry[3], entry[4]; entry[0]=0; entry[1]=0; entry[2]=0; entry[3]=0; entry[4]=0}'" > ${mapDir}/master.map
		
	rm 	$tmpDir/*
	
	(printf ":) Finished merging map files: " && date) >&1 | tee -a ${logfile}
		
	# Once main mapping is finished merge indels
	(printf "\n:) Starting to merge indels: " && date) >&1 | tee -a ${logfile}
	
	seq -f "%03g" 1 $jobcount | xargs -i bash -c "test -f ${splitDir}/job_{}.indel.txt && cat ${splitDir}/job_{}.indel.txt" | sort --parallel=16 -k 1,1n -k 2,2 | awk 'NR==1{prev1=$1; prev2=$2; prev3=$3; next}$1==prev1&&$2==prev2{prev3+=$3;next}{print prev1, prev2, prev3; prev1=$1; prev2=$2; prev3=$3}END{print prev1, prev2, prev3}' > ${mapDir}/master.indel.txt
	
	(printf ":) Finished merging indels: " && date) >&1  | tee -a ${logfile}
	
	# Once main merging of indels is finished merge clips
	(printf "\n:) Starting to merge clips: " && date) >&1 | tee -a ${logfile}
	
	seq -f "%03g" 1 $jobcount | xargs -i bash -c "test -f ${splitDir}/job_{}.clip.txt && cat ${splitDir}/job_{}.clip.txt" | sort --parallel=16 -k 1,1n | awk 'NR==1{prev1=$1; prev2=$2; next}$1==prev1{prev2+=$2;next}{print prev1, prev2; prev1=$1; prev2=$2}END{print prev1, prev2}' > ${mapDir}/master.clip.txt
	
	(printf ":) Finished merging clips: " && date) >&1  | tee -a ${logfile}

fi

######################### BASIC REDUCE STEP ###############################

if [ "$skip_reduce" = false ]
then
	# check that stuff is in from previous steps
	if [ ! -f $mapDir/ref.map ] || [ ! -f $mapDir/master.map ] || [ ! -d $splitDir ] || [ ! -d $mapDir ] || [ ! -d $tmpDir ]
	then
		echo >&2 ":( Can't find necessary files associated with previous steps! Please make sure you have preprocessed data or start from scratch. Exiting!" && exit 1
	fi 

	(printf "\n:) Starting to filter clips: " && date) >&1 | tee -a ${logfile}

	(paste ${mapDir}/master.map ${mapDir}/ref.map | parallel --will-cite --pipe -j ${jobcount} -k -N 100000 --block 1G "var=`echo {#}` && var=\$(((var-1)*100000)) && awk -v consensus_threshold=${consensus_threshold} -v shift=\${var} -v ignore_reference=${ignore_reference} -f ${pipeline}/indel-reduce-map.awk ${mapDir}/master.indel.txt ${mapDir}/master.clip.txt - " | awk -v seq=${restriction_sequence} -f ${pipeline}/filter-clips.awk ${mapDir}/master.clip.txt - > ${mapDir}/master.clip.filtered.txt
	(printf "\n:) Finished filtering clips: " && date) >&1 | tee -a ${logfile}

	(printf "\n:) Starting to reduce map into fasta: " && date) >&1 | tee -a ${logfile}

  	(paste ${mapDir}/master.map ${mapDir}/ref.map | parallel --will-cite --pipe -j ${jobcount} -k -N 100000 --block 1G "var=`echo {#}` && var=\$(((var-1)*100000)) && awk -v consensus_threshold=${consensus_threshold} -v shift=\${var} -v ignore_reference=${ignore_reference} -f ${pipeline}/indel-reduce-map.awk ${mapDir}/master.indel.txt ${mapDir}/master.clip.txt - " | awk -f ${pipeline}/utils/wrap-fasta-sequence.awk - > out.fasta) 2> >(tee -a ${logfile} 2>&1)

		
	(printf ":) Finished reducing map into fasta: " && date) >&1 | tee -a ${logfile}

	
	## temp
	#echo ">output" out.fasta > out.fasta.tmp && mv out.fasta.tmp out.fasta
fi


	
exit
	
	
	
	




















	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
#	seq -f "%03g" 0 $(( jobcount )) | xargs -P8 -I % touch "$splitDir/job_"%".sam"

	for (( i=1; i <= $jobcount; i++ ))
	do
		filename=${splitDir}"/job_"$(printf "%03g" $i)".sam"

		cmd="echo \"(: Starting to map ${filename}.\" &&
			date &&
			awk -v mapq_threshold=${mapq_threshold} -f ${pipeline}/map-sam-to-pos-array.awk > ${filename}.map 2> ${filename}.feature &&
			echo \"(: Finished mapping ${filename}.\" &&
			date;"

		case $cluster in
			uger)
				qsub -o ${debugDir}/uger.out -j y -q $long_queue -r y -N "${groupname}_sort" -l m_mem_free=2g ${uger_wait} <<-MAP
					$eval $cmd
				MAP
            ;;
			slurm)
				jid=`sbatch <<- MAP | egrep -o -e "\b[0-9]+$"
				#!/bin/bash -l
				#SBATCH -p $long_queue
				#SBATCH -t 1440
				#SBATCH -c 8
				#SBATCH --ntasks=1
				#SBATCH --mem-per-cpu=25G
				#SBATCH -o $debugDir/map-%j.out
				#SBATCH -e $debugDir/map-%j.err
				#SBATCH -J "${groupname}_map_${jname}"
				${sbatch_wait}
				$cmd
				MAP`
				dependmap="afterok:$jid"
			;;
			local)
            	eval ${cmd} | tee -a ${logfile}
		esac

	done
fi




####################   MAP REFERENCE   #####################
# awk -f ${pipeline}/utils/generate-cprops-file.awk ${reference} | LC_ALL=C sort -k 1,1 > ${reference}.cprops
# 
# bash ${pipeline}/finalize/construct-fasta-from-asm.sh ${reference}.cprops <(LC_ALL=C sort -k 1,1 ${reference}.cprops | awk '{print $2}') ${reference} | awk -f ${pipeline}/utils/parse-reference.awk | awk 'FILENAME==ARGV[1]{name[FNR]=$1;next}$0~/>/{counter++; $0=name[counter]}1' <(LC_ALL=C sort -k 1,1 ${reference}.cprops) - > map.txt


####################   MAP SAM FILES   #####################
addmapfrag(){
	LC_ALL=C sort -k3,3 -k 4,4n -s | awk -v mapq_threshold=${mapq_threshold} -f ${pipeline}/map-sam-to-pos-array.awk 
}
export -f addmapfrag

parallel --will-cite -a $2 -j 80% --pipepart "LC_ALL=C sort -k3,3 -k 4,4n -s | awk -v mapq_threshold=${mapq_threshold} -f ${pipeline}/map-sam-to-pos-array-new.awk | LC_ALL=C sort -k1,1 -k 2,2" > temp.map 2> ${filename}.feature

exit


#******




exit 0

cd ${topDir}


################### ADD REFERENCE DATA (optional?) ##########################

if [ "$skip_addref" = false ]; then

	echo "(: Launching jobs to add reference data to the map (optional step)."

	for (( i=1; i <= $jobcount; i++ ))
	do
		if [ "$skip_map" = false ]; then
			sbatch_wait="#SBATCH -d ${dependmap[$i]}"
		else
			sbatch_wait=""
		fi

		jid=`sbatch <<- ADDREF | egrep -o -e "\b[0-9]+$"
		#!/bin/bash -l
		#SBATCH -p $queue
		#SBATCH -o $debugDir/addref-%j.out
		#SBATCH -e $debugDir/addref-%j.err
		#SBATCH -t 1200
		#SBATCH -c 1
		#SBATCH --ntasks=1
		#SBATCH -J "${groupname}_addref_${jname}"
		${sbatch_wait}
		srun bash ${add_ref_script} ${debugDir}"/cprops_partition.$i.txt" ${reference}
		ADDREF`
		dependaddref[$i]="afterok:$jid"
	done
fi

################## CONSOLIDATE PRELIMINARY FASTA ##########################

cd ${topDir}
if [ "$skip_reduce" = false ]; then

	echo "(: Launching jobs to make preliminary fasta."

	for (( i=1; i <= $jobcount; i++ ))
	do
#		if [ "$skip_addref" = false ]; then
#			sbatch_wait="#SBATCH -d ${dependaddref[$i]}"
#		else
#			sbatch_wait=""
#		fi

		if [ "$skip_addref" = true ] || [ -z ${dependaddref[$i]} ]; then
			sbatch_wait=""
		else
			sbatch_wait="#SBATCH -d ${dependaddref[$i]}"
		fi

#		cmd=""
#		while read chrname skip
#		do
#			if [ -f ${mapDir}/${chrname}.map.txt ]; then
#				cmd="bash ${reduce_script} ${chrname} ${consensus_threshold}; $cmd"
#			fi
#		done <${debugDir}"/cprops_partition.$i.txt"

#		if [ "$cmd" = "" ]; then
#			continue
#		fi

		jid=`sbatch <<- reduce | egrep -o -e "\b[0-9]+$"
		#!/bin/bash -l
		#SBATCH -p $queue
		#SBATCH -o $debugDir/reduce-%j.out
		#SBATCH -e $debugDir/reduce-%j.err
		#SBATCH -t 1440
		#SBATCH -n 1
		#SBATCH --ntasks=1
		#SBATCH --mem-per-cpu=32G
		#SBATCH -J "${groupname}_reduce_${jname}"
		${sbatch_wait}
		srun bash ${reduce_script} ${debugDir}"/cprops_partition.$i.txt" ${consensus_threshold}
		reduce`
		dependreduce[$i]="afterok:$jid"
	done
fi

exit 0

#			awk -v consensus_threshold="$consensus_threshold" -f ${reduce_script} ${mapDir}/{$chrname}.indel.txt ${mapDir}/${chrname}.map.txt | awk -f ${break_script} > ${mapDir}/${chrname}.prelim.fasta



### PRELIMINARY WORK - extract relevant reference fragment and map gaps in it
awk -v chr="$chr" -v start="$region_start" -v end="$region_end" -f ${local_fasta} ${reference} > local_ref.fa
awk -f ${map_gaps} local_ref.fa > ref_gaps.bed
ref_gaps="ref_gaps.bed"

### MAP SAM FILE TO POSITIONS
echo "(: Building a library of observed nucleotides."

if [ -z "$chr" ]; then
    filter="!(\$0~/@/)"
elif [ -z "$region_start" ]; then
    filter="!(\$0~/@/)&&\$3==\"$chr\""
else
    filter="!(\$0~/@/)&&\$3==\"$chr\"&&\$4>=$region_start&&\$4<=$region_end"
fi

cmd="awk '{if (int(NR/2)!= NR/2){chr=substr(\$1,2); start=\$2; len=\$3; end=start+len-1}else{print chr\"_\"start\"_\"end, 0, chr, start, 60, len\"M\", \"*\", 0, 0, \$0}}' local_ref.fa | cat - $@ | awk '${filter}{print}'"
eval $cmd | awk -v mapq_threshold="$mapq_threshold" -f ${map_script} > map.txt
sort -k1,1 -k2,2n map.txt > map.txt.sort
mv map.txt.sort map.txt

sort -k1,1 -k2,2n -k3,3n -k4,4 indels.txt > indels.txt.sort
mv indels.txt.sort indels.txt

### DUMP PRELIMINARY RECONSTRUCTION
echo "(: Reconstructing preliminary fasta."
awk -v consensus_threshold="$consensus_threshold" -f ${reduce_script} indels.txt map.txt > prelim.reconstructed.fasta
awk -f ${break_script} prelim.reconstructed.fasta > prelim.reconstructed.contigs.fasta

################ supp for visualization ################
bash ${prep_for_igv} prelim.reconstructed.contigs.fasta ${panTro4}
bash ${prep_for_igv} prelim.reconstructed.contigs.fasta ${hg38}
awk -f ${print_contig_N50} prelim.reconstructed.contigs.fasta
################ supp for visualization ################

### POLISHING: try to assemble through reference gaps
awk '{print $1, $2-100, $3-$2+200}' ${ref_gaps} > active_regions.txt

sort -k1,1 -k2,2n -k3,3n -k4,4 active_regions.txt > active_regions.txt.sort
mv active_regions.txt.sort active_regions.txt

### POLISHING: flag indels and merge suspicious regions


################ supp for visualization ################
#awk -f ${txt_to_bed_script} indels.txt > indel.bed
awk -f ${txt_to_bed_script} active_regions.txt > active_regions.bed
################ supp for visualization ################


### REASSEMBLE SUSPICIOUS REGIONS
echo "(: Reassembling suspicious regions."
if [ -f reassembled.indels.txt ]; then
    rm -f reassembled.indels.txt
fi
cat active_regions.bed | parallel -j10 --will-cite "awk -v region={.} -f ${extract_script} $@ | awk -v region={.} -v anchored="1" -f ${reassemble_script} >> reassembled.indels.txt"

sort -k1,1 -k2,2n -k3,3n -k4,4 reassembled.indels.txt > reassembled.indels.txt.sort
mv reassembled.indels.txt.sort reassembled.indels.txt

### RECONSTRUCTING FASTA
echo "(: Reconstructing polished fasta."
awk -v consensus_threshold="$consensus_threshold" -f ${reduce_script} reassembled.indels.txt map.txt > polished.reconstructed.fasta
awk -f ${break_script} polished.reconstructed.fasta > polished.reconstructed.contigs.fasta

################ supp for visualization ################
hg38="/Users/olga/WORK/CAMBRIDGE/Genetics/ASSEMBLY/Assisted_assembly/anonamos/hg38.fa"
panTro4=${reference}
bash ${prep_for_igv} polished.reconstructed.contigs.fasta ${hg38}
bash ${prep_for_igv} polished.reconstructed.contigs.fasta ${panTro4}
awk -f ${print_contig_N50} polished.reconstructed.contigs.fasta
################ supp for visualization ################
exit 0

### FLAG MORE SUSPICIOUS REGIONS
echo "(: Listing candidate regions for reassembly."
awk -v sig_threshold=2 -f ${flag_script} indels.txt > flagged_regions.bed
awk -v sig_threshold=1 -f ${flag_script} sig_indels.txt > sig_flagged_regions.bed





#### MAKE A PRELIMINARY RECONSTRUCTION
#echo "(: Building a library of observed nucleotides."
#cmd="awk '!(\$0~/@/)${add_filter}{print}' ${1}"
#eval $cmd | awk -v mapq_threshold="$mapq_threshold" -f ${map_script} > map.txt
#sort -k1,1 -k2,2n map.txt > map.txt.sort
#mv map.txt.sort map.txt
#sort -k1,1 -k2,2n -k3,3n -k4,4 indels.txt > indels.txt.sort
#mv indels.txt.sort indels.txt

################temp - supp for visualization
awk -f ${txt_to_bed_script} indels.txt > indel.bed
#################temp - supp for visualization

echo "(: Dumping reconstructed fasta."
awk -v consensus_threshold="$consensus_threshold" -f ${reduce_script} indels.bed map.txt > prelim.reconstructed.fasta

#### POLISH RECONSTRUCTION BY REALIGNING READS AND REASSEMBLING SUSPICIOUS SECTIONS
echo "(: Polishing the new genome: realigning reads to the new sequence. NOTE: It is better if this step is done genome-wide."
bwa index prelim.reconstructed.fasta
cmd="awk '!(\$0~/@/)${add_filter}{if (\$10!~/GATCGATC/){print \">\"\$1; print \$10;} else {n=split(\$10,tmp,\"GATCGATC\"); for (i=2; i<=n-1; i++) {print \">\"\$1\"_\"i; print \"GATC\"tmp[i]\"GATC\";}; print \">\"\$1\"_\"1; print tmp[1]\"GATC\"; print \">\"\$1\"_\"n; print \"GATC\"tmp[n];}}' ${1}" #break reads - should not be needing to do this if the reads are broken in the first place.
#echo "$cmd"
bwa mem prelim.reconstructed.fasta <(eval ${cmd})> reads.to.prelim.reconstructed.fasta.sam

echo "(: Creating a map of observed nucleotides for the preliminary genome assembly."
awk '($0!~/@/)&&($2==0){print}' reads.to.prelim.reconstructed.fasta.sam | awk -v mapq_threshold="$mapq_threshold" -f ${map_script} > map.txt    #additional filtering is needed to reduce the amount of spurious alignments when doing a single segment: maybe primary only, forward strand only or both? Do liftover instead of realignment? If ends up like that should do it on chimp, not on preliminary assembly.
sort -k1,1 -k2,2n map.txt > map.txt.sort
mv map.txt.sort map.txt
sort -k1,1 -k2,2n -k3,3n -k4,4 indels.bed > indels.bed.sort
mv indels.bed.sort indels.bed

echo "(: Listing candidate regions for reassembly."
awk -v sig_threshold=2 -f ../../Anonamos-bin-160305/list-active-regions.awk indels.bed > regions.for.reassembly.txt

echo "(: Reassembling suspicious regions."

rm -f reassembled.indels.bed
cat regions.for.reassembly.txt | parallel -j10 --will-cite "awk -v region={.} -f ${extract_script} reads.to.prelim.reconstructed.fasta.sam | awk -v region={.} -f ${reassemble_script} >> reassembled.indels.bed"

echo "(: Reconstructing polished fasta."
sort -k1,1 -k2,2n -k3,3n -k4,4 reassembled.indels.bed > reassembled.indels.bed.sort
mv reassembled.indels.bed.sort reassembled.indels.bed
awk -v consensus_threshold="$consensus_threshold" -f ${reduce_script} reassembled.indels.bed map.txt > reconstructed.fasta

##
##
##
##
##
##
##
