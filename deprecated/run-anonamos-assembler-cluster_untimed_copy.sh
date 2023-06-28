#!/bin/bash -x

USAGE="
*****************************************************
This is a wrapper for a the Anonamos comparative genome assembler which analyzes the alignment data of the dna sequences derived from species A to a genome of closely related species B to create ref-assisted contigs for species A. (This is a short workflow version of Anonamos with no provisional genome creation and a single map-reduce procedure.) The input to the assembler are single-end SAM alignment files (note that performance might increase if the files are preprocessed: split at chimeric junctions, dedupped, filtered for overlaps).

	The first step (-S sort) the sam file(s) are merge-sorted into a master sorted SAM file. The master SAM is then split into manageable chunks. The split size is adjusted dynamically to result in a predefined number of files determined by jobcount parameter (unprompted).
	In the second step (-S map) a per base nucleotide map is constructed for each scaffold in the reference assembly. Indels and read breaks positions with respect to the referene are also recorded. The step includes merging of the chunk map results. The result is a master map file.
	The third (optional?) step (-S addref) is adding reference nucleotides to the map and reference gap data.
	The forth step (-S reduce)is consolidation of preliminary fasta. Should probably include some reassembly but right now is just heuristics.

Version date: July 19, 2016.

Usage: ./run-anonamos-assembler.sh [-h] [-c min_coverage] [-q min_mapq] [-r chr:start-end] [-R reference_fasta] path_to_sam_file(s)

ARGUMENTS:
path_to_sam_file        Path to sam file describing read alignment to sequence used as assisting reference

OPTIONS:
-h                      Shows this help
-c min_coverage         Minimum read coverage required for the reference fragment to participate in genome reconstruction [default is 3]
-q min_mapq             Minimum mapping quality of the read to be considered in the recontruction process [default is 1]
-R reference_fasta		Path to genome fasta used to generate the sam files.

-S stage				Start from a particular stage. Can be split, map, addref, consol etc.
-r chr:start-end        Reconstruct fasta from a limited genomic region of a reference. The option can be an explicit listing of the chromosome with genomic start and end location in 1-base coordinate system or just the name of a chromosome for whole-chromosome reconstruction. This option is used for filtering out the sam file [default is empty for no filtering]. Not functional at the moment.
-x cluster				Cluster name (uger or slurm).

Uses map-sam-to-pos-array.awk and reduce-map-into-fasta.awk that should be in the same folder as this wrapper script.
*****************************************************
"
###################### SET DEFAULTS #############################

shopt -s extglob # enable pathname expansions

# problem-specific
min_coverage=3
min_mapq=1
skip_prep=false
skip_sort=false
skip_map=false
skip_addref=false
skip_reduce=false

# organizational
topDir=$(pwd)
debugDir=${topDir}"/debug"
samDir=${topDir}"/sam"
sortDir=${topDir}"/sam"
tmpDir=${topDir}"/tmp"
splitDir=${topDir}"/map"
mapDir=${topDir}"/map"

#prelimDir=${topDir}"/prelim"
#declare -A dependmap
#declare -A dependaddref
#declare -A dependreduce

###################### SCRIPT SHORTCUTS ############################

##	TODO: review and purge

basepath=`cd "$( dirname $0)" && pwd`
# preparatory
index_fasta="$basepath""/index-fasta.awk"
# sort & index
#index_sam="$basepath""/index-sam.awk"
index_parallel="$basepath""/index-parallel.awk"
# map nucleotides
map_script="$basepath""/map-script.sh"
map_awk_script="$basepath""/map-sam-to-pos-array-po.awk"
add_ref_script="$basepath""/add-ref-script.sh"
add_ref_awk_script="$basepath""/supplement-with-reference-po.awk"

# didn't get here yet
reduce_script="$basepath""/reduce-script.sh"
reduce_awk_script="$basepath""/reduce-map-into-fasta-stringent-po.awk"
flag_script="$basepath""/list-active-regions.awk"
extract_script="$basepath""/cut-reads-from-sam.awk"
reassemble_script="$basepath""/de-bruijn-assemble-DFS-both-strands.awk"
break_script="$basepath""/break-scaffold-into-contigs.awk"
local_fasta="$basepath""/cut-seq-from-fasta.awk"
map_gaps="$basepath""/make-gap-bed-for-local-ref.awk"
merge_maps="$basepath""/merge-maps.awk"
merge_indels="$basepath""/merge-indels.awk"
merge_skips="$basepath""/merge-skips.awk"

print_contig_N50="$basepath""/print-contig-n50.awk"
txt_to_bed_script="$basepath""/indel-txt-to-bed.awk"
prep_for_igv="$basepath""/prep-for-igv.sh"

# deprecated
#partition_cprops="$basepath""/partition-cprops.sh"


# TODO: Check that expected scripts are in the basepath folder
#if [ ! -f $map_script ] || [ ! -f $reduce_script ] ; then
#    echo "): Relevant scripts not found in bin folder, exiting!"
#    exit 1
#else
#    echo "(: Check succeeded."
#fi

###################### HANDLE OPTIONS ###############################

while getopts "c:q:S:r:x:h" opt; do
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
	S)	stage=$OPTARG
	;;
	r)  re1='^[A-Za-z0-9]+$'; re2='^[A-Za-z0-9]+\:[0-9]+\-[0-9]+$'
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
	x)	cluster=$OPTARG
	;;
	*) echo "$USAGE"
		exit 1
	;;
esac
done

shift $(( OPTIND-1 ))

if [ -z $cluster ]; then
	echo "Please specify cluster. Exiting!"
	exit 1
fi

# set cluster-specific defaults

case $cluster in
	uger)
		reference="/broad/aidenlab2/olga/panTro4.fa"	# temporary
		usePath=/broad/software/scripts/useuse
		load_cluster="use UGER"
		load_coreutils="use Coreutils"
		load_parallel="use .parallel-20140722"
		short_queue="short"
		long_queue="long"
		jobcount=50
		source $usePath
		$load_cluster
	;;
	slurm)
		reference="/work/ea14/juicer/references/genome_collection/temp/panTro4.fa"	# temporary
		short_queue="commons"
		long_queue="commons"
#		splitsize=22500000
		jobcount=10
		export PATH=/projects/ea14/parallel/bin:$PATH
	;;
	*)
		echo ":( Unrecognized cluster name. Exiting!"
		exit 1
esac

###################### TEMPORARY: FOR QC ###########################
panTro4=${reference}
#hg38="/poscratch/aidenlab/olga/assisted_asm/references/hg38.fa"
####################################################################

cov_threshold=$((min_coverage-1)) # ge rather than gt criteria
mapq_threshold=$((min_mapq-1)) # ge rather than gt criteria

if [ ! -z "$stage" ]; then
	case $stage in
		reduce)	skip_addref=true
		;&
		addref)		skip_map=true
		;&
		map)		skip_sort=true
		;&
		sort)		skip_prep=true
		;;
		*)			echo "$USAGE"
					exit 1
	esac
fi

####################### CHECK INPUT ################################

echo "(: Checking input."

if [ "$skip_prep" = false ]; then
	if [ $# -lt 1 ]; then
		echo "): The set of arguments provided is incomplete. Please double-check your input, exiting!"
		echo "$USAGE"
		exit 1
	fi
elif [ "$skip_prep" = true ]; then
	if [ ! -d "$samDir" ] || [ ! -d "$sortDir" ] || [ ! -d "$mapDir" ] || [ ! -d "$splitDir" ]; then
		echo "): Can't find expected folders $sortDir and $mapDir! Please make sure you have preprocessed data or start from scratch."
		echo "$USAGE"
		exit 1
	fi
fi

groupname="a$(date +%s)"

############### ORGANIZE FOLDER IF RUNNING FROM SCRATCH ##############

if [ "$skip_prep" = false ]; then

	echo "(: Organizing workspace."

	if [ ! -d $debugDir ] && [ ! -d $samDir ] && [ ! -d $sortDir ] && [ ! -d $splitDir ] && [ ! -d $prelimDir ]; then
		mkdir "$debugDir"; chmod 777 "$debugDir";
		mkdir "$samDir"; chmod 777 "$samDir";
		mkdir "$sortDir"; chmod 777 "$sortDir";
		mkdir "$splitDir"; chmod 777 "$splitDir";
		mkdir "$mapDir"; chmod 777 "$mapDir";
		mkdir "$prelimDir"; chmod 777 "$prelimDir";
		mkdir "$tmpDir"; chmod 777 "$tmpDir";
	else
		echo "Warning: folders with clashing names found in working directory. Note that some files might get overwritten!"
	fi

	echo $@ | xargs -n 1 echo | xargs -I % ln -sf % $samDir
	touch reference.cprops

#	touch "$splitDir/reference.split.cprops"
#	seq 1 $jobcount | xargs -P8 -I % touch "$splitDir/job."%".cprops"

#	cmd="awk -f ${index_fasta} ${reference} > ${topDir}/reference.cprops; cd ${splitDir}; bash ${partition_cprops} ${topDir}/reference.cprops ${jobcount} ${reference} && echo \":) Done with preparatory reference analysis.\""

	cmd="echo \":) Starting preparatory reference analysis.\" && date; awk -f ${index_fasta} ${reference} > ${topDir}/reference.cprops && echo \":) Done with preparatory reference analysis.\" && date;"

	case $cluster in
		uger)
			qsub -o ${debugDir}/uger.out -j y -q $short_queue -r y -N "${groupname}_prep" <<-PREP
			eval $cmd
			PREP
		;;
		slurm)
			jid=`sbatch <<- PREP | egrep -o -e "\b[0-9]+$"
			#!/bin/bash -l
			#SBATCH -p $short_queue
			#SBATCH -t 100
			#SBATCH -c 1
			#SBATCH --ntasks=1
			#SBATCH --mem-per-cpu=2G
			#SBATCH -o $debugDir/prep-%j.out
			#SBATCH -e $debugDir/prep-%j.err
			#SBATCH -J "${groupname}_prep"
			$cmd
			PREP`
			dependprep="afterok:$jid"
	esac

fi

##################### SORT AND SPLIT SAM ####################

cd "$topDir"

if [ "$skip_sort" = false ]; then

	echo "(: Launching jobs to split, sort and index sam file(s)."

	cmd="echo \":) Starting sorting input sam files into a master sam file.\" && date; linecount=\$(LC_ALL=C cat ${samDir}/* | awk '!/@/' | sort -T ${tmpDir} -k3,3 -k 4,4n --parallel=8 -S 20G -s | tee ${sortDir}/master.sam | wc -l) && split -a 3 -l \$(( \${linecount}/${jobcount} + 1 )) -d --additional-suffix=.sam ${sortDir}/master.sam ${splitDir}/job_ && echo \":) Finished sorting input sam files into a master sam file.\" && date;"

	case $cluster in
		uger)
			qsub -o ${debugDir}/uger.out -j y -q $long_queue -r y -N "${groupname}_sort" -l m_mem_free=25g <<-SORT
			source $usePath; $load_coreutils; eval $cmd
			SORT
		;;
		slurm)
			jid=`sbatch <<- SORT | egrep -o -e "\b[0-9]+$"
			#!/bin/bash -l
			#SBATCH -p $long_queue
			#SBATCH -t 1440
			#SBATCH -c 2
			#SBATCH --ntasks=1
			#SBATCH --mem-per-cpu=25G
			#SBATCH -o $debugDir/sort-%j.out
			#SBATCH -e $debugDir/sort-%j.err
			#SBATCH -J "${groupname}_sort"
			$cmd
			SORT`
			dependsort="afterok:$jid"
	esac

fi

######################### MAPPING STEP ################################

cd ${topDir}

echo "(: Launching jobs to make a per base map of alignments."

if [ "$skip_map" = false ]; then
	echo "(: Launching jobs to map sam file(s)."

	if [ "$skip_sort" = false ]; then
		sbatch_wait="#SBATCH -d $dependsort"
		uger_wait="-hold_jid ${groupname}_sort"
	else
		sbatch_wait=""
		uger_wait=""
	fi

	seq -f "%03g" 0 $(( jobcount - 1 )) | xargs -P8 -I % touch "$splitDir/job_"%".sam"
	doawk() {
		FILE=$1
		awk 'NR==1{start=$3" "$4}END{print start, $3, $4}' "$FILE" > "$FILE".index
	}
	export -f doawk

	for (( i=1; i <= $jobcount; i++ ))
	do
		filename=${splitDir}"/job_"$(printf "%03g" $(( $i -1)))".sam"

#		cmd="echo \"(: Starting mapping $i file.\" && date; bash ${map_script} $i ${mapq_threshold} && echo \"(: Finished mapping $i file.\" && date;"

#		cmd="echo \"(: Starting mapping $filename file.\" && date; parallel -a ${filename} --pipepart --block 10M awk \'NR==1{chr=\\\$3\; pos=\\\$4}END{print chr, pos, \\\$3, \\\$4}\' > ${filename}.index && echo \"(: Finished mapping ${filename} file.\" && date;"

		cmd="echo \"(: Starting mapping $filename file.\" && date; parallel -a ${filename} --pipepart -j8 -k --block 500M doawk && parallel -a ${filename} --pipepart -k --block 500M awk -f ${map_awk_script} && echo \"(: Finished mapping ${filename} file.\" && date;"

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

#				jid=`sbatch <<- MAP | egrep -o -e "\b[0-9]+$"
#				#!/bin/bash -l
#				#SBATCH -p $queue
#				#SBATCH -o $debugDir/map-%j.out
#				#SBATCH -e $debugDir/map-%j.err
#				#SBATCH -t 1440
#				#SBATCH -c 1
#				#SBATCH --ntasks=1
#				#SBATCH -J "${groupname}_map_${jname}"
#				${sbatch_wait}
#				srun --ntasks=1 ${cmd}
#				MAP`
#				dependmap[$i]="afterok:$jid"
		esac

	done
fi

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
#				cmd="bash ${reduce_script} ${chrname} ${cov_threshold}; $cmd"
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
		srun bash ${reduce_script} ${debugDir}"/cprops_partition.$i.txt" ${cov_threshold}
		reduce`
		dependreduce[$i]="afterok:$jid"
	done
fi

exit 0

#			awk -v cov_threshold="$cov_threshold" -f ${reduce_script} ${mapDir}/{$chrname}.indel.txt ${mapDir}/${chrname}.map.txt | awk -f ${break_script} > ${mapDir}/${chrname}.prelim.fasta



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
awk -v cov_threshold="$cov_threshold" -f ${reduce_script} indels.txt map.txt > prelim.reconstructed.fasta
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
awk -v cov_threshold="$cov_threshold" -f ${reduce_script} reassembled.indels.txt map.txt > polished.reconstructed.fasta
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
awk -v cov_threshold="$cov_threshold" -f ${reduce_script} indels.bed map.txt > prelim.reconstructed.fasta

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
awk -v cov_threshold="$cov_threshold" -f ${reduce_script} reassembled.indels.bed map.txt > reconstructed.fasta

##
##
##
##
##
##
##
