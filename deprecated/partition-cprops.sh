#!/bin/bash
#### Description: Script to break the reference into more or less equally sized chunks to be processed in parallel. May involve breaking larger reference contigs or scaffolds into pieces. Uses greedy approximation to do multiway paritition.
#### Input: Original reference file cprops and (optional) number of chunks to break the reference into. The default number of chunks is 25.
#### Output: A set of files (N<=jobcount) in the working directory with scaffold names and start and end positions.
#### Dependency: ?
#### Written by: Olga Dudchenko - olga.dudchenko@bcm.edu on 06/11/2016

cprops=$1		# Reference cprops file
jobcount=$2		# Expected job count
reference=$3

if [ -z $jobcount ]; then
	jobcount=25	# Default job count
fi

cmd="awk '{count+=\$3}END{print count}' ""$cprops"
totlen=`eval $cmd`	# Scaffold reference length
approx=$((totlen / jobcount))	# Approximate max length of an individual contiguous piece

if [ -f reference.split.cprops ]; then
	rm reference.split.cprops	# Cleanup before starting
fi

while read chrname bytecount chrlen
do
#	split=$((chrlen  / (chrlen / approx + 1) + 1))
	split=$approx
	start=1
	p=0
	while true; do
		((p++))
		end=$((start + split-1))
		if [ "$end" -ge "$chrlen" ]; then
			end=$chrlen
			echo $chrname"_"$start"_"$end" "$bytecount" "$((end-start+1))" "$p >> reference.split.cprops
			break
        else
			echo $chrname"_"$start"_"$end" "$bytecount" "$((end-start+1))" "$p >> reference.split.cprops
			byteadd=$(tail -c +${bytecount} $reference | awk -v splitpos="$((end-start+1))" '{counter+=length; if(counter>splitpos){print NR-1; exit}}')
			let bytecount=$bytecount+$((end-start+1))+$byteadd
			let start=$start+$split
		fi
	done
done<$cprops

# Perform greedy heuristic multi-way number partitioning
find . -type f -name 'job.*.cprops' -delete 	# Cleanup before starting
declare -a len=( $(for i in {1..$jobcount}; do echo 0; done) )
while read fragname fragbyte fraglen fragno
do
	minset=1
	for (( i=2; i <= $jobcount; i++ ))
    do
		if [[ ${len[i]} -lt ${len[minset]} ]]; then
			minset=$i
		fi
    done
	echo $fragname" "$fragbyte" "$fraglen" "$fragno >> "job.""$minset"".cprops"
	let len[$minset]=${len[minset]}+$fraglen
done< <(sort -k3,3nr reference.split.cprops)
