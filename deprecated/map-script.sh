#!/bin/bash

chrnamefile=$1
mapq_threshold=$2
topDir=$(pwd)
unset index
declare -A index

basepath=`cd "$( dirname $0)" && pwd`
map_awk_script="$basepath""/map-sam-to-pos-array-po.awk"

sam_file=${topDir}"/sort/master.sam"
index_file=${topDir}"/sort/master.index"

while read tmpchr mbpos bytecount
do
	index["$tmpchr $mbpos"]=$bytecount
done <"$index_file"

awk '{print $1}' ${chrnamefile} | xargs -L1 -P 8 tail -c +${index[$chrname]} | awk -v filesuffix="${topDir}/map/$chrname" -v mapq_threshold=$mapq_threshold -f ${map_awk_script}


while read chrname skip
do
#	if [ -z ${index["$chrname"]} ]; then
#		index["$chrname"]=${index["EOF"]}
#	fi
	if [ ! -z ${index["$chrname"]} ]; then
		tail -c +${index[$chrname]} ${sam_file} | awk -v filesuffix="${topDir}/map/$chrname" -v mapq_threshold=$mapq_threshold -f ${map_awk_script}
		if [ -f ${out_suf}.indel.txt ]; then
			sort -k1,1n -k3,3n ${out_suf}.indel.txt > ${out_suf}.indel.sorted && mv ${out_suf}.indel.sorted ${out_suf}.indel.txt
		fi
		if [ -f ${out_suf}.skip.txt ]; then
			sort -k1,1n ${out_suf}.skip.txt > ${out_suf}.skip.sorted && mv ${out_suf}.skip.sorted ${out_suf}.skip.txt
		fi
	fi
done <${chrnamefile}


#cmd="tail -c +${index[$chrname]} ${sam_file} | awk -v filesuffix=\"$out_suf\" -v mapq_threshold=$mapq_threshold -f ${map_awk_script}"
#echo $cmd
#eval $cmd
#tail -c +${index[$chrname]} ${sam_file} | awk -v filesuffix="\"${out_suf}\"" -f ${map_awk_script}
