#!/bin/bash
fasta_file=$1
chr=$2

# DEPENDENCIES:
pipeline=`cd "$( dirname $0)" && cd .. && pwd`
basename=`basename $fasta_file .fna`
basename=`basename $basename .fasta`
basename=`basename $basename .fa`

awk -v chr=${chr} '$0~/>/{if(test){exit};test=0}($1==">"chr){test=1}test{print}' ${fasta_file} > ${basename}.${chr}.fasta

