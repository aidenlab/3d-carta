#!/bin/bash

fasta_file=$1
block=$2
pipeline=`cd "$( dirname $0)" && cd .. && pwd`

[ -z ${block} ] && block=100000

basename=`basename $fasta_file .fna`
basename=`basename $basename .fasta`
basename=`basename $basename .fa`

awk -v block=${block} -f ${pipeline}/utils/wrap-fasta-sequence.awk ${fasta_file} | awk '$0~/>/{seq=$1; counter=0;; next}{counter++; print seq"_"counter; print}' > $basename.split_for_alignment.fasta
