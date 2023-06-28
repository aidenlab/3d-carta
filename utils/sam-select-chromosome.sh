#!/bin/bash
sam_file=$1
chr=$2
parallel -a ${sam_file} --will-cite --pipe-part --block-size 1G -k "awk -v chr=${chr} '\$3==chr'" > `basename $1 .sam`.${chr}.sam
