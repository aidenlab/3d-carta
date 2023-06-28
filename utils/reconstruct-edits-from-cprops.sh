#!/bin/bash
# NOTE: Relies on standard labels :::fragment_ and :::debris. 

pipeline=`cd "$( dirname $0)" && cd .. && pwd`

current_cprops=$1

# 1) reconstruct original cprops file
awk '{split($1,a,":::fragment_||:::debris")}(a[1]!=a_prev){if(NR!=1){counter++; print a_prev, counter, len}; a_prev=a[1]; len=$3; next}{len+=$3}END{counter++; print a_prev, counter, len}' ${current_cprops} > original.cprops

# 2) reconstruct all the edits
awk 'BEGIN{OFS="\t"; print "chr1", "x1", "x2", "chr2", "y1", "y2", "color", "id", "X1", "X2", "Y1", "Y2"}$1~/:::debris/{print $1, 0, $3, $1, 0, $3, "0,0,0", "debris", 0, $3, 0, $3}' ${current_cprops} | awk -f ${pipeline}/lift/lift-input-annotations-to-asm-annotations.awk ${current_cprops} <(awk '{print $2}' ${current_cprops}) - | awk -f ${pipeline}/lift/lift-asm-annotations-to-input-annotations.awk original.cprops <(awk '{print $2}' original.cprops) - 

rm original.cprops
# TODO: redo everything in terms of .assembly files

