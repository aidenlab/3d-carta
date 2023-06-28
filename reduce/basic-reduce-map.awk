## Primitive reduce procedure which reports nucleotide corresponding to top value in the map array.
## awk -v c=consensus_threshold -f basic-reduce-map.awk <path-to-map-file>
## INPUT: map array file.
## OUTPUT: single-letter stdout in global coordinate system 
## OPTIONS: consensus_threshold. Default: 1.
## NOTE: systematic bias towards 'smaller' nucleotide characters (A<C<T<G) due to how max is defined.
## Written by: OD
## Version: 1605**?
BEGIN{
        nucl[1]="A"; nucl[2]="C"; nucl[3]="G"; nucl[4]="T"; nucl[5]="N"
        if (!length(c)){c=1}
}
{
	mmax=0;
	for(i=1 ;i<=4; i++){
		if ($i>mmax){
			mmax=$i
			emax=i
		}
	}
	if (mmax>c){print nucl[emax]}else{print "N"}
}
