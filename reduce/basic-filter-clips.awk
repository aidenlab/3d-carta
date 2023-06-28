## A version of filtering script to remove clip specific to Hi-C data, i.e. that are nearby a restriction site.
## awk -v t=tolerance -f remove-clips-from-restriction.awk <path-clip-list-in-global-coordinates> <path-to-unwrapped-fasta>
## Input: List of clips in global-position" "coverage formate and unwrapped fasta against which to check (typically used with original but can be applied to assisted reference). 
## Output: Filtered clips stdout in global-position" "coverage format.
## NOTE: unclear if using with primitive reduce is much better than using with vanilla assisted reference which could be done easier and faster. Another problem is that doing it via sequence precludes careful handling of heterozygous restriction sites. It is possible that we will switch to filtering during mapping.
## Written by: OD
## Version: 171120
BEGIN{
	len=length(seq)
	if(!length(t)){tolerance=1}else{tolerance=t}
}
# clips
FILENAME==ARGV[1]{
	cov[$1]=$2
	next
}
{
	if (length(tmp)==len+2*tolerance){
		tmp=substr(tmp,2)""$0
	}
	else
	{
		tmp=tmp""$0
	}
}
(tmp~seq){
	for(i=FNR-len-2*tolerance+1; i<=FNR; i++){
		delete cov[i]
	}
}
END{
	for (i in cov){
		print i, cov[i]
	}
}