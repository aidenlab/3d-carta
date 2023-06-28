## Tentatively we list features from unaccounted for list such that density in a given window is larger than?

## default window is 50bp with 25bp shift? 

## this version requires unaccounted for list to be sorted!

# read in unaccounted
FILENAME==ARGV[1]{

	# test with just indels
#	if(NF==2){next}

#	unaccounted_pos[FNR]=$1
	unaccounted_cov[$1]+=$NF

# 	special handling for long deletions	? Or perhaps explicitly check when list of active regions is available and widen if overlaps - maybe both necessary
	if ( ($2~/-/) && int($1/25)!=int(($1+length($2))/25))
	{
		unaccounter_cov[$1+length($2)]+=$NF
	}
	counter=FNR
	next
}
# read in master map ## TODO: put in condition not to bother counting unless necessary
{
	cov=$1+$2+$3+$4+$5-1
	if (FNR in unaccounted_cov)
	{
		cov < prevcov ? usecov=prevcov : usecov=cov
		if(usecov<=1){next}	# avoid regions where we can't do much?
		first+=unaccounted_cov[FNR]/usecov
		second+=unaccounted_cov[FNR]/usecov
	}
	prevcov=cov
}
(FNR+25)%50==0 && first{
	if (first>0.5){
		print FNR-50, FNR
	}
	first=0
}
FNR%50==0 && second{
	if (second>0.5){
		print FNR-50, FNR
	}
	second=0
}
END{
	##TODO: what to do about the last bit?
}