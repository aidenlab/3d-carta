## Tentatively we list features from unaccounted for list such that density in a given window is larger than?

## default window is 50bp with 25bp shift? 

## this version requires unaccounted for list to be sorted!

# read in unaccounted
FILENAME==ARGV[1]{

	# test with just indels
#	if(NF==2){next}

	unaccounted_pos[FNR]=$1
	unaccounted_cov[FNR]=$NF
	counter=FNR
	next
}
# read in master map
{
	first+=$1+$2+$3+$4+$5-1
	second+=$1+$2+$3+$4+$5-1
}
(FNR+25)%50==0{
	tmp=0
	while(unaccounted_pos[first_counter]<=FNR && first_counter<=counter){
		tmp+=unaccounted_cov[first_counter]
		first_counter++
	}
	if (tmp>first/100 && first>50){
		print FNR-50, FNR
		print first>"coverage_distrib.txt"
	}
#	if(FNR%100000==0){print "I am here: ", FNR, tmp, first/50 > "/dev/stderr"}
	first=0
}
FNR%50==0{
	tmp=0
	while(unaccounted_pos[second_counter]<=FNR && second_counter<=counter){
		tmp+=unaccounted_cov[second_counter]
		second_counter++
	}
	if (tmp>second/100 && second>50){
		print FNR-50, FNR
		print first>"coverage_distrib.txt"

	}
#	if(FNR%100000==0){print "I am here: ", FNR, tmp, second/50 > "/dev/stderr"}
	second=0
}
END{
	##TODO: what to do about the last bit?
}