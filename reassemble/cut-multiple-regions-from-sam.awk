#note that includes secondary alignments - after cutting at RE sites should probably ignore secondary alignments?
BEGIN{
#	split(region, tmp)
#	chr=tmp[1]
#	start=tmp[2]
#	end=tmp[3]
#	split(region,tmp,":")
#	chr=tmp[1]
#	split(tmp[2],tmpp,"-")
#	start=tmpp[1]
#	end=tmpp[2]
	mapq_threshold=0
	
#	readlength=150 - to empirically skip regions
	
	active_region_counter=1

	chromshifter=0
} 
function shift(ARRAY,	i, n)
{
	for (i in ARRAY)
		n++
	for (i=1; i<=n-1; i++)
	{
		ARRAY[i] = ARRAY[i+1]
	}
	delete ARRAY[i]
}
# read in sorted cprops for global positioning
FILENAME==ARGV[1]{
	chromshift[$1]=chromshifter
	chromshifter+=$3
	next
}
# read in active regions in global positioning
## TODO: what if none? this might break, handle.
FILENAME==ARGV[2]{
	start_active_region[FNR]=$1
	end_active_region[FNR]=$2
	
	start=start_active_region[active_region_counter]
	end=end_active_region[active_region_counter]

	next
}
# print first element projected patch position
FNR==1{
	print "--- "
	print start, end
}
# don't bother if mapping quality below threshold
($5<mapq_threshold){next}
# convert alignment coordinate to global coordinate system
{
	$4+=chromshift[$3]
## TODO: skip reads that empirically won't contribute
}
{

	# Ugly workaround for old gawk, split does not spit out separators, have to split separately and shift
	arrlength = split($6, cigar, "[MIDNSHPX=]") # Only deal with M, I, D, H and S.
	split($6, seps, "[0-9]+")		
	shift(seps)
	
	s=1 #string position of first match
	pos=$4 #reference position of first match

	dels=0
	ins=0
	if (seps[1]=="S")
		s+=cigar[1]

	for (k in seps)	# count ins and dels to accurately estimate read end position
	{
		if (seps[k]=="D")
			dels+=cigar[k]
		else if (seps[k]=="I")
			ins+=cigar[k]
	}
	
	while (pos>end){	# this does not take all if multiple regions overlap within the same read but keeps things sorted		
		if (active_region_counter==length(start_active_region)){exit}
		active_region_counter++
		start=start_active_region[active_region_counter]
		end=end_active_region[active_region_counter]
		print "--- "
		print start, end
	}
	
#	if (200<FNR&&FNR<250){print pos+length($10)-s+dels-ins, start, pos-s+1, end >"/dev/stderr"}


	if ((pos+length($10)-s+dels-ins >= start) && (pos-s+1 <= end))	# if there is a chance that you overlap with the region of interest
	{
	
#	print "I am here: "pos+length($10)-s+dels-ins, start, pos-s+1, end >"/dev/stderr"

		overlap=""
		is_source=0
		is_sink=0
		if (start < pos-s+1)
			substrstart = 1
		else
		{
			is_source=1
			substrstart = s-pos+start
		}

		if (end > pos+length($10)-s+dels-ins)
			substrlength = length($10)-substrstart+1		
		else
		{
			is_sink=1

			for (k = 1; k <= arrlength; k++)
			{
				if (pos >= end)
					break
				
				if (seps[k] == "S")
				{
					if (k!=1)
					{		
						s += cigar[k]
						pos += cigar[k]
					}
				}
				else if (seps[k] == "H")
				{
					continue
				}			
				else if (seps[k] == "D")
				{
					pos = pos + cigar[k]
				}
				else if (seps[k] == "I")
				{
					s += cigar[k]
				}
				else
				{
					pos += cigar[k]
					s += cigar[k]
				}
			}
			substrlength = s-substrstart+end-pos+1
		}
					
		overlap = substr($10, substrstart, substrlength)
		if (overlap!="")
		{	
			if (is_source)
				overlap="^"overlap
			if (is_sink)
				overlap=overlap"$"
			print overlap
		}
	}

}
