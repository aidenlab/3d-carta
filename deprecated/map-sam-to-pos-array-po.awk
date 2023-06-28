# note that this assumes file sorted both by chr and by position
BEGIN{
	split("", nset)
	sig_threshold=2 # TODO: deside if we want to filter out some of the indels and skips here
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
($5>mapq_threshold){

#	if (chrname=="")
#	{
#		chrname=$3
#		poshift=$4
#	}	
	if ($3!=chrname)	# start of file or next chromosome
	{
		print ">"$3
		chrname=$3
		poshift=$4
	}
		
	if ($4 > poshift)	# dump data when moving to next interval
	{
		n = asorti(nset, sortpos, "@ind_num_asc")

		for (i=1; i<=n; i++)
		{
			if (sortpos[i]+poshift>=$4)
				break
			print sortpos[i]+poshift, nset[sortpos[i]]
			delete nset[sortpos[i]]
		}
		for (j=i; j<=n; j++)
		{
			nset[sortpos[j]+poshift-$4]=nset[sortpos[j]]
			delete nset[sortpos[j]]
		}
		poshift=$4	
	}
	pos=0
	
	# Ugly workaround for old gawk, split does not spit out separators, have to split separately and shift
	arrlength = split($6, cigar, "[MIDNSHPX=]") # Only deal with M, I, D, H and S.
	split($6, seps, "[0-9]+")		
	shift(seps)

	# Extract_data_from_read()
	s=1 #string position tracker

	for (k = 1; k <= arrlength; k++) #process each cigar block appropriately
		{
			if (seps[k] == "S")
			{
				skip[$3" "(pos+poshift)]+=1
				read_ids[$3" "(pos+poshift)]=read_ids[$3" "(pos+poshift)]" "$1
				s += cigar[k]
			}
			else if (seps[k] == "H")
			{
				skip[$3" "(pos+poshift)]+=1
				read_ids[$3" "(pos+poshift)]=read_ids[$3" "(pos+poshift)]" H_"$1
				continue
			}			
			else if (seps[k] == "D")
			{
				indel[$3" "(pos+poshift)" "cigar[k]]+=1
				pos = pos + cigar[k]
			}
			else if (seps[k] == "I")
			{
				indel[$3" "(pos+poshift)" "substr($10, s, cigar[k])]+=1
				s += cigar[k]
			}	
			else # mostly seps[k] = "M" TODO: handle other seps in case occur
			{
				for (i = 1; i <= cigar[k]; i++)
				{
					if (pos in nset)
						split(nset[pos], nucl)
					else
						split("0 0 0 0 0", nucl)
					call=toupper(substr($10,s,1))
					if (call == "A")
						nucl[1]+=1
					else if (call == "C")
						nucl[2]+=1
					else if (call == "G")
						nucl[3]+=1
					else if (call == "T")
						nucl[4]+=1
					else
						nucl[5]+=1
					nset[pos]=nucl[1]" "nucl[2]" "nucl[3]" "nucl[4]" "nucl[5]
					pos++
					s++
				}
			}
			
		}
}
END{

	n = asorti(nset, sortpos, "@ind_num_asc")
	for (i=1; i<=n; i++)
	{
		print sortpos[i]+poshift, nset[sortpos[i]]
	}

	for (i in skip)
	{
		print i, 1, skip[i], read_ids[i] >> "skips.txt"
	}
	for (i in indel)
	{
		split(i, tmp)
		if (tmp[3]~/^[0-9]+$/)
			print tmp[1], tmp[2], tmp[3], indel[i] >> "indels.txt"
		else
			print tmp[1], tmp[2], 1, indel[i], tmp[3] >> "indels.txt"
	}

}

