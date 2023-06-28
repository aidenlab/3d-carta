#### Description: Script to map sam file or sam file fragments.
#### Input: *cprops and *.sam file. Note that the files are assumed to be sorted!
#### Parameters: mapq_threshold (default is 0) and filename. Reads with mapping quality lesser than mapq_threshold will be discarded.
#### Output: map data in format "A C G T N". The first line of the file is @(position of the first match in file - 1) which is used during merge. 
#### Written by: Olga Dudchenko - olga.dudchenko@bcm.edu, version dating 10/18/2017
#### Version with 0.5 weight to mapq0..?
BEGIN{
	split("", nset)
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
# parse cprops file for global positioning
FILENAME==ARGV[1]{
	chromshift[$1] = total
	total+=$3
	next
}
($5<mapq_threshold){next}
($5==mapq_threshold){w=1}
($5>mapq_threshold){w=1}
{
	if ($3!=chrname)	# start of file or next chromosome
	{
				
		if (globalshift==""){
			# print header for global positioning during merge
			globalshift = chromshift[$3] + $4 - 1
			print "@"(globalshift) > filename".map.txt"
		} else {
			# print leftovers from previous chromosome
			n = asorti(nset, sortpos, "@ind_num_asc")
			for (i=1; i<=n; i++)
			{
				while (globalshift != chromshift[chrname] + sortpos[i]+poshift-1){
					print "" > filename".map.txt"
					globalshift++
				}
				print nset[sortpos[i]] > filename".map.txt"
				delete nset[sortpos[i]]
				globalshift++
			}
			
		
			# print leftovers of chromosome
			while (globalshift!=chromshift[$3])
			{
				print "" > filename".map.txt"
				globalshift++
			}
		}

		chrname=$3
		poshift=$4
		leftmost=$4
		rightmost=0
		

	}
		
	if ($4 > poshift)	# dump data when moving to next interval
	{
		n = asorti(nset, sortpos, "@ind_num_asc")

		for (i=1; i<=n; i++)
		{			
			if (sortpos[i]+poshift>=$4)
				break	

			while (globalshift != chromshift[chrname] + sortpos[i] + poshift - 1){
				print "" > filename".map.txt"
				globalshift++
			}

			print nset[sortpos[i]] > filename".map.txt"
			delete nset[sortpos[i]]
			globalshift++
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
				clip[chromshift[$3]+pos+poshift]+=1*w
#				read_id[chromshift[$3]+pos+poshift]=read_id[chromshift[$3]+pos+poshift]" "$1

				s += cigar[k]
			}
			else if (seps[k] == "H")
			{
				clip[chromshift[$3]+pos+poshift]+=1*w
#				read_id[chromshift[$3]+pos+poshift]=read_id[chromshift[$3]+pos+poshift]" H_"$1
			
				continue
			}			
			else if (seps[k] == "D")
			{
				tmp=0
				tmpstr=""
				while (tmp<cigar[k])
				{
					tmpstr=tmpstr"-"
					tmp++
				}
				indel[(chromshift[$3]+pos+poshift)" "tmpstr]+=1*w
				pos = pos + cigar[k]
			}
			else if (seps[k] == "I")
			{
				indel[(chromshift[$3]+pos+poshift)" "substr($10, s, cigar[k])]+=1*w
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
						nucl[1]+=1*w
					else if (call == "C")
						nucl[2]+=1*w
					else if (call == "G")
						nucl[3]+=1*w
					else if (call == "T")
						nucl[4]+=1*w
					else
						nucl[5]+=1*w
					nset[pos]=nucl[1]" "nucl[2]" "nucl[3]" "nucl[4]" "nucl[5]
					pos++
					s++
				}
				if ( pos-1 > rightmost)
					rightmost=poshift+pos-1
			}
			
		}
}
END{

	n = asorti(nset, sortpos, "@ind_num_asc")
			
	for (i=1; i<=n; i++)
	{
		while (globalshift != chromshift[chrname] + sortpos[i]+poshift -1){
			print "" > filename".map.txt"
			globalshift++
		}
		print nset[sortpos[i]] > filename".map.txt"
		globalshift++
	}

	for (i in clip)
	{
#		print i, 1, skip[i], read_ids[i] > filename".feature.txt"
		print i, clip[i] > filename".clip.txt"
	}
	for (i in indel)
	{
		print i, indel[i] > filename".indel.txt"
		
# 		split(i, tmp)
# 		if (tmp[3]~/^[0-9]+$/)
# 		{
# 			print tmp[1], tmp[2], tmp[3], indel[i] > filename".feature.txt"
# 			#>> "indels.txt"
# 		}
# 		else
# 		{
# 			print tmp[1], tmp[2], 1, indel[i], tmp[3] > filename".feature.txt"
# 			#>> "indels.txt"
# 		}
	}
}

