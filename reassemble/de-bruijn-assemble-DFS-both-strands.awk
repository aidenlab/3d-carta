BEGIN{
	if(!k){k=50}	# default kmer size
	min_support=1	# minimal weight required to accept edge
	RE="GATC"
	
	nblock_threshold=25  # mock string for no traversal or circular graph
	for (i=1; i<=nblock_threshold; i++)
		nstring=nstring"N"
}

# recursive depth-first search through kmer graph
function DFS(path){
	if (substr(path,length(path)-k+1) == end)
		traversal[path]=weight[path]
	else
	{
		for (i in counter)
		{
			if (substr(i,1,k)==substr(path,length(path)-k+1))
			{
					weight[path""substr(i,length(i))]=weight[path]+counter[i]

				if (path ~ substr(i,2))
				{
					message="Cycle warning!"
					cycle_flag[substr(i,2)]=1
########????					delete traversal
					continue
				}
				else
				{
#					weight[path""substr(i,length(i))]=weight[path]+counter[i]
#print path, weight[path] > "/dev/stderr"
					DFS(path""substr(i,length(i)))
				}
				
			}
			
		}
	}
}
function complement(nt)
{
	if (nt=="A")
		return "T"
	else if (nt=="T")
		return "A"
	else if (nt=="G")
		return "C"
	else if (nt=="C")
		return "G"
	else
		return nt
}
function reverse(sequence,		revstr, i, n){
	n=split(sequence,nt,"")
	revstr=complement(nt[n])
	for (i=n-1; i>=1; i--)
	{
		revstr=revstr""complement(nt[i])
	}
	
	return revstr
}
	


# read in the reads
{
	if ($0~/\^/)	# read overlapping start of interval
	{
		$0=substr($0,2)
		source_read[FNR]=1
	}
	if ($0~/\$/)		# read overlapping end of interval
	{
		$0=substr($0,1,length($0)-1)
		sink_read[FNR]=1
	}

	read[FNR]=$0
}
END{

# break reads into kmers given k

	for (readid in read)
	{
		if (length(read[readid]) < k)
			continue
	
		if (readid in source_read)	# read overlapping start of interval
		{
			kmer=substr(read[readid],1,k)
			source[kmer]+=1
			if (source[kmer]>=maxsource)
			{
				start=kmer
				maxsource=source[kmer]
			}
		}
		if (readid in sink_read)		# read overlapping end of interval
		{
			kmer=substr(read[readid],length(read[readid])-k+1)
			sink[kmer]+=1
			if (sink[kmer]>=maxsink)
			{
				end=kmer
				maxsink=sink[kmer]
			}
		}

		for (i=2; i<=length(read[readid])-k+1; i++)
		{
			edge=substr(read[readid],i-1,k+1)
			counter[edge]+=1	# count edges
			if (! anchored)		# don't bother keeping track of reverse strand if reassembling from aligned reads, reduce complexity
			{
				counter[reverse(edge)]+=1
			}
		}
	}


	for (i in counter)	# remove poorly supported edges
	{
		teststr="."RE"."
		if ((counter[i] < min_support)&&(i!~teststr))	# leave precious edges spanning RE sites (alt - increase weight of reads upfront)
		{
			delete counter[i]
		}
	}

	str=start

	DFS(str)

	if ((length(traversal)==0)|| message)
	{
		if (!message)
			message="No contiguous paths found!"
		
		max=0
		maxpath=""
		for (p in weight)
		{
			if (weight[p] > max)
			{
				maxpath=p
				max=weight[p]
			}
		}

		if (substr(maxpath, length(maxpath)-k+1)!=end) # most likely not a traversing path - only if cycle is part of buble
		{
			outstring = maxpath
			delete weight
			max=0
			maxpath=""

			str=reverse(end)
			end=reverse(start)

			if (anchored)	# explicitely convert read orientation for anchored processing
			{
				for (edge in counter)
				{
					rev_counter[reverse(edge)]=counter[edge]
				}
				delete counter
				for (edge in rev_counter)
				{
					counter[edge]=rev_counter[edge]
				}
				delete rev_counter
			}

			DFS(str)
			for (p in weight)
			{
				if (weight[p] > max)
				{
					maxpath=p
					max=weight[p]
				}
			}
			outstring = outstring""nstring""reverse(maxpath)
			maxpath=outstring
		}

	}
	else
	{
		for (i in traversal)
		{
			if (traversal[i]>max)
			{
				max=traversal[i]
				maxpath=i
			}
		}
	}

	split(region, tmp)
	print maxpath, message

}
