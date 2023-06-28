#### Description: Script to index the master sam file for quick access downstream scripts.
#### Input: Path to reference.split.cprops and master.sam files.
#### Output: master sam index stdout in format "split_reference_frag_name start_byte".
#### Written by: Olga Dudchenko - olga.dudchenko@bcm.edu on 06/14/2016
{
	# read in the refrence.split.cprops file
	if (FILENAME==ARGV[1])
	{
		# split name to get reference data
		n=split($1,a,"_")
		chrname=a[1]
		for (i=2; i<=n-2; i++)
			chrname=chrname"_"a[i]
		fragcount[chrname]+=1		
		start[chrname" "fragcount[chrname]]=a[n-1]
		end[chrname" "fragcount[chrname]]=a[n]
		next
	}
}
{
	c+=length+1

	if ($3!=prev)
	{
		if (FNR>1)
		{
			print prev"_"start[prev" "counter]"_"end[prev" "counter], start_byte
		}
		start_byte = sprintf("%4i", c-length)
		counter = 1
	}

	prev=$3

	if (fragcount[prev] == 1)
	{
		next		
	}	
	
	if ((counter <= fragcount[prev] - 1) && ($4 >= start[prev" "(counter+1)]))
	{
		print prev"_"start[prev" "counter]"_"end[prev" "counter], start_byte
		start_byte = sprintf("%4i", c-length)
		counter += 1
	}
		
}
END{print prev"_"start[prev" "counter]"_"end[prev" "counter], start_byte; print "EOF", sprintf("%4i", c+1)}
