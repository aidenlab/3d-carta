#### Description: Script to index the master sam file for quick access of downstream scripts. The reference is annotated with a given resolution (by default every 1Mb of every reference fragment is annotated).
#### Input: master.sam file
#### Parameters: bin zise (bin, optional).
#### Output: master sam index stdout in format "ref_frag_name bin_start_pos bin_start_byte".
#### Written by: Olga Dudchenko - olga.dudchenko@bcm.edu on 06/14/2016
BEGIN{
	if (bin==0)
		bin = 1000000	# annotate every megabase by default
}
{
#	print $0
	c+=length+1	# byte track

	if ($3!=chrname || $4 >= compar)
	{
		tmp = bin*int(($4-1)/bin) + 1
		print $3, tmp, c-length
		chrname = $3
		compar = tmp + bin
	}
}
END{
	print "EOF", "-", c+1
}
