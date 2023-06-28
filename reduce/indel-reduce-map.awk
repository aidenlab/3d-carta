## Reduce procedure that reports nucleotide corresponding to top value in the map array and accepts some of the indel calls, based on some empirical criteria.
## awk -v consensus_threshold=${consensus_threshold} -v shift=${shift} -v ignore_reference=${ignore_reference} -f indel-reduce-map.awk <indel_file> <clip_file> <extended_map_file>
## INPUT: indel and clip call files (w coverage, in global positioning system) and extended map file (reference file passed on together with master map to help distinguish reference calls from read calls).
## OUTPUT: single-letter stdout in global coordinate system and unaccounted indels (i.e. ones that were not incorporated based on empirical criteria) as stderr.
## OPTIONS: consensus_threshold (default: 1); shift (default: 0); ignore_reference (default: 0).
## NOTE: systematic bias towards 'smaller' nucleotide characters (A<C<T<G) due to how max is defined.
## Written by: OD
## Version: 171129
BEGIN{
        nucl[1]="A"; nucl[2]="C"; nucl[3]="G"; nucl[4]="T"; nucl[5]="N"
        
        # just in case make sure that options passed are meaningful and compatible. If not change to defaults
        if (!length(consensus_threshold)){
        	print ":( Wrong value passed for consensus_threshold. Setting parameter to default consensus_threshold=1" > "/dev/stderr"
        	consensus_threshold = 1
        }
        if (consensus_threshold==1 && !ignore_reference){
         	print ":( Cannot use reference when consensus_threshold is set to 1. Setting parameter to default ignore_reference=1" > "/dev/stderr"
        	ignore_reference=1
        }	   	
}
# read in the sorted indel file. Per position only one indel with highest number of calls is considered. Not full-proof, should reassemble?
FILENAME==ARGV[1]{
	if ($3>=cov[$1])  # note bias
	{	
		if (cov[$1]!=0){
			counter++
			alt[counter]=$1" "seq[$1]" "cov[$1]
		}
		seq[$1]=$2
		cov[$1]=$3
	} else {
		counter++
		alt[counter]=$1" "$2" "$3
	}
	next
}
# let's try this: read in the sorted clip file
FILENAME==ARGV[2]{
	clip[$1]=$2
	next
}
# shift records if running in parallel
shift{FNR+=shift; start=shift; shift=0}
# skip any deletions if decided to report
deletioncounter{print ""; deletioncounter--; next}
# read in the sorted map
{	
	mmax=0;
	for(i=1 ;i<=4; i++) # note bias
	{
		if ($i>mmax)
		{
			mmax=$i
			emax=i
		}
	}
	
	report=nucl[emax]

	
	if (FNR in seq){
		if (seq[FNR]~/-/)
		{ # deletion
			if (cov[FNR] + clip[FNR] + clip[FNR+length(seq[FNR])] >= mmax)
			{
				report=""
				mmax=cov[FNR]
			}
		}
		else
		{	# insertion, consider clips contributing
			#if (cov[FNR] >= mmax/2) # mmax is at least 2 for insertions
			if (cov[FNR] + clip[FNR] >= mmax/2) # mmax is at least 2 for insertions
			{
				report=seq[FNR]""report
				mmax=cov[FNR]
			}		
		}
	}

	# if too few calls report N else make sure that boundary calls come from reads rather than reference
	if (mmax<consensus_threshold)
	{
		report="N"
	} else if (mmax==consensus_threshold && ignore_reference)
	{
		# bother only if not indel as latter cannot come from reference
		if (report==nucl[emax])
		{
			if (NF<=5)
			{
				report="N" # shouldn't happen in complete pipeline
			} else
			{
				# if not indel choose call different from that of reference if ref map is provided
				while (emax<=4)
				{
					if ($emax==consensus_threshold && !$(emax+5))
						break
					emax++
				}
				report=nucl[emax]
			}
		}
	} 
	
	print report
	
	if (report=="")
	{
		deletioncounter=length(seq[FNR])-1
		# remove accounted clips and indels from list
		delete clip[FNR]; delete clip[FNR+length(seq[FNR])]; delete seq[FNR]
	} else if (report!=nucl[emax]){ # what to do about Ns? &&report!="N"
		#print "I am here "report", "nucl[emax] > "/dev/stderr"
		# remove clips from list that are accounted for
		delete clip[FNR]; delete seq[FNR]	
	}
# 	prevmmax=mmax
}
END{
	start+=0
	for (i in seq){
		i+=0
		if (i>start && i<=FNR && seq[i]!=""){
			print i, seq[i], cov[i] > "/dev/stderr"
		}
	}
	for (i in clip){
		i+=0
		if (i>start && i<=FNR && clip[i]!=""){
			print i, clip[i] > "/dev/stderr"
		}
	}
	for (i in alt){
		split(alt[i],a," ")
		a[1]+=0
		if (a[1]>start && a[1]<=FNR){
			print alt[i] > "/dev/stderr"
		}
	}
}
