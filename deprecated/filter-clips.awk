BEGIN{
	len=length(seq)
	if(!tolerance){tolerance=1}
}
# clips
FILENAME==ARGV[1]{
	cov[$1]=$2
	next
}
# keep track of last len+2*tolerance lines (tolerance on either side)
{
	str[(FNR-1) % (len+2*tolerance)]=$0
}
{
	tmp=""
	for (i=1; i<=len+2*tolerance;i++)
	{
		tmp=tmp""str[(FNR-1+i) % (len+2*tolerance)]
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