FILENAME==ARGV[1]{
	len[$1]=$2-$1
	seq[$1]=$3
	next
}
counter{
	print ""
	counter--
	next
}
(FNR in len){
	print seq[FNR]
	counter=len[FNR]
	next
}
1