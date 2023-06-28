BEGIN{
	# defaults
	if(!width){width=100}
}
($1<prev+width){
	prev=$1
	counter+=$NF
	next
}
start{
	print start, prev, counter
}
{
	start=$1
	prev=$1
	counter=$NF
}
END{
	print start, prev, counter	
}