#### Description: Convenience awk script to calculate global coordinate in the genome from conventional chr:start-end.
#### Usage: awk -v query="chr:start-end" -f compute-global-coordinate.awk <path_to_sorted_cprops>
#### Input: "Sorted" cprops (LC_ALL=C).
#### Output: stdout "global_start_pos global_end_pos".
#### Written by: Olga Dudchenko - olga.dudchenko@bcm.edu. Version date 11/02/2017.

# parse query
BEGIN{
	gsub(",","",query)
	n=split(query,a,":||-")
}
$1!=a[1]{
	shift+=$3
	next
}
{
	if (a[2])
		a[2]+=shift
	else
		a[2]=shift+1
	if (a[3])
		a[3]+=shift
	else
		a[3]=shift+$3
	exit
}
END{
	if (a[2]&&a[3])
		print a[2], a[3]
	else
		exit 1
}