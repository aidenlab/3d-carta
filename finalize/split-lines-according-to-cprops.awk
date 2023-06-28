#### Description: Script to interleave globally-positioned map/reduce file or stdout with fasta markers listed in cprops file.
#### Usage: awk -f split-map-according-to-cprops.awk <cprops-file> <line-to-position-file>
#### Input: Cprops file, globally-positioned-file.
#### Output: Up to wrapping standard fasta stdout.
#### NOTE: Ordering in cprops is key.
#### TODO (maybe): do wrapping internally?
#### Written by: Olga Dudchenko - olga.dudchenko@bcm.edu, version 11/29/2017

# read in cprops
BEGIN{shift=0}
FILENAME==ARGV[1]{
	name[shift]=$1
	shift+=$3
	next
}
# read in globally-positioned-file
((FNR-1) in name){print ">"name[FNR-1]}
1