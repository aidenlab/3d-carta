## Reference to create a map file from the reference fasta.
## awk -f parse-reference.awk <path-to-reference-fasta>
## Output: map file formatted stdout.
## Written by: OD
## Version: 171120 
$0~/>/{next}
{
	n=split($0,nucl,"");
	for(i=1;i<=n;i++){
		if (nucl[i]=="A"||nucl[i]=="a")
			print "1 0 0 0 0"
		else if (nucl[i]=="C" || nucl[i]=="c")
			print "0 1 0 0 0"
		else if (nucl[i]=="G" || nucl[i]=="g")
			print "0 0 1 0 0"
		else if (nucl[i]=="T" || nucl[i]=="t")
			print "0 0 0 1 0" 
		else print "0 0 0 0 1"
	}
}
