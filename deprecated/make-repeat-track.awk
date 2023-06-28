BEGIN{
	# defaults
	if(!k){k=100}
	if(!kk){kk=15}
	if(!threshold){threshold=5}
}
{
	if (FNR>=k){
		str=substr(str,2)""$0
	} else {
		str=str""$0
		next
	}
	
	for(i=kk;i<=k;i++){
		tmp[i]=substr(str,i-kk+1,kk)
	}
	
	max=0
	tent_counter=1
	
	n=asort(tmp)
	for(i=2; i<=n+1; i++){
		if(tmp[i-1]==tmp[i]){
			tent_counter++
		}else{
			if(tent_counter>max){
				max=tent_counter
				tent_counter=1
			}
		}
	}
	print max
	if(max>threshold){print FNR}
}