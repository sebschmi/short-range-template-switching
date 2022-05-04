{OFS="\t";
if($6~/S/ && $6!~/[DI]/){
	bwd=0; if(and(16,$2)==16){bwd=1}
	sms=match($6,"[0-9]+S[0-9]+M[0-9]+S",asms);
	sm=match($6,"[0-9]+S[0-9]+M",asm);
	ms=match($6,"[0-9]+M[0-9]+S",ams);
	if(sms==1 && length(asms[0])==length($6)) {
		match($6,"([0-9]+).([0-9]+).([0-9]+)",len);
#		print $2,bwd,$4,$6,"sms",len[1],len[2],len[3]
		print $4
		print $4+len[2]
#		print substr($10,len[1]+1,len[2])
	} else if(sm==1 && length(asm[0])==length($6)) {
		match($6,"([0-9]+).([0-9]+)",len);
#		print $2,bwd,$4,$6,"sm",len[1],len[2]		
		print $4
#		print substr($10,len[1]+1)
	} else if(ms==1 && length(ams[0])==length($6)) {
		match($6,"([0-9]+).([0-9]+)",len);
#		print $2,bwd,$4,$6,"ms",len[1],len[2] 
		print $4+len[1]
#		print substr($10,1,len[1])
	}
}
}
