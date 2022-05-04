BEGIN{pp=0;OFS="\t";}
{
mut=1
ldf=length($4)-length($5)
if(ldf>1) mut=ldf
if(ldf*-1>1) mut=ldf*-1
if(beg==0){sum=mut;beg=$2;end=$2}
if(pp+20>$2)
{sum=sum+mut;end=$2}
else
{
if(beg>0 && end>0 && sum>1){print $1,beg,end,end-beg+1,sum}
sum=mut;beg=$2;end=$2
}
pp=$2
}
