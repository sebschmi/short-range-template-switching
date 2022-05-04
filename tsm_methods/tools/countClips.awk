BEGIN{pp=0;OFS="\t";}
{
if(beg==0){sum=$1;beg=$2;end=$2}
if(pp+20>$2)
{sum=sum+$1;end=$2}
else
{
if(beg>0 && end>0 && sum>10){print chr,beg,end,end-beg+1,sum}
sum=$1;beg=$2;end=$2
}
pp=$2
}
