#!/bin/bash

bam=$1
ref=$2
outputname0=${bam%%.bam}
outputname=${outputname0##*/}

zero=$(bedtools genomecov -ibam $bam -bga | awk '$4==0 {bpCountZero+=($3-$2)} {print bpCountZero}' | tail -1)
nonzero=$(bedtools genomecov -ibam $bam -bga | awk '$4>0 {bpCountNonZero+=($3-$2)} {print bpCountNonZero}' | tail -1)

if [ -z "$nonzero" ]
then
	nonzero=0
	percent=$(bc <<< "scale=6; ($nonzero / ($zero + $nonzero))*100")
elif [ -z "$zero" ]
then
	zero=0
	percent=$(bc <<< "scale=6; ($nonzero / ($zero + $nonzero))*100")
else
	percent=$(bc <<< "scale=6; ($nonzero / ($zero + $nonzero))*100")
fi

echo -e "id\tpercent_mapped\tplasmid\n$outputname\t$percent\t$ref"
