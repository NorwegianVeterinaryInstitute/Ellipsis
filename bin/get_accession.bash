#!/bin/bash
file=$1
filename0=${file/_mobtyper}
filename=${filename0%%.fasta_report.txt}.acc
tail --lines=+2 $file | awk -F " " '{print $14}' > ./${filename}
