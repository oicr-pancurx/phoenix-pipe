#!/bin/sh

FILE=$1;

module load tabix;

(grep start $FILE | sed 's/^/#/'; grep -v start $FILE | awk 'BEGIN {OFS="\t"} { $2=$2-1; print $0 }' | sed 's/^chr//' | sed 's/^X/23/' | sed 's/^Y/24/' | sed 's/^M/25/' | sort -k1,1n -k2,2n | sed 's/^25/M/' | sed 's/^24/Y/' | sed 's/^23/X/' | sed 's/^/chr/') | bgzip > ${FILE}.sorted.bed.gz
tabix -p bed ${FILE}.sorted.bed.gz

