#!/usr/bin/env bash

cut -f1,2,12 $1 \
    | awk '$1>$2{print $2,$1,$3} $1<=$2{print $0}' \
    | sort -k1,1 -k2,2 \
    | awk -v OFS='\t' \
        '
        NR==1{pair[1]=$1; pair[2]=$2; tally=$3}
        NR>1{if (pair[1]!=$1 || pair[2]!=$2)
                {if (pair[1]==pair[2]) {sim = tally} else {sim = tally / 2};
                 print pair[1],pair[2],sim;
                 tally=$3;pair[1]=$1;pair[2]=$2
            }
             else
                {tally+=$3}}
        '
