#!/bin/bash

DIM="7"
CONVERTER=../../../converter/convert
LABELS=("o" "t" "a" "b" "c" "d" "e")
FILES=("o/testp1.lp" "o/testp2.lp" "a/testp2_a.lp" "b/testp2_b.lp" "c/testp2_c.lp" "d/testp2_d.lp" "e/testp2_e.lp")

for a in $( ls );
do
    echo $a
    cd $a
    for ((i = 0; i < 7; i++));
    do
        rm -R ${LABELS[$i]}
        rm *.lp
        rm *.dat
    done 
    cd ../
done  
