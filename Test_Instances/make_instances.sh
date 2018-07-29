#!/bin/bash

DIM="7"
CONVERTER=../../../converter/convert
LABELS=("o" "t" "a" "b" "c" "d" "e")
FILES=("o/testp1.lp" "o/testp2.lp" "a/testp2_a.lp" "b/testp2_b.lp" "c/testp2_c.lp" "d/testp2_d.lp" "e/testp2_e.lp")

for a in $( ls );
do
    echo $a
    cd $a
    rm -R $DIM
    mkdir $DIM
    for ((i = 0; i < 7; i++));
    do
        for ((j = i+1; j < 7; j++));
        do
            for ((k = j+1; k < 7; k++));
            do
                for ((l = k+1; l < 7; l++));
                do
                    for ((m = l+1; m < 7; m++));
                    do
                        for ((n = m+1; n < 7; n++));
                        do
                            for ((o = n+1; o < 7; o++));
                            do
                                $CONVERTER $DIM ${FILES[$i]} ${FILES[$j]} ${FILES[$k]} ${FILES[$l]} ${FILES[$m]} ${FILES[$n]} ${FILES[$o]}
                                NAME=${a}_${DIM}_${LABELS[$i]}${LABELS[$j]}${LABELS[$k]}${LABELS[$l]}${LABELS[$m]}${LABELS[$n]}${LABELS[$o]}
                                echo $NAME
                                mv multiobj_prob_MIP.mps $DIM/${NAME}.mps
                            done
                        done
                    done
                done
            done
        done
    done 
    rm *.mps
    cd ../
    #break
done  
