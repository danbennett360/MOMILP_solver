#!/bin/bash

SEED=4001
DIM=10
CONVERTER=../../converter/convert
GENERATOR=generator/generate_instances
FILES=("testp1.lp" "testp2.lp" "testp3.lp" "testp4.lp" "testp5.lp" "testp6.lp" "testp7.lp" "testp8.lp" "testp9.lp" "testp10.lp")
    
for ((VAR = 25; VAR < 800; VAR *= 2));
do
    echo $VAR
    for ((CON = 25; CON < 800; CON *= 2));
    do
        N=${VAR}_${CON}
        mkdir $DIM/$N
        for ((a = 0; a < 5; a++));
        do
            $GENERATOR $SEED $DIM $VAR $CON $(($VAR/2)) 0 $(($VAR/2))
            $CONVERTER $DIM ${FILES[0]} ${FILES[1]} ${FILES[2]} ${FILES[3]} ${FILES[4]} ${FILES[5]} ${FILES[6]} ${FILES[7]} ${FILES[8]} ${FILES[9]}
            NAME=belottiStyle${SEED}_${DIM}_${VAR}_${CON}_$(($VAR/2))_0_$(($VAR/2))
            mv multiobj_prob_MIP.mps $DIM/$N/${NAME}.mps
            SEED=$(($SEED+1))
            $GENERATOR $SEED $DIM $VAR $CON $(($VAR/2)) $(($VAR/2)) 0
            $CONVERTER $DIM ${FILES[0]} ${FILES[1]} ${FILES[2]} ${FILES[3]} ${FILES[4]} ${FILES[5]} ${FILES[6]} ${FILES[7]} ${FILES[8]} ${FILES[9]}
            NAME=belottiStyle${SEED}_${DIM}_${VAR}_${CON}_$(($VAR/2))_$(($VAR/2))_0
            mv multiobj_prob_MIP.mps $DIM/$N/${NAME}.mps
            SEED=$(($SEED+1))
            $GENERATOR $SEED $DIM $VAR $CON $(($VAR/2)) $(($VAR/4)) $(($VAR/4))
            $CONVERTER $DIM ${FILES[0]} ${FILES[1]} ${FILES[2]} ${FILES[3]} ${FILES[4]} ${FILES[5]} ${FILES[6]} ${FILES[7]} ${FILES[8]} ${FILES[9]}
            NAME=belottiStyle${SEED}_${DIM}_${VAR}_${CON}_$(($VAR/2))_$(($VAR/4))_$(($VAR/4))
            mv multiobj_prob_MIP.mps $DIM/$N/${NAME}.mps
            SEED=$(($SEED+1))
            $GENERATOR $SEED $DIM $VAR $CON $(($VAR/4)) $(($VAR/4)) $(($VAR/2))
            $CONVERTER $DIM ${FILES[0]} ${FILES[1]} ${FILES[2]} ${FILES[3]} ${FILES[4]} ${FILES[5]} ${FILES[6]} ${FILES[7]} ${FILES[8]} ${FILES[9]}
            NAME=belottiStyle${SEED}_${DIM}_${VAR}_${CON}_$(($VAR/4))_$(($VAR/4))_$(($VAR/2))
            mv multiobj_prob_MIP.mps $DIM/$N/${NAME}.mps
            SEED=$(($SEED+1))
        done
    done
done

echo $SEED

#for ((i = 0; i < 7; i++));
#do
#    for ((j = i+1; j < 7; j++));
#    do
#        for ((k = j+1; k < 7; k++));
#        do
#            for ((l = k+1; l < 7; l++));
#            do
#                for ((m = l+1; m < 7; m++));
#                do
#                    for ((n = m+1; n < 7; n++));
#                    do
#                        for ((o = n+1; o < 7; o++));
#                        do
#                            $CONVERTER $DIM ${FILES[$i]} ${FILES[$j]} ${FILES[$k]} ${FILES[$l]} ${FILES[$m]} ${FILES[$n]} ${FILES[$o]}
#                            NAME=${a}_${DIM}_${LABELS[$i]}${LABELS[$j]}${LABELS[$k]}${LABELS[$l]}${LABELS[$m]}${LABELS[$n]}${LABELS[$o]}
#                            echo $NAME
#                            mv multiobj_prob_MIP.mps $DIM/${NAME}.mps
#                        done
#                    done
#                done
#            done
#        done
#    done
#done 
