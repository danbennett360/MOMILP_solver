#!/bin/bash

DIM="5"
FILES="DansStuff/rnd_5_225_30.mop"
FLAGS="-store T -lexopt F -epsilon .0001 -showprogress F -progressvalue 20 -reldist T -maxiter 1000000 -debug 0"


echo "GLPK time using simplex method:"
/usr/bin/time -p ./gMOS $DIM $FILES $FLAGS > gout1.txt
mv points.txt points3.txt
echo
#echo "GLPK time using interior point method:"
#/usr/bin/time -p ./gMOS $DIM $FILES $FLAGS -interior T > gout2.txt
#mv points.txt points2.txt
#echo
echo "CPLEX time:"
/usr/bin/time -p ./cMOS $DIM $FILES $FLAGS > cout.txt
