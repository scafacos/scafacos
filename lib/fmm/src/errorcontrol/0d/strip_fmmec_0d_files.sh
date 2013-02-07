#!/bin/bash

grep " = 0" fmmmn.f | colrm 1 15 | colrm 26 1000 > fmmec_0d_consts.h

COUNTER=0
CUT_S=1
CUT_E=$((2*969))
CUT_L=$CUT_E
while [ $COUNTER -lt 51 ];
do
 echo Extracting multipole order $COUNTER  \(from line $CUT_S to line $CUT_E\)
 FILECOUNTER=`printf "%02d" "$COUNTER"`
 # split files | add comma | remove comma from last line's end
 sed -n ${CUT_S},${CUT_E}p fmmec_0d_consts.h | sed '''s/$/,/' | sed '''$s/,//' > include/fmmec_0d_p${FILECOUNTER}.h
 CUT_S=$((CUT_E+1))
 CUT_L=$((CUT_L+969))
 CUT_E=$((CUT_E+CUT_L))
 let COUNTER=COUNTER+1
done
echo "Remove temporary generated file"
rm -v fmmec_0d_consts.h

