#/bin/bash

COUNTER=0

BBB=-1
AAA=0
while [ $COUNTER -lt 51 ];
do
 echo Building from template for order $COUNTER 
 COUNTER2=$((COUNTER-1))
 FILECOUNTER=`printf "%02d" "$COUNTER"`
 cmdline="sed 's/AAA/$FILECOUNTER/' fmmec_3d_pXX.tmpl | sed '''s/BBB/$COUNTER2/' > fmmec_3d_p${FILECOUNTER}.c"
 eval $cmdline
 let COUNTER=COUNTER+1
done


