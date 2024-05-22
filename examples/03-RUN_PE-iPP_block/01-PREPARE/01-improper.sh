#!/bin/bash
  
nkind=1

nimprA=23
iatA=50
jatA=49
katA=52
latA=51

nimprB=0
iatB=2
jatB=1
katB=4
latB=3

delta=3

echo "nkinds: ${nkind}" >impropers.ndx
echo "# ikind number of impropers" >>impropers.ndx
echo "0 ${nimprA}" >>impropers.ndx
for (( i=0; i<${nimprA}; i++ )); do
    echo "${iatA}       ${jatA}       ${katA}       ${latA}   2   -30.500 0.5178E+03" >>impropers.ndx
    iatA=`echo $iatA+$delta|bc -l`
    jatA=`echo $jatA+$delta|bc -l`
    katA=`echo $katA+$delta|bc -l`
    latA=`echo $latA+$delta|bc -l`
done

#echo "# ikind number of impropers" >>impropers.ndx
#echo "1 ${nimprB}" >>impropers.ndx
#for (( i=0; i<${nimprB}; i++ )); do
#    echo "${iatB}       ${jatB}       ${katB}       ${latB}   2   -30.500 0.5178E+03" >>impropers.ndx
#    iatB=`echo $iatB+$delta|bc -l`
#    jatB=`echo $jatB+$delta|bc -l`
#    katB=`echo $katB+$delta|bc -l`
#    latB=`echo $latB+$delta|bc -l`
#done
