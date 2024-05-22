#!bin/bash
nchains=3
nchA=3
nchB=0
natchA=365
# First kind chain
idxs=1
idxe=357

echo "nmols: $nchains" >PE24iPP24_model_headtail.dat
echo "#ich idx-head idx-tail (indexes start at 0)" >>PE24iPP24_model_headtail.dat
for (( i=0; i<$nchA; i++ )) do
    echo "$i $idxs $idxe" >>PE24iPP24_model_headtail.dat
    idxs=`echo $idxs+$natchA|bc -l`
    idxe=`echo $idxe+$natchA|bc -l`
done

# Second kind chain
natchB=464
idxsB=6161
idxeB=6620
for (( i=$nchA; i<$nchains; i++ )) do
    echo "$i $idxsB $idxeB" >>PE24iPP24_model_headtail.dat
    idxsB=`echo $idxsB+$natchB|bc -l`
    idxeB=`echo $idxeB+$natchB|bc -l`
done

