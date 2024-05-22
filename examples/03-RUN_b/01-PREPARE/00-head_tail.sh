#!bin/bash
nchains=2
nchA=2
nchB=0
natchA=80
# First kind chain
idxs=1
idxe=72

echo "nmols: $nchains" >PE10-iPP10m_block_2Ch_headtail.dat
echo "#ich idx-head idx-tail (indexes start at 0)" >>PE10-iPP10m_block_2Ch_headtail.dat
for (( i=0; i<$nchA; i++ )) do
    echo "$i $idxs $idxe" >>PE10-iPP10m_block_2Ch_headtail.dat
    idxs=`echo $idxs+$natchA|bc -l`
    idxe=`echo $idxe+$natchA|bc -l`
done

# Second kind chain
natchB=464
idxsB=6161
idxeB=6620
for (( i=$nchA; i<$nchains; i++ )) do
    echo "$i $idxsB $idxeB" >>PE10-iPP10m_block_2Ch_headtail.dat
    idxsB=`echo $idxsB+$natchB|bc -l`
    idxeB=`echo $idxeB+$natchB|bc -l`
done

