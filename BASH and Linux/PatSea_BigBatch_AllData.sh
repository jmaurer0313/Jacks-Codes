#!/bin/bash

#cur input flags list in order: numberOfModels modelSetName modelSetIdx 
#ModelSetIdx pings the modelSet name in the underlying matlab code
#order is : 1=model_lin 2=model_loop1 3=model_loopN 4=model_brch

#loop over the conditions array within the GenAlg code itself
for ((k=1; k<= $4 ; k++))
do

#create the directories for each models err and out files
for ((i=1 ; i<=$1 ; i++))

do
mkdir -p errOutFiles_cond$k/myDir$i$2_$k
done

#move copies of the parent level codes into the individual err and out file folders
for ((j=1; j<=$1 ; j++))

do
cp PatSea_Nstate_polz_wGActrl_FUNC_multiOut_v2_weightCalc_ModGen.m errOutFiles_cond$k/myDir$j$2_$k/
cp PatSea_BigBatch_AllData.srun errOutFiles_cond$k/myDir$i$2_$k/
done

#now enter each output folder, start a unique bash job in each one, return to starting dir each time
for ((z=1; z<=$1; z++))

do
cd errOutFiles_cond$k/myDir$i$2_$k
#submit the batch file here
sbatch PatSea_BigBatch_AllData.srun $z $3 $k;
cd ..
cd .. 
done
done