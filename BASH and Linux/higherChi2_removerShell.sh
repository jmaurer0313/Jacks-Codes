module load matlab

#for ((z=1; z<=$1; z++))

#do
#flag 1 is model set 

matlab -nodisplay -r  "higherchisq_remover($1)";
 

#done
