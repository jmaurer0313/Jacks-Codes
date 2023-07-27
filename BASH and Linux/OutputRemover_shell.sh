module load matlab

#for ((z=1; z<=$1; z++))

#do
#flag 1 is output number 

matlab -nodisplay -r  "OutputRemover_v1($1)";
 

#done