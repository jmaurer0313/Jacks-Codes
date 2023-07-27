% quick script to copy and move select files based on scanID list


allFiles=dir('*.mat');
destinationFold ='C:\Users\Ryzen 5\Dropbox\chosenAPD2mat_output\E5E6\gp32_0p0uM\FPGA\100mM_NaCl\SortedSubset';

for j=1:length(scanIDList)
    for i=1:length(allFiles)
        if contains(allFiles(i).name, scanIDList(j))
           curName=allFiles(i).name; 
           copyfile(curName, destinationFold)  
        end
    
   end
end