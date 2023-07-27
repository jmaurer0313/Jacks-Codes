% read in 2D bulk data 

function [dmat, smat, t21ax, t43ax] = bulk2Dread(base_folder)
 

        t21ax = dlmread(fullfile(base_folder, 't21base'));
        t21ax= t21ax(:,1); 
        t43ax = dlmread(fullfile(base_folder, 't43base'));
        dmatR = dlmread(fullfile(base_folder, 't21XQuad'));
        dmatI = dlmread(fullfile(base_folder, 't21YQuad'));
        smatR = dlmread(fullfile(base_folder, 't43XQuad'));
        smatI = dlmread(fullfile(base_folder, 't43YQuad'));
        
        dmat= dmatR +1i.*dmatI; 
        smat= smatR +1i*smatI; 
       

end 