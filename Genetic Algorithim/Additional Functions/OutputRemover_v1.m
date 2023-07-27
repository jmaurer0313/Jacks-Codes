% given a base folder, look through all the Pat sea outputs and delete the
% outputN folders based on the input index

function [] = OutputRemover_v1(outputNumber)

% should be the folder level containign the outputs, so i.e. xxxx/Outputs/
baseFold='C:\Users\Jack\Dropbox\TempHolder_modGenPS\TestFolder\Outputs\';

programName= 'GenAlg_Nstate_updated2';
patSeaVersion=6; 

% the output folder number to remove
% outputNumber=2; 


cd(baseFold); 
modelFolders= dir([programName '*']);


for i=1:length(modelFolders)
   cd(modelFolders(i).name); 
   
if exist([baseFold modelFolders(i).name filesep() 'PatSea_Nstate_polz_outputs_v' num2str(patSeaVersion)])    
    cd([baseFold modelFolders(i).name filesep() 'PatSea_Nstate_polz_outputs_v' num2str(patSeaVersion)]);    
    if exist([baseFold modelFolders(i).name filesep() 'PatSea_Nstate_polz_outputs_v' num2str(patSeaVersion) filesep() 'output' num2str(outputNumber)])
        
        rmdir([baseFold modelFolders(i).name filesep() 'PatSea_Nstate_polz_outputs_v' num2str(patSeaVersion) filesep() 'output' num2str(outputNumber)],'s')
        cd(baseFold); 
    
    end
    
end

end

end