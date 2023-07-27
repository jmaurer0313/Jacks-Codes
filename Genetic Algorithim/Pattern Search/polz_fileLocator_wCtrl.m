function [GenAlgFilesLocationPrefix, genAlgDataLocationPrefix, histogram_FilePath, C2_FilePath,C4_FilePath] = polz_fileLocator(terminalID, constructName, modelName, programName, protein_str, saltConditions, fileSuffix)

global verboseMode model_name normalizeMode testTimesMode 
global computer_terminal_str

% PURPOSE:  (1) finds GenAlg files
%           (2) loads targetdata
%
% INPUT:    terminalID from computerMode (don't need computer_terminal_str)
%           constructName, modelName, protein_str (in GenAlg script)

% verboseMode = 1;

% Part 1: Find the folder for the GenAlg you are running
%--------------------------------------------------------------------------
% Navigate to the correct folder with the target data
%--------------------------------------------------------------------------
    switch computer_terminal_str
        case 'computer_baimi_mode' %Work Desktop
            if normalizeMode
            genAlgDataLocationPrefix =  ['C:\Users\baimi\Documents\MATLAB\GenAlg\' constructName filesep() saltConditions '\genAlgData\targetData'];
            else
                if testTimesMode
                    genAlgDataLocationPrefix =  ['C:\Users\baimi\Documents\MATLAB\GenAlg\' constructName filesep() saltConditions '\genAlgData_UNorm_testwCtrl_' fileSuffix '\targetData'];
                else
                    genAlgDataLocationPrefix =  ['C:\Users\baimi\Documents\MATLAB\GenAlg\' constructName filesep() saltConditions '\genAlgData_UNorm\targetData'];
                end
%             genAlgDataLocationPrefix =  ['C:\Users\Jack\Dropbox\chosenAPD2mat_output\' constructName '\gp32_0p0uM\FPGA\' saltConditions '\genAlgData_UNorm\targetData'];
            end

        case 'computer_bisraels_mode' %Macbook pro
            genAlgDataLocationPrefix =  ['/Users/bisraels/Dropbox/MarcusLab/Data/smData_Processed/' constructName '/gp32_0p0uM/ChosenMolecules/genAlgFiles/targetData'];
        case 'computer_claire_mode'
            %             genAlgDataLocationPrefix =  ['/Users/clairealbrecht/Dropbox/MarcusLab/Data/smData_Processed/' constructName '/gp32_0p0uM/ChosenMolecules/genAlgFiles/targetData'];
%             if BrettData_mode == 1
%                 genAlgDataLocationPrefix =  ['/Users/clairealbrecht/Dropbox/MarcusLab/Data/smData_Processed/' constructName '/gp32_0p5uM/ChosenMolecules/genAlgFiles/TargetData'];
%             else
                genAlgDataLocationPrefix =  ['/Users/clairealbrecht/Dropbox/MATLAB_programs/claire_programs/from_Jack/updatedgenalgcodes/targetData'];
%             end
        case 'computer_calbrecht_mode'
%             if BrettData_mode == 1
%                 genAlgDataLocationPrefix =  ['/Users/calbrecht/Dropbox/MarcusLab/Data/smData_Processed/' constructName '/gp32_0p5uM/ChosenMolecules/genAlgFiles/TargetData'];
%             else
            genAlgDataLocationPrefix =  ['/Users/calbrecht/Dropbox/MATLAB_programs/claire_programs/from_Jack/updatedgenalgcodes/targetData'];
%             end
        case 'computer_JackLaptop_mode' %Jacks Laptop
            if normalizeMode
            genAlgDataLocationPrefix =  ['C:\Users\Jack\Dropbox\chosenAPD2mat_output\' constructName '\gp32_0p0uM\FPGA\' saltConditions '\genAlgData\targetData'];
            else
                if testTimesMode                    
                    genAlgDataLocationPrefix =  ['C:\Users\Jack\Dropbox\chosenAPD2mat_output\' constructName '\gp32_0p0uM\FPGA\' saltConditions '\genAlgData_UNorm_testwCtrl_' fileSuffix '\targetData'];                   
                else                    
                    genAlgDataLocationPrefix =  ['C:\Users\Jack\Dropbox\chosenAPD2mat_output\' constructName '\gp32_0p0uM\FPGA\' saltConditions '\genAlgData_UNorm\targetData'];
                end
%             genAlgDataLocationPrefix =  ['C:\Users\Jack\Dropbox\chosenAPD2mat_output\' constructName '\gp32_0p0uM\FPGA\' saltConditions '\genAlgData_UNorm\targetData'];
            end  
         case 'computer_talapasDesktop_mode' %Jacks Laptop
            if normalizeMode
            genAlgDataLocationPrefix =  ['C:\Users\Jack\Dropbox\chosenAPD2mat_output\' constructName '\gp32_0p0uM\FPGA\' saltConditions '\genAlgData\targetData'];
            else
                if testTimesMode                    
                    genAlgDataLocationPrefix =  ['/gpfs/home/jmaurer3/Documents/MATLAB/GenAlg_Data/' constructName '/' saltConditions '/genAlgData_UNorm_testwCtrl_' fileSuffix '/targetData'];                   
                else                    
                    genAlgDataLocationPrefix =  ['C:\Users\Jack\Dropbox\chosenAPD2mat_output\' constructName '\gp32_0p0uM\FPGA\' saltConditions '\genAlgData_UNorm\targetData'];
                end
%             genAlgDataLocationPrefix =  ['C:\Users\Jack\Dropbox\chosenAPD2mat_output\' constructName '\gp32_0p0uM\FPGA\' saltConditions '\genAlgData_UNorm\targetData'];
            end  
            
        case 'computer_labTower_mode' %Jacks Laptop
            if normalizeMode
            genAlgDataLocationPrefix =  ['C:\Users\jmaurer3\Documents\KineticModeling\Data\' constructName filesep() saltConditions '\genAlgData\targetData'];
            else
                if testTimesMode                    
                    genAlgDataLocationPrefix =  ['C:\Users\jmaurer3\Documents\KineticModeling\Data\' constructName filesep() saltConditions '\genAlgData_UNorm_testwCtrl_' fileSuffix '\targetData'];                    
                else                    
                    genAlgDataLocationPrefix =  ['C:\Users\jmaurer3\Documents\KineticModeling\Data\' constructName filesep() saltConditions '\genAlgData_UNorm\targetData'];
                end
            end
        otherwise
            error('Cannot detect computer');
    end
    
    
    %%
    
    %--------------------------------------------------------------------------
    % Navigate to the correct folder with the GenAlg output
    %--------------------------------------------------------------------------
    switch computer_terminal_str
        case 'computer_baimi_mode' %Work Desktop
            if normalizeMode
            GenAlgFilesLocationPrefix =  ['C:\Users\baimi\Documents\MATLAB\GenAlg\' constructName filesep() saltConditions '\genAlgData\Outputs\' programName '_' model_name];
            else
                if testTimesMode
                    GenAlgFilesLocationPrefix =  ['C:\Users\baimi\Documents\MATLAB\GenAlg\' constructName filesep() saltConditions '\genAlgData_UNorm_testwCtrl_' fileSuffix '\Outputs\' programName '_' model_name];
%                GenAlgFilesLocationPrefix =  ['C:\Users\jmaurer3\Documents\KineticModeling\Data\' constructName filesep() saltConditions '\genAlgData_UNorm_test\Outputs\' programName '_' model_name];
                else
               GenAlgFilesLocationPrefix =  ['C:\Users\baimi\Documents\MATLAB\GenAlg\' constructName filesep() saltConditions '\genAlgData_UNorm\Outputs\' programName '_' model_name];
                end
%             GenAlgFilesLocationPrefix =  ['C:\Users\Jack\Dropbox\chosenAPD2mat_output\' constructName '\gp32_0p0uM\FPGA\' saltConditions '\genAlgData_UNorm\Outputs\' programName '_' model_name];
            end
        case 'computer_bisraels_mode' %Macbook pro
            GenAlgFilesLocationPrefix =  ['/Users/bisraels/Dropbox/MarcusLab/Data/smData_Processed/' constructName '/gp32_0p0uM/ChosenMolecules/genAlgFiles/3state123_cyclical'];
        
        case 'computer_claire_mode' %Macbook pro
%             GenAlgFilesLocationPrefix =  ['/Users/clairealbrecht/Dropbox/MarcusLab/Data/smData_Processed/' constructName '/gp32_0p0uM/ChosenMolecules/genAlgFiles/3state123_cyclical'];
%            if BrettData_mode == 1
%               GenAlgFilesLocationPrefix = ['/Users/clairealbrecht/Dropbox/MarcusLab/Data/smData_Processed/' constructName '/gp32_0p5uM/ChosenMolecules/genAlgFiles/outputs/'  programName '_' model_name];
%            else
              GenAlgFilesLocationPrefix = ['/Users/clairealbrecht/Dropbox/MATLAB_programs/claire_programs/from_Jack/updatedgenalgcodes/outputs/' programName];% '_' model_name];
%            end
           
        case 'computer_calbrecht_mode' %Macbook pro
%             GenAlgFilesLocationPrefix =  ['/Users/clairealbrecht/Dropbox/MarcusLab/Data/smData_Processed/' constructName '/gp32_0p0uM/ChosenMolecules/genAlgFiles/3state123_cyclical'];
%             if polzMode == 1
%                 GenAlgFilesLocationPrefix = ['/Users/calbrecht/Dropbox/MATLAB_programs/claire_programs/from_Jack/updatedgenalgcodes/outputs/' programName '_' model_name];                
%             else
                GenAlgFilesLocationPrefix = ['/Users/calbrecht/Dropbox/MarcusLab/Data/smData_Processed/' constructName '/gp32_0p5uM/ChosenMolecules/genAlgFiles/outputs/'  programName ];%'_' model_name];
%             end
            
        case 'computer_JackLaptop_mode' %Jacks Laptop
            if normalizeMode
            GenAlgFilesLocationPrefix =  ['C:\Users\Jack\Dropbox\chosenAPD2mat_output\' constructName '\gp32_0p0uM\FPGA\' saltConditions '\genAlgData\Outputs\' programName '_' model_name];
            else
                if testTimesMode                   
                    GenAlgFilesLocationPrefix =  ['C:\Users\Jack\Dropbox\chosenAPD2mat_output\' constructName '\gp32_0p0uM\FPGA\' saltConditions '\genAlgData_UNorm_testwCtrl_' fileSuffix '\Outputs\' programName '_' model_name];
                   
%                GenAlgFilesLocationPrefix =  ['C:\Users\jmaurer3\Documents\KineticModeling\Data\' constructName filesep() saltConditions '\genAlgData_UNorm_test\Outputs\' programName '_' model_name];
                else
               GenAlgFilesLocationPrefix =  ['C:\Users\Jack\Dropbox\chosenAPD2mat_output\' constructName '\gp32_0p0uM\FPGA\' saltConditions '\genAlgData_UNorm\Outputs\' programName '_' model_name];
                end
%             GenAlgFilesLocationPrefix =  ['C:\Users\Jack\Dropbox\chosenAPD2mat_output\' constructName '\gp32_0p0uM\FPGA\' saltConditions '\genAlgData_UNorm\Outputs\' programName '_' model_name];
            end
        case 'computer_talapasDesktop_mode' %Jacks Laptop
            if normalizeMode
            GenAlgFilesLocationPrefix =  ['C:\Users\Jack\Dropbox\chosenAPD2mat_output\' constructName '\gp32_0p0uM\FPGA\' saltConditions '\genAlgData\Outputs\' programName '_' model_name];
            else
                if testTimesMode                   
                    GenAlgFilesLocationPrefix =  ['/gpfs/home/jmaurer3/Documents/MATLAB/GenAlg_Data/' constructName '/' saltConditions '/genAlgData_UNorm_testwCtrl_' fileSuffix '/Outputs/' programName '_' model_name];
                   
%                GenAlgFilesLocationPrefix =  ['C:\Users\jmaurer3\Documents\KineticModeling\Data\' constructName filesep() saltConditions '\genAlgData_UNorm_test\Outputs\' programName '_' model_name];
                else
               GenAlgFilesLocationPrefix =  ['C:\Users\Jack\Dropbox\chosenAPD2mat_output\' constructName '\gp32_0p0uM\FPGA\' saltConditions '\genAlgData_UNorm\Outputs\' programName '_' model_name];
                end
%             GenAlgFilesLocationPrefix =  ['C:\Users\Jack\Dropbox\chosenAPD2mat_output\' constructName '\gp32_0p0uM\FPGA\' saltConditions '\genAlgData_UNorm\Outputs\' programName '_' model_name];
            end
            
        case 'computer_labTower_mode' %Jacks Laptop
             if normalizeMode
            GenAlgFilesLocationPrefix =  ['C:\Users\jmaurer3\Documents\KineticModeling\Data\' constructName filesep() saltConditions '\genAlgData\Outputs\' programName '_' model_name];
            else
                if testTimesMode                    
                    GenAlgFilesLocationPrefix =  ['C:\Users\jmaurer3\Documents\KineticModeling\Data\' constructName filesep() saltConditions '\genAlgData_UNorm_testwCtrl_' fileSuffix '\Outputs\' programName '_' model_name];                    
%                GenAlgFilesLocationPrefix =  ['C:\Users\jmaurer3\Documents\KineticModeling\Data\' constructName filesep() saltConditions '\genAlgData_UNorm_test\Outputs\' programName '_' model_name];
                else
               GenAlgFilesLocationPrefix =  ['C:\Users\jmaurer3\Documents\KineticModeling\Data\' constructName filesep() saltConditions '\genAlgData_UNorm\Outputs\' programName '_' model_name];
                end
            end
    end
    if exist(GenAlgFilesLocationPrefix,'dir') ~= 7
        disp('Making a directory to hold the outputs of the model');
        mkdir(GenAlgFilesLocationPrefix);
    end
    cd(GenAlgFilesLocationPrefix);
    wd = pwd;
    if verboseMode == 1
        disp(['You are now in the folder: ' wd]);
    end



if exist(GenAlgFilesLocationPrefix,'dir') ~= 7
    disp('Making a directory to hold the outputs of the model');
    mkdir(GenAlgFilesLocationPrefix);
end




wd = pwd;
if verboseMode == 1
    disp(['You are now in the folder: ' wd]);
end


if exist(genAlgDataLocationPrefix,'dir') ~= 7
    error(['Cannot locate data folder ' genAlgDataLocationPrefix]);
end

if verboseMode == 1
    disp(['Looking for data in ' genAlgDataLocationPrefix]);
end


% Define file paths for the data surfaces
HistfNameKeyWord = '*HistData.mat';
fileNames = dir([genAlgDataLocationPrefix filesep() HistfNameKeyWord]);
if isempty(fileNames) == 1
    error(['Could not Find any files with the name '  fNameKeyWord ' in ' genAlgDataLocationPrefix]);
end
histogram_FileName = fileNames(1).name;
histogram_FilePath = [genAlgDataLocationPrefix filesep() histogram_FileName];


C2fNameKeyWord = '*StitchTCF__wAVG.mat';
fileNames = dir([genAlgDataLocationPrefix filesep() C2fNameKeyWord]);
if isempty(fileNames) == 1
    error(['Could not Find any files with the name '  C2fNameKeyWord]);
end
C2_FileName = fileNames(1).name;
C2_FilePath = [genAlgDataLocationPrefix filesep() C2_FileName];


C4fNameKeyWord = 'C4*.mat';
fileNames = dir([genAlgDataLocationPrefix filesep() C4fNameKeyWord]);
if isempty(fileNames) == 1
    error(['Could not Find any files with the name '  C4fNameKeyWord]);
end
C4_FileName = fileNames(1).name;
C4_FilePath = [genAlgDataLocationPrefix filesep() C4_FileName];

end
