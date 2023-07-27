% GOAL: loop over all XYquadCalc folders (10, 100,1000usec) to calcualte
% the 4pt SL at each resolution.
%%
% SET VALUES FOR THE 4-PT RUN TO BE CALCUALTED -- MUST SPECIFIY THE
% CONSTRCUT NAME AND VALUE OF TAU2 (ZERO BY DEFAULT FOR NOW)
% THE BINSIZE LIST DOENS'T NEED ADJUSTMENT - SHOULD REMAIN 10-100-1000
% AS THE INPUT/OUTPUTS OF THE STITCH 4PT RELY ON THESE RESOLUTIONS BEING
% USED - CAN BE ADJUSTED IF NEEDED....
% *****************************************************************
constructName='E5E6';
quantityName='P1AmpAvg';
crossCorr=0; 
% tau2 list must be in usec format, i.e. 1000=1ms
% tau2usecList=[0 1000 2000 5000 10000 50000 100000 250000]; 
% binSizeList=[1000];
tau2usecList=[0 1000]; 
binSizeList=[1000];
% determine wheter you want to caluclate tau1=tau3=0
calcZeroTau=1;
% turn on/off a subset calc, specify the name ('Uright' 'Total' 'boxSet'
% etc.) and then provide the file location of the Bivar_wAVG you would like
% to determine the trace files from
subsetOpt=0;
subsetName='Uright'; 
subsetLoc='C:\Users\jmaurer3\Documents\4ptTCFs\FPGA_Data\100mM_NaCl\(+15)Dimer\CombineData\XYQuadData_(+15)Dimer10usec\10usec DenseTau BF 2Pt_TCFs_(+15)Dimer\LD TCFs_v2\P1PhaseFact\Bivar_wAVG.mat';
% ******************************************************************
% BE IN THE COMBINED DATA FOLDER (OR WHICHEVER FOLDER HOLDS ALL THE XYQUAD
% FOLDERS AT EVERY RESOLUTION) -- THIS WILL BE THE "STARTFOLD"
startFold=pwd;
homeFold = dir(pwd);
% This will be the list of directories containing (*usec) in the name
dirList=[];
resList=[];
if strcmp(quantityName,'FRET')
keyWord= 'FRET';   
else
keyWord= 'XYQuadData';
end

%%

% Loops over the high level directory  to find where the sub folders
% containing the XYQUADCALE_#usec resolution folders are, ignoring the data files 
for i=1:numel(homeFold)
    
    if homeFold(i).isdir
        dirList= [dirList, homeFold(i)];
    end
    
end

% With the list of directories in the home folder, filter them to only the
% folders cotaining the string 'usec' in order to elimate the '.' and '..'
% directories. This is then the list of folders which contain the proper
% XYQuadCalc folders. 
for j=1:numel(dirList)
    if contains(dirList(j).name, keyWord)
        resList = [resList, dirList(j)];
    end
end 


% Now with only the res folders, The resulting 'resList' will be sorted
% automaticallty from the longest name to shortest name. So Assuming that
% 1-10000 usec data is available, the list will be 1usec at the last index,
% 10usec at the second to last index, and so on. 
numFolds= numel(resList);
% TenUsecFold= [startFold filesep() resList(numFolds).name]; 
% HundUsecFold= [startFold filesep() resList(numFolds-1).name]; 
% ThouUsecFold= [startFold filesep() resList(numFolds-2).name]; 
runTime10us=0;
fileTime10us=0;
runTime100us=0;
fileTime100us=0;
runTime1000us=0;
fileTime1000us=0;
runTime10000us=0;
fileTime10000us=0;

if size(binSizeList)>1
if binSizeList(1)>binSizeList(2)
    disp('Please generate binSizeList in ascending order - i.e. [10 100 1000]');
    return;
end 
end

%%

for m=1:numel(tau2usecList)
%%
for i=1:numel(binSizeList)
    if binSizeList(i)==10
        for j=1:numFolds
            if contains(resList(j).name,'10usec')
                target=resList(j).name; 
            end
        end
    cd([startFold filesep() target]);
    
    elseif binSizeList(i)==50
       for j=1:numFolds
            if contains(resList(j).name,'50usec')
                target=resList(j).name; 
            end
        end
    cd([startFold filesep() target]);
    
    elseif binSizeList(i)==100
       for j=1:numFolds
            if contains(resList(j).name,'100usec')
                target=resList(j).name; 
            end
        end
    cd([startFold filesep() target]);
    
    elseif binSizeList(i)==250
       for j=1:numFolds
            if contains(resList(j).name,'250usec')
                target=resList(j).name; 
            end
        end
    cd([startFold filesep() target]);
   
    elseif binSizeList(i)==1000
       for j=1:numFolds
            if contains(resList(j).name,'1000usec')
                target=resList(j).name; 
            end
        end
    cd([startFold filesep() target]);
    
    elseif binSizeList(i)==10000
       for j=1:numFolds
            if contains(resList(j).name,'10000usec')
                target=resList(j).name; 
            end
        end
    cd([startFold filesep() target]);
    end
    
    if ~crossCorr
    [totalRunTime,avgFileRunTime,tau1L,tau3L,NfilesToIterate,filesList] = Par_trace2C4_FPGA_FUNC_sub_v3(binSizeList(i),constructName,tau2usecList(m),subsetOpt,subsetLoc,subsetName,quantityName,calcZeroTau); 
    elseif crossCorr
    [totalRunTime,avgFileRunTime,tau1L,tau3L,NfilesToIterate,filesList] = Par_trace2C4_FPGA_FUNC_sub_v3_crossCorr(binSizeList(i),constructName,tau2usecList(m),subsetOpt,subsetLoc,subsetName,quantityName,calcZeroTau,crossCorr); 
    end
    
    if binSizeList(i)==10
        runTime10us=totalRunTime;
        fileTime10us=avgFileRunTime;  
    elseif binSizeList(i)==100
         runTime100us=totalRunTime;
         fileTime100us=avgFileRunTime;
    elseif binSizeList(i)==1000
         runTime1000us=totalRunTime;
         fileTime1000us=avgFileRunTime;
    elseif binSizeList(i)==10000
         runTime10000us=totalRunTime;
         fileTime10000us=avgFileRunTime;
    end
    
end
cd(startFold);
end 