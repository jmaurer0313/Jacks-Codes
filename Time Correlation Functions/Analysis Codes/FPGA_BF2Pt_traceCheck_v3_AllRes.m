%AUTHOR: Jack Maurer 
% Goal: perform brute force 2pt TCF calculations on the FPGA data arrays.
% Store outputs for both P1 and P2 assuming they both carry some
% information 

disp(['Your current directory is' string(pwd)]);

tic
tracesCalcd=0;
% traceFolder= pwd;

% ------------SET CONSTRUCT AND RESOLUTION PARAMETERS----
ConstructName='(+15)Dimer'; 
enableP2=0; 
% useCounts=0;
binSizeList= [250 1000 10000]; 
% binSizeList= [2000]; 


% set the indices to calculate and where the Autocorr to crosscorr switch
% occurs
startInd=1;
endInd=11; 
AutoToCrossInd=10; 

% % set a resolution in terms of seconds
longTau=0;
denseTau=1;

% start the process of organizing directory for AllRes loop
startFold=pwd;
homeFold = dir(pwd);
% This will be the list of directories containing (*usec) in the name
dirList=[];
resList=[];
usec= 'usec';

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
    if contains(dirList(j).name, usec)
        resList = [resList, dirList(j)];
    end
end 


for z=1:numel(binSizeList)

   folderName=(['XYQuadData_' ConstructName num2str(binSizeList(z)) 'usec']); 

   res=binSizeList(z)*1e-6; 
   resName= res/10^-6;

    cd(folderName);

% % % set the options for file types
% % appendFile=0;  
% % mapped=0;  
% % simulated=0;  

if ~longTau
% CREATE TAU ARRAY TO STEP THROUGH
NptsTCFDesired = 160;
NptsToAdd = NptsTCFDesired - 10;
tauArrayUsec = [1:9, round(logspace(1, 6.4771212,NptsToAdd))]';
tauArrayUsec = tauArrayUsec(tauArrayUsec>=res*1e6);
% tauArrayUsec = [usres:usres:10*usres]';
end

if ~longTau && ~denseTau 
% CREATE TAU ARRAY TO STEP THROUGH
NptsTCFDesired = 250;
NptsToAdd = NptsTCFDesired - 10;
tauArrayUsec = [1:9, round(logspace(1, 8,NptsToAdd))]';
tauArrayUsec = tauArrayUsec(tauArrayUsec>=res*1e6);
end 

if denseTau
NptsTCFDesired = 450;
NptsToAdd = NptsTCFDesired - 10;
tauArrayUsec = [1:9, round(logspace(1, 8,NptsToAdd))]';
tauArrayUsec = tauArrayUsec(tauArrayUsec>=res*1e6);   
end 
% This resulting Tau array will be equal at each index to the number of
% bins which seperate the two points in the TCF. 
% IN OTHER WORDS: THIS IS THE BINS TAU ARRAY, DENOTING THE SPACING OF BINS
% AT THE CURRENT RESOLTUTION TCF
tauArray = unique(floor(tauArrayUsec/(res*1e6)));

% The tauArray will be in sec, running from 0 to 3 seconds.
tauArray = [0; tauArray];

tauArraySec= tauArray*res;

% Specify an output folder for the TCFs to be saved by construct, function
% and time resolution

% NOTE: This while loop will always be true, simply in place so that the large
% body of code to create the proper directories is collapsable
l=1;

while l==1 
   
    if ~longTau && ~denseTau

outFolder1= ([num2str(resName) 'usec BF 2Pt_TCFs_' ConstructName filesep() 'LD TCFs_v2' filesep() 'P1Amp']);
outFolder2= ([num2str(resName) 'usec BF 2Pt_TCFs_' ConstructName filesep() 'LD TCFs_v2' filesep() 'P1Phase']);
outFolder3= ([num2str(resName) 'usec BF 2Pt_TCFs_' ConstructName filesep() 'LD TCFs_v2' filesep() 'Rate']);
outFolder4= ([num2str(resName) 'usec BF 2Pt_TCFs_' ConstructName filesep() 'LD TCFs_v2' filesep() 'Xquad']);
outFolder5= ([num2str(resName) 'usec BF 2Pt_TCFs_' ConstructName filesep() 'LD TCFs_v2' filesep() 'Yquad']);
outFolder6= ([num2str(resName) 'usec BF 2Pt_TCFs_' ConstructName filesep() 'LD TCFs_v2' filesep() 'P1PhaseFact']);
outFolder7= ([num2str(resName) 'usec BF 2Pt_TCFs_' ConstructName filesep() 'LD TCFs_v2' filesep() 'XcrossY']);
outFolder8= ([num2str(resName) 'usec BF 2Pt_TCFs_' ConstructName filesep() 'LD TCFs_v2' filesep() 'AmpCrossPhs']);

    end
  if longTau

outFolder1= ([num2str(resName) 'usec LongTau BF 2Pt_TCFs_' ConstructName filesep() 'LD TCFs_v2' filesep() 'P1Amp']);
outFolder2= ([num2str(resName) 'usec LongTau BF 2Pt_TCFs_' ConstructName filesep() 'LD TCFs_v2' filesep() 'P1Phase']);
outFolder3= ([num2str(resName) 'usec LongTau BF 2Pt_TCFs_' ConstructName filesep() 'LD TCFs_v2' filesep() 'Rate']);
outFolder4= ([num2str(resName) 'usec LongTau BF 2Pt_TCFs_' ConstructName filesep() 'LD TCFs_v2' filesep() 'Xquad']);
outFolder5= ([num2str(resName) 'usec LongTau BF 2Pt_TCFs_' ConstructName filesep() 'LD TCFs_v2' filesep() 'Yquad']);
outFolder6= ([num2str(resName) 'usec LongTau BF 2Pt_TCFs_' ConstructName filesep() 'LD TCFs_v2' filesep() 'P1PhaseFact']);  
outFolder7= ([num2str(resName) 'usec LongTau BF 2Pt_TCFs_' ConstructName filesep() 'LD TCFs_v2' filesep() 'XcrossY']);        
outFolder8= ([num2str(resName) 'usec LongTau BF 2Pt_TCFs_' ConstructName filesep() 'LD TCFs_v2' filesep() 'AmpCrossPhs']);        
  end 
    
    if denseTau

outFolder1= ([num2str(resName) 'usec DenseTau BF 2Pt_TCFs_' ConstructName filesep() 'LD TCFs_v3' filesep() 'P1Amp']);
outFolder2= ([num2str(resName) 'usec DenseTau BF 2Pt_TCFs_' ConstructName filesep() 'LD TCFs_v3' filesep() 'P1Phase']);
outFolder3= ([num2str(resName) 'usec DenseTau BF 2Pt_TCFs_' ConstructName filesep() 'LD TCFs_v3' filesep() 'Rate']);
outFolder4= ([num2str(resName) 'usec DenseTau BF 2Pt_TCFs_' ConstructName filesep() 'LD TCFs_v3' filesep() 'P1AmpAvg']);
% outFolder5= ([num2str(resName) 'usec DenseTau BF 2Pt_TCFs_' ConstructName filesep() 'LD TCFs_v3' filesep() 'P1PhaseAvg']);
outFolder5= ([num2str(resName) 'usec DenseTau BF 2Pt_TCFs_' ConstructName filesep() 'LD TCFs_v3' filesep() 'P1PhaseFact']);
% outFolder7= ([num2str(resName) 'usec DenseTau BF 2Pt_TCFs_' ConstructName filesep() 'LD TCFs_v3' filesep() 'XcrossY']);        
% outFolder8= ([num2str(resName) 'usec DenseTau BF 2Pt_TCFs_' ConstructName filesep() 'LD TCFs_v3' filesep() 'AmpCrossPhs']);
outFolder6= ([num2str(resName) 'usec DenseTau BF 2Pt_TCFs_' ConstructName filesep() 'LD TCFs_v3' filesep() 'P1PhaseFactAvg']);
outFolder7= ([num2str(resName) 'usec DenseTau BF 2Pt_TCFs_' ConstructName filesep() 'LD TCFs_v3' filesep() 'AvgAmpCrossRate']);          
% outFolder8= ([num2str(resName) 'usec DenseTau BF 2Pt_TCFs_' ConstructName filesep() 'LD TCFs_v3' filesep() 'AvgAmpCrossRawAmp']);
outFolder8= ([num2str(resName) 'usec DenseTau BF 2Pt_TCFs_' ConstructName filesep() 'LD TCFs_v3' filesep() 'P1AmpAvg2f']);
outFolder9= ([num2str(resName) 'usec DenseTau BF 2Pt_TCFs_' ConstructName filesep() 'LD TCFs_v3' filesep() 'P1AmpAvg3f']);
outFolder10= ([num2str(resName) 'usec DenseTau BF 2Pt_TCFs_' ConstructName filesep() 'LD TCFs_v3' filesep() 'P1AmpAvg4f']);
outFolder11= ([num2str(resName) 'usec DenseTau BF 2Pt_TCFs_' ConstructName filesep() 'LD TCFs_v3' filesep() 'P1AmpAvg5f']);

    end 

% 
% 
% 

if exist(outFolder1, 'dir')~= 7
    mkdir(outFolder1);
end


if exist(outFolder2, 'dir')~= 7
    mkdir(outFolder2);
end


if exist(outFolder3, 'dir')~= 7
    mkdir(outFolder3);
end

if exist(outFolder4, 'dir')~= 7
    mkdir(outFolder4);
end


if exist(outFolder5, 'dir')~= 7
    mkdir(outFolder5);
end


if exist(outFolder6, 'dir')~= 7
    mkdir(outFolder6);
end

if exist(outFolder7, 'dir')~= 7
    mkdir(outFolder7);
end

if exist(outFolder8, 'dir')~= 7
    mkdir(outFolder8);
end

if exist(outFolder9, 'dir')~= 7
    mkdir(outFolder9);
end

if exist(outFolder10, 'dir')~= 7
    mkdir(outFolder10);
end

if exist(outFolder11, 'dir')~= 7
    mkdir(outFolder11);
end

break
end

% Create one list of upper level files (results of XYQuad Calc). Then make
% a list for every output folders contents, when the 'case' structure below
% chooses the current quantity, cross reference the destination folder
% scanIDs with the upper level folder scanIDs
start= pwd;
% homeFold = dir(pwd);
% % This will be the list of directories containing (*usec) in the name
% dirList=[];
% resList=[];
% usec= 'usec';
% 
% % Loops over the high level directory  to find where the sub folders
% % containing the XYQUADCALE_#usec resolution folders are, ignoring the data files 
% for i=1:numel(homeFold)
%     
%     if homeFold(i).isdir
%         dirList= [dirList, homeFold(i)];
%     end
% end
% 
% % With the list of directories in the home folder, filter them to only the
% % folders cotaining the string 'usec' in order to elimate the '.' and '..'
% % directories. This is then the list of folders which contain the proper
% % XYQuadCalc folders. 
% 
% 
%     
% for j=1:numel(dirList)
%     if contains(dirList(j).name, usec)
%         resList = [resList, dirList(j)];
%     end
% end  

TrajFiles= dir('*.mat');

for k=1:numel(TrajFiles)
%     1:numel(TrajFiles)
    
    savedTrace=0; 
    fileName = TrajFiles(k).name;
    scanID = fileName([1:15]);
    disp(['Currently on File ' num2str(k) ' of ' num2str(numel(TrajFiles))]);
    load(fileName);
    
    
    
%     ****Control loop iterations here with AppendFile on to selectivley
%     alter saved contents of file. 
for i=startInd:endInd
%     for i=9:9
      tcfArr= zeros(1,numel(tauArray));
      tcfNpairs= zeros(1,numel(tauArray));
%     
%     Case Switch for the various quantities which need TCFs calculated.
%     Sets the array, name and output folder for each case prior to
%     calculating the TCF.
%  

    if i==1
        
    F=(abs(P1PhaseFact/res));
    name= 'P1Amp';
    outFolder=outFolder1;
    
    %     GENERAL BLOCK OF CODE TO CHECK DESINATION FOLDERS FOR SAVED TRACES
    cd(outFolder1);
    targetFiles=dir('*usec.mat');
    namesList=['blankblank00000'];
%     generate a list of strings of the currently saved BF2PT Files of the
%     current Quantitty F
        for l= 1:numel(targetFiles)
            namesList= [namesList; targetFiles(l).name(1:15)];
        end      
         
            for n=1:size(namesList,1)
                if contains(namesList(n, 1:15), scanID)
                    savedTrace=1;
                    break                    
                else 
                    savedTrace=0;
                end 
            end
            
            if savedTrace==1 
            cd(start);
            continue
            else
                cd(start);
            end
            
    end
    
    
    if i==2 
        
    F=(angle(P1PhaseFactN));
    name= 'P1Phase';
    outFolder=outFolder2;
    
    %     GENERAL BLOCK OF CODE TO CHECK DESINATION FOLDERS FOR SAVED TRACES
    cd(outFolder2);
    targetFiles=dir('*usec.mat');
    namesList=['blankblank00000'];
%     generate a list of strings of the currently saved BF2PT Files of the
%     current Quantitty F
        for l= 1:numel(targetFiles)
            namesList= [namesList; targetFiles(l).name(1:15)];
        end      
         
            for n=1:size(namesList,1)
                if contains(namesList(n, 1:15), scanID)
                    savedTrace=1;
                    break                    
                else 
                    savedTrace=0;
                end 
            end
            
            if savedTrace==1 
            cd(start);
            continue
            else
                cd(start);
            end
    end
    
%    
    
    if i==3
    F=(counts/res);
    name= 'Rate';
    outFolder=outFolder3;
    %     GENERAL BLOCK OF CODE TO CHECK DESINATION FOLDERS FOR SAVED TRACES
    cd(outFolder3);
    targetFiles=dir('*usec.mat');
    namesList=['blankblank00000'];
%     generate a list of strings of the currently saved BF2PT Files of the
%     current Quantitty F
        for l= 1:numel(targetFiles)
            namesList= [namesList; targetFiles(l).name(1:15)];
        end      
         
            for n=1:size(namesList,1)
                if contains(namesList(n, 1:15), scanID)
                    savedTrace=1;
                    break                    
                else 
                    savedTrace=0;
                end 
            end
            
            if savedTrace==1 
            cd(start);
            continue
            else
                cd(start);
            end
    end
    
      if i==4
         
    F=abs(P1PhaseFactN);
    name= 'P1AmpAvg';
    outFolder=outFolder4;
    %     GENERAL BLOCK OF CODE TO CHECK DESINATION FOLDERS FOR SAVED TRACES
    cd(outFolder4);
    targetFiles=dir('*usec.mat');
    namesList=['blankblank00000'];
%     generate a list of strings of the currently saved BF2PT Files of the
%     current Quantitty F
        for l= 1:numel(targetFiles)
            namesList= [namesList; targetFiles(l).name(1:15)];
        end      
         
            for n=1:size(namesList,1)
                if contains(namesList(n, 1:15), scanID)
                    savedTrace=1;
                    break                    
                else 
                    savedTrace=0;
                end 
            end
            
            if savedTrace==1 
            cd(start);
            continue
            else
                cd(start);
            end
      end
      
       
    
       if i==5
           
    F=(P1PhaseFact/res);
    name= 'P1PhaseFact';
    outFolder=outFolder5;
    %     GENERAL BLOCK OF CODE TO CHECK DESINATION FOLDERS FOR SAVED TRACES
    cd(outFolder5);
    targetFiles=dir('*usec.mat');
    namesList=['blankblank00000'];
%     generate a list of strings of the currently saved BF2PT Files of the
%     current Quantitty F
        for l= 1:numel(targetFiles)
            namesList= [namesList; targetFiles(l).name(1:15)];
        end      
         
            for n=1:size(namesList,1)
                if contains(namesList(n, 1:15), scanID)
                    savedTrace=1;
                    break                    
                else 
                    savedTrace=0;
                end 
            end
            
            if savedTrace==1 
            cd(start);
            continue
            else
                cd(start);
            end
       end
             
           
    if i==6
           
    F=(P1PhaseFactN);
    name= 'P1PhaseFactAvg';
    outFolder=outFolder6;
    %     GENERAL BLOCK OF CODE TO CHECK DESINATION FOLDERS FOR SAVED TRACES
    cd(outFolder6);
    targetFiles=dir('*usec.mat');
    namesList=['blankblank00000'];
%     generate a list of strings of the currently saved BF2PT Files of the
%     current Quantitty F
        for l= 1:numel(targetFiles)
            namesList= [namesList; targetFiles(l).name(1:15)];
        end      
         
            for n=1:size(namesList,1)
                if contains(namesList(n, 1:15), scanID)
                    savedTrace=1;
                    break                    
                else 
                    savedTrace=0;
                end 
            end
            
            if savedTrace==1 
            cd(start);
            continue
            else
                cd(start);
            end
    end
    
     if i==7
       
    F=abs(P1PhaseFactN2f);
    name= 'P1AmpAvg2f';
    outFolder=outFolder8;
    %     GENERAL BLOCK OF CODE TO CHECK DESINATION FOLDERS FOR SAVED TRACES
    cd(outFolder8);
    targetFiles=dir('*usec.mat');
    namesList=['blankblank00000'];
%     generate a list of strings of the currently saved BF2PT Files of the
%     current Quantitty F
        for l= 1:numel(targetFiles)
            namesList= [namesList; targetFiles(l).name(1:15)];
        end      
         
            for n=1:size(namesList,1)
                if contains(namesList(n, 1:15), scanID)
                    savedTrace=1;
                    break                    
                else 
                    savedTrace=0;
                end 
            end
            
            if savedTrace==1 
            cd(start);
            continue
            else
                cd(start);
            end
     end
    
     if i==8
       
    F=abs(P1PhaseFactN3f);
    name= 'P1AmpAvg3f';
    outFolder=outFolder9;
    %     GENERAL BLOCK OF CODE TO CHECK DESINATION FOLDERS FOR SAVED TRACES
    cd(outFolder9);
    targetFiles=dir('*usec.mat');
    namesList=['blankblank00000'];
%     generate a list of strings of the currently saved BF2PT Files of the
%     current Quantitty F
        for l= 1:numel(targetFiles)
            namesList= [namesList; targetFiles(l).name(1:15)];
        end      
         
            for n=1:size(namesList,1)
                if contains(namesList(n, 1:15), scanID)
                    savedTrace=1;
                    break                    
                else 
                    savedTrace=0;
                end 
            end
            
            if savedTrace==1 
            cd(start);
            continue
            else
                cd(start);
            end
     end
    
  if i==9
       
    F=abs(P1PhaseFactN4f);
    name= 'P1AmpAvg4f';
    outFolder=outFolder10;
    %     GENERAL BLOCK OF CODE TO CHECK DESINATION FOLDERS FOR SAVED TRACES
    cd(outFolder10);
    targetFiles=dir('*usec.mat');
    namesList=['blankblank00000'];
%     generate a list of strings of the currently saved BF2PT Files of the
%     current Quantitty F
        for l= 1:numel(targetFiles)
            namesList= [namesList; targetFiles(l).name(1:15)];
        end      
         
            for n=1:size(namesList,1)
                if contains(namesList(n, 1:15), scanID)
                    savedTrace=1;
                    break                    
                else 
                    savedTrace=0;
                end 
            end
            
            if savedTrace==1 
            cd(start);
            continue
            else
                cd(start);
            end
            
  end
    
  
   if i==10
       
    F=abs(P1PhaseFactN5f);
    name= 'P1AmpAvg5f';
    outFolder=outFolder11;
    %     GENERAL BLOCK OF CODE TO CHECK DESINATION FOLDERS FOR SAVED TRACES
    cd(outFolder11);
    targetFiles=dir('*usec.mat');
    namesList=['blankblank00000'];
%     generate a list of strings of the currently saved BF2PT Files of the
%     current Quantitty F
        for l= 1:numel(targetFiles)
            namesList= [namesList; targetFiles(l).name(1:15)];
        end      
         
            for n=1:size(namesList,1)
                if contains(namesList(n, 1:15), scanID)
                    savedTrace=1;
                    break                    
                else 
                    savedTrace=0;
                end 
            end
            
            if savedTrace==1 
            cd(start);
            continue
            else
                cd(start);
            end
    end
       
    if i==11
       
    F=abs(P1PhaseFactN);
    G=(counts/res);
    name= 'AvgAmpCrossRate';
    outFolder=outFolder7;
    %     GENERAL BLOCK OF CODE TO CHECK DESINATION FOLDERS FOR SAVED TRACES
    cd(outFolder7);
    targetFiles=dir('*usec.mat');
    namesList=['blankblank00000'];
%     generate a list of strings of the currently saved BF2PT Files of the
%     current Quantitty F
        for l= 1:numel(targetFiles)
            namesList= [namesList; targetFiles(l).name(1:15)];
        end      
         
            for n=1:size(namesList,1)
                if contains(namesList(n, 1:15), scanID)
                    savedTrace=1;
                    break                    
                else 
                    savedTrace=0;
                end 
            end
            
            if savedTrace==1 
            cd(start);
            continue
            else
                cd(start);
            end
    end
  
  if i==12
       
    F=abs(P1PhaseFactN);
    G=abs(P1PhaseFact);
    name= 'AvgAmpCrossRawAmp';
    outFolder=outFolder8;
    %     GENERAL BLOCK OF CODE TO CHECK DESINATION FOLDERS FOR SAVED TRACES
    cd(outFolder8);
    targetFiles=dir('*usec.mat');
    namesList=['blankblank00000'];
%     generate a list of strings of the currently saved BF2PT Files of the
%     current Quantitty F
        for l= 1:numel(targetFiles)
            namesList= [namesList; targetFiles(l).name(1:15)];
        end      
         
            for n=1:size(namesList,1)
                if contains(namesList(n, 1:15), scanID)
                    savedTrace=1;
                    break                    
                else 
                    savedTrace=0;
                end 
            end
            
            if savedTrace==1 
            cd(start);
            continue
            else
                cd(start);
            end
  end
    
%        create an if statement for the case of the auto corr and cross corr
%        functions
if i<=AutoToCrossInd 
%     For each file in the list, will need to calculated a TCF. so must
%     loop over the array of Tau values and then step the starting time
%     over the full range of the trajectory. 
    avg= nanmean(F);
    varOfArray= nanvar(F);
 if strcmp(name, 'P1Amp') || strcmp(name, 'P1PhaseFact')
    for j=1:numel(tauArray)
%         this is the current value of Tau in seconds
       curTau=tauArray(j);
       
       for t=1:numel(F)-curTau
%          if strcmp(name, 'P1Amp') || strcmp(name, 'P1PhaseFact')

             
             if (counts(t)==0)||(counts(t+curTau)==0)
               continue
             else
             TCFtemp = (conj(F(t))-conj(avg))*(F(t+curTau)-avg);
             tcfArr(j) = tcfArr(j) + TCFtemp; 
             tcfNpairs(j) = tcfNpairs(j)+1; 
             end
           
       end
       
    end
    
  else
             
    for j=1:numel(tauArray)
%         this is the current value of Tau in seconds
       curTau=tauArray(j);
       
       for t=1:numel(F)-curTau
           if isnan(F(t)) || isnan(F(t+curTau))
%            if (counts(t)<2)||(counts(t+curTau)<2)
           continue
           else
             TCFtemp = (conj(F(t))-conj(avg))*(F(t+curTau)-avg);
             tcfArr(j) = tcfArr(j) + TCFtemp; 
             tcfNpairs(j) = tcfNpairs(j)+1; 
           end
           
       end
     end
%         
    
    
    
   end
    

    
   elseif i>AutoToCrossInd 
          
    avgF= nanmean(F);
    avgG= nanmean(G);
    varOfArray= nanvar(F);
  if strcmp(name, 'AvgAmpCrossRawAmp')
    for j=1:numel(tauArray)
%         this is the current value of Tau in seconds
       curTau=tauArray(j);      
          
       for t=1:numel(G)-curTau
           if (counts(t)==0)||(counts(t+curTau)==0)
           continue
           else
             TCFtemp = (conj(F(t))-conj(avgF))*(G(t+curTau)-avgG);
             tcfArr(j) = tcfArr(j) + TCFtemp; 
             tcfNpairs(j) = tcfNpairs(j)+1; 
           end
%          
       end
    end
    
  else
      
     for j=1:numel(tauArray)
%         this is the current value of Tau in seconds
       curTau=tauArray(j);  
       for t=1:numel(G)-curTau
           if isnan(F(t))||isnan(G(t+curTau))
           continue
           else
             TCFtemp = (conj(F(t))-conj(avgF))*(G(t+curTau)-avgG);
             tcfArr(j) = tcfArr(j) + TCFtemp; 
             tcfNpairs(j) = tcfNpairs(j)+1; 
           end
%          
       end
       
      end
%            
    end
    
   end
% Now the tcfArr should be filled and ready to save to its appropiate
% folder. weight the TCF by the number of pairs which contributed to that
% value of tau

tcfArr= tcfArr./tcfNpairs;


% *****BLOCK TO APPEND EXISTIN FILES WITH THE VARIANCE OF THE UNDERLYING
% ARRAY BEING OPERATED ON ********
% foutName = [outFolder filesep() scanID name '_2PtTCFs_' num2str(resName) '_usec.mat'];
% save(foutName,'varOfArray', '-append'); 

    
foutName = [outFolder filesep() scanID name '_2PtTCFs_' num2str(resName) '_usec.mat'];
save(foutName,'tcfArr', 'tcfNpairs','tauArraySec','tauArray', 'trajLength','varOfArray');


% if appendFile==1
%     save(foutName,'tcfAmpTest','-append'); 
%     
% else 
%     save(foutName,'tcfAmpTest','tcfArr','tcfNpairs','tauArraySec','tauArray', 'trajLength','varOfArray');
% end 

tracesCalcd = tracesCalcd + 1;

    end

end

cd(startFold);

toc
end 
%  clear all
