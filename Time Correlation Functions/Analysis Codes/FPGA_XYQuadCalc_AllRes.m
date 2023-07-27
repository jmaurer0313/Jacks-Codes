% GOAL: generate data traces at the specified resolutions for FPGA data
% streams, in the same fashion as 'XYQuadCalc_AllResAtOnce'


tic

% disp(['Your current directory is' string(pwd)]);
traceFolder= pwd;
% set a bin size in uS for the phase histograms to be looped over
% binSizeList=[500];
binSizeList=[5000 8000 12000 15000];
% binSizeList=[1000000];
ConstructName= 'E3E4';
enableP2=0;
unwrapPhase=1; 

% temporary bin size list to test out factor of 2 in complex exponentials 
% binSizeList=[10000 1000 100 10];

numbRes= numel(binSizeList);
folderSize= numel(traceFolder);

% bin= 0:63;
TrajFiles= dir('*cut.mat');


start = pwd;


for d=1:numel(binSizeList)
    if unwrapPhase==0
folderName=(['XYQuadData_NuW_' ConstructName num2str(binSizeList(d)) 'usec']);
    else 
folderName=(['XYQuadData_' ConstructName num2str(binSizeList(d)) 'usec']); 
    end

% Make a new folder to store histogram counts if it doesn't exist already
if exist(folderName, 'dir')~= 7
    mkdir(folderName);
end
end 
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


% Now with only the res folders, The resulting 'resList' will be sorted
% automaticallty from the longest name to shortest name. So Assuming that
% 1-10000 usec data is available, the list will be 1usec at the last index,
% 10usec at the second to last index, and so on. 

numFolds= numel(resList);

% oneUsecFold = resList(numFolds).name;
% tenUsecFold = resList(numFolds-1).name;
% hundUsecFold= resList(numFolds-2).name;
% thouUsecFold= resList(numFolds-3).name;




for z=1:numel(binSizeList)
%     numel(binSizeList)
    binSize=binSizeList(z); 
    
    if unwrapPhase==0
folderName=(['XYQuadData_NuW_' ConstructName num2str(binSizeList(z)) 'usec']);
    else 
folderName=(['XYQuadData_' ConstructName num2str(binSizeList(z)) 'usec']); 
    end
%     folderName=(['XYQuadData_' ConstructName num2str(binSizeList(z)) 'usec']);
%     This will start with the 1000usec data folder in this case (works
%     in oppisite direction as above) 
    curXYFolder=resList(z).name; 
%     cd(curXYFolder);
%     targetFiles=dir('*.mat');
%     
%     namesList=[];
% %     generate a list of strings of the currently saved XYQuadCalc Files
%         for i= 1:numel(targetFiles)
%             namesList= [namesList; targetFiles(i).name(1:12)];
%         end 
%         
% %         change back to raw data directory from the curXYFolder 
%         cd ..

%     change dir to the curXY Folder and generate a list of all the SCAN
%     ids to check against the TrajFiles List 

    



for f=1:(numel(TrajFiles))
%     (numel(TrajFiles))
%     
    binNum1=1;
  
%    
    disp(['Currently on File ' num2str(f) ' of ' num2str(numel(TrajFiles))]);
%     should be k=1:numel(ch1Files) in order to go over all the files in
%     the folder. Made it 1-2 to restrict the run time and plot density
    fileName = TrajFiles(f).name;
    scanID = fileName(1:15);
    
%     run over the namesList of the saved data and check to see that this
%     file is not already contained in that list for operation
        savedTrace=0; 
%     for n=1:size(namesList,1)
%         if contains(namesList(n, 1:12), scanID)
%         savedTrace=1;
%         break
%         else
%         savedTrace=0;
%         end 
%     end 
    
    if savedTrace==0
        
    load(fileName);
   
    if isempty(p1) || isempty(p2) 
        continue
    end 
%     Subtract the minimum time from all elements of the hitList in order
%     to have the intialized arrays from (1,Bins) be matched in terms of
%     the starting times and the relevant bin ranges which will be covered
%     by determination of the trajLength. 
    timesUsec= timesUsec - (min(timesUsec));
  
    
%      *****THIS IS THE NEW REPHASE METHOD*******
% This method takes the X and Y quadrature for the ENTIRE TRACE, which
% should reflect the overall average phase as a value away from zero phase.
% By subtracting off this angle from all phases within the list the
% distribution of phases will be centered at zero on the interval of -pi to
% pi whicch atan2 returns. This minimizes "toggling" erros across -Theta to
% Theta in neighboring bins. 
% 
% ONLY NEEDED FOR P1 P2 SINCE FPGA P1MAP AND P2MAP ARE RECENTERED PRIOR
  XquadP1 = mean(cos(p1.*2*pi));
  YquadP1 = mean(sin(p1.*2*pi));
   
  if enableP2
   XquadP2 = mean(cos(p2.*2*pi));
   YquadP2 = mean(sin(p2.*2*pi));
   
   shiftAngle2=atan2(YquadP2,XquadP2);
   
   p2 = (p2.*2*pi)-(shiftAngle2); 
    curTrajP2= [timesUsec, p2];
  end
%    could be something about going from 0->2pi to -angle -> angle, if
%    issues arise could be from this subtraction 

if unwrapPhase==1
   shiftAngle1=atan2(YquadP1,XquadP1);
   
   p1 = (p1.*2*pi)-(shiftAngle1);
   
else
    p1 = (p1.*2*pi);
end 
   
   curTrajP1= [timesUsec, p1];
      
       
if size(p1Map,1)<size(p1Map,2)
       p1Map=p1Map.';
       if enableP2
       p2Map=p2Map.';
       end
end 

       curTrajP1Map= [timesUsec, p1Map];
       if enableP2
       curTrajP2Map= [timesUsec, p2Map];
       end 
               
       
% create an empty array for X and Y total by dividing the length of the current
% traj by the time width of the window setting. WindSize will be in usec.
% Set the intial X/Y value to zero and the counts contributing to zero. 
trajLength=(max(timesUsec)-min(timesUsec));
Bins=floor(trajLength/binSize);
    
%     Ch1X= zeros(1,Bins); 
%     Ch1Y= zeros(1,Bins);
%     Ch2X= zeros(1,Bins); 
%     Ch2Y= zeros(1,Bins); 
    P1PhaseFact= zeros(1,Bins);
    P1PhaseFact2f= zeros(1,Bins);
    P1PhaseFact3f= zeros(1,Bins);
    P1PhaseFact4f= zeros(1,Bins);
    P1PhaseFact5f= zeros(1,Bins);
    P1PhaseFactMap= zeros(1,Bins);
    
    if enableP2
    P2PhaseFact= zeros(1,Bins);    
    P2PhaseFactMap= zeros(1,Bins);
    tempP2PhFct=0;
    tempP2PhFctMap=0;
    end 
    
    counts=zeros(1,Bins);
  
% set the temp holders to zero for the start of a new file/trajectory
%   
    tempCounts=0;
    tempP1PhFct=0;
    tempP1PhFct2f=0;
    tempP1PhFct3f=0;
    tempP1PhFct4f=0;
    tempP1PhFct5f=0;
    tempP1PhFctMap=0;
    
    
 

    for j=1:(numel(timesUsec))
%         disp(['Current Ch1 iteration is ' int2str(j)]);
        if(curTrajP1(j,1)<(binNum1*binSize))
        tempCounts=tempCounts+1;
        tempP1PhFct=tempP1PhFct + exp(1i*(curTrajP1(j,2)));
        tempP1PhFct2f = tempP1PhFct2f + exp(2i*(curTrajP1(j,2)));
        tempP1PhFct3f = tempP1PhFct3f + exp(3i*(curTrajP1(j,2)));
        tempP1PhFct4f = tempP1PhFct4f + exp(4i*(curTrajP1(j,2)));
        tempP1PhFct5f = tempP1PhFct5f + exp(5i*(curTrajP1(j,2)));
        tempP1PhFctMap=tempP1PhFctMap + exp(1i*(curTrajP1Map(j,2)));
         
        if enableP2
        tempP2PhFct=tempP2PhFct + exp(1i*(curTrajP2(j,2)));       
        tempP2PhFctMap=tempP2PhFctMap + exp(1i*(curTrajP2Map(j,2)));
        end
           
        else                   
            
%           
           counts(binNum1) = tempCounts;
           P1PhaseFact(binNum1)= tempP1PhFct;   
           P1PhaseFact2f(binNum1)= tempP1PhFct2f;  
           P1PhaseFact3f(binNum1)= tempP1PhFct3f; 
           P1PhaseFact4f(binNum1)= tempP1PhFct4f; 
           P1PhaseFact5f(binNum1)= tempP1PhFct5f; 
           P1PhaseFactMap(binNum1)= tempP1PhFctMap;
           
           if enableP2
                P2PhaseFact(binNum1)= tempP2PhFct;
                P2PhaseFactMap(binNum1)= tempP2PhFctMap;
                tempP2PhFct=exp(1i*(curTrajP2(j,2)));
                tempP2PhFctMap=exp(1i*(curTrajP2Map(j,2)));
           end 
           
%            take the bin num to the proper bin which the current time value falls 
%            in rather than update the counter by 1
            binNum1= floor(curTrajP1(j,1)/binSize)+1;
        
        tempCounts=1;
        tempP1PhFct= exp(1i*(curTrajP1(j,2))); 
        tempP1PhFct2f= exp(2i*(curTrajP1(j,2)));
        tempP1PhFct3f= exp(3i*(curTrajP1(j,2)));
        tempP1PhFct4f= exp(4i*(curTrajP1(j,2)));
        tempP1PhFct5f= exp(5i*(curTrajP1(j,2)));
        tempP1PhFctMap= exp(1i*(curTrajP1Map(j,2)));
        
                                
        end
        
    end 
    
 
    
%  create the total arrays for all quantities and save them to a folder


% FIX 2/19/2019: Given that the total counts can be zero at the lowest time
% resolutions, the instances of divide by zero must be corrected for to
% eliminate the NaNs. the command XtotalN(isnan(XtotalN))=0; will replace
% all instanceso of NaN in the array being operated on with zero. Could in
% principal also replace Nans with the average of the array. But with the X
% and Y resultant components, the component should be zero if no counts
% were recorded. 
if enableP2
    P2PhaseFactN= P2PhaseFact./counts; 
     P2PhaseFactMapN= P2PhaseFactMap./counts; 
end

 P1PhaseFactN= P1PhaseFact./counts;  
 P1PhaseFactN2f= P1PhaseFact2f./counts; 
 P1PhaseFactN3f= P1PhaseFact3f./counts; 
 P1PhaseFactN4f= P1PhaseFact4f./counts; 
 P1PhaseFactN5f= P1PhaseFact5f./counts; 
 P1PhaseFactMapN= P1PhaseFactMap./counts;  

    
    foutName = [folderName filesep() scanID '_XYQuadData_ ' num2str(binSize) '_usec.mat'];   
    if enableP2
    save(foutName,'trajLength' ,'counts',  'P1PhaseFactN', 'P2PhaseFactN','P1PhaseFactMapN','P2PhaseFactMapN','P1PhaseFact','P2PhaseFact','P1PhaseFactMap','P2PhaseFactMap');     
    else
     save(foutName,'trajLength' ,'counts',  'P1PhaseFactN', 'P1PhaseFactMapN','P1PhaseFact','P1PhaseFactMap','P1PhaseFactN2f', 'P1PhaseFactN3f',  'P1PhaseFactN4f', 'P1PhaseFactN5f');     
    end 
    binNum1=1;
  
    
     
    else 
        continue
   
    end
    
 end


end 
toc