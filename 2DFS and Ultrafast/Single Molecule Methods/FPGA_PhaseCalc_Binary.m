% GOAL: process data files from the FPGA using the frameowrk set up y XY
% quad calc taking in arbitrary timesadn phases. 

% **********BE IN THE DEFRAG FOLDER FOR THE STREAM DATA***********

tic

% ********INPUT PARAMETERS************
 binSizeList=[1000000 100000 10000 1000 300];
% binSizeList=[300];
ConstructName= 'E5E6';
% ************************************




disp(['Your current directory is' string(pwd)]);
traceFolder= pwd;
filesList= dir('*-Stream');
mainDir=dir; 


for d=1:numel(binSizeList)
%     switch to the upper level directory to make the output folder
%     cd(mainDir(1).folder)
folderName=(['FPGAPhaseCalc_' ConstructName '_' num2str(binSizeList(d)) 'usec']);

% Make a new folder to store histogram counts if it doesn't exist already
if exist (folderName, 'dir')~= 7
    mkdir(folderName);
end
% switch back to data file folder 
%     cd(tempDir(1).folder)
end 



% Create a loop to search over the directory and load in a new file name
% with each iteration as the full path to the TDMS file
for z=1:numel(binSizeList)
    binSize=binSizeList(z);
    % ***********INILIZATION FOR THE TCF CALC**********
res=binSizeList (z)/1e6;  
resName= res/10^-6; 
l=1;  

% CREATE TAU ARRAY TO STEP THROUGH
NptsTCFDesired = 160;
NptsToAdd = NptsTCFDesired - 10;
tauArrayUsec = [1:9, round(logspace(1, 6.4771212,NptsToAdd))]';
tauArrayUsec = tauArrayUsec(tauArrayUsec>=res*1e6);
% tauArrayUsec = [usres:usres:10*usres]';

% This resulting Tau array will be equal at each index to the number of
% bins which seperate the two points in the TCF. 
% IN OTHER WORDS: THIS IS THE BINS TAU ARRAY, DENOTING THE SPACING OF BINS
% AT THE CURRENT RESOLTUTION TCF
tauArray = unique(floor(tauArrayUsec/(res*1e6)));

% The tauArray will be in sec, running from 0 to 3 seconds.
tauArray = [0; tauArray];
tauArraySec= tauArray*res;

folderName=(['FPGAPhaseCalc_' ConstructName '_' num2str(binSizeList(z)) 'usec']);



    for f=1:numel(filesList)
%       for f=3:3
binNum1=1; 
 disp(['Currently on File ' num2str(f) ' of ' num2str(numel(filesList))]);
% mainDir and Fileslisr will ne offset by 2 in the loop to account for the
% '.' and '..' directories in mainDir
targetFolder=[traceFolder,'\',filesList(f).name,'\'];
scanID = filesList(f).name(1:15);

timeID = fopen([targetFolder 'time.bin']);
p1ID = fopen([targetFolder 'p1.bin']);
p2ID = fopen([targetFolder 'p2.bin']);

time = fread(timeID,Inf,'uint64=>uint64',0,'s');
p1 = fread(p1ID,Inf,'float64=>double',0,'s');
p2 = fread(p2ID,Inf,'float64=>double',0,'s');


% % just for a bad run with bubbles (set up for 2nd run on 10/17 on E5E6)
if f==2
    IbadTimes = time>=26612215345+1063636;
    time(IbadTimes)=[];
    p1(IbadTimes)=[];
    p2(IbadTimes)=[];
end
    

% Since the sorting script starts at binNum=1, all times must be 1 greater
% to avoid indexing error
times = ((double (time-time(1))./(8e7))+1);
%  if f==2
%      times=times(1:round(0.8*length(time))
p1r = (2.*pi).*p1;
p2r = (2.*pi).*p2;
% histogram(p1,40)
if numel (times) ~= numel(p1r)
           disp(['Error, number of times and phases dont match in file ' scanID]); 
           cd(mainDir(1).folder);
       continue 
       end 
% plot(times, p1r)

 usecHitList_ch1= times*(1e6);  
%  curTrajCh1= [usecHitList_ch1.', p1r.'];
 curTrajCh1= [usecHitList_ch1, p1r];
 trajLength= times(end); 
 

    
    
     
    
disp(['Currently resetting counters for file ' num2str(f)]);
 Bins = floor(trajLength/binSize);
 Ch1PhaseFact= zeros(1,Bins);
 Ch1Counts=zeros(1,Bins);

    
%     Subtract the minimum time from all elements of the hitList in order
%     to have the intialized arrays from (1,Bins) be matched in terms of
%     the starting times and the relevant bin ranges which will be covered
%     by determination of the trajLength. 
   
       
% create an empty array for X and Y total by dividing the length of the current
% traj by the time width of the window setting. WindSize will be in usec.
% Set the intial X/Y value to zero and the counts contributing to zero. 

   

       
%       Take all the events in the first window set by windSize and create
%       the above quanties for each and every bin

% set the temp holders to zero for the start of a new file/trajectory
   
    tempCh1Counts=0;
  
    tempCh1PhFct=0;

    for j=1:(numel(usecHitList_ch1))
%         disp(['Current Ch1 iteration is ' int2str(j)]);
        if(curTrajCh1(j,1)<(binNum1*binSize)) 
            tempCh1PhFct = tempCh1PhFct + exp(1i*(curTrajCh1(j,2)));
            tempCh1Counts = tempCh1Counts + 1; 

        else                   
            
           
           Ch1Counts(binNum1) = tempCh1Counts;
           Ch1PhaseFact(binNum1)= tempCh1PhFct;
%         
            
%            take the bin num to the proper bin which the current time value falls 
%            in rather than update the counter by 1
            binNum1= floor(curTrajCh1 (j,1)/binSize)+1;
%             j=j+1;
          
            tempCh1PhFct= exp(1i*(curTrajCh1(j,2)));
            tempCh1Counts=1;          
                                
        end
        
    end
    disp(['Currently finishing Data Sorting for file' num2str(f)]);
    
%     ******* BEGIN TCF CALC********
for r=1:3
    if r==1
    F=Ch1PhaseFact; 
    end
    if r==2
    F=abs(Ch1PhaseFact);
    end 
    if r==3
    F=angle(Ch1PhaseFact); 
    end 
    
    avg= nanmean(F);
    varOfArray= nanvar(F);
    tcfArr = zeros(1,numel(tauArray));
    tcfNpairs = zeros(1,numel(tauArray));
    
    for j=1:numel(tauArray)
%         this is the current value of Tau in seconds
       curTau=tauArray(j);
       
       for t=1:numel(F)-curTau
           if isnan(F(t))||isnan(F(t+curTau))
           continue
           else
             TCFtemp = (conj(F(t))-conj(avg))*(F(t+curTau)-avg);
             tcfArr(j) = tcfArr(j) + TCFtemp; 
             tcfNpairs(j) = tcfNpairs(j)+1; 
           end
%          NOTE: learned that for loops need no iteration update! will run
%          through the values 1 at a time. While loops require iteration
%          update.
%          t=t+1;
       end
%         j=j+1;
    
    
    
    end
    
    
    Ch1tcfArr= tcfArr./tcfNpairs;
    if r==1
    PhaseFactTCF=Ch1tcfArr; 
    PhaseFactNpairs= tcfNpairs; 
    end
    if r==2
    AbsTCF=Ch1tcfArr;
    AbsNpairs= tcfNpairs; 
    end 
    if r==3
    AngleTCF=Ch1tcfArr;
    AngleNpairs= tcfNpairs; 
    end
    
    disp(['Currently Finishing TCF Calc for file ' num2str(f)]);
    
end 
timeAtRes = resName.*(0:(numel(Ch1Counts)-1));

cd(mainDir(1).folder);

foutName = [folderName filesep() scanID '_FPGAPhaseCalc _ ' num2str(binSize) '_usec.mat'];   
    save(foutName,'trajLength' , 'AngleNpairs', 'AbsNpairs','PhaseFactNpairs','Ch1Counts', 'Ch1PhaseFact','AngleTCF','AbsTCF','PhaseFactTCF','times','p1r','tauArray','tauArraySec','tauArrayUsec', 'timeAtRes');
  
%     if  z < numel(binSizeList)  
% 	cd(tempDir(1).folder); 
%     else 
%       cd(mainDir(1).folder);
%     end
    
end

end

        
    