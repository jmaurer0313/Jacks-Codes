% c4 stitch plotter based on set of input resolutions

constructName='(+15)Dimer';
quantityName='P1AmpAvg';
% quantityName='Rate';
saveMode=0;
normalizeMode=0;
% make the tau2 label obey this conevtion to match the folders within each
% directory of 4PtTCFs -- this way the stitchplotting can happen for an
% abritray tau2  --always 6 numbers long in the format "000100us"
tau2usec='000000us'; 
calcZeroTau=1;
removePlotDetails=1;

% User Options
% SET OF RESOLUTIONS TO STITCH - BE IN ASCENDING ORDER
binSizeList=[250]; 
% *****************************************************************
% BE IN THE COMBINED DATA FOLDER (OR WHICHEVER FOLDER HOLDS ALL THE XYQUAD
% FOLDERS AT EVERY RESOLUTION) -- THIS WILL BE THE "STARTFOLD"


generateTaus=0; 
% startRes=250; % in usec, choices of 10, 100 (for now)
subsetOpt=0; % used to look at Uright for now, others to be added. turning off will look at all traces


startFold=pwd;
homeFold = dir(pwd);
% This will be the list of directories containing (*usec) in the name
dirList=[];
resList=[];
targetList={};
if strcmp(quantityName,'FRET')
keyWord= 'FRET';   
else
keyWord= 'XYQuadData';
end

if generateTaus
% IF NEEDED, THE TOTALTAUARRAY CAN BE REGENRATED HERE TO MATCH THE DATA
% FILES... FIRST RUN OF THE PAR CODE WASN'T SAVING THE TOTALTAUARRAY
  NptsTCFDesired = 325;
    NptsToAdd = NptsTCFDesired - 200;
    tauArray_msec = ([1:9, round(logspace(1, 3,NptsToAdd))]')*1000;
    tauArray_msec = tauArray_msec(tauArray_msec>=1000);
    tauArray_msec = (unique(tauArray_msec/1000))*1000;
    tauArrayUsecTotal = [linspace(10,100,10), linspace(100,1000,10),(tauArray_msec)']';
%     all based on starting at 10usec resolution for the stitching...
    tauArrayTotal = unique(floor(tauArrayUsecTotal/(10)));
    tauArrayUsecTotal=tauArrayTotal*10;
    tauArraySecTotal= tauArrayUsecTotal*1e-6;
end
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

for m=1:numel(binSizeList)
    for i=1:numFolds    
    if contains(resList(i).name, [num2str(binSizeList(m)) 'usec'])
        targetList=[targetList; resList(i).name];
    end
    end
    
end
%%


for j=1:numel(binSizeList)
   
    cd(targetList{j});
    
    curDir=dir(pwd);
    FourPtFolder=[];
    for m=1:numel(curDir)
        
        if  ~calcZeroTau
            if contains(curDir(m).name, '3ptTCFs') && contains(curDir(m).name, quantityName) && ~contains(curDir(m).name,'calcZeroTau')
            FourPtFolder = curDir(m).name;
            end
        elseif calcZeroTau
            if contains(curDir(m).name, '3ptTCFs_SL') && contains(curDir(m).name, quantityName) && contains(curDir(m).name,'calcZeroTau')
            FourPtFolder = curDir(m).name;
            end
        end
    end 
    cd(FourPtFolder);
    
    
    tau2Fold=[];
    curDir2=dir(pwd);
    for m=1:numel(curDir2)
        if contains(curDir2(m).name, tau2usec)
        tau2Fold = curDir2(m).name;
        end
    end 
    
    cd(tau2Fold);
    
    tcfFiles=dir('*fourptTCF.mat'); 
%     load in the first file to set the size of the array
    load(tcfFiles(1).name);
    total4Pt=zeros(numel(tau1L),numel(tau3L));
    totalPairs=zeros(numel(tau1L),numel(tau3L));
    
    for i=1:numel(tcfFiles)

    load(tcfFiles(i).name);
    total4Pt=total4Pt+(fourptTCF.*NpairsTau_fourptTCF); 
    totalPairs=totalPairs+NpairsTau_fourptTCF; 
    
    end
    
    if j==1
     firstSurf=total4Pt./totalPairs;  
     if numel(binSizeList)>1
     firstTauArray=tauArrayUsecTotal; 
     else
      firstTauArray=tau1L;   
     end
    elseif j==2
     secondSurf=total4Pt./totalPairs;  
     secondTauArray=tau1L*1e-6;
   elseif j==3
     thirdSurf=total4Pt./totalPairs;  
    end
    
    cd(startFold);
end

% setup the index sets to match the surfaces
NpointsRatio=15;
% this will be 1:N in the hund but x:N in the ten (to match at 100usec tau)
if numel(binSizeList)==2
    
%     first find the inde at which the first surface meets up with the 2nd.
%     stitch from here on. 
firstStitchIndices = find(firstTauArray<=binSizeList(2));
if calcZeroTau
firstStitchIndices=firstStitchIndices(2:end);
end
firstSurfIndex = firstStitchIndices(end);


NRatios=0;
finalRatio=0;

for i=0:NpointsRatio-1
    if secondTauArray(1)~=0
    temp1=firstSurf(firstSurfIndex+i,firstSurfIndex)/secondSurf(1+i,1);
    temp2=firstSurf(firstSurfIndex,firstSurfIndex+i)/secondSurf(1,1+i);
    elseif secondTauArray(1)==0
    temp1=firstSurf(firstSurfIndex+i,firstSurfIndex)/secondSurf(2+i,2);
    temp2=firstSurf(firstSurfIndex,firstSurfIndex+i)/secondSurf(2,2+i);   
    end
    finalRatio=finalRatio+temp1+temp2;
    NRatios=NRatios+2; 
    
    
end

if calcZeroTau
    if secondTauArray(1)==0
   secondSurf=secondSurf((2:end),(2:end)); 
    end
end

scaleFactor=finalRatio/NRatios; 
secondSurf=secondSurf.*scaleFactor;
firstSurf(firstSurfIndex:end,firstSurfIndex:end)=secondSurf;
C4finalArray=firstSurf;

end

if numel(binSizeList)==1    
    C4finalArray=firstSurf;    
end
    

if calcZeroTau
finalPlot=C4finalArray((2:end),(2:end));
C4finalArray=C4finalArray((2:end),(2:end)); 
else
finalPlot=C4finalArray;   
end
if normalizeMode
finalPlot=finalPlot./(abs(finalPlot(1,1)));
end
C4finalArrayNorm=C4finalArray./(C4finalArray(1,1));

%% Set up the figure for plotting
figure();
clf
set(gcf,'Color','w');
set(gcf,'Name','3-point Time Correlation Function');

%% Assign variables
% load in the first file to get the time axes
% load(tcfFiles(1).name);
if ~generateTaus
    tauArraySecTotal=tauArrayUsecTotal*1e-6;
end 

time= firstTauArray*1e-6; 
if time(1)==0
    time=time(2:end); 
end
% time = tau1arrayUsec *1e-6
% tau2 = tauArraySecTotal;
% tau2 = tau2ValUsec*1e-6;
% C4 = FourPtTCF_avg;

% finalPlot

%% Plot the 4point TCF (birds eye view)
% NEED TO DETERMINE HOW TO HANDLE THE COMPLEXITY IN 4PT CALCULATION
% surf(time(2:end),time(2:end),real(C4(2:end,2:end)));
surf(time,time,real(finalPlot));
shading interp;
colormap('jet');
% caxis([0,1]);
caxis([min(min(real(finalPlot))),max(max(real(finalPlot)))]);

hold on;
contour3(time,time,real(finalPlot),5,'k','LineWidth',2);



title_str = [constructName ' 3-point TCF: \tau_2 = ' tau2usec ' N=' num2str(numel(tcfFiles))];
title(title_str,'FontSize',12);
xlabel('\tau_1 (sec)','fontsize',12);
ylabel('\tau_3 (sec)','fontsize',12);
zlabel('C^{(3)}(\tau_1,\tau_2,\tau_3)','fontsize',12);
set(gca,'yscale','log');
set(gca,'xscale','log');
% set(gca,'zscale','log');
set(gca,'FontSize',12);
grid on
axis image
axis tight;
colorbar();
xlim([-inf,time(end)]);
ylim([-inf,time(end)]);
zlim([min(min(real(finalPlot))),max(max(real(finalPlot)))]);
view(0,90);


%% Set up the figure for plotting
figure();
clf
set(gcf,'Color','w');
set(gcf,'Name','3-point Time Correlation Function');

%% Plot the 4point TCF (angled view)
% NEED TO DETERMINE HOW TO HANDLE THE COMPLEXITY IN 4PT CALCULATION
% surf(time(2:end),time(2:end),real(C4(2:end,2:end)));
surf(time,time,real(finalPlot));
shading interp;
colormap('jet');
% caxis([0,1]);
caxis([min(min(real(finalPlot))),max(max(real(finalPlot)))]);

hold on;
contour3(time,time,real(finalPlot),5,'k','LineWidth',2);


if removePlotDetails
title_str = [constructName ' 3-point TCF: \tau_2 = ' tau2usec ' N=' num2str(numel(tcfFiles))];
title(title_str,'FontSize',12);
xlabel('\tau_1 (sec)','fontsize',12);
ylabel('\tau_3 (sec)','fontsize',12);
zlabel('C^{(3)}(\tau_1,\tau_2,\tau_3)','fontsize',12);
set(gca,'yscale','log');
set(gca,'xscale','log');
% set(gca,'zscale','log');
set(gca,'FontSize',12);
grid off
% axis image
% axis tight;
xlim([-inf,time(end)]);
ylim([-inf,time(end)]);
zlim([min(min(real(finalPlot))),max(max(real(finalPlot)))]);
% view(75,25);
view([2,-1,1]); 

% set(gca,'YTickLabel',[]);
% set(gca,'XTickLabel',[]);
% set(gca,'ZTickLabel',[]);
% set(gca,'XMinorTick','off','YMinorTick','off')
else
title_str = [constructName ' 3-point TCF: \tau_2 = ' tau2usec ' N=' num2str(numel(tcfFiles))];
title(title_str,'FontSize',12);
xlabel('\tau_1 (sec)','fontsize',12);
ylabel('\tau_3 (sec)','fontsize',12);
zlabel('C^{(3)}(\tau_1,\tau_2,\tau_3)','fontsize',12);
set(gca,'yscale','log');
set(gca,'xscale','log');
% set(gca,'zscale','log');
set(gca,'FontSize',12);
grid on
% axis image
% axis tight;
xlim([-inf,time(end)]);
ylim([-inf,time(end)]);
zlim([min(min(real(finalPlot))),max(max(real(finalPlot)))]);
% view(75,25);
view([2,-1,1]);
end
tauArrayUsecTotal=firstTauArray; 

if saveMode
    outputFolder= [constructName '_C3wAVG_' quantityName '_' num2str(binSizeList(1)) 'usec']; 
    
    if exist(outputFolder,'dir') ~= 7        
        mkdir(outputFolder);
    end
    
    save([outputFolder  filesep() 'C4_tau2=' tau2usec],'C4finalArray', 'time', 'C4finalArrayNorm', 'tauArrayUsecTotal', 'tau2usec'); 
   
end
        
