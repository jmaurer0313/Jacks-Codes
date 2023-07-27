% Goal: create bivariates of the mean sq versus the square of the mean for
% the raw phasefactor to se if there is any clustering of the data in a
% particular location. 

% ******SET PARAMETERS WITH EACH RUN*********
constructName='(+15)Dimer';
% with the "Avg" PF carrying the same string to start, there could be a
% redundancy issue -- Luckily it pulls the first comparison so will still
% take the raw PF
quantity='P1AmpAvg'; 
% quantity='AvgAmpCrossRate';
% quantity='Rate'; 
TCFname='LD TCFs_v3';
res=1000e-6;
saveOpt=0;
% ***********************************************

% *************OTHER ADJUSTMENTS FOR BIVAR SORTING AND DIVISON****
v2TCFs=1; 
resUsec=res*1e6;
TrajFiles= dir('*.mat');
% define a fraction of the interval between the min and max vlaue fo the
% two quantities which will serve a the cutoff for the bivariate 
% horizontal line
fracMaxMeanSq=0.1;
% vertical line
fracMaxSqMean=0.1;
% horz upper bound
horzUp=0.8;
% vertical upper bound
vertUp=0.8; 

% interval boundaries for the <e^2>-<e>^2 array as percents of the total
% range of the distribution. In theory large values are better and would
% display large decays in the TCf, small values are likely noise. 
upBound=1.0;
lowBound=0.3;

% *************************************************

% intialize an array that has two columns (to start), 1. square of the
% mean 2. mean of the square  and second arrya with 3. scanID

dataArray=zeros(numel(TrajFiles),2); 
scanIDarray= [];

for i=1:numel(TrajFiles)
    
    load(TrajFiles(i).name); 
    
    if strcmp('P1PhaseFact',quantity)
        
        if v2TCFs
            meanSq=mean((P1PhaseFact/res).*(conj((P1PhaseFact/res))));
            sqMean= mean((P1PhaseFact/res))*mean(conj((P1PhaseFact/res)));
        else
    
     meanSq=mean(P1PhaseFact.*(conj(P1PhaseFact)));
     sqMean= mean(P1PhaseFact)*mean(conj(P1PhaseFact));
        end
        
    elseif  strcmp('Counts',quantity)
    
     meanSq=mean(counts.*(conj(counts)));
     sqMean= mean(counts)*mean(conj(counts));
     
     elseif  strcmp('Rate',quantity)
    
     meanSq=mean((counts/res).*(conj((counts/res))));
     sqMean= mean((counts/res))*mean(conj((counts/res)));
     
    elseif  strcmp('P1Phase',quantity)
    
     meanSq=mean(angle(P1PhaseFact).*(conj(angle(P1PhaseFact))));
     sqMean= abs(mean(angle(P1PhaseFact))*mean(conj(angle(P1PhaseFact))));
     
     elseif  strcmp('P1AmpAvg',quantity)
    
     meanSq=nanmean(abs(P1PhaseFactN).*(conj(abs(P1PhaseFactN))));
     sqMean= abs(nanmean(abs(P1PhaseFactN))*nanmean(conj(abs(P1PhaseFactN))));
     
      elseif  strcmp('P1Amp',quantity)
    
     meanSq=nanmean(abs(P1PhaseFact).*(conj(abs(P1PhaseFact))));
     sqMean= abs(nanmean(abs(P1PhaseFact))*nanmean(conj(abs(P1PhaseFact))));
     
     elseif  strcmp('P1PhaseFactAvg',quantity)
    
     meanSq=nanmean((P1PhaseFactN).*(conj((P1PhaseFactN))));
     sqMean= abs(nanmean((P1PhaseFactN))*nanmean(conj((P1PhaseFactN))));
     
     elseif  strcmp('AvgAmpCrossRate',quantity)
    
     meanSq=nanmean(abs(P1PhaseFactN).*(conj((counts/res))));
     sqMean= abs(nanmean(abs(P1PhaseFactN))*nanmean(conj((counts/res))));
         
    end 
    
    fileName = TrajFiles(i).name;
    scanID = convertCharsToStrings(fileName([1:15]));
    
    scanIDarray=[scanIDarray ; scanID]; 
    dataArray(i,1)= meanSq;
    dataArray(i,2)=sqMean;    
    
    
end 

figure(5)
    
     scatter(dataArray(:,2), dataArray(:,1));
     title([constructName ' Bivariate - MeanSq vs SqMean ' num2str(resUsec) 'usec ' 'X line=' num2str(fracMaxSqMean*100) '% ' 'Y line=' num2str(fracMaxMeanSq*100) '%' ],'FontSize',14); 
%     bivar=[dataArray(:,2), dataArray(:,1)]; 
%     colormap('hot');
% %     hist3(bivar,'CDataMode','auto','FaceColor','interp','Nbins',[25 25]);
%     edges_meanSq = linspace(min(dataArray(:,1)), max(dataArray(:,1)), 50);
%     edges_sqMean = linspace(min(dataArray(:,2)), max(dataArray(:,2)),50);
%     hist3(bivar,'CDataMode','auto','FaceColor','interp','EdgeColor','interp','edges',{edges_meanSq, edges_sqMean});
% %     title(scanID,'FontSize',7);
%     view(2); 
    xlim([min(dataArray(:,2)), max(dataArray(:,2))]);
    ylim([min(dataArray(:,1)), max(dataArray(:,1))]);
    xlabel(['Square Mean']);
    ylabel(['Mean Square']);
    meanSqInterval=max(dataArray(:,1))-min(dataArray(:,1));
    sqMeanInterval=max(dataArray(:,2))-min(dataArray(:,2));
    
    meanSqCutoff=(fracMaxMeanSq*meanSqInterval)+min(dataArray(:,1));
    meanSqUp=(horzUp*meanSqInterval)+min(dataArray(:,1));
    sqMeanCutoff= (fracMaxSqMean*sqMeanInterval)+min(dataArray(:,2));
    SqMeanUp= (vertUp*sqMeanInterval)+min(dataArray(:,2));
    
    hold on
    line([SqMeanUp,SqMeanUp],[0,1e10], 'LineStyle', '--', 'Color','r');
    line([sqMeanCutoff,sqMeanCutoff],[0,1e10], 'LineStyle', '--', 'Color','r');
    line([0,1e10],[meanSqCutoff,meanSqCutoff], 'LineStyle', '--', 'Color','r');
    line([0,1e10],[meanSqUp,meanSqUp], 'LineStyle', '--', 'Color','r');
    hold off
%   generate the scanID list corresponding to each subdivision, 
% which we'll call Uleft Lleft Uright Lright (where U and L designate upper and lower)  
Uleft=[];
Lleft=[];
Uright=[];
Lright=[];
boxSet=[];

deltaArr=dataArray(:,1)-dataArray(:,2);
ratioArr=dataArray(:,1)./dataArray(:,2);
deltaInterval=max(deltaArr)-min(deltaArr);
ratioInterval=max(ratioArr)-min(ratioArr);
deltaUpper=(upBound*deltaInterval)+min(deltaArr);
deltaLower=(lowBound*deltaInterval)+min(deltaArr);

ratioUpper=(upBound*ratioInterval)+min(ratioArr);
ratioLower=(lowBound*ratioInterval)+min(ratioArr);

deltaTraces=[]; 
ratioTraces=[];


%   dataArray(i,1)= meanSq;
%     dataArray(i,2)=sqMean;  
for j=1:size(dataArray,1)
    
    if dataArray(j,1)>=meanSqCutoff && dataArray(j,2)>=sqMeanCutoff && dataArray(j,1)<=meanSqUp && dataArray(j,2)<=SqMeanUp
        boxSet=[boxSet,scanIDarray(j)];
    end 
    
    if deltaArr(j)>=deltaLower && deltaArr(j)<=deltaUpper
        deltaTraces=[deltaTraces, scanIDarray(j)]; 
    end 
    
    if ratioArr(j)>=ratioLower && ratioArr(j)<=ratioUpper
        ratioTraces=[ratioTraces, scanIDarray(j)]; 
    end
    
    if dataArray(j,1)>=meanSqCutoff && dataArray(j,2)<=sqMeanCutoff
        Uleft=[Uleft,scanIDarray(j)];
        
    elseif dataArray(j,1)<=meanSqCutoff && dataArray(j,2)<=sqMeanCutoff
         Lleft=[Lleft,scanIDarray(j)];
    
    elseif dataArray(j,1)>=meanSqCutoff && dataArray(j,2)>=sqMeanCutoff
         Uright=[Uright,scanIDarray(j)];
         
    elseif dataArray(j,1)<=meanSqCutoff && dataArray(j,2)>=sqMeanCutoff
         Lright=[Lright,scanIDarray(j)];
         
    end 
    
end 

% now with the scanID lists for each subdivisoion, move into the
% PhaseFactor folder and make a tcf for each subdivision 
homeFold = dir(pwd);
startFolder=pwd;
% This will be the list of directories containing (*usec) in the name
dirList=[];
tcfList=[];
LDList=[];
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
    if contains(dirList(j).name, 'DenseTau')
        tcfList = [tcfList, dirList(j)];
    end
end 

% now naviagte to the TCFs foler themselves
cd(tcfList(1).name);

% get the LD TCFs folder
LDfolder=dir(pwd);
% since only the LD folder is present at this level (for now) it will
% always be the third folder in the list with the first two being '.' and
% '..'
% UPDATE: since there is now LD_v2, will need ot find Dir with that name    
    
    for j=1:numel(LDfolder)
    if strcmp(LDfolder(j).name, TCFname)
        LDList = [LDList, LDfolder(j)];
    end
    end
    
    if v2TCFs && numel(LDList)>1
        
    cd(LDList(2).name);
    else
    
    cd(LDList(1).name);
    
    end
    
 

% make a new list of all TCF folders
tcfFolders=dir(pwd);
targetFolders=[];

for j=1:numel(tcfFolders)
    if strcmp(tcfFolders(j).name, quantity)
        targetFolders = [targetFolders, tcfFolders(j)];
    end
end 

% based on alphabetical order, P1PhaseFact will always be the first in the
% list, with  FactMapN and  FactN coming 2nd and 3rd respectively
cd(targetFolders(1).name);

wAVGFiles=dir('*wAVG.mat');

% now inside the P1PhaseFact folder, generate four TCFs based on the
% scanIDs that were found from the bivariate subdivsions 
curWindowFiles= dir('*.mat');
% load in the first file to get the size correct to intialize the resulting
% TCFs
load(curWindowFiles(1).name); 

UleftTCF= zeros(1,numel(tauArray)); 
LleftTCF= zeros(1,numel(tauArray)); 
UrightTCF= zeros(1,numel(tauArray));
LrightTCF= zeros(1,numel(tauArray));
boxSetTCF=zeros(1,numel(tauArray));

deltaTCF=zeros(1,numel(tauArray));
ratioTCF=zeros(1,numel(tauArray));

totalTCF= zeros(1,numel(tauArray));


UleftPairs= zeros(1,numel(tauArray)); 
LleftPairs= zeros(1,numel(tauArray)); 
UrightPairs= zeros(1,numel(tauArray));
LrightPairs= zeros(1,numel(tauArray));
boxSetPairs=zeros(1,numel(tauArray));

deltaPairs=zeros(1,numel(tauArray));
ratioPairs=zeros(1,numel(tauArray));

totalPairs=zeros(1,numel(tauArray));

if numel(wAVGFiles)==0
    loopIndex=numel(curWindowFiles);
else
    loopIndex=numel(curWindowFiles)-1;
end 

for k=1:loopIndex
    
    curID=curWindowFiles(k).name;
    curID=curID(1:15); 
    load(curWindowFiles(k).name); 
    
    totalTCF=totalTCF+ (tcfArr.*tcfNpairs);
    totalPairs = totalPairs + tcfNpairs;
    
    if ismember(curID, deltaTraces)
        deltaTCF= deltaTCF + (tcfArr.*tcfNpairs);
        deltaPairs= deltaPairs + tcfNpairs; 
    end 
    
      if ismember(curID, ratioTraces)
        ratioTCF= ratioTCF + (tcfArr.*tcfNpairs);
        ratioPairs= ratioPairs + tcfNpairs; 
    end 
    
     if ismember(curID, boxSet)
        boxSetTCF= boxSetTCF + (tcfArr.*tcfNpairs);
        boxSetPairs= boxSetPairs + tcfNpairs; 
    end 
     
    if ismember(curID, Uleft)
        UleftTCF = UleftTCF + (tcfArr.*tcfNpairs); 
        UleftPairs = UleftPairs + tcfNpairs;
        
    elseif ismember(curID, Lleft)
%         load(curWindowFiles(k).name);
        LleftTCF = LleftTCF + (tcfArr.*tcfNpairs); 
        LleftPairs = LleftPairs + tcfNpairs;
        
     elseif ismember(curID, Uright)
%          load(curWindowFiles(k).name);
        UrightTCF = UrightTCF + (tcfArr.*tcfNpairs); 
        UrightPairs = UrightPairs + tcfNpairs;
        
     elseif ismember(curID, Lright)
%          load(curWindowFiles(k).name);
        LrightTCF = LrightTCF + (tcfArr.*tcfNpairs); 
        LrightPairs = LrightPairs + tcfNpairs;
        
    end    
  
end 

% average and normalize the TCFs
UleftTCF=UleftTCF./UleftPairs;
LleftTCF=LleftTCF./LleftPairs;
UrightTCF=UrightTCF./UrightPairs;
LrightTCF=LrightTCF./LrightPairs;
totalTCF=totalTCF./totalPairs;
deltaTCF=deltaTCF./deltaPairs; 
ratioTCF=ratioTCF./ratioPairs; 
boxSetTCF=boxSetTCF./boxSetPairs; 

UrightTCFraw=UrightTCF;
totalTCFraw=totalTCF; 

UleftTCF=UleftTCF/UleftTCF(2);
LleftTCF=LleftTCF/LleftTCF(2);
UrightTCF=UrightTCF/UrightTCF(2);
LrightTCF=LrightTCF/LrightTCF(2);
totalTCF=totalTCF./totalTCF(2);
deltaTCF=deltaTCF/deltaTCF(2); 
ratioTCF=ratioTCF/ratioTCF(2); 
boxSetTCF=boxSetTCF/boxSetTCF(2); 

numMolecules=numel(Uright); 

if saveOpt
save('Bivar_wAVG','tauArraySec','totalTCF','UrightTCF','totalTCFraw','UrightTCFraw','boxSetTCF','deltaTCF','numMolecules','Uright','boxSet','deltaTraces','scanIDarray'); 
end

figure(1);
myTitle= ([constructName ' Uleft ' quantity ' BF 2pt TCF ' num2str(resUsec) 'usec wAVG N=' num2str(numel(Uleft))]);
    plot(tauArraySec' , real(UleftTCF), 'LineWidth', 1.5);
    hold on
%     plot(tauArraySec' , imag(pairsTcfFinal), 'LineWidth', 1.5);
%     xlim([binSize/1e6 3]); 
    xlim([tauArraySec(2) inf]); 
    ylim([0.01 1.5]); 
    title(myTitle,'color', 'red','fontsize', 18);
    xlabel('Tau (sec)','color','blue','fontsize', 14);
    ylabel('TCF(Tau) by Pairs','color','blue','fontsize', 14);
%     text(time(5),tcf(5),scanID);
    set(gca,'xscale','log');
    set(gca,'yscale','log');
    hold off
    
    figure(2);
myTitle= ([constructName ' Lleft ' quantity ' BF 2pt TCF ' num2str(resUsec) 'usec wAVG N=' num2str(numel(Lleft))]);
    plot(tauArraySec' , real(LleftTCF), 'LineWidth', 1.5);
    hold on
%     plot(tauArraySec' , imag(pairsTcfFinal), 'LineWidth', 1.5);
%     xlim([binSize/1e6 3]); 
    xlim([tauArraySec(2) inf]); 
    ylim([0.01 1.5]); 
    title(myTitle,'color', 'red','fontsize', 18);
    xlabel('Tau (sec)','color','blue','fontsize', 14);
    ylabel('TCF(Tau) by Pairs','color','blue','fontsize', 14);
%     text(time(5),tcf(5),scanID);
    set(gca,'xscale','log');
    set(gca,'yscale','log');
    hold off

    figure(3);
myTitle= ([constructName ' Uright ' quantity ' BF 2pt TCF ' num2str(resUsec) 'usec wAVG N=' num2str(numel(Uright))]);
    plot(tauArraySec' , real(UrightTCF), 'LineWidth', 1.5);
    hold on
%     plot(tauArraySec' , imag(pairsTcfFinal), 'LineWidth', 1.5);
%     xlim([binSize/1e6 3]); 
    xlim([tauArraySec(2) inf]); 
    ylim([0.01 1.5]); 
    title(myTitle,'color', 'red','fontsize', 18);
    xlabel('Tau (sec)','color','blue','fontsize', 14);
    ylabel('TCF(Tau) by Pairs','color','blue','fontsize', 14);
%     text(time(5),tcf(5),scanID);
    set(gca,'xscale','log');
    set(gca,'yscale','log');
    hold off
    
    
    figure(4);
myTitle= ([constructName ' Lright ' quantity ' BF 2pt TCF ' num2str(resUsec) 'usec wAVG N=' num2str(numel(Lright))]);
    plot(tauArraySec' , real(LrightTCF), 'LineWidth', 1.5);
    hold on
%     plot(tauArraySec' , imag(pairsTcfFinal), 'LineWidth', 1.5);
%     xlim([binSize/1e6 3]); 
    xlim([tauArraySec(2) inf]); 
    ylim([0.01 1.5]); 
    title(myTitle,'color', 'red','fontsize', 18);
    xlabel('Tau (sec)','color','blue','fontsize', 14);
    ylabel('TCF(Tau) by Pairs','color','blue','fontsize', 14);
%     text(time(5),tcf(5),scanID);
    set(gca,'xscale','log');
    set(gca,'yscale','log');
    hold off
    
    figure(6);
myTitle= ([constructName ' All ' quantity ' BF 2pt TCF ' num2str(resUsec) 'usec wAVG ']);
    plot(tauArraySec' , real(UleftTCF), 'LineWidth', 1.5);
    hold on
    plot(tauArraySec' , real(LleftTCF), 'LineWidth', 1.5);
    plot(tauArraySec' , real(UrightTCF), 'LineWidth', 1.5);
    plot(tauArraySec' , real(LrightTCF), 'LineWidth', 1.5);
%     plot(tauArraySec' , imag(pairsTcfFinal), 'LineWidth', 1.5);
%     xlim([binSize/1e6 3]); 
    xlim([tauArraySec(2) inf]); 
    ylim([0.01 1.2]); 
    title(myTitle,'color', 'red','fontsize', 18);
    xlabel('Tau (sec)','color','blue','fontsize', 14);
    ylabel('TCF(Tau) by Pairs','color','blue','fontsize', 14);
    legend(['Upper left N=' num2str(numel(Uleft))],['Lower left N=' num2str(numel(Lleft))],['Upper right N=' num2str(numel(Uright))],['Lower right N=' num2str(numel(Lright))]); 
%     text(time(5),tcf(5),scanID);
    set(gca,'xscale','log');
    set(gca,'yscale','log');
      hold off
      
      figure(7);
myTitle= ([constructName ' Total ' quantity ' BF 2pt TCF ' num2str(resUsec) 'usec wAVG N=' num2str(numel(curWindowFiles))]);
    plot(tauArraySec' , real(totalTCF), 'LineWidth', 1.5);
   
%     plot(tauArraySec' , imag(pairsTcfFinal), 'LineWidth', 1.5);
%     xlim([binSize/1e6 3]); 
    xlim([tauArraySec(2) inf]); 
    ylim([0.01 1.2]); 
    title(myTitle,'color', 'red','fontsize', 18);
    xlabel('Tau (sec)','color','blue','fontsize', 14);
    ylabel('TCF(Tau) by Pairs','color','blue','fontsize', 14);
%     legend(['Upper left N=' num2str(numel(Uleft))],['Lower left N=' num2str(numel(Lleft))],['Upper right N=' num2str(numel(Uright))],['Lower right N=' num2str(numel(Lright))]); 
%     text(time(5),tcf(5),scanID);
    set(gca,'xscale','log');
    set(gca,'yscale','log');
     
    
    figure(8);
myTitle= ([constructName ' Delta ' quantity ' BF 2pt TCF ' num2str(resUsec) 'usec wAVG N=' num2str(numel(deltaTraces))]);
    plot(tauArraySec' , real(deltaTCF), 'LineWidth', 1.5);
    hold on
%     plot(tauArraySec' , imag(pairsTcfFinal), 'LineWidth', 1.5);
%     xlim([binSize/1e6 3]); 
    xlim([tauArraySec(2) inf]); 
    ylim([0.01 1.5]); 
    title(myTitle,'color', 'red','fontsize', 18);
    xlabel('Tau (sec)','color','blue','fontsize', 14);
    ylabel('TCF(Tau) by Pairs','color','blue','fontsize', 14);
%     text(time(5),tcf(5),scanID);
    set(gca,'xscale','log');
    set(gca,'yscale','log');
    hold off
    
    figure(11);
myTitle= ([constructName ' Ratio ' quantity ' BF 2pt TCF ' num2str(resUsec) 'usec wAVG N=' num2str(numel(ratioTraces))]);
    plot(tauArraySec' , real(ratioTCF), 'LineWidth', 1.5);
    hold on
%     plot(tauArraySec' , imag(pairsTcfFinal), 'LineWidth', 1.5);
%     xlim([binSize/1e6 3]); 
    xlim([tauArraySec(2) inf]); 
    ylim([0.01 1.5]); 
    title(myTitle,'color', 'red','fontsize', 18);
    xlabel('Tau (sec)','color','blue','fontsize', 14);
    ylabel('TCF(Tau) by Pairs','color','blue','fontsize', 14);
%     text(time(5),tcf(5),scanID);
    set(gca,'xscale','log');
    set(gca,'yscale','log');
    hold off
    
    figure(9)
    histogram(deltaArr,30);
    hold on
    line([deltaLower,deltaLower],[0,20],'LineStyle', '--', 'Color','r');
    line([deltaUpper,deltaUpper],[0,20],'LineStyle', '--', 'Color','r');
    title('Histogram - <E^2>-<E>^2');
    xlabel('<E^2>-<E>^2 - delta '); 
    ylabel('Occurences'); 
    hold off
    
     figure(12)
    histogram(ratioArr,30);
    hold on
    line([ratioLower,ratioLower],[0,20],'LineStyle', '--', 'Color','r');
    line([ratioUpper,ratioUpper],[0,20],'LineStyle', '--', 'Color','r');
    title('Histogram - <E^2>/<E>^2');
    xlabel('<E^2>/<E>^2 - ratio '); 
    ylabel('Occurences'); 
    hold off
    
    
     figure(10);
myTitle= ([constructName ' BoxSet ' quantity ' BF 2pt TCF ' num2str(resUsec) 'usec wAVG N=' num2str(numel(boxSet))]);
    plot(tauArraySec' , real(boxSetTCF), 'LineWidth', 1.5);
    hold on
%     plot(tauArraySec' , imag(pairsTcfFinal), 'LineWidth', 1.5);
%     xlim([binSize/1e6 3]); 
    xlim([tauArraySec(2) inf]); 
    ylim([0.01 1.5]); 
    title(myTitle,'color', 'red','fontsize', 18);
    xlabel('Tau (sec)','color','blue','fontsize', 14);
    ylabel('TCF(Tau) by Pairs','color','blue','fontsize', 14);
%     text(time(5),tcf(5),scanID);
    set(gca,'xscale','log');
    set(gca,'yscale','log');
    hold off
    
    
    
    cd(startFolder);