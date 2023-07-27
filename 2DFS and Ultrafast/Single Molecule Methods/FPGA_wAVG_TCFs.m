% GOAL: Plot all the TCFs within a Particular Folder and save the Weighted
% Average plot to that same folder.


traceFolder= pwd;

folderSize= numel(traceFolder);

% ----MUST SET THESE VALUES FOR EACH FOLDER UNTIL UPDATED---
ConstructName= 'E3E4';
binSize= 1000;
% percentToKeep= .65;
% tcfQuant= ('PhaseFactTotalN');

% **** IF OPERATING ON CHOSEN TRACES, SET ITERATIONS MANUALLY*******
% ******************************************************************

% maxBinNum= (30000000/binSize);

% % **** NOTE: FOR Rtotal NEED TO REMOVE *usec WILDCARD (since Rtot TCFs are ScanID.mat) ****
curWindowFiles= dir('*.mat');
wAVGFiles= dir('*wAVG.mat');
load(curWindowFiles(1).name);



% LengfigTCF=figure; 
PairsfigTCF=figure; 
totalLength=0;
% ----NOTE-----:The number of zeros here is set by the numbef of Taus being sampled for
% the BF TCF. Should be reset if the number of Taus sampled is ever
% changed (in this case nTau= 160 currently). 

% intialize the empty arrays
AbsPairsTotal= zeros(1,numel(tauArray));
AbsTotalTcf= zeros(1,numel(tauArray));
     
     
AnglePairsTotal= zeros(1,numel(tauArray)); 
AngleTotalTcf= zeros(1,numel(tauArray)); 

PhaseFactPairsTotal= zeros(1,numel(tauArray)); 
PhaseFactTotalTcf= zeros(1,numel(tauArray)); 


% FIX 3/27/2019: rather than check the number of files in upper level
% directories, a new check is now performed to see if the wAVG file already
% exists. If so, the number of Plots is 1 less than the total .mat files. else it
% is equal. 
if numel(wAVGFiles)~=0
    numScans=numel(curWindowFiles)-1;
else
    numScans=numel(curWindowFiles);
end


fOutName=([ ConstructName '_'  'TCFs_' num2str(binSize) '_usec_wAVG']); 

% Create a list of the indivudula scanIDs and the variance of each TCF
% stored as varOfArray. Sort the list (descending) by var and then make a
% new list of the top N traces with the highest variance. Then construct
% the wAVG from this list. 
% varList=[]; 

% for j=1:numel(curWindowFiles)-1
%     
%     fileName = curWindowFiles(j).name;
%     scanID = fileName([1:12]);
%     load(fileName);
%     firstTau = real(tcfArr(3));
%     secondTau = real(tcfArr(5));
%     thirdTau = real(tcfArr(7));
%     fourthTau = real(tcfArr(10));
%     newEntry = [varOfArray,firstTau,secondTau,thirdTau,fourthTau,j];
%     varList = [varList;newEntry];
%     
%     
% end 

% CURRENTLY SORTS BASED ON THE VALUES OF VARIOUS TAUS AND NOT THE VARIANCE,
% WHICH PRODUCES A FAR SMOOTHER AND BETTER DEFINED wAVG CURVE
% varListSortedFirst = sortrows(varList,[2 1],{'descend' 'descend'});
% varListSortedSecond = sortrows(varList,[3 1],{'descend' 'descend'});
% varListSortedThird = sortrows(varList,[4 1],{'descend' 'descend'});
% varListSortedFourth = sortrows(varList,[5 1],{'descend' 'descend'});
% varListSortedFifth = sortrows(varList,[1 2],{'descend' 'descend'});
% numToKeep= floor(percentToKeep*numel(curWindowFiles));
% upperListFirst=varListSortedFirst((1:numToKeep),6); 
% upperListSecond=varListSortedSecond((1:numToKeep),6);
% upperListThird=varListSortedThird((1:numToKeep),6);
% upperListFourth=varListSortedFourth((1:numToKeep),6);
% upperListFifth=varListSortedFifth((1:numToKeep),6);
% 
% 
% totalInList=0; 
%     
% 
% 
% scanIDList=[];
% IF OPERATING ON CHOSEN TRACES, SET ITERATIONS MANUALLY
for k= 1:numScans
%     1:numScans
%     numPlots
%     numel(AbsWindowFiles)

       
%     given the number of plots that need to be made, designate the number
%     of columns and rows for subplot by finding the truncated root of the
%     number of plots and then add one to the number of rows so the remainder 
%     spills over onto the next row, so that the grid created will
%     be at least large enough to fit all plots created from subplot. 
   

% -----ONLY NEEDED FOR SUB PLOT IMPLEMENTATION--------
%     DimRow= floor(sqrt(numPlots)); 
%     DimCol= floor(sqrt(numPlots)); 
%     
%     if ((DimRow*DimCol)< numPlots)
%         DimRow=DimRow+1;
%         if (DimRow*DimCol)< (numPlots)
%             DimCol=DimCol+1;
%         end
%     end
%     if ismember(k,upperListFirst)&& ismember(k,upperListSecond) && ismember(k,upperListThird)&& ismember(k,upperListFourth)
%         ismember(k,upperListFourth)
%         && ismember(k,upperListThird)
%         totalInList= totalInList+1; 

    fileName = curWindowFiles(k).name;
    scanID = convertCharsToStrings(fileName([1:15]));
%     scanIDList=[scanIDList;scanID];
        
    load(fileName);
    
%     The total TCFs averaged and weighted by pairs

     AbsPairsTcf= AbsTCF.*AbsNpairs;  
     AbsPairsTotal= AbsPairsTotal + AbsNpairs; 
     AbsTotalTcf= AbsTotalTcf + AbsPairsTcf; 
     
     AnglePairsTcf= AngleTCF.*AngleNpairs;  
     AnglePairsTotal= AnglePairsTotal + AngleNpairs; 
     AngleTotalTcf= AngleTotalTcf + AnglePairsTcf;
     
     PhaseFactPairsTcf= PhaseFactTCF.*PhaseFactNpairs;  
     PhaseFactPairsTotal= PhaseFactPairsTotal + PhaseFactNpairs; 
     PhaseFactTotalTcf= PhaseFactTotalTcf + PhaseFactPairsTcf;
        
    
    
    
end 
    
AbsWeightTcf= AbsTotalTcf./AbsPairsTotal; 
AngleWeightTcf= AngleTotalTcf./AnglePairsTotal; 
PhaseFactWeightTcf= PhaseFactTotalTcf./PhaseFactPairsTotal; 

AbsMaxPairs= maxk(AbsWeightTcf, 3);
AngleMaxPairs= maxk(AngleWeightTcf, 3);
PhaseFactMaxPairs= maxk(PhaseFactWeightTcf, 3);

% As of 4/16/19, the normalization max is now the first entry. Since the Tau=0
% is no longer being included in the TCF.
% As of 06/18/2019 the normalization factor is back to the Tau=1st bin
% since Tau=0 is now being calculated for variance comparison 
% normFactorPairs= maxPairs(1,1); 

AbsNormFactor= AbsMaxPairs(1,2);
AngleNormFactor= AngleMaxPairs(1,2);
PhaseFactNormFactor= PhaseFactMaxPairs(1,2);

AbsTcfFinal= AbsWeightTcf/AbsNormFactor; 
AngleTcfFinal= AngleWeightTcf/AngleNormFactor; 
PhaseFactTcfFinal= PhaseFactWeightTcf/PhaseFactNormFactor; 
% if binSize==1
%     pairsTcfFinal= pairsWeightTcf/pairsWeightTcf(3); 
%     
% else 
% pairsTcfFinal = pairsWeightTcf/normFactorPairs;
% end 

numMolecules=k; 

% ------CODE FOR PLOTTING THE AVERAGE TCF WEIGHTED BY LENGTH--------
%    figure(LengfigTCF);
%     plot(tauArraySec' , lengthTcfFinal, 'LineWidth', 2);
%     xlim([-inf 3]); 
%     title(myTitle1,'color', 'red', 'fontsize', 16);
%     xlabel('Tau (sec)','color', 'blue','fontsize', 14);
%     ylabel('TCF(Tau) by Length','color', 'blue','fontsize', 14);
% %     text(time(5),tcf(5),scanID);
%     set(gca,'xscale','log');
%     set(gca,'yscale','log');
%    
%     savefig(myTitle1);

% myTitle1= ([ConstructName ' ' tcfQuant ' BF 2pt TCF ' num2str(binSize) 'usec (weighted by Length)']);
myTitle2= ([ConstructName ' 2pt TCF ' num2str(binSize) 'usec wAVG N=' num2str(numMolecules)]);
% ------- CODE FOR PLOTTING THE AVG TCF BY NUMBER OF PAIRS -------
    figure(PairsfigTCF);
    plot(tauArraySec' , real(PhaseFactTcfFinal), 'LineWidth', 1.5);
    hold on
     plot(tauArraySec' , real(AbsTcfFinal), 'LineWidth', 1.5);
    hold on
     plot(tauArraySec' , real(AngleTcfFinal), 'LineWidth', 1.5);
    hold on
    legend('PhaseFact TCF','Amp TCF','Phase TCF');
%     plot(tauArraySec' , imag(pairsTcfFinal), 'LineWidth', 1.5);
    xlim([0 tauArraySec(end)]); 
    ylim([0.01 1.5]); 
    title(myTitle2,'color', 'red','fontsize', 18);
    xlabel('Tau (sec)','color','blue','fontsize', 14);
    ylabel('TCF(Tau) by Pairs','color','blue','fontsize', 14);
%     text(time(5),tcf(5),scanID);
    set(gca,'xscale','log');
    set(gca,'yscale','log');
    
usres=binSize;  
% savefig(myTitle2);
% save(fOutName, 'usres','tauArraySec', 'tauArray', 'pairsTcfFinal', 'pairsTotal','numMolecules','totalInList','pairsTotalTcf','scanIDList');

% clear all