%AUTHOR: Jack Maurer 
% GOAL: Plot all the TCFs within a Particular Folder and save the Weighted
% Average plot to that same folder.


traceFolder= pwd;

folderSize= numel(traceFolder);
targetFolder='C:\Users\Ryzen 5\Dropbox\chosenAPD2mat_output\E1E2\gp32_0p0uM\FPGA\100mM_NaCl\SortedSubset';

% ----MUST SET THESE VALUES FOR EACH FOLDER UNTIL UPDATED---
ConstructName= 'E3E4';
binSize= 1000;
percentToKeep= 0.99999;
testSet=0;
tcfQuant= ('P1AmpAvg');
saveOpt=0; 
moveOpt=0; 
TCF2D=0;

% is this a cross corr TCF?
crossCorr=0;

% **** IF OPERATING ON CHOSEN TRACES, SET ITERATIONS MANUALLY*******
% ******************************************************************

maxBinNum= (30000000/binSize);

% % **** NOTE: FOR Rtotal NEED TO REMOVE *usec WILDCARD (since Rtot TCFs are ScanID.mat) ****
curWindowFiles= dir('*.mat');
wAVGFiles= dir('*wAVG.mat');
load(curWindowFiles(5).name);
numMolecules= numel(curWindowFiles)-1; 


% LengfigTCF=figure; 
PairsfigTCF=figure; 
totalLength=0;
% ----NOTE-----:The number of zeros here is set by the numbef of Taus being sampled for
% the BF TCF. Should be reset if the number of Taus sampled is ever
% changed (in this case nTau= 160 currently). 
pairsTotalTcf= zeros(1,numel(tauArray)); 
testTotalTCF= zeros(1,numel(tauArray)); 
lengthTotalTcf= zeros(1,numel(tauArray));
pairsTotal= zeros(1,numel(tauArray));
timeArr=[0:(binSize/1000000):(30-1*(binSize/1000000))];

% FIX 3/27/2019: rather than check the number of files in upper level
% directories, a new check is now performed to see if the wAVG file already
% exists. If so, the number of Plots is 1 less than the total .mat files. else it
% is equal. 
if numel(wAVGFiles)~=0
    numPlots=numel(curWindowFiles)-1;
else
    numPlots=numel(curWindowFiles);
end
% homeDir=pwd; 
% cd ..
% cd ..
% cd ..
% AbsWindowFiles= dir('*.mat');
% numPlots=numel(AbsWindowFiles);
% 
% cd(homeDir);

fOutName=([ ConstructName '_' tcfQuant  'TCF_' num2str(binSize) '_usec_wAVG']); 

% Create a list of the indivudula scanIDs and the variance of each TCF
% stored as varOfArray. Sort the list (descending) by var and then make a
% new list of the top N traces with the highest variance. Then construct
% the wAVG from this list. 
varList=[]; 

for j=1:numel(curWindowFiles)-1
    
    fileName = curWindowFiles(j).name;
    scanID = fileName([1:15]);
    load(fileName);
    firstTau = real(tcfArr(3));
    secondTau = real(tcfArr(4));
    thirdTau = real(tcfArr(6));
    fourthTau = real(tcfArr(8));
    if crossCorr
    newEntry = [varOfArrayF,firstTau,secondTau,thirdTau,fourthTau,j];
    else
    newEntry = [varOfArray,firstTau,secondTau,thirdTau,fourthTau,j]; 
    end
    varList = [varList;newEntry];
    
    if testSet==1
        firstTau = real(tcfAmpTest(2));
        secondTau = real(tcfAmpTest(4));
        thirdTau = real(tcfAmpTest(6));
        fourthTau = real(tcfAmpTest(8));
        newEntry = [varOfArray,firstTau,secondTau,thirdTau,fourthTau,j];
        varList = [varList;newEntry];
    end 
    
    
end 

% CURRENTLY SORTS BASED ON THE VALUES OF VARIOUS TAUS AND NOT THE VARIANCE,
% WHICH PRODUCES A FAR SMOOTHER AND BETTER DEFINED wAVG CURVE
varListSortedFirst = sortrows(varList,[2 1],{'descend' 'descend'});
varListSortedSecond = sortrows(varList,[3 1],{'descend' 'descend'});
varListSortedThird = sortrows(varList,[4 1],{'descend' 'descend'});
varListSortedFourth = sortrows(varList,[5 1],{'descend' 'descend'});
varListSortedFifth = sortrows(varList,[1 2],{'descend' 'descend'});
numToKeep= floor(percentToKeep*numel(curWindowFiles));
% numToKeep= 187;
upperListFirst=varListSortedFirst((1:numToKeep),6); 
upperListSecond=varListSortedSecond((1:numToKeep),6);
upperListThird=varListSortedThird((1:numToKeep),6);
upperListFourth=varListSortedFourth((1:numToKeep),6);
upperListFifth=varListSortedFifth((1:numToKeep),6);


totalInList=0; 
    


scanIDList=[];
% IF OPERATING ON CHOSEN TRACES, SET ITERATIONS MANUALLY
for k=1:numPlots
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
    if ismember(k,upperListFirst)&& ismember(k,upperListSecond) && ismember(k,upperListThird)&& ismember(k,upperListFourth)
%         ismember(k,upperListFourth)
%         && ismember(k,upperListThird)
        totalInList= totalInList+1; 
    fileName = curWindowFiles(k).name;
    scanID = convertCharsToStrings(fileName([1:15]));
    scanIDList=[scanIDList;scanID];
        
    load(fileName);
    
   if TCF2D
    tauArraySec=tauArrayUsec/1e6;
   end
    
    
    if moveOpt
        save([targetFolder filesep() fileName], 'tauArraySec','tauArray','tcfArr','tcfNpairs','trajLength','varOfArray'); 
    end 
        %     The total TCFs averaged by length and by pairs
     totalLength=totalLength+trajLength;     
     lengthTcf= tcfArr*trajLength;      
     lengthTotalTcf=lengthTotalTcf + lengthTcf; 
     
     pairsTcf= tcfArr.*tcfNpairs;  
     pairsTotal= pairsTotal + tcfNpairs; 
     pairsTotalTcf= pairsTotalTcf + pairsTcf; 
     
     if testSet==1
        testTCF = tcfAmpTest;
        testTotalTCF = testTotalTCF + tcfAmpTest;
     end
        
    else
        continue;
    
    end
end 

testTotalTCF=testTotalTCF/totalInList; 
testTotalTCF=testTotalTCF./testTotalTCF(3);
    
LengthWeightTcf= lengthTotalTcf/totalLength; 
pairsWeightTcf= pairsTotalTcf./pairsTotal; 

maxLength= maxk(LengthWeightTcf, 3);
maxPairs= maxk(pairsWeightTcf, 3);

% As of 4/16/19, the normalization max is now the first entry. Since the Tau=0
% is no longer being included in the TCF.
% As of 06/18/2019 the normalization factor is back to the Tau=1st bin
% since Tau=0 is now being calculated for variance comparison 
normFactorLength= maxLength(1,2); 
normFactorPairs= maxPairs(1,2); 
if binSize==1
    pairsTcfFinal= pairsWeightTcf/pairsWeightTcf(3); 
    
else 
lengthTcfFinal= LengthWeightTcf/normFactorLength;
pairsTcfFinal = pairsWeightTcf/normFactorPairs;
end 

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

tauArraySec=tauArraySec(1:length(pairsTcfFinal)); 

myTitle1= ([ConstructName ' ' tcfQuant ' BF 2pt TCF ' num2str(binSize) 'usec (weighted by Length)']);
myTitle2= ([ConstructName ' ' tcfQuant ' BF 2pt TCF ' num2str(binSize) 'usec wAVG N=' num2str(totalInList)]);
% ------- CODE FOR PLOTTING THE AVG TCF BY NUMBER OF PAIRS -------
    figure(PairsfigTCF);
    plot(tauArraySec' , abs(real(pairsTcfFinal)), 'LineWidth', 1.5);
    hold on
%     plot(tauArraySec' , imag(pairsTcfFinal), 'LineWidth', 1.5);
%     xlim([binSize/1e6 3]); 
    xlim([tauArraySec(2) inf]); 
    ylim([0.01 abs(1.1*real(pairsTcfFinal(2)))]); 
    title(myTitle2,'color', 'red','fontsize', 18);
    xlabel('Tau (sec)','color','blue','fontsize', 14);
    ylabel('TCF(Tau) by Pairs','color','blue','fontsize', 14);
%     text(time(5),tcf(5),scanID);
    set(gca,'xscale','log');
    set(gca,'yscale','log');
    
usres=binSize;  

if testSet==1
figure(2)
plot(tauArraySec', testTotalTCF);
set(gca,'xscale','log');
set(gca,'yscale','log');

end 



if saveOpt==1
savefig(myTitle2);
save(fOutName, 'usres','tauArraySec', 'tauArray', 'pairsTcfFinal', 'pairsTotal','numMolecules','totalInList','pairsTotalTcf','scanIDList','firstTau','secondTau','thirdTau','fourthTau');
end
% clear all
