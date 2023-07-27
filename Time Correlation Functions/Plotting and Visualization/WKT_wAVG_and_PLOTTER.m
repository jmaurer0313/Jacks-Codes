% GOAL: Plot all the TCFs within a Particular Folder and save the Weighted
% Average plot to that same folder.


traceFolder= pwd;

folderSize= numel(traceFolder);

% ----MUST SET THESE VALUES FOR EACH FOLDER UNTIL UPDATED---
ConstructName= 'E3E4';
binSize= 10000;
tcfQuant= ('Ytotal'); 

maxBinNum= (30000000/binSize); 

curWindowFiles= dir('*.mat');

myTitle= ([ConstructName ' ' tcfQuant ' Avg 2pt TCF ' num2str(binSize) 'usec']);
figTCF=figure; 
totalLength=0;
totalTcf= zeros(1,maxBinNum);  
timeArr=[0:(binSize/1000000):(30-1*(binSize/1000000))];

for k=1:numel(curWindowFiles)

    numPlots=numel(curWindowFiles);
    
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
    
    fileName = curWindowFiles(k).name;
    scanID = fileName([1:12]);
        
    load(fileName);
    
     totalLength=totalLength+trajLength;
     
     if length(tcf)<maxBinNum
         numZeros= maxBinNum-length(tcf);
         padTcf= [tcf';zeros(numZeros,1)];
     else
         padTcf= tcf'; 
     end
     
     nanCheck= isnan(padTcf);
     
     if sum(nanCheck)>0
         continue 
     else
        
     
     curTcf= (padTcf')*trajLength; 
     totalTcf= totalTcf + curTcf; 
     end
    
    
    
    
end
weightedTcf= totalTcf/totalLength; 
maxList= maxk(weightedTcf, 2);
normFactor= maxList(1,2); 
weightedTcf= weightedTcf/normFactor; 

   figure(figTCF);
    plot(timeArr , weightedTcf);
    xlim([-inf 3]);    
    title(myTitle);
    xlabel('Tau (sec)');
    ylabel('Avg TCF(Tau)');
%     text(time(5),tcf(5),scanID);
    set(gca,'xscale','log');
    

savefig(myTitle);

clear all