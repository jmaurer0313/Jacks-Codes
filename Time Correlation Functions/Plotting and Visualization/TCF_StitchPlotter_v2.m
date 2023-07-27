% 02/07/2022 New version of Stitch Plotter - should take in a set of resolution
% values, for a specified TCF, then stitch together the TCF at each
% specified resolution. 

% ****** SET THE NAME OF THE TCF TO BE STICHED ******
ConstructName='(+15)Dimer';
% TCFname= 'P1AmpAvg'; 
% TCFname= 'P1PhaseFactAvg'; 
TCFname= 'Rate';
TCFfolder='LD TCFs_v3';
Normalize=1;
saveMode=0;
saveNewStitch=0;
overlayWithNext=1;

%choose to smooth curve, if so, select to use spline or not. Otherwise turn
%smooth off and use a reduced density. If none on, standard plot. 
smoothCurve=0; 
plotSpline=0;
reduceDensity=0;

pointsToSmooth=5;
simpleAnchor=1;

% set of resolutions to be stitched
resSet= [250];

% ----------SECTION TO SET STITCH OPTIONS------------
% this sets the number of points which will be used to form the simpleAnchor
% scaling between curve segments. too many points will ensure noisy
% sections of the preceeding curve are used. Too few will be insufficent
% sampling to form a reasonable anchor. something in the 3-5 range to
% start. 
numPointsForStitch=4;
indexToJoin=1;
useSpline=0;
finalScaleFactorAjd=1.0;
ScaleFactorAjd2ndSeg=1.0;

% set the following options true to perform a qudratic and cubic fit to the
% data segments where they overlap, take the overage of the best fit
% options to form a scaling value for one curve to the other
fitAnchor=0;


startFolder=pwd; 
stitchTaus=[];
stitchTCF=[]; 

for j=1:numel(resSet)
  
    folderName=(['XYQuadData_' ConstructName num2str(resSet(j)) 'usec']); 
    
    cd(folderName)
    
    TCFlocation=([num2str(resSet(j)) 'usec' ' ' 'DenseTau BF 2Pt_TCFs_' ConstructName filesep() TCFfolder ...
        filesep() TCFname]); 
    
    cd(TCFlocation); 
    
%     now load the first file for tauArraySec to compare against the next
%     set of Taus
    TCFfiles=dir('*2PtTCF*'); 
    load(TCFfiles(1).name); 
    
%     if this is the beginning/first resolution, simply hold the tauArray
%     as a temp vairable and add together the files in the wAVG fashion.
%     Else, find an overlap region between the lower and upper resolutions
%     to stitch across, with a few methods.
    
     if j==1
        pairsTotalTcf= zeros(1,numel(tauArray)); 
        pairsTotal= zeros(1,numel(tauArray));
        tauArrayTemp=tauArraySec; 
        
        for n=1:numel(TCFfiles)           
           load(TCFfiles(n).name); 
           pairsTcf= tcfArr.*tcfNpairs;  
           pairsTotal= pairsTotal + tcfNpairs; 
           pairsTotalTcf= pairsTotalTcf + pairsTcf; 
        end
        if tauArraySec(1)==0
          pairsWeightTcf= pairsTotalTcf(2:end)./pairsTotal(2:end);
        else
          pairsWeightTcf= pairsTotalTcf./pairsTotal;
        end
           cd(startFolder);
         
     elseif j>1
%          first constrain the taus to be shared for stitching, hold the
%          previous segment as a temp, attach it to the next segment
%            load(TCFfiles(1).name)
            if j==2
           previousTaus=tauArrayTemp((tauArrayTemp>0)); 
            elseif j>2
%           dont overwrite previousTaus in this case
            end
     
           
     if simpleAnchor
%            first method to establish anchoring of one curve to another -
%            sample the preceeding tauArray for points where it is equal to
%            the succeeding tauArray. This will be integer mulltiples of
%            the current res in resSet
            
%             Preceeding Tau Indices         
        remainderArr=mod(previousTaus,resSet(j)/1e6);
        tausToStitch=previousTaus(remainderArr==0);
        indicesToStitchPrev=(find(remainderArr==0));
        
%         current tau indices
%%          
       CurTaus=tauArraySec((tauArraySec>0));
        indicesToStitchCur=[];
        indicesToKeepPrev=[]; 
        for m=1:numPointsForStitch
           indicesToStitchCur=[indicesToStitchCur , find(CurTaus==tausToStitch(m))]; 
            if ~isempty(find(CurTaus==tausToStitch(m)))
              indicesToKeepPrev=[indicesToKeepPrev,indicesToStitchPrev(m)]; 
            end
        end
        %%
        
%         now preceeding and suceeding indices have been establish for
%         which the two curves share an exact tau value. The index used to
%         conjoin the two curves is the earliest available for now. May be
%         useful to change it and delay the cross over point between
%         resolutions should a critical time scale fall at the specified cross over point 
        if j==2
        previousCurve=pairsWeightTcf; 
%         stitchTaus = [stitchTaus,previousTaus(1:indicesToKeepPrev(1)-1)];
        
        pairsTotalTcf= zeros(1,numel(tauArray)); 
        pairsTotal= zeros(1,numel(tauArray));
        tauArrayTemp=tauArraySec; 
        
        for n=1:numel(TCFfiles)           
           load(TCFfiles(n).name); 
           pairsTcf= tcfArr.*tcfNpairs;  
           pairsTotal= pairsTotal + tcfNpairs; 
           pairsTotalTcf= pairsTotalTcf + pairsTcf; 
        end
        
        if tauArraySec(1)==0
          pairsWeightTcf= pairsTotalTcf(2:end)./pairsTotal(2:end);
        else
          pairsWeightTcf= pairsTotalTcf./pairsTotal;
        end
        
%         now to compare the 2 curves at the indices which they match
            scaleFactor=0;
          for z=1:numel(indicesToKeepPrev)
%              take the points in the suceeding res and divide them by the value at the previous
            scaleFactor= scaleFactor + pairsWeightTcf(indicesToStitchCur(z))/previousCurve(indicesToKeepPrev(z));
          end
%           average of the point wise comparisons
          scaleFactor=scaleFactor/(numel(indicesToKeepPrev))*finalScaleFactorAjd;
          oldTaus=previousTaus;
          oldCurve=previousCurve;
%           now the previous taus and TCF curve are the stiched taus and
%           stithced TCF awaiting a final segment (or not) 
          previousTaus = [previousTaus(1:indicesToKeepPrev(indexToJoin)-1) ; CurTaus(indicesToStitchCur(indexToJoin):end)]; 
          previousCurve = [previousCurve(1:indicesToKeepPrev(indexToJoin)-1) , pairsWeightTcf(indicesToStitchCur(indexToJoin):end)./scaleFactor]; 
           
          splineCurve = spline(previousTaus',previousCurve,oldTaus'); 
          splineFactor = nanmean(oldCurve(indicesToKeepPrev(indexToJoin)+10:end-30)./splineCurve(indicesToKeepPrev(indexToJoin)+10:end-30));
          splineScaledCurve =[previousCurve(1:indicesToKeepPrev(indexToJoin)-1) , previousCurve(indicesToKeepPrev(indexToJoin):end)*splineFactor];
          
          if useSpline
          spline2stitchAvg = (splineScaledCurve(indicesToKeepPrev(indexToJoin))+previousCurve(indicesToKeepPrev(indexToJoin)))/2;
          finalFactor = previousCurve(indicesToKeepPrev(indexToJoin))/spline2stitchAvg; 
          previousCurve = [previousCurve(1:indicesToKeepPrev(indexToJoin)-1), previousCurve(indicesToKeepPrev(indexToJoin):end)/finalFactor];
          end
          cd(startFolder);
           
        
        
        
        elseif j>2
            
%         previousCurve=pairsWeightTcf; 
%         stitchTaus = [stitchTaus,previousTaus(1:indicesToKeepPrev(1)-1)];

%             Preceeding Tau Indices         
        remainderArr=mod(previousTaus,resSet(j)/1e6);
        tausToStitch=previousTaus(remainderArr==0);
        indicesToStitchPrev=(find(remainderArr==0));
        
%         current tau indices
%%          
       CurTaus=tauArraySec((tauArraySec>0));
        indicesToStitchCur=[];
        indicesToKeepPrev=[]; 
        for m=1:numPointsForStitch
           indicesToStitchCur=[indicesToStitchCur , find(CurTaus==tausToStitch(m))]; 
            if ~isempty(find(CurTaus==tausToStitch(m)))
              indicesToKeepPrev=[indicesToKeepPrev,indicesToStitchPrev(m)]; 
            end
        end
        
        pairsTotalTcf= zeros(1,numel(tauArray)); 
        pairsTotal= zeros(1,numel(tauArray));
        tauArrayTemp=tauArraySec; 
        
        for n=1:numel(TCFfiles)           
           load(TCFfiles(n).name); 
           pairsTcf= tcfArr.*tcfNpairs;  
           pairsTotal= pairsTotal + tcfNpairs; 
           pairsTotalTcf= pairsTotalTcf + pairsTcf; 
        end
        
        if tauArraySec(1)==0
          pairsWeightTcf= pairsTotalTcf(2:end)./pairsTotal(2:end);
        else
          pairsWeightTcf= pairsTotalTcf./pairsTotal;
        end
        
%         now to compare the 2 curves at the indices which they match
            scaleFactor=0;
          for z=1:numel(indicesToKeepPrev)
%              take the points in the suceeding res and divide them by the value at the previous
            scaleFactor= scaleFactor + pairsWeightTcf(indicesToStitchCur(z))/previousCurve(indicesToKeepPrev(z));
          end
%           average of the point wise comparisons
          scaleFactor=scaleFactor/(numel(indicesToKeepPrev))*ScaleFactorAjd2ndSeg;
          oldTaus=previousTaus;
          oldCurve=previousCurve;
%           now the previous taus and TCF curve are the stiched taus and
%           stithced TCF awaiting a final segment (or not) 

          previousTaus = [previousTaus(1:indicesToKeepPrev(indexToJoin)-1) ; CurTaus(indicesToStitchCur(indexToJoin):end)]; 
          previousCurve = [previousCurve(1:indicesToKeepPrev(indexToJoin)-1) , pairsWeightTcf(indicesToStitchCur(indexToJoin):end)./scaleFactor];
%           previousTaus=[previousTaus(1:indicesToKeepPrev(1)-1) ; CurTaus]; 
%           previousCurve=[previousCurve(1:indicesToKeepPrev(1)-1) , pairsWeightTcf./scaleFactor]; 
                    
          splineCurve = spline(previousTaus',previousCurve,oldTaus'); 
          splineFactor = nanmean(oldCurve(indicesToKeepPrev(indexToJoin)+10:end-30)./splineCurve(indicesToKeepPrev(indexToJoin)+10:end-30));
          splineScaledCurve =[previousCurve(1:indicesToKeepPrev(indexToJoin)-1) , previousCurve(indicesToKeepPrev(indexToJoin):end)*splineFactor];
          
          if useSpline
          spline2stitchAvg = (splineScaledCurve(indicesToKeepPrev(indexToJoin))+previousCurve(indicesToKeepPrev(indexToJoin)))/2;
          finalFactor = previousCurve(indicesToKeepPrev(indexToJoin))/spline2stitchAvg; 
          previousCurve = [previousCurve(1:indicesToKeepPrev(indexToJoin)-1), previousCurve(indicesToKeepPrev(indexToJoin):end)/finalFactor];
          end
          cd(startFolder);  
            
            
            
        end
        
        
     end
     
     end

        
 %%
end
 if numel(resSet)>1
stitchTCF=previousCurve;
stitchTaus=previousTaus;
 else
stitchTCF=pairsWeightTcf; 
stitchTaus=tauArraySec(tauArraySec>0); 
 end

if Normalize
    stitchTCF=stitchTCF./(stitchTCF(1)); 
end

if plotSpline
    splineTaus=[];
   for i=1:numel(stitchTaus)
       if i==numel(stitchTaus)
           splineTaus=[splineTaus, stitchTaus(i)];
       else
       splineTaus=[splineTaus, stitchTaus(i), (stitchTaus(i)+stitchTaus(i+1))/2 ]; 
       end
   end
end
% splineTauAxis = linspace(stitchTaus(1),stitchTaus(end),2000);
if reduceDensity
    reducedTaus=stitchTaus(1:2:end);
    reducedTCF=stitchTCF(1:2:end); 
end

figure(1)
if smoothCurve
    if plotSpline
        splineTCF=spline(stitchTaus,stitchTCF,splineTaus);
        plot(splineTaus,splineTCF,'LineWidth',2);
    else
        filtCurve=medfilt1(abs(real(stitchTCF)),3);
        plot(stitchTaus',[abs(real(stitchTCF(1))),filtCurve(2:end)],'LineWidth',2); 
    end
elseif reduceDensity
        plot(reducedTaus,reducedTCF,'LineWidth',2);
else
  plot(stitchTaus',(real(stitchTCF)),'LineWidth',2);   
end
set(gca,'xscale','log');
set(gca,'yscale','log');
xlim([stitchTaus(1), 3]); 
xlabel('Tau (sec)');
ylabel('C2 Mag (arb)');
legend('stitchedPlot','250usec curve');
title(['Overlay of stitch and original data - ' ConstructName]);
if overlayWithNext
    hold on;
end

totalMatches=length(TCFfiles);

if saveMode
  
    if saveNewStitch   
      folderName=([ConstructName num2str(resSet(1)) 'usec' TCFname '_StitchWAVG_v2']); 
    else
      folderName=([ConstructName num2str(resSet(1)) 'usec' TCFname '_StitchWAVG']); 
    end
     
    if exist(folderName, 'dir')~= 7
    mkdir(folderName);
    end
%     
%     if ~Normalize
%     if exist(folderNameRaw, 'dir')~= 7
%     mkdir(folderNameRaw);
%     end
%     end
   
    if saveNewStitch
    fOutName=([ folderName filesep() ConstructName TCFname ' StitchTCF_' '_wAVG']); 
    save(fOutName, 'stitchTCF', 'stitchTaus','totalMatches','resSet','finalScaleFactorAjd','indexToJoin');    
    else
    fOutName=([ folderName filesep() ConstructName TCFname ' StitchTCF_' '_wAVG']); 
    save(fOutName, 'stitchTCF', 'stitchTaus','totalMatches','resSet');    
    end
   
end 
        
    

   
    
    
    
