
function [weightC4Func,weightC2func,c2Scalar,c4Scalar,histScalar] = weightsCalculator_v4_chi2(noiseFloorHist,noiseFloorC2,noiseFloorC4,noisePercentDiffC2,noisePercentDiffC4,C2_exp_x, C2_exp_y, C4_tau1range, C4_tau2eq0_exp, targetHistogram,c2Mag,c4Mag,histMag,segmentsDesired,c4Slider,c4WeightTime,Amp_bins,histUpLim,c2WeightTime,C2LateMag,C3LateMag)
% write a short routine to find the optimal weight functions for the C2 and
% C4 givne some arbitrary data. The contribution from each decade of Tau
% should contribute equally given both the differences in the magnitudes,
% but also the number of points in each decade. Once the optimal
% 1/nthroot() is determined for the C2 and C4 at 1% deviation, determine the
% scalars that would match a 1% deviation from the Hist, C2 and C4.

%V2 version will try to segment the decades in N pieces and then perform a
%best fit of the constructed weight surfaces to a 1/nthroot() function for
%the value of 'n'

noiseC2perc = noisePercentDiffC2; 
noiseC4perc = noisePercentDiffC4;
noiseFloorValC2=noiseFloorC2; 
noiseFloorValC4=noiseFloorC4; 
noiseHist=noiseFloorHist;
% start with just the c2 to see if an optimal 1/nthroot() can be found
c2Data = C2_exp_y;
c2Tau = C2_exp_x;

c4Data = C4_tau2eq0_exp;
c4Tau = C4_tau1range;
% % short long slope will dictate the the ratio between the first decade and
% % the subsquent decade - value here dictates the 1st decade scalar, with
% % the last taking 1x scalar. 
% shortLongSlope=1.0;

lowExp=floor(log10(c2Tau(1)));
highExp=ceil(log10(c2Tau(end))); 
lowLim= 1*10^(lowExp); 
highLim= 1*10^(highExp); 
numDecades=log10(highLim/lowLim); 
idxPairList=[];

    histLB=0.0025;
    histUB=histUpLim;
    [I_LB,val_LB]=find(Amp_bins>histLB);
    [I_UB,val_UB]=find(Amp_bins<histUB);
    startIdx=val_LB(1);
    endIdx=val_UB(end);

%%
decadeSegments=segmentsDesired;

% find the points in the tau array where each decade begins and ends (keep
% the indexes)

for j=1:numDecades
%     curScaledTau= c2Tau*1*10^(abs(lowExp)-j);
    [val,idx]=min(abs((c2Tau*1*10^(abs(lowExp)-j))-1));
    
    if j>1
    idxPairList=[idxPairList ; idxPairList(j-1,2)+1 , idx];   
    else
    idxPairList=[1 , idx];
    end
    
end

idxPairListSegments=[];
% Now that the decades have been found, split them into the number of
% segments desired (maximum of 4 given the first decade has only 4 points)
for j=1:length(idxPairList)
    mySpan = idxPairList(j,2)-idxPairList(j,1);
    idxSeg=ceil(mySpan/decadeSegments); 
    for m=1:decadeSegments
        if m==1
        idxPairListSegments=[idxPairListSegments ; idxPairList(j,1),idxPairList(j,1)+ idxSeg];    
        else
            if (idxSeg + (idxPairListSegments(end,2)+1))<idxPairList(j,2) 
             idxPairListSegments=[idxPairListSegments ; (idxPairListSegments(end,2)+1),(idxSeg + idxPairListSegments(end,2)+1)];
            elseif idxPairListSegments(end,2)<idxPairList(j,2) 
             idxPairListSegments=[idxPairListSegments ; (idxPairListSegments(end,2)+1), idxPairList(j,2)]; 
             break;
            end
       
        end
    end
end

% if idxPairListSegments(end,2)>length(c2Tau)
%     idxPairListSegments(end,2)=length(c2Tau);
% end
    
    


% so now, sum the data (c2Data) at 1% deviation to get the raw (unaveraged)
% chi2 across each of these segments. Then establish the multiplicative
% values needed to balance the segments in overall contribution to the
% chi-2
dataChi2Segs=zeros(1,length(idxPairListSegments)); 
dataChi2SegsScaled=zeros(1,length(idxPairListSegments)); 
perErrRampC2=linspace(0,noiseC2perc,length(idxPairListSegments));

[idxs]=find(c2Tau>c2WeightTime); 
dotDivLateC2=mean(c2Data((idxs(1)-1):idxs(1)));
dotDivC2=c2Data;
dotDivC2(idxs(1):end)=dotDivLateC2;

for i=1:length(idxPairListSegments)
    dataChi2Segs(i)=sum(((c2Data(idxPairListSegments(i,1):idxPairListSegments(i,2)) - ...
        (1+(noiseFloorValC2)+perErrRampC2(i))*c2Data(idxPairListSegments(i,1):idxPairListSegments(i,2))).^2) ...
        ./dotDivC2(idxPairListSegments(i,1):idxPairListSegments(i,2)));    
    
end 

% INTERESTING NOTE: the chi2 max contribution is both dependent on the data
% set (in terms of where the peak contribution at x% deviation occurs) but
% is also not strcitly increasing (i.e. the final decade is not the most
% significantly weighted, it tends to the middle 2-3 decades). 

% POSSIBLE SOLUTION: create a set of scales that causes all N decades to
% become equal in terms of chi2 contribution, then apply a variable linear
% ramp to the reuslting weights within each decade, that allows for a
% slight overweighting of early times compared to late time in Tau

scalesC2=1./(dataChi2Segs./max(dataChi2Segs));
test= dataChi2Segs.*(scalesC2);
% plot(test.*fliplr(linspace(1,shortLongSlope,length(test))))
% hold on
% plot(dataChi2Segs)

% check that if this scaling is applied to the ctual chi2 that the sum
% obeys unity across each decade
% for i=1:length(idxPairList)
%     dataChi2SegsScaled(i)=sum((c2Data(idxPairList(i,1):idxPairList(i,2)) - ...
%         (1.1)*c2Data(idxPairList(i,1):idxPairList(i,2))).^2).*scalesC2(i);     
% end 

% scalesWramp = scalesC2.*fliplr(linspace(1,shortLongSlope,length(test)));
weightC2func = [];
noiseMaskC2=[];

for j=1:length(idxPairListSegments)
    weightC2func = [weightC2func , ones(1,idxPairListSegments(j,2)-idxPairListSegments(j,1)+1)*scalesC2(j)];
    noiseMaskC2 = [noiseMaskC2, ones(1,idxPairListSegments(j,2)-idxPairListSegments(j,1)+1)*(1+(noiseFloorValC2)+perErrRampC2(j))];
end 
%%
% appears the c2 weight func is properly constrcuted, now to prepare the
% equivalent for the C4

%************************************************************************
%*************START THE C4 WEIGHT FUNCTION **************************

% short long slope will dictate the the ratio between the first decade and
% the subsquent decade - value here dictates the 1st decade scalar, with
% the last taking 1x scalar. 
%shortLongSlope=1.2;

lowExp=floor(log10(c4Tau(1)));
highExp=ceil(log10(c4Tau(end))); 
lowLim= 1*10^(lowExp); 
highLim= 1*10^(highExp); 
numDecades=log10(highLim/lowLim); 
idxPairListC4=[];

% find the points in the tau array where each decade begins and ends (keep
% the indexes)

for j=1:numDecades
%     curScaledTau= c2Tau*1*10^(abs(lowExp)-j);
    [val,idx]=min(abs((c4Tau*1*10^(abs(lowExp)-j))-1));
    
    if j>1
    idxPairListC4=[idxPairListC4 ; idxPairListC4(j-1,2)+1 , idx];   
    else
    idxPairListC4=[1 , idx];
    end
    
end

% now segment the C4 tau range into the same number of segments per decade
% as the c2
idxPairListSegmentsC4=[];
% Now that the decades have been found, split them into the number of
% segments desired (maximum of 4 given the first decade has only 4 points)
for j=1:length(idxPairListC4)
    mySpan = idxPairListC4(j,2)-idxPairListC4(j,1);
    idxSeg=ceil(mySpan/decadeSegments); 
    for m=1:decadeSegments
        if m==1
        idxPairListSegmentsC4=[idxPairListSegmentsC4 ; idxPairListC4(j,1),idxPairListC4(j,1)+ idxSeg];    
        else
            if (idxSeg + (idxPairListSegmentsC4(end,2)+1))<idxPairListC4(j,2) 
             idxPairListSegmentsC4=[idxPairListSegmentsC4 ; (idxPairListSegmentsC4(end,2)+1),(idxSeg + idxPairListSegmentsC4(end,2)+1)];
            elseif idxPairListSegmentsC4(end,2)<idxPairListC4(j,2) 
             idxPairListSegmentsC4=[idxPairListSegmentsC4 ; (idxPairListSegmentsC4(end,2)+1), idxPairListC4(j,2)]; 
             break;
            end
       
        end
    end
end

% in order to isolate the relevant sections of the 4point, construct a
% boolean mask for each decade and apply to the data matrix
dataChi2SegsC4=zeros(1,length(idxPairListSegmentsC4)); 
perErrRampC4=linspace(0,noiseC4perc,length(idxPairListSegmentsC4));

baseMask=zeros(length(c4Tau),length(c4Tau)); 

[idxs,bool]=find(c4Tau>c4WeightTime); 
dotDivLate=mean(c4Data(idxs(1):idxs(1),1:idxs(1)));
dotDivSurf=c4Data;
dotDivSurf(idxs(1):end,1:end)=dotDivLate;
dotDivSurf(1:end,idxs(1):end)=dotDivLate;

for i=1:length(idxPairListSegmentsC4)
    curMaskTau=baseMask;
%     curMaskTau3=baseMask;
    
    curMaskTau(idxPairListSegmentsC4(i,1):idxPairListSegmentsC4(i,2),(1:idxPairListSegmentsC4(i,2)))=1;  
    curMaskTau((1:idxPairListSegmentsC4(i,2)),idxPairListSegmentsC4(i,1):idxPairListSegmentsC4(i,2))=1;
    
        dataChi2SegsC4(i) = sum(sum(abs((((c4Data*curMaskTau)-(c4Data*curMaskTau)*(1 + (noiseFloorValC4) + perErrRampC4(i))).^2)./(dotDivSurf*curMaskTau)),'omitnan'),'omitnan');
%     dataChi2SegsC4(i) = sum(sum(((c4Data*curMaskTau)-(c4Data*curMaskTau)*(1 + (noiseFloorValC4) + perErrRampC4(i))).^2));

       
end 

scalesC4 = 1./(dataChi2SegsC4./max(dataChi2SegsC4));
test = dataChi2SegsC4.*(scalesC4);
% scalesC4=scalesC4.*fliplr(linspace(1,shortLongSlope,length(test)));
weightC4Func = zeros(length(c4Tau),length(c4Tau));
noiseMaskC4 = zeros(length(c4Tau),length(c4Tau));


for j=1:length(scalesC4)
   
    weightC4Func(idxPairListSegmentsC4(j,1):idxPairListSegmentsC4(j,2),(1:idxPairListSegmentsC4(j,2)))=scalesC4(j);  
    weightC4Func((1:idxPairListSegmentsC4(j,2)),idxPairListSegmentsC4(j,1):idxPairListSegmentsC4(j,2))=scalesC4(j);
    
    noiseMaskC4(idxPairListSegmentsC4(j,1):idxPairListSegmentsC4(j,2),(1:idxPairListSegmentsC4(j,2)))=(1 + (noiseFloorValC4) + perErrRampC4(j));
    noiseMaskC4((1:idxPairListSegmentsC4(j,2)),idxPairListSegmentsC4(j,1):idxPairListSegmentsC4(j,2))=(1 + (noiseFloorValC4) + perErrRampC4(j));
end

% due to the low density of points in the C2 at the very end, there is an
% uptick in the magnitude of the weight func, when those should be the
% least signifcant points. apply the last "valid" point to all points
% beyond 3 seconds
[val,idxs]=find(c2Tau>c2WeightTime); 
weightC2func(idxs(1):idxs(end))=C2LateMag*mean(weightC2func(idxs(1)-5:idxs(1)));

% set the C4 "steepness" by alterign the maxima in the first grid of points
% (which is always too high a value for a shallow surface)
% [val,idxs]=(max(max(weightC4Func))); 
[list]=find(weightC4Func>=weightC4Func(1,1));
weightC4Func(list)=weightC4Func(list)*c4Slider;
% weightC2func(idxs(1):idxs(end))=mean(weightC2func(idxs(1)-5:idxs(1)));

%try to set the final segment of the c4 surface equal to 1/2 the avergae of
%the preceding decade, start from 500ms and flatten surface past that mark 
[idxs,bool]=find(c4Tau>c4WeightTime); 
lateTimeWeight=C3LateMag*weightC4Func(idxs(1),idxs(1));
weightC4Func(idxs(1):end,1:end)=lateTimeWeight;
weightC4Func(1:end,idxs(1):end)=lateTimeWeight;
% now with the weighting surfaces determined, try to find the scalars that
% make the 3 surfaces match in magnitude 

% first get the overall Chi2 of the C2 with weightFunc applied, C4 with
% weight func applied, and then hist with no weight func applied -all at
% 10% deviation

% c2Chi2 = mean(((c2Data-c2Data.*noiseMaskC2).^2).*weightC2func); 
% c4Chi2 = mean(mean(((c4Data-c4Data.*noiseMaskC4).^2).*weightC4Func));
% histChi2 = mean((targetHistogram-targetHistogram*(1+(noiseHist))).^2); 

% Need to implement a toggle for alterante data sets on histogrtma bounds
% like in chiSqCalc_v4
c2Chi2 = mean((((c2Data-c2Data.*noiseMaskC2).^2)./dotDivC2).*weightC2func); 
c4Chi2 = mean(mean(abs((((c4Data-c4Data.*noiseMaskC4).^2)./dotDivSurf).*weightC4Func)));
histChi2 = mean((( (targetHistogram(startIdx:endIdx)-targetHistogram(startIdx:endIdx)*(1+(noiseHist))).^2)./targetHistogram(startIdx:endIdx))); 

chi2Set=[ c2Chi2 , c4Chi2, histChi2 ];
chi2Scalars = 1./(chi2Set./max(chi2Set));

% c2Scalar = c2Mag*chi2Scalars(1);
% c4Scalar = c4Mag*chi2Scalars(2);
% histScalar = histMag*chi2Scalars(3);

c2Scalar = c2Mag;
c4Scalar = c4Mag;
histScalar = histMag;

% figure()
% plot(c2Tau, weightC2func);
% set(gca,'xscale','log');
% 
% c4Test = 1./(nthroot(c4Tau,2)).*(1./(nthroot(c4Tau,2)))';
% figure()
% surf(c4Tau,c4Tau,weightC4Func./(weightC4Func(1,1)));
% hold on
% surf(c4Tau,c4Tau,c4Test./(c4Test(1,1)));
% view([45 0]); 
% 
% set(gca,'xscale','log');
% set(gca,'yscale','log');


end
% *****************************************************************************
% OLD VERSION OF THE FUNCTION - TALAPAS VERSION OF THE FUNCTION UP TOP^^^^^
% *************************************************************************
% function [weightC4Func,weightC2func,c2Scalar,c4Scalar,histScalar] = weightsCalculator_v4_chi2(noiseFloorHist,noiseFloorC2,noiseFloorC4,noisePercentDiffC2,noisePercentDiffC4,C2_exp_x, C2_exp_y, C4_tau1range, C4_tau2eq0_exp, targetHistogram,c2Mag,c4Mag,histMag,segmentsDesired,c4Slider,c4WeightTime,Amp_bins,histUpLim,c2WeightTime,C2LateMag,C3LateMag)
% % write a short routine to find the optimal weight functions for the C2 and
% % C4 givne some arbitrary data. The contribution from each decade of Tau
% % should contribute equally given both the differences in the magnitudes,
% % but also the number of points in each decade. Once the optimal
% % 1/nthroot() is determined for the C2 and C4 at 1% deviation, determine the
% % scalars that would match a 1% deviation from the Hist, C2 and C4.
% 
% %V2 version will try to segment the decades in N pieces and then perform a
% %best fit of the constructed weight surfaces to a 1/nthroot() function for
% %the value of 'n'
% 
% noiseC2perc = noisePercentDiffC2; 
% noiseC4perc = noisePercentDiffC4;
% noiseFloorValC2=noiseFloorC2; 
% noiseFloorValC4=noiseFloorC4; 
% noiseHist=noiseFloorHist;
% % start with just the c2 to see if an optimal 1/nthroot() can be found
% c2Data = C2_exp_y;
% c2Tau = C2_exp_x;
% 
% c4Data = C4_tau2eq0_exp;
% c4Tau = C4_tau1range;
% % % short long slope will dictate the the ratio between the first decade and
% % % the subsquent decade - value here dictates the 1st decade scalar, with
% % % the last taking 1x scalar. 
% % shortLongSlope=1.0;
% 
% lowExp=floor(log10(c2Tau(1)));
% highExp=ceil(log10(c2Tau(end))); 
% lowLim= 1*10^(lowExp); 
% highLim= 1*10^(highExp); 
% numDecades=log10(highLim/lowLim); 
% idxPairList=[];
% 
%     histLB=0.0025;
%     histUB=histUpLim;
%     [I_LB,val_LB]=find(Amp_bins>histLB);
%     [I_UB,val_UB]=find(Amp_bins<histUB);
%     startIdx=val_LB(1);
%     endIdx=val_UB(end);
% 
% %%
% decadeSegments=segmentsDesired;
% 
% % find the points in the tau array where each decade begins and ends (keep
% % the indexes)
% 
% for j=1:numDecades
% %     curScaledTau= c2Tau*1*10^(abs(lowExp)-j);
%     [val,idx]=min(abs((c2Tau*1*10^(abs(lowExp)-j))-1));
%     
%     if j>1
%     idxPairList=[idxPairList ; idxPairList(j-1,2)+1 , idx];   
%     else
%     idxPairList=[1 , idx];
%     end
%     
% end
% 
% idxPairListSegments=[];
% % Now that the decades have been found, split them into the number of
% % segments desired (maximum of 4 given the first decade has only 4 points)
% for j=1:length(idxPairList)
%     mySpan = idxPairList(j,2)-idxPairList(j,1);
%     idxSeg=ceil(mySpan/decadeSegments); 
%     for m=1:decadeSegments
%         if m==1
%         idxPairListSegments=[idxPairListSegments ; idxPairList(j,1),idxPairList(j,1)+ idxSeg];    
%         else
%             if (idxSeg + (idxPairListSegments(end,2)+1))<idxPairList(j,2) 
%              idxPairListSegments=[idxPairListSegments ; (idxPairListSegments(end,2)+1),(idxSeg + idxPairListSegments(end,2)+1)];
%             elseif idxPairListSegments(end,2)<idxPairList(j,2) 
%              idxPairListSegments=[idxPairListSegments ; (idxPairListSegments(end,2)+1), idxPairList(j,2)]; 
%              break;
%             end
%        
%         end
%     end
% end
% 
% % if idxPairListSegments(end,2)>length(c2Tau)
% %     idxPairListSegments(end,2)=length(c2Tau);
% % end
%     
%     
% 
% 
% % so now, sum the data (c2Data) at 1% deviation to get the raw (unaveraged)
% % chi2 across each of these segments. Then establish the multiplicative
% % values needed to balance the segments in overall contribution to the
% % chi-2
% dataChi2Segs=zeros(1,length(idxPairListSegments)); 
% dataChi2SegsScaled=zeros(1,length(idxPairListSegments)); 
% perErrRampC2=linspace(0,noiseC2perc,length(idxPairListSegments));
% 
% [idxs]=find(c2Tau>c2WeightTime); 
% dotDivLateC2=mean(c2Data((idxs(1)-1):idxs(1)));
% dotDivC2=c2Data;
% dotDivC2(idxs(1):end)=dotDivLateC2;
% 
% for i=1:length(idxPairListSegments)
%     dataChi2Segs(i)=sum(((c2Data(idxPairListSegments(i,1):idxPairListSegments(i,2)) - ...
%         (1+(noiseFloorValC2)+perErrRampC2(i))*c2Data(idxPairListSegments(i,1):idxPairListSegments(i,2))).^2) ...
%         ./dotDivC2(idxPairListSegments(i,1):idxPairListSegments(i,2)));    
%     
% end 
% 
% % INTERESTING NOTE: the chi2 max contribution is both dependent on the data
% % set (in terms of where the peak contribution at x% deviation occurs) but
% % is also not strcitly increasing (i.e. the final decade is not the most
% % significantly weighted, it tends to the middle 2-3 decades). 
% 
% % POSSIBLE SOLUTION: create a set of scales that causes all N decades to
% % become equal in terms of chi2 contribution, then apply a variable linear
% % ramp to the reuslting weights within each decade, that allows for a
% % slight overweighting of early times compared to late time in Tau
% 
% scalesC2=1./(dataChi2Segs./max(dataChi2Segs));
% test= dataChi2Segs.*(scalesC2);
% % plot(test.*fliplr(linspace(1,shortLongSlope,length(test))))
% % hold on
% % plot(dataChi2Segs)
% 
% % check that if this scaling is applied to the ctual chi2 that the sum
% % obeys unity across each decade
% % for i=1:length(idxPairList)
% %     dataChi2SegsScaled(i)=sum((c2Data(idxPairList(i,1):idxPairList(i,2)) - ...
% %         (1.1)*c2Data(idxPairList(i,1):idxPairList(i,2))).^2).*scalesC2(i);     
% % end 
% 
% % scalesWramp = scalesC2.*fliplr(linspace(1,shortLongSlope,length(test)));
% weightC2func = [];
% noiseMaskC2=[];
% 
% for j=1:length(idxPairListSegments)
%     weightC2func = [weightC2func , ones(1,idxPairListSegments(j,2)-idxPairListSegments(j,1)+1)*scalesC2(j)];
%     noiseMaskC2 = [noiseMaskC2, ones(1,idxPairListSegments(j,2)-idxPairListSegments(j,1)+1)*(1+(noiseFloorValC2)+perErrRampC2(j))];
% end 
% %%
% % appears the c2 weight func is properly constrcuted, now to prepare the
% % equivalent for the C4
% 
% %************************************************************************
% %*************START THE C4 WEIGHT FUNCTION **************************
% 
% % short long slope will dictate the the ratio between the first decade and
% % the subsquent decade - value here dictates the 1st decade scalar, with
% % the last taking 1x scalar. 
% %shortLongSlope=1.2;
% 
% lowExp=floor(log10(c4Tau(1)));
% highExp=ceil(log10(c4Tau(end))); 
% lowLim= 1*10^(lowExp); 
% highLim= 1*10^(highExp); 
% numDecades=log10(highLim/lowLim); 
% idxPairListC4=[];
% 
% % find the points in the tau array where each decade begins and ends (keep
% % the indexes)
% 
% for j=1:numDecades
% %     curScaledTau= c2Tau*1*10^(abs(lowExp)-j);
%     [val,idx]=min(abs((c4Tau*1*10^(abs(lowExp)-j))-1));
%     
%     if j>1
%     idxPairListC4=[idxPairListC4 ; idxPairListC4(j-1,2)+1 , idx];   
%     else
%     idxPairListC4=[1 , idx];
%     end
%     
% end
% 
% % now segment the C4 tau range into the same number of segments per decade
% % as the c2
% idxPairListSegmentsC4=[];
% % Now that the decades have been found, split them into the number of
% % segments desired (maximum of 4 given the first decade has only 4 points)
% for j=1:length(idxPairListC4)
%     mySpan = idxPairListC4(j,2)-idxPairListC4(j,1);
%     idxSeg=ceil(mySpan/decadeSegments); 
%     for m=1:decadeSegments
%         if m==1
%         idxPairListSegmentsC4=[idxPairListSegmentsC4 ; idxPairListC4(j,1),idxPairListC4(j,1)+ idxSeg];    
%         else
%             if (idxSeg + (idxPairListSegmentsC4(end,2)+1))<idxPairListC4(j,2) 
%              idxPairListSegmentsC4=[idxPairListSegmentsC4 ; (idxPairListSegmentsC4(end,2)+1),(idxSeg + idxPairListSegmentsC4(end,2)+1)];
%             elseif idxPairListSegmentsC4(end,2)<idxPairListC4(j,2) 
%              idxPairListSegmentsC4=[idxPairListSegmentsC4 ; (idxPairListSegmentsC4(end,2)+1), idxPairListC4(j,2)]; 
%              break;
%             end
%        
%         end
%     end
% end
% 
% % in order to isolate the relevant sections of the 4point, construct a
% % boolean mask for each decade and apply to the data matrix
% dataChi2SegsC4=zeros(1,length(idxPairListSegmentsC4)); 
% perErrRampC4=linspace(0,noiseC4perc,length(idxPairListSegmentsC4));
% 
% baseMask=zeros(length(c4Tau),length(c4Tau)); 
% 
% [idxs,bool]=find(c4Tau>c4WeightTime); 
% dotDivLate=mean(c4Data(idxs(1):idxs(1),1:idxs(1)));
% dotDivSurf=c4Data;
% dotDivSurf(idxs(1):end,1:end)=dotDivLate;
% dotDivSurf(1:end,idxs(1):end)=dotDivLate;
% 
% for i=1:length(idxPairListSegmentsC4)
%     curMaskTau=baseMask;
% %     curMaskTau3=baseMask;
%     
%     curMaskTau(idxPairListSegmentsC4(i,1):idxPairListSegmentsC4(i,2),(1:idxPairListSegmentsC4(i,2)))=1;  
%     curMaskTau((1:idxPairListSegmentsC4(i,2)),idxPairListSegmentsC4(i,1):idxPairListSegmentsC4(i,2))=1;
%     
%     dataChi2SegsC4(i) = sum(sum(abs((((c4Data*curMaskTau)-(c4Data*curMaskTau)*(1 + (noiseFloorValC4) + perErrRampC4(i))).^2)./(dotDivSurf*curMaskTau)),'omitnan'),'omitnan');
% %     dataChi2SegsC4(i) = sum(sum(((c4Data*curMaskTau)-(c4Data*curMaskTau)*(1 + (noiseFloorValC4) + perErrRampC4(i))).^2));
% 
%        
% end 
% 
% scalesC4 = 1./(dataChi2SegsC4./max(dataChi2SegsC4));
% test = dataChi2SegsC4.*(scalesC4);
% % scalesC4=scalesC4.*fliplr(linspace(1,shortLongSlope,length(test)));
% weightC4Func = zeros(length(c4Tau),length(c4Tau));
% noiseMaskC4 = zeros(length(c4Tau),length(c4Tau));
% 
% 
% for j=1:length(scalesC4)
%    
%     weightC4Func(idxPairListSegmentsC4(j,1):idxPairListSegmentsC4(j,2),(1:idxPairListSegmentsC4(j,2)))=scalesC4(j);  
%     weightC4Func((1:idxPairListSegmentsC4(j,2)),idxPairListSegmentsC4(j,1):idxPairListSegmentsC4(j,2))=scalesC4(j);
%     
%     noiseMaskC4(idxPairListSegmentsC4(j,1):idxPairListSegmentsC4(j,2),(1:idxPairListSegmentsC4(j,2)))=(1 + (noiseFloorValC4) + perErrRampC4(j));
%     noiseMaskC4((1:idxPairListSegmentsC4(j,2)),idxPairListSegmentsC4(j,1):idxPairListSegmentsC4(j,2))=(1 + (noiseFloorValC4) + perErrRampC4(j));
% end
% 
% % due to the low density of points in the C2 at the very end, there is an
% % uptick in the magnitude of the weight func, when those should be the
% % least signifcant points. apply the last "valid" point to all points
% % beyond 3 seconds
% [val,idxs]=find(c2Tau>c2WeightTime); 
% weightC2func(idxs(1):idxs(end))=C2LateMag*mean(weightC2func(idxs(1)-5:idxs(1)));
% 
% % set the C4 "steepness" by alterign the maxima in the first grid of points
% % (which is always too high a value for a shallow surface)
% % [val,idxs]=(max(max(weightC4Func))); 
% [list]=find(weightC4Func>=weightC4Func(1,1));
% weightC4Func(list)=weightC4Func(list)*c4Slider;
% % weightC2func(idxs(1):idxs(end))=mean(weightC2func(idxs(1)-5:idxs(1)));
% 
% %try to set the final segment of the c4 surface equal to 1/2 the avergae of
% %the preceding decade, start from 500ms and flatten surface past that mark 
% [idxs,bool]=find(c4Tau>c4WeightTime); 
% lateTimeWeight=C3LateMag*weightC4Func(idxs(1),idxs(1));
% weightC4Func(idxs(1):end,1:end)=lateTimeWeight;
% weightC4Func(1:end,idxs(1):end)=lateTimeWeight;
% % now with the weighting surfaces determined, try to find the scalars that
% % make the 3 surfaces match in magnitude 
% 
% % first get the overall Chi2 of the C2 with weightFunc applied, C4 with
% % weight func applied, and then hist with no weight func applied -all at
% % 10% deviation
% 
% % c2Chi2 = mean(((c2Data-c2Data.*noiseMaskC2).^2).*weightC2func); 
% % c4Chi2 = mean(mean(((c4Data-c4Data.*noiseMaskC4).^2).*weightC4Func));
% % histChi2 = mean((targetHistogram-targetHistogram*(1+(noiseHist))).^2); 
% 
% % Need to implement a toggle for alterante data sets on histogrtma bounds
% % like in chiSqCalc_v4
% c2Chi2 = mean((((c2Data-c2Data.*noiseMaskC2).^2)./dotDivC2).*weightC2func); 
% c4Chi2 = mean(mean((((c4Data-c4Data.*noiseMaskC4).^2)./dotDivSurf).*weightC4Func));
% histChi2 = mean((( (targetHistogram(startIdx:endIdx)-targetHistogram(startIdx:endIdx)*(1+(noiseHist))).^2)./targetHistogram(startIdx:endIdx))); 
% 
% chi2Set=[ c2Chi2 , c4Chi2, histChi2 ];
% chi2Scalars = 1./(chi2Set./max(chi2Set));
% 
% c2Scalar = c2Mag*chi2Scalars(1);
% c4Scalar = c4Mag*chi2Scalars(2);
% histScalar = histMag*chi2Scalars(3);
% 
% % figure()
% % plot(c2Tau, weightC2func);
% % set(gca,'xscale','log');
% % 
% % c4Test = 1./(nthroot(c4Tau,2)).*(1./(nthroot(c4Tau,2)))';
% % figure()
% % surf(c4Tau,c4Tau,weightC4Func./(weightC4Func(1,1)));
% % hold on
% % surf(c4Tau,c4Tau,c4Test./(c4Test(1,1)));
% % view([45 0]); 
% % 
% % set(gca,'xscale','log');
% % set(gca,'yscale','log');
% 
% 
% end