% Author: Jack Maurer

% Goal: Generate an average picture across the top N models for the
% relative contributions to the chi2 from the current weighjt surface.
% rescale and recalculate given new surface params, then recalc the
% expected weighted average terms from the new params for a scalar
% adjustment prior to next refinement. 
global fitHistMode fitC2Mode fitC4Mode Amp_bins
fitHistMode=1;
fitC2Mode=1;
fitC4Mode=1;
% supply a baseFolder for the outputs to be looped over at a certain output
% number in the case of PatSea
% **************SET THESE INPUTS FOR YOUR USE CASE**********************
% startFold= 'C:\Users\Ryzen 5\Dropbox\TempHolder_C3_250usecC2\GenAlg_Data\E3E4\300mM_NaCl\genAlgData_UNorm_wCtrlwWeights_10ms_c2_250usec_v4\Outputs';
startFold=  'C:\Users\Ryzen 5\Dropbox\TempHolder_C3_250usecC2\GenAlg_Data\E5E6\20mM_NaCl_0mM_MgCl2\genAlgData_UNorm_wCtrlwWeights_10ms_c2_250usec_v4\Outputs';
outputNumber=12;
PatSea_version_Num=6; 
numModels2Keep=5; 
genAlgMode=0;
PatSeaMode=1;
showWeightSurfs=1;
programName= 'GenAlg_Nstate_updated2';
% ***********************************************************************

% BLOCK FOR THE DEFAULT PARAMETERS THAT EVERY GENALG OR PATSEA RUN TAKE (UNLESS OVERWRITTEN LATER ON)
%***********************************************************************
%**************current values for PatSea 03-07-2023 *******************
noisePercentDiffC2=0.15;
noisePercentDiffC4=0.15;
decadeSegments=10;
histUpLim=0.37;
c4Slider=0.5;
%in Seconds (500ms for now)
c4WeightTime=0.8;
c2WeightTime=1.0;
C2LateMag=0.5;
C3LateMag=1.5; 
% set the base line noise on the c2/c4 as well as hist. less baseline noise
% will steepen the weight function for the c2/c4
noiseFloorHist=0.02;
noiseFloorC2=0.01;
noiseFloorC4=0.4;
% set the relative contributions of the 3 surfaces 
c2Mag=8000*4;
c4Mag=200e3;
histMag=2e5*10;
%********************************************************************

% BLOCK FOR ANY CHANGES/UPDATES/TRIAL AND ERROR PERMUTATIONS OF THE ABOVE PARAMETERS 
% TO SEE HOW THE RELATIVE CHI2 CONTRIBUTIUONS ARE AFFECTED 
        c2Mag=c2Mag*2*5*5*2*2*2*3*10*10*20*10*8*5*8*4*2;
        c4Mag=c4Mag*3*10*5*5*6*4*3*100*15*5*8*4*1.25;
        histMag=histMag*10*10*1e4*10*3*1.25*1.5*1.2*1.5*1.5;
        loadPatSeaMode=1;
        outputToLoad=9;
        histUpLim=0.41;
        %c4WeightTime=0.6;
        c4WeightTime=0.9;
        c4Slider=1.0;
        
        noisePercentDiffC4=0.01;
        noiseFloorC4=0.4;
        noisePercentDiffC2=0.17;
        noiseFloorC2=0.02;
        
        C3LateMag=6.0; 
        C2LateMag=0.65;
        
        %altBoundsID=4;
        decadeSegments=40;
        
        %************secondary set to try out for c2 vs c3 early/late tau issues
        c4WeightTime=1.0;
        c4Slider=1.0;
        
        noisePercentDiffC4=0.05;
        noiseFloorC4=0.1;
        noisePercentDiffC2=0.17;
        noiseFloorC2=0.02;
        
        C3LateMag=1.2; 
        C2LateMag=0.65;
        
         altBoundsID=2; 
        decadeSegments=40;
% ***********************************************************************************
% *********************************************************************************

    histLB=0.0;
    histUB=histUpLim;
    [I_LB,val_LB]=find(Amp_bins>histLB);
    [I_UB,val_UB]=find(Amp_bins<histUB);
    startIdx=val_LB(1);
    endIdx=val_UB(end); 

% input a table of sorted model outcomes from either the GenAlg or PatSea
% (myTableSorted vs myTableSortedPS)then loop over the top N models and add
% the relative contributions from each model to the runnign average.
% Consider also a subplot (for N<15) for each model, to try and rule out
% outliers    

% first call and define the weight functions presently being used (adjusted
% for changed by the above block)
[weightC4Func,weightC2func,c2Scalar,c4Scalar,histScalar] = weightsCalculator_v4_chi2(noiseFloorHist,noiseFloorC2,noiseFloorC4,noisePercentDiffC2,noisePercentDiffC4,C2_exp_x, C2_exp_y, C4_tau1range, C4_tau2eq0_exp, targetHistogram,c2Mag,c4Mag,histMag,decadeSegments,c4Slider,c4WeightTime,Amp_bins,histUpLim,c2WeightTime,C2LateMag,C3LateMag);
c4Tau=C4_tau1range;
c2Tau=C2_exp_x;


if showWeightSurfs
figure(1)
plot(c2Tau, weightC2func);
set(gca,'xscale','log');
c4Test = 1./(nthroot(c4Tau,2)).*(1./(nthroot(c4Tau,2)))';
figure(2)
hold off;
surf(c4Tau,c4Tau,weightC4Func./(weightC4Func(1,1)));
hold on
surf(c4Tau,c4Tau,c4Test./(c4Test(1,1)));
view([45 0]); 
set(gca,'xscale','log');
set(gca,'yscale','log');
end

c3Chi2SurfTot=[];
c2Chi2CurveTot=[];
c3Chi2SegmentsTot=[];
c2Chi2SegmentsTot=[];
chisquared_Weighted_arrayTot=zeros(numModels2Keep,3);

for i=1:numModels2Keep
    if PatSeaMode        
    curModel=myTableSortedPS.(1)(i);
    curModel=curModel{1};
    targetFolderFit=[startFold filesep() curModel filesep() 'PatSea_Nstate_polz_outputs_v' num2str(PatSea_version_Num) filesep() ...
           'output' num2str(outputNumber) filesep() 'BestFitResults_PatSea_Nstate_polz.mat'];
    targetFolderData=[startFold filesep() curModel filesep() 'PatSea_Nstate_polz_outputs_v' num2str(PatSea_version_Num) filesep() ...
           'output' num2str(outputNumber) filesep() 'fitInputData.mat'];
    elseif genAlgMode        
    curModel=myTableSorted.(1)(i);
    curModel=[programName curModel{1}];
    targetFolderFit=[startFold filesep() curModel filesep() 'lowestChiSquare' filesep() 'BestFitResults.mat']
    targetFolderData=[startFold filesep() curModel filesep() 'lowestChiSquare' filesep() 'fitInputData.mat'];
    end
    
   if exist(targetFolderFit)
        varlist = {'weightC4Func','wC4Func','weightC2func','weightingFactor_C2','weightingFactor_C4_t0','weightingFactor_Amphist'}; %Find the variables that already exist
        varlist =strjoin(varlist','$|'); %Join into string, separating vars by '|'
        load(targetFolderFit,'-regexp', ['^(?!' varlist ')\w']);   
        load(targetFolderData,'-regexp', ['^(?!' varlist ')\w']); 
%        load(targetFolderFit,);
%        load(targetFolderData);
   else
       continue;      
   end
% define the controls, and simulated data surfaces for the presenet model
% in the loop
if genAlgMode
    c2CtrlAmp=c2ctrlAmp;
    C2C4_ctrlTime=c2ctrlTau;
    c4CtrlAmp=c4ctrlAmp; 
    sigmas=sigma_A; 
end

addControlMode=1;
tau2 = 0;
c2Exp = @ (x, A, tau) A*exp(-x/tau);
c4Exp = @ (x, A, tau) 2*((A*exp(-x/tau)).*(A*exp(-x/tau)).');
c2Control = c2Exp(C2_exp_x,c2CtrlAmp,C2C4_ctrlTime)*C2_exp_y(1);
c4Control = c4Exp(C4_tau1range,c4CtrlAmp,C2C4_ctrlTime)*C4_tau2eq0_exp(1,1); 
[C2_sim,C2time] = PA2C2(P,A,C2_exp_x,yoff,addControlMode,c2Control);
[C4_sim,C4time] = PAK2C4(P,A,K,C4_tau1range,tau2,zoff,addControlMode,c4Control);
[Peq, denom_hist_sim, hist_sim, Peq_names] = histMaker_Nstate_v3normCorrect(P, A, sigmas, Amp_bins);


if i==1
[idxs]=find(C2time>c2WeightTime); 
dotDivLateC2=mean(C2_exp_y((idxs(1)-1):idxs(1)));
dotDivC2=C2_exp_y;
dotDivC2(idxs(1):end)=dotDivLateC2;
    
[idxs,bool]=find(C4_tau1range>c4WeightTime); 
dotDivLate=mean(C4_tau2eq0_exp(idxs(1):idxs(1),1:idxs(1)));
dotDivSurf=C4_tau2eq0_exp;
dotDivSurf(idxs(1):end,1:end)=dotDivLate;
dotDivSurf(1:end,idxs(1):end)=dotDivLate; 
end

% run the chi2 stats calc on the current model, tabulate total and store
% individual outcomes in z-stack matrix 
[c3Chi2Surf, c2Chi2Curve, c3Chi2Segments, c2Chi2Segments,c2IdxTimes,c3IdxTimes] = chi2StatsCalc(weightC2func, weightC4Func, C2_exp_y, C4_tau2eq0_exp, C2_sim, C4_sim, C2_exp_x, C4_tau1range, c2WeightTime, c4WeightTime, decadeSegments);

% run each model for the chi2 weighted array that corresponds with the
% current choice of weightFunction parameters
% [chisquared,chisquared_array, chisquared_unweighted, chisquared_unweighted_array,chisquared_Weighted_array] = chiSqCalc_v4_polz(K, A, C2time, C4time,addControlMode,c2Control,c4Control,sigmas,yoff,zoff,c4WeightTime);
chisquared_Weighted_arrayTot(i,1) = mean(abs((((hist_sim(startIdx:endIdx)'-targetHistogram(startIdx:endIdx)).^2)./(targetHistogram(startIdx:endIdx)))))*histScalar;
chisquared_Weighted_arrayTot(i,2) = mean( abs((((C2_sim - reshape(C2_exp_y,1,length(C2_exp_y))).^2)./(reshape(dotDivC2,1,length(dotDivC2)))).*reshape(weightC2func,1,length(weightC2func))))*c2Scalar;
chisquared_Weighted_arrayTot(i,3) = mean(mean(abs((((C4_sim - C4_tau2eq0_exp).^2)./(dotDivSurf)).*wC4func)))*c4Scalar;

% chisquared_Weighted_arrayTot=[chisquared_Weighted_arrayTot;chisquared_Weighted_array'];
c3Chi2SurfTot=cat(3,c3Chi2SurfTot,c3Chi2Surf);
c2Chi2CurveTot=[c2Chi2CurveTot;c2Chi2Curve];
c3Chi2SegmentsTot=[c3Chi2SegmentsTot;c3Chi2Segments];
c2Chi2SegmentsTot=[c2Chi2SegmentsTot;c2Chi2Segments];
end

c3Chi2SurfTotAvg=sum(c3Chi2SurfTot,3)./(numModels2Keep);
c2Chi2CurveTotAvg=sum(c2Chi2CurveTot,1)./(numModels2Keep); 
c3Chi2SegmentsTotAvg=sum(c3Chi2SegmentsTot,1)./(numModels2Keep);
c2Chi2SegmentsTotAvg=sum(c2Chi2SegmentsTot,1)./(numModels2Keep);
chisquared_Weighted_arrayAvg=sum(chisquared_Weighted_arrayTot,1)./(numModels2Keep); 

% Make the plots to capture the essential behavior of the surfaces
figure(3)
plot(c2Tau, c2Chi2CurveTotAvg); 
set(gca,'xscale','log');
title('Average Chi2 of C2 across entire curve')
xlabel('Tau (sec)');
ylabel('Mag chi2'); 

figure(4)
surf(c4Tau,c4Tau,c3Chi2SurfTotAvg);
colormap('jet'); 
set(gca,'xscale','log');
set(gca,'yscale','log');
title('Average Chi2 of C3 across entire surface')
xlabel('Tau1 (sec)');
ylabel('Tau2 (sec)'); 
zlabel('Mag chi2'); 

figure(5)
plot(c2IdxTimes,c2Chi2SegmentsTotAvg)
set(gca,'xscale','log');
title('Average Chi2 of C2 across discretized segments from WeightCalc')
xlabel('Tau (sec)');
ylabel('Mag chi2'); 

figure(6)
plot(c3IdxTimes,c3Chi2SegmentsTotAvg) 
set(gca,'xscale','log');
title('Average Chi2 of C3 across discretized segments from WeightCalc')
xlabel('Tau (sec)');
ylabel('Mag chi2'); 






