% NOTE: the c2 can be steepened by letting the noiseFloor go to small
% values, while the percentDiff ramp is taken to large values. In the case
% of the c4, the steepness is signifcant nearly always due to the sheer
% number of points increasingly nonlinearly in the successive decacdes of
% tau. Taking the c4 to the "flattened" limit (i.e. noisefloor small and
% noise ramp almost nothing (same or less than floor) causes the c4 weight
% func to saturate in its steepness. Thus a slider is introdcued to take
% the maximal grid at the start down by X% in order to "flatten" the c4
% weight func beyond the limits allowed for by the noise floor and ramp.


%***********************************************************************
%**************current values for PatSea 03-07-2023 *******************
noisePercentDiffC2=0.15;
noisePercentDiffC4=0.01;
% noisePercentDiffC4=0.15;
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
%noiseFloorC4=0.15;
% noiseFloorC4=0.2;
noiseFloorC4=0.4;


% set the relative contributions of the 3 surfaces 
c2Mag=8000*4;
c4Mag=200e3;
histMag=2e5*4;
%********************************************************************

%drop values here for specific data set testing
    statsThisModel=1
        
      c4WeightTime=0.9;
            
            noisePercentDiffC4=0.6;
            noiseFloorC4=0.1;
            noisePercentDiffC2=0.2;
            noiseFloorC2=0.03;
            
            c4Slider=2.5;
            
            C2LateMag=0.45;
            C3LateMag=0.5;
            
            decadeSegments=20;
        
[weightC4Func,weightC2func,c2Scalar,c4Scalar,histScalar] = weightsCalculator_v4_chi2(noiseFloorHist,noiseFloorC2,noiseFloorC4,noisePercentDiffC2,noisePercentDiffC4,C2_exp_x, C2_exp_y, C4_tau1range, C4_tau2eq0_exp, targetHistogram,c2Mag,c4Mag,histMag,decadeSegments,c4Slider,c4WeightTime,Amp_bins,histUpLim,c2WeightTime,C2LateMag,C3LateMag);
% [weightC4Func,weightC2func,c2Scalar,c4Scalar,histScalar] = weightsCalculator_v3_chi2(noiseFloorHist,noiseFloorC2,noiseFloorC4,noisePercentDiffC2,noisePercentDiffC4,C2_exp_x, C2_exp_y, C4_tau1range, C4_tau2eq0_exp, targetHistogram,c2Mag,c4Mag,histMag,decadeSegments,c4Slider)

if statsThisModel
% Do a quick check on the chi2 contribution from a partocualr model fit, to
% check for the regions that are really driving the calculation 
addControlMode=1;
tau2 = 0;
c2Exp = @ (x, A, tau) A*exp(-x/tau);
c4Exp = @ (x, A, tau) 2*((A*exp(-x/tau)).*(A*exp(-x/tau)).');
c2Control = c2Exp(C2_exp_x,c2CtrlAmp,C2C4_ctrlTime)*C2_exp_y(1);
c4Control = c4Exp(C4_tau1range,c4CtrlAmp,C2C4_ctrlTime)*C4_tau2eq0_exp(1,1); 

[C2_sim,C2time] = PA2C2(P,A,C2_exp_x,yoff,addControlMode,c2Control);
[C4_sim,C4time] = PAK2C4(P,A,K,C4_tau1range,tau2,zoff,addControlMode,c4Control);

[c3Chi2Surf, c2Chi2Curve, c3Chi2Segments, c2Chi2Segments] = chi2StatsCalc(weightC2func, weightC4Func, C2_exp_y, C4_tau2eq0_exp, C2_sim, C4_sim, C2_exp_x, C4_tau1range, c2WeightTime, c4WeightTime, decadeSegments);
end

% set(gca,'zscale','log');
%some quick checks on use of an altered C3/C2 surface to "bubble out"
%missed sections
% surf(C4_tau1range,C4_tau1range,((C4_tau2eq0_exp-C4_sim).^2)./(max(max((C4_tau2eq0_exp-C4_sim).^2))));
% surf(C4_tau1range,C4_tau1range,((C4_tau2eq0_exp-C4_sim).^2)./(C4_tau2eq0_exp));
% [idxs,bool]=find(C4_tau1range>0.5); 
% dotDivLate=mean(C4_tau2eq0_exp(idxs(1):idxs(1),1:idxs(1)));
% dotDivSurf=C4_tau2eq0_exp;
% dotDivSurf(idxs(1):end,1:end)=dotDivLate;
% dotDivSurf(1:end,idxs(1):end)=dotDivLate;
%  set(gca,'zscale','log');
% set(gca,'xscale','log');
% set(gca,'yscale','log');
% myNewC3surf=weightC4Func.*(1+((C4_tau2eq0_exp-C4_sim).^2)./(max(max((C4_tau2eq0_exp-C4_sim).^2))));
% myNewC3surf=weightC4Func.*(1+abs((C4_tau2eq0_exp-C4_sim)./C4_tau2eq0_exp)./(max(max((C4_tau2eq0_exp-C4_sim)./C4_tau2eq0_exp))));