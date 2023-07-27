% AUTHOR: Claire Albrecht
% CREATED: September 2021
% PURPOSE: general minimization function for model constructed from model_builder()
%

function [chisquared, chisquared_array, chisquared_unweighted, chisquared_unweighted_array] = multigoaltcf_Nstate_polz_v4_ctrl(x0)
global normalizeMode diagnoseMode
global targetHistogram weightingFactor_FREThist
global C2_exp_x C2_exp_y weightingFactor_C2 weightC2func
global C4_tau1range C4_tau2eq0_exp weightingFactor_C4_t0 wC4func
global fitHistMode fitC2Mode fitC4Mode
global FRET_bins
global Nparam Nstates
global model_name  modelName
global verboseMode ctrlScaling ctrlTime
global c2Exp c4Exp NctrlParams
global c2ctrlAmp_GenAlg c4ctrlAmp_GenAlg c2ctrlTau_GenAlg 


global SelectedTimeRange  paramAdjs addControlMode param_strings weightHistFunc


verbose_mode = 0;

params = x0;
Nparam_ctrlAdj = Nparam - NctrlParams;
% tijs = x(1:(Nparam_ctrlAdj - 2*Nstates));
% A = x((Nparam_ctrlAdj - 2*Nstates+1):Nparam_ctrlAdj-Nstates);
% sigmas = x(((Nparam_ctrlAdj-Nstates)+1:Nparam_ctrlAdj));
% 
% % define control parameters from guesses
% ctrlParams = x(Nparam_ctrlAdj+1:Nparam);
% c2ctrlAmp = ctrlParams(1);
% c4ctrlAmp = ctrlParams(2);
% c2ctrlTau = ctrlParams(3);

Ntijs = Nparam_ctrlAdj - 2*Nstates;

tijs_temp = params(1:(Nparam_ctrlAdj - 2*Nstates));
tijs_new = tijs_temp ;

A_temp = x0((Nparam_ctrlAdj - 2*Nstates+1):Nparam_ctrlAdj-Nstates);
A_new = A_temp;

sigma_A_temp = x0((Nparam_ctrlAdj-Nstates)+1:Nparam_ctrlAdj);
sigma_A_new = sigma_A_temp;

ctrlParams = x0(Nparam_ctrlAdj+1:Nparam);
c2Amp_temp = ctrlParams(1);
c2Amp_new = c2Amp_temp;

c4Amp_temp = ctrlParams(2);
c4Amp_new = c4Amp_temp;

c2c4Time = ctrlParams(3);
c2c4Time_new = c2c4Time;

c2ControlCur = c2Exp(C2_exp_x,c2Amp_new,c2c4Time_new)*C2_exp_y(1);
c2CorrectedSurf = C2_exp_y - c2ControlCur; 
yoff = abs(mean(c2CorrectedSurf(end-6:end)));

c4ControlCur = c4Exp(C4_tau1range,c4Amp_new,c2c4Time_new)*C4_tau2eq0_exp(1,1);        
c4CorrectedSurf = C4_tau2eq0_exp - c4ControlCur; 
zoff = abs(mean(mean(c4CorrectedSurf(end-6:end,end-6:end))));


% [K, A, rates, tijs, ~, ~, ~, ~, ~, ~,~, ~] = model_builder(tijs_new, A_new, modelName);
% [K, A, rates, tijs, Nstates, Nparam, sigma_A, param_strings, boundsArray, maxMutationArray,maxMutationCountsArray, minMutationCountsArray] = model_builderBI(tijs, A,sigma_A,modelName);
% [K, A, rates, tijs, Nstates, Nparam, sigma_A, param_strings, boundsArray, maxMutationArray,maxMutationCountsArray, minMutationCountsArray, atEdges, boundsArray_noAdj] = model_builderBI_fullRange(tijs_new, A_new, sigma_A_new, modelName,SelectedTimeRange, paramAdjs);
[K, A, rates, tijs, Nstates, Nparam, param_strings, boundsArray, maxMutationArray,maxMutationCountsArray, minMutationCountsArray, boundsArray_noAdj] = model_builder_Polz_v4_10ms_PatSea(tijs_new, A_new, model_name,paramAdjs,c2ctrlAmp_GenAlg,c4ctrlAmp_GenAlg,c2ctrlTau_GenAlg);
% boundsArray = [ boundsArray; ctrlScaling ; ctrlTime]; 
            
 
C2time = C2_exp_x;
C4time = C4_tau1range;
% [chisquared,chisquared_array, chisquared_unweighted, chisquared_unweighted_array] = chiSqCalc_v2_BI_v2(K, A, sigma_A, C2_time, C4_time, addControlMode, c2Control, c4Control);
% [chisquared,chisquared_array, chisquared_unweighted, chisquared_unweighted_array,chisquared_Weighted_array] = chiSqCalc_v3_polz(K, A, C2time, C4time,addControlMode,c2ControlCur,c4ControlCur,sigma_A_new);
[chisquared,chisquared_array, chisquared_unweighted, chisquared_unweighted_array,chisquared_Weighted_array] = chiSqCalc_v4_polz(K, A, C2time, C4time,addControlMode,c2ControlCur,c4ControlCur,sigma_A_new,yoff,zoff);
disp(['control scales are:' 'C2Amp=' num2str(c2Amp_new) ' C4Amp=' num2str(c4Amp_new) ' Time=' num2str(c2c4Time_new)]);

if verbose_mode == 1
    disp(['multigoal Nstates: ', num2str(Nstates)])
    disp(['multigoal x0: ', num2str(x0)])
    disp(['multigoal params: ', num2str(params)])
    disp(['multigoal Nparam: ',num2str(Nparam)])
    disp(['multigoal Ntijs: ',num2str(Ntijs)])
    
    disp(['multigoal len tijs: ',num2str(length(tijs_new))])
    disp(['multigoal tijs: ',num2str(tijs_new)])
    
    disp(['multigoal len A: ',num2str(length(A_new))])
    disp(['multigoal A: ',num2str(A_new)])
end


end

% model_name = x0{1};
% disp(['multigoal Nstates: ', num2str(Nstates)])
% disp(['multigoal Ntijs: ',num2str(Ntijs)])

% params = x0{2};
% disp(['params: ',num2str(params)])
% disp(['len params: ',num2str(length(params))])
%
% [t12, t13, t21, t23, t31, A1, A2, A3] = feval(@(x) x{:}, num2cell(params));
%
% tijs = [t12, t13, t21, t23, t31];
% A = [A1, A2, A3];
%
% % tijs = params(1:(Nparam - Nstates));

%
% % A = params((Nparam - Nstates + 1):Nparam);

% tijs_new= [];
% for i = 1:Ntijs
%     tijs_temp = feval(@(x) x{i}, num2cell(params));
%     disp(['tij: ',num2str(tijs_temp)])
%     tijs_new = [tijs_new, tijs_temp];
% end

% A_new= [];
% for i = (Ntijs+1):Nparam
%     A_temp = feval(@(x) x{i}, num2cell(params));
%     disp(['A: ',num2str(A_temp)])
%     A_new = [A_new, A_temp];
% end