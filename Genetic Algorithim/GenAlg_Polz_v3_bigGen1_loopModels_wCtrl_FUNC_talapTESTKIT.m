%__________________________________________________________________________
% AUTHOR:  Jack Maurer + Claire Albrecht  (September 2019 Creation)
%
% NAME: GenAlg_3state123_cyclical_Polz.m
%
% FUNCTION: % Fits (1) the Amp Histogram, (2) the 2pt TCF and (3) the 4-point TCF
%
%--------------------------------------------------------------------------
% PROCEDURE:
%--------------------------------------------------------------------------
% (1) Load the experimental histogram
% (2) Load the experimental 2pt time correlation functions to fit
% (3) Load the experimental 4pt TCFs for various tau2 values
% (4) Make an initial population of guesses for the model paramaters
% (5) Construct Histograms, C2, and {C4} for each guess
% (6) Compare each simulated plot to the experimental one to calculate rms
% (7) Choose the best elements of each and mix together.
% (8) Repeat steps 5-7 until the simulation converges below the threshold
%
%--------------------------------------------------------------------------
% INPUT:
%--------------------------------------------------------------------------
% (1) FRET Histogram : '3p15mer_0p0uMgp32_GaussianFit_N3.mat'
% (2) Two-Point TCF : '3p15mer_0p0uMgp32_tcfavg_2exp_fitresult.mat'
% (3) Four-Point TCF : '000010us_tau2-000000us_Neq035_fourptTCFavg.mat'
% (4) Solutions to Diff Eqns: 'symCondProb_3state123_cyclical.mat'
% OUTPUT:
% (1) BestFitResults.mat         %All the fitting paramaters
% (2) BestFitRestults_hist.fig   %
% (3) fitInputData.mat
% (4) genAlgParamaters.mat
% (5) ModelResultsFigure.fig
% (6) plottingParamaters.mat
%
%--------------------------------------------------------------------------
% EXTERNAL PROGRAMS CALLS: (needs these codes in MATLAB PATH to function)
%--------------------------------------------------------------------------
% ANALYTICAL Algorithms (fast)
% (1) histMaker_3state123_cyclical_analytical    % Calculates Peq
% (2) C2maker_3state123_cyclical_analytical      % Calculates C2
% (3) C4maker_3state123_cyclical_analytical      % Calculates C4

% NUMERICAL Algorithms (slow)
% (0) ODEsolver_3state123_cyclical.m    % Calculates the conditional probabilities (symCondProb_3state123_cyclical.mat)
% (1) histMaker_3state123_cyclical      % Calculates Peq
% (2) C2Maker_3state123_cyclical        % Calculates C2
% (3) C4Maker_3state123_cyclical        % Calculates C4
%
%--------------------------------------------------------------------------
% MODIFICATION LOG:
%--------------------------------------------------------------------------
% BI 20190924 Added Claire and I's Updated scripts to calculate Eq Pop.
% BI 20190926 *** Added the yoffset to the simulations of C2 (for better chisquared comparison)
% BI 20191007 Replaced all numerical scripts with updated analytical options
% BI 20191008 Adding the yoffset^2 to the 4 point TCF
% BI 20191016 Added ability to save best three fits instead of just the best fit
% BI 20191018 Saves construct information and adds titles to all plots
% BI 20200130 Saving additional information
% BI 20200130 New definition of yoffset based on data
% JM 20210901: converted to my PC, eliminated all bu MatrixMode cases 
% JM 20210907: converted to the polarization case 
%__________________________________________________________________________
% function [] = GenAlg_Polz_v3_bigGen1_loopModels_wCtrl_FUNC_talapTESTKIT(model_idx,setIdx)
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% addpath statments for talapas use case - 06012022
%   addpath(genpath('/gpfs/home/jmaurer3/Documents/MATLAB/MatlabCodes/KineticModeling/KineticNetwork'));

% [computer_terminal_str, terminalID,~] = computerMode_JM();
computer_terminal_str ='computer_JackLaptop_mode';
terminalID = 'C:\Users\Ryzen 5\';
% if computer_terminal_str(10:18) == 'calbrecht'
%     cd('/Users/calbrecht/Dropbox/MATLAB_programs/claire_programs/from_Jack/updatedgenalgcodes')
%     cc;
% elseif computer_terminal_str(10:15) == 'claire'
%     cc;
% end

programName = 'GenAlg_Nstate_updated2';
save_prefix = 'SaveName';
% disp(['Now Running ' programName '.m']);
close all
%//////////////////////////////////////////////////////////////////////////
% PART 1: Set the program up with various options
%\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

%constructFolderNames = {'E5E6'};
%constructFolderNames = {'E1E2'};
constructFolderNames = {'E3E4'};
% constructFolderNames = {'(+15)Dimer'};
saltConditions='100mM_NaCl';
% saltConditions='TestSalt';

% addpath statments for talapas use case - 06012022
%   addpath(genpath('/gpfs/home/jmaurer3/Documents/MATLAB/MatlabCodes/KineticModeling/KineticNetwork'));

% -------------------------------------------------------------------------
% Choose which model you???d like to use
% -------------------------------------------------------------------------
% model_name = ???ssDNA_3state_cyclical???; % This is called by model_builder()
% model_name = ???ssDNA_3state_cyclical_strictBounds???;
% Nstates = 3;
% ------------------------------------------------------------------------
mod_startIdx=1;
mod_endIdx=1;

global Nstates Nparam
% Nstates = 3;
Nstates = 4;

model_list = { 'model_lin', 'model_loop1','model_loopN','model_brch'}; 
model_set = model_list(setIdx); 
modelNum=model_idx; 
% loadMode=1;
% noise increase on the c2 and c4 at the latest tau compared to th first
% tau (across decades) - in terms of % noise change, i.e. last decade 15%
% noisier than first decade . Base on 1/sqrt(NtotalPairs) of TCF wAVG
% % noisePercentDiffC2=0.1;
% % noisePercentDiffC4=0.1; 
% % 
% % decadeSegments=1;
% % % c4 slider less than 1 is flatter, greate than 1 is steeper (at a fixed
% % % noise floor and ramp) - lower limit of ~0.15 to avoid going below
% % % secondary maxima on c4 weight func surface.
% % c4Slider=1.0;
% % 
% % % set the base line noise on the c2/c4 as well as hist. less baseline noise
% % % will steepen the weight function for the c2/c4
% % noiseFloorHist=0.02;
% % noiseFloorC2=0.1;
% % noiseFloorC4=0.15;

noisePercentDiffC2=0.15;
noisePercentDiffC4=0.15;
decadeSegments=10;
% c4 slider less than 1 is flatter, greate than 1 is steeper (at a fixed
% noise floor and ramp) - lower limit of ~0.15 to avoid going below
% secondary maxima on c4 weight func surface.
%c4Slider=1.5;

c4Slider=1.0;
%noisePercentDiffC4=0.05; 
c4WeightTime=0.5;

% set the base line noise on the c2/c4 as well as hist. less baseline noise
% will steepen the weight function for the c2/c4
noiseFloorHist=0.02;
noiseFloorC2=0.01;
%noiseFloorC4=0.15;
noiseFloorC4=0.2;
% set the relative contributions of the 3 surfaces 
c2Mag=8000;
c4Mag=200e3;
histMag=2e5;

% conditions_arr = {[{'S1S2'}, {'100mM NaCl'}]; [{'S4S5'}, {'100mM NaCl'}]; [{'S18S2'},{'100mM NaCl'}]; [{'S19S5'},{'100mM NaCl'}]};
% constructFolderNames =conditions_arr{condition}{1};
% saltConditions = conditions_arr{condition}{2};
% model_set={ 'linear1324' 'linear1342' 'linear1432' ...
%           'linear1234' 'linear1423' 'linear1243' ...
%           'linear2134' 'linear2143' 'linear2314' ...
%           'linear2413' 'linear3124' 'linear3214' ...
%           'pyramid12' 'pyramid13' 'pyramid14' ...
%           'pyramid21' 'pyramid23' 'pyramid24' ...
%            'pyramid31' 'pyramid32' 'pyramid34' ...
%            'pyramid41' 'pyramid42' 'pyramid43' ...
%            'outer_loop1324' 'outer_loop1234' 'outer_loop1243' ...
%            'fully_connected_minus43' 'fully_connected_minus41' 'fully_connected_minus23' 'fully_connected_minus24' ...
%            'fully_connected_minus21' 'fully_connected_minus31' ...
%              'linear123' 'linear231' 'linear312' 'outer_loop123'};


% cutoff time in Seconds 
cutOffTime=2.5;
startTime=250e-6;
% startTime=1000e-6;
addC2resSuffix=1;
weightVersionNum=2;

% toggle for the histogram integration time
hist500ms=0;
hist10ms=1;
hist100ms=0;

% create a suffix for the output folder
if hist10ms
    fileSuffix='10ms';
elseif hist100ms
    fileSuffix='100ms';
elseif hist500ms
    fileSuffix='500ms';
end 

if addC2resSuffix
    fileSuffix=[fileSuffix '_c2_' num2str(startTime*1e6) 'usec_v' num2str(weightVersionNum)];
end
% values for artificial extension of data (special case - see model builder
% and below code)
endRibbon=20;
ribbonPoints=35;

% *******parameters and functions for addition of control timescales*******
global addControlMode varyCtrlMode NctrlParams ctrlMutateTau ctrlMutateAmpc2 ctrlMutateAmpc4 

% IF VARY CTRL MODE - TURN YOFF AND ZOFF TO NOT GLOBAL
% global yoff zoff

addControlMode=1;
varyCtrlMode=1;
NctrlParams=3;
ctrlMutateTau=0.2;
ctrlMutateAmpc2=0.2; 
ctrlMutateAmpc4=0.2; 



global c2Exp 
c2Exp = @ (x, A, tau) A*exp(-x/tau);
global c4Exp 
% 02102022 - eliminated the 2* since it seemed to be overshooting the c4
% surface
% c4Exp = @ (x, A, tau) ((A*exp(-x/tau)).*(A*exp(-x/tau)).');
c4Exp = @ (x, A, tau) 2*((A*exp(-x/tau)).*(A*exp(-x/tau)).');

% this is the amplitude of the control decay component from a noramlized C2
% fit. The c2Exp then has to be scaled bythe true amplitude of the overall
% C2 at the first tau, the c4Exp has to also be scaled by the first C4
% value with an addtional factor of 2 (for now).
global c2AmpFitVal 
global c2TauFitVal

if startTime==1e-3
c2AmpFitVal = 0.509;
c2TauFitVal = 3; 

elseif startTime==250e-6 || startTime==10000e-6
%     ***** (+1) Dimer values *****
%     based on 4 ExP fit
% c2AmpFitVal = 0.30837;
% c2TauFitVal = 2.5; 

% arbitrary... but closer? +1 Dimer
%c2AmpFitVal = 0.3;
%c2TauFitVal = 3; 

% ***** (+15) Dimer Values ******
% c2AmpFitVal = 0.16;
% c2TauFitVal = 3.6; 

% ***** (-1) Dimer Values ******
%c2AmpFitVal = 0.23;
%c2TauFitVal = 6; 

% ***** (-2) Dimer Values ******
c2AmpFitVal = 0.1;
c2TauFitVal = 5; 

end
% ***********************************************************************

% enter a value for the C2(0) for the data set to be run (the variance of
% the overall data set)
timeZeroExp=0.095;

% Choose GenX method
BrettGenX=0;
JackGenX=1; 

% set FRET mode on to work on old FRET data
global FRETmode
FRETmode=0;
% create an empty parameter boolean array for the intial pass to model
% builder
paramAdjs=[];
recallAdjustmentsMode=1; 
maxWindIters=3;
% still relevant
pointsPerDecade=50; 


% model_name = 'ssDNA_3state_cyclical'; % This is called by model_builder()
% Nstates = 3;


%--------------------------------------------------------------------------
% Genetic Algorithm Patamaters
%--------------------------------------------------------------------------
NmembersInitPop = 200;
gen1Nmembers= 500;
randomizeModeMode = 0;
if randomizeModeMode == 0
    percentReproduce = 25;
    percentToKeep = 50;
    %-------------------DERIVED PARAMATERS------------------
    Nreproduce = round(NmembersInitPop*percentReproduce/100);% Number of individuals who reproduce per generation.
    membersToKeep = NmembersInitPop*percentToKeep/100;
end
Nmutations = 1*NmembersInitPop; % Number of mutations per generation.

maxGenerations = 1000;
maxRepeats = 35;
threshold = 0.01;   % Minimal allowable percentage difference before program quits.
maxIterations = 2500;
forceMoreIterationsMode = 1;
NumberToExtend = 1;
%NOTES: lowering the percent that reproduces but raising the percent to
%keep as high seems to help.

%--------------------------------------------------------------------------
% Declare global variables
%--------------------------------------------------------------------------
global normalizeMode verboseMode guessUpdateMode diagnoseMode useAnalyticalAlgorithmsMode useMatrixMethodMode testTimesMode
global genNum fitHistMode fitC2Mode fitC4Mode
global targetHistogram weightingFactor_Amphist Amp_bins weightHistFunc % Histogram optimization
global C2_exp_x C2_exp_y weightingFactor_C2  weightC2func 
global C4_tau1range C4_tau2eq0_exp weightingFactor_C4_t0 wC4func 
global sigma_A 

%--------------------------------------------------------------------------
% User Options                                                  (Part 1)
%--------------------------------------------------------------------------
normalizeMode = 0;%Normalizes the 2pt and the 4pt correlation functions
verboseMode = 0;
guessUpdateMode = 0; %Very detailed: shows changes to any paramaters
diagnoseMode = 0; %Shows relative contribution of the 2pt TCF to the histogram at the end of each chi square calculation
testTimesMode=1; %will adjust the C2 and C4 time range so it can be matched

% if testTimesMode
%     startTime=1e-3;
% end

clockMode = 0; %Times various features of the algorithm
saveModeDefault = 1;
plotMode = 1;                               %Makes plots.
useFigurePosnMode = 1;

plot_Gen1_memberGuessesMode = 0; % figure(1-2-3) Plots Histograms, C2, C4 of gen1 (SLOW)
plot_Gen1ChiSquaredMode = 0;%figure(4)  Adds a point to the chisquare vs member plot(plot4) SLOW %Plot statistics from the first generation (informative)

plot_GenerationalProgressMode = 1;%figure(4) Gray circles on figure(4) %Plot a generational update with best fit from the current generation
plot_paramater_histogram_mode = 1;%figure(5)

plot_model_GenerationalUpdatesMode = 0; %Figure(6) Makes a model of the system with rates between states EACH GENERATION
plot_GenerationalChiSquaredHistogramMode =0; %figure(7) Histogram of chisquare values
pauseBetweenGenerationMode = 0; %Requires manual button press to continue

%At the end of the program
showBestFitMode = 0; %Replaces the current guess with the best guess (end of code)

fitHistMode = 1; fitHistData_mode = 1;%If 0 you will fit to the histogram fit
fitC2Mode = 1;
fitC4Mode = 1;

% UPDATE: 3-14-2022 : changed the weights to maintain a prodcut of ~15 but
% cause the hist to be weight 10x less, the c2 5x more and the c4 5x less -
% did this based on runaway case in fixedDecade sliding where the c4 fit
% dominates the windows movement

%Choose Weighting Amount
% weightingFactor_Amphist = 5e-1;% Weighting for Amp hist comparison
% weightingFactor_C2 = 1.502;
% weightingFactor_C4_t0 = 1.99e4;

% weightingFactor_Amphist = 5e-1;% Weighting for Amp hist comparison
% weightingFactor_C2 = 0.3;
% weightingFactor_C4_t0 = 9.95e4;

if sum([fitHistMode,fitC2Mode,fitC4Mode]) ==0
    error('No surfaces to optimize to! Pick at least one.');
end


%**************************************************************************
useAnalyticalAlgorithmsMode = 	2; % 0: analytical, 1: numerical, 2: for matrix method
useMatrixMethodMode = 1;


if  useMatrixMethodMode == 1
    if verboseMode == 1
        disp('     ***using Matrix method algorithms');
    end
end
%**************************************************************************

if verboseMode == 1
    fprintf('     Each Generation will have %d/%d of its members reproduce(%d%%) \r',Nreproduce,NmembersInitPop,percentReproduce);
    fprintf('     Each Generation will see a total of %d mutations. (%d per individual on average).\r',Nmutations,Nmutations/NmembersInitPop);
    fprintf('     Each Generation will have the top %d percent unmolested (%d/%d individuals).\r',percentToKeep,membersToKeep,NmembersInitPop);
end

% start the loop that goes over all the Models

for model_Idx = mod_startIdx:mod_endIdx
% model_name= model_set{model_Idx}; 
%//////////////////////////////////////////////////////////////////////////
% PART 1: Model Specific Paramaters
%\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
A = zeros([Nstates, 1]);
tijs = ones([5*Nstates, 1]);
if hist500ms
  [~, ~, ~, ~, Nstates, Nparam, sigma_A, param_strings, boundsArray, maxMutationArray,maxMutationCountsArray, minMutationCountsArray,atEdges, boundsArray_noAdj] = model_builder_Polz_v3_500ms(tijs, A, model_name,paramAdjs);            
elseif hist10ms
    if varyCtrlMode == 1
        c2ctrlAmp_0= c2AmpFitVal;
        c4ctrlAmp_0= c2AmpFitVal;
        c2ctrlTau_0= c2TauFitVal; 
        c2ctrlAmp= c2AmpFitVal;
        c4ctrlAmp= c2AmpFitVal;
        c2ctrlTau= c2TauFitVal; 
        
   [model_lin,model_loop1,model_loopN,model_brch] = model_generator_v5(Nstates,model_set);
   
   % initialize dummy variables
    loadMode = 1;
    K = zeros(Nstates);
    model_name = 'model_name_temp';
    dbCell={};
    rates_str_woDB = 'rates_str_temp';
    rates_idx_woDB = [];
    rates_idx_DB = 0;
    Nparam = 3*Nstates;
    param_strings = {};
    boundsArray = [];
    maxMutationArray = [];
    maxMutationCountsArray = [];
    minMutationCountsArray = [];
    atEdges = [];
    boundsArray_noAdj = [];
[K, rates, tijs, Nstates, Nparam, param_strings, boundsArray, maxMutationArray,maxMutationCountsArray, minMutationCountsArray, atEdges, boundsArray_noAdj, model_name,dbCell,rates_str_woDB,rates_idx_woDB,rates_idx_DB] = libraryLoader(loadMode, model_set, modelNum , Nstates, tijs, c2ctrlAmp, c2ctrlTau, c4ctrlAmp, K, model_name,dbCell,rates_str_woDB,rates_idx_woDB,rates_idx_DB,Nparam,param_strings, boundsArray, maxMutationArray,maxMutationCountsArray, minMutationCountsArray, atEdges, boundsArray_noAdj,model_lin,model_loop1,model_loopN,model_brch) ;   
    loadMode = 0;
        
%         [K, A, rates, tijs, Nstates, Nparam, param_strings, boundsArray, maxMutationArray,maxMutationCountsArray, minMutationCountsArray,atEdges, boundsArray_noAdj] = model_builder_Polz_v4_10ms(tijs, A, model_name,paramAdjs,c2ctrlAmp_0,c2ctrlTau_0,c4ctrlAmp_0);            
%     [K, rates, tijs, Nstates, Nparam, param_strings, boundsArray, maxMutationArray,maxMutationCountsArray, minMutationCountsArray, atEdges, boundsArray_noAdj, model_name] = libraryLoader(loadMode, model_set, modelNum , Nstates, tijs, c2ctrlAmp, c2ctrlTau, c4ctrlAmp);
    else
        [~, ~, ~, ~, Nstates, Nparam, param_strings, boundsArray, maxMutationArray,maxMutationCountsArray, minMutationCountsArray,atEdges, boundsArray_noAdj] = model_builder_Polz_v3_10ms(tijs, A, model_name,paramAdjs);            
    end
    end
% [~, ~, ~, ~, Nstates, Nparam, sigma_A, param_strings, boundsArray, maxMutationArray,maxMutationCountsArray, minMutationCountsArray,atEdges, boundsArray_noAdj] = model_builder_Polz_v3(tijs, A, model_name,paramAdjs);

% block to set and intialize the relevant parameters and arrays for a
% moving window in the parameter boundaries - must come after the first
% instance of model builder to get the dimensions correct (except for
% paramBool)

% just going to adjust the rates for now (amps leftoff)
paramAdjs=zeros(1,Nparam-2*Nstates); 
histCell=cell(1,Nparam-2*Nstates);
% interaton tracker to hold values within histCell until every Nth
% iteration (5 to start with)
windowIter=0; 

%Start in the single molecule folder (smData_Processed): comp specific
wd = pwd;
if contains(wd,'C:\Users\baimi\')%Work computer
    computer_baimi_mode = 1;
    computer_terminal_str = 'computer_baimi_mode';
    if verboseMode == 1
        disp(['Turning on on ' computer_terminal_str]);
    end
elseif contains(wd,'/Users/bisraels')%macbook
    computer_bisraels_mode = 1;
    computer_terminal_str = 'computer_bisraels_mode';
    if verboseMode == 1
        disp(['Turning on on ' computer_terminal_str]);
    end
elseif contains(wd,'/Users/clairealbrecht')%claire
    computer_claire_mode = 1;
    computer_terminal_str = 'computer_claire_mode';
    if verboseMode == 1
        disp(['Turning on on ' computer_terminal_str]);
    end
else
    general_computer_mode = 1;
    computer_terminal_str = 'computer_general_mode';
    if verboseMode == 1
        disp(['Turning on on ' computer_terminal_str]);
        useFigurePosnMode = 0;
        disp('Setting useFigurePosnMode to zero');
    end
end

% [computer_terminal_str, terminalID] = computerMode_JM();
computer_terminal_str ='computer_JackLaptop_mode';
terminalID = 'C:\Users\Ryzen 5\';

NconstructFolderNames = numel(constructFolderNames);
for construct_idx = 1:NconstructFolderNames
    %         for construct_idx = NconstructFolderNames:-1:1
%     constructName = char(constructFolderNames(construct_idx));
    constructName = constructFolderNames;

    disp(['       Construct ' num2str(construct_idx) '/' num2str(NconstructFolderNames) ': ' constructName]);
    
    if exist('constructName','var') ~= 1
        error('Not sure what construct to analyze');
    end
    
    %--------------------------------------------------------------------------
    % Navigate to the correct folder
    %--------------------------------------------------------------------------
    switch computer_terminal_str
        case 'computer_baimi_mode' %Work Desktop
            GenAlgFilesLocationPrefix =  ['C:\Users\baimi\Dropbox\MarcusLab\Data\smData_Processed\' constructName '\gp32_0p0uM\ChosenMolecules\genAlgFiles\3state123_cyclical'];
        case 'computer_bisraels_mode' %Macbook pro
            GenAlgFilesLocationPrefix =  ['/Users/bisraels/Dropbox/MarcusLab/Data/smData_Processed/' constructName '/gp32_0p0uM/ChosenMolecules/genAlgFiles/3state123_cyclical'];
        case 'computer_claire_mode' %Macbook pro
%             GenAlgFilesLocationPrefix =  ['/Users/clairealbrecht/Dropbox/MarcusLab/Data/smData_Processed/' constructName '/gp32_0p0uM/ChosenMolecules/genAlgFiles/3state123_cyclical'];
           GenAlgFilesLocationPrefix = ['/Users/clairealbrecht/Dropbox/MATLAB_programs/claire_programs/from_Jack/updatedgenalgcodes/outputs/' programName '_' model_name];
        case 'computer_calbrecht_mode' %Macbook pro
%             GenAlgFilesLocationPrefix =  ['/Users/clairealbrecht/Dropbox/MarcusLab/Data/smData_Processed/' constructName '/gp32_0p0uM/ChosenMolecules/genAlgFiles/3state123_cyclical'];
           GenAlgFilesLocationPrefix = ['/Users/calbrecht/Dropbox/MATLAB_programs/claire_programs/from_Jack/updatedgenalgcodes/outputs/' programName '_' model_name];
        case 'computer_JackLaptop_mode' %Jacks Laptop
            if normalizeMode
            GenAlgFilesLocationPrefix =  ['C:\Users\Ryzen 5\Dropbox\chosenAPD2mat_output\' constructName '\gp32_0p0uM\FPGA\' saltConditions '\genAlgData\Outputs\' programName '_' model_name];
            else
                if testTimesMode                   
                    GenAlgFilesLocationPrefix =  ['C:\Users\Ryzen 5\Dropbox\chosenAPD2mat_output\' constructName '\gp32_0p0uM\FPGA\' saltConditions '\genAlgData_UNorm_wCtrlwWeights_' fileSuffix '\Outputs\' programName '_' model_name];
                   
%                GenAlgFilesLocationPrefix =  ['C:\Users\jmaurer3\Documents\KineticModeling\Data\' constructName filesep() saltConditions '\genAlgData_UNorm_test\Outputs\' programName '_' model_name];
                else
               GenAlgFilesLocationPrefix =  ['C:\Users\Ryzen 5\Dropbox\chosenAPD2mat_output\' constructName '\gp32_0p0uM\FPGA\' saltConditions '\genAlgData_UNorm\Outputs\' programName '_' model_name];
                end
%             GenAlgFilesLocationPrefix =  ['C:\Users\Jack\Dropbox\chosenAPD2mat_output\' constructName '\gp32_0p0uM\FPGA\' saltConditions '\genAlgData_UNorm\Outputs\' programName '_' model_name];
            end
        case 'computer_labTower_mode' %lab desktop
            if normalizeMode
            GenAlgFilesLocationPrefix =  ['C:\Users\jmaurer3\Documents\KineticModeling\Data\' constructName filesep() saltConditions '\genAlgData\Outputs\' programName '_' model_name];
            else
                if testTimesMode                    
                    GenAlgFilesLocationPrefix =  ['C:\Users\jmaurer3\Documents\KineticModeling\Data\' constructName filesep() saltConditions '\genAlgData_UNorm_wCtrlwWeights_' fileSuffix '\Outputs\' programName '_' model_name];                    
%                GenAlgFilesLocationPrefix =  ['C:\Users\jmaurer3\Documents\KineticModeling\Data\' constructName filesep() saltConditions '\genAlgData_UNorm_test\Outputs\' programName '_' model_name];
                else
               GenAlgFilesLocationPrefix =  ['C:\Users\jmaurer3\Documents\KineticModeling\Data\' constructName filesep() saltConditions '\genAlgData_UNorm\Outputs\' programName '_' model_name];
                end
            end
            
        case 'computer_vostro1_mode'             
            if normalizeMode
            GenAlgFilesLocationPrefix =  ['C:\Users\Public\Documents\MATLAB\GenAlg\Data\' constructName '\' saltConditions '\genAlgData\Outputs\' programName '_' model_name];
            else
                if testTimesMode                   
                    GenAlgFilesLocationPrefix =  ['C:\Users\Public\Documents\MATLAB\GenAlg\Data\' constructName '\' saltConditions '\genAlgData_UNorm_wCtrlwWeights_' fileSuffix '\Outputs\' programName '_' model_name];
                   
%                GenAlgFilesLocationPrefix =  ['C:\Users\jmaurer3\Documents\KineticModeling\Data\' constructName filesep() saltConditions '\genAlgData_UNorm_test\Outputs\' programName '_' model_name];
                else
               GenAlgFilesLocationPrefix =  ['C:\Users\Public\Documents\MATLAB\GenAlg\Data\' constructName '\' saltConditions '\genAlgData_UNorm\Outputs\' programName '_' model_name];
                end
%             GenAlgFilesLocationPrefix =  ['C:\Users\Jack\Dropbox\chosenAPD2mat_output\' constructName '\gp32_0p0uM\FPGA\' saltConditions '\genAlgData_UNorm\Outputs\' programName '_' model_name];
            end
            
        case 'computer_labMAC_mode'             
            if normalizeMode
            GenAlgFilesLocationPrefix =  ['/Users/ahmarcus/Documents/MATLAB/GenAlg/Data/' constructName '/' saltConditions '/genAlgData/Outputs/' programName '_' model_name];
            else
                if testTimesMode                   
                    GenAlgFilesLocationPrefix =  ['/Users/ahmarcus/Documents/MATLAB/GenAlg/Data/' constructName '/' saltConditions '/genAlgData_UNorm_wCtrlwWeights_' fileSuffix '/Outputs/' programName '_' model_name];
                   
%                GenAlgFilesLocationPrefix =  ['C:\Users\jmaurer3\Documents\KineticModeling\Data\' constructName filesep() saltConditions '\genAlgData_UNorm_test\Outputs\' programName '_' model_name];
                else
               GenAlgFilesLocationPrefix =  ['/Users/ahmarcus/Documents/MATLAB/GenAlg/Data/' constructName '/' saltConditions '/genAlgData_UNorm/Outputs/' programName '_' model_name];
                end
%             GenAlgFilesLocationPrefix =  ['C:\Users\Jack\Dropbox\chosenAPD2mat_output\' constructName '\gp32_0p0uM\FPGA\' saltConditions '\genAlgData_UNorm\Outputs\' programName '_' model_name];
            end
            
        case 'computer_talapasDesktop_mode' %Jacks Laptop
            if normalizeMode
            GenAlgFilesLocationPrefix =  ['C:\Users\Jack\Dropbox\chosenAPD2mat_output\' constructName '\gp32_0p0uM\FPGA\' saltConditions '\genAlgData\Outputs\' programName '_' model_name];
            else
                if testTimesMode                   
                    GenAlgFilesLocationPrefix =  ['/gpfs/home/jmaurer3/Documents/MATLAB/GenAlg_Data/' constructName '/' saltConditions '/genAlgData_UNorm_wCtrlwWeights_' fileSuffix '/Outputs/' programName '_' model_name];
                   
%                GenAlgFilesLocationPrefix =  ['C:\Users\jmaurer3\Documents\KineticModeling\Data\' constructName filesep() saltConditions '\genAlgData_UNorm_test\Outputs\' programName '_' model_name];
                else
               GenAlgFilesLocationPrefix =  ['C:\Users\Jack\Dropbox\chosenAPD2mat_output\' constructName '\gp32_0p0uM\FPGA\' saltConditions '\genAlgData_UNorm\Outputs\' programName '_' model_name];
                end
%             GenAlgFilesLocationPrefix =  ['C:\Users\Jack\Dropbox\chosenAPD2mat_output\' constructName '\gp32_0p0uM\FPGA\' saltConditions '\genAlgData_UNorm\Outputs\' programName '_' model_name];
            end
    end
    if exist(GenAlgFilesLocationPrefix,'dir') ~= 7
        disp('Making a directory to hold the outputs of the model');
        mkdir(GenAlgFilesLocationPrefix);
    end
    cd(GenAlgFilesLocationPrefix);
    wd = pwd;
    if verboseMode == 1
        disp(['You are now in the folder: ' wd]);
    end
    
    adjustList=dir('*ParameterAdjustments');
    if recallAdjustmentsMode && isempty(adjustList)~=1
        load(adjustList(1).name);
    end
    
    if hist500ms
       [K, A, rates, tijs, Nstates, Nparam, sigma_A, param_strings, boundsArray, maxMutationArray,maxMutationCountsArray, minMutationCountsArray,atEdges, boundsArray_noAdj] = model_builder_Polz_v3_500ms(tijs, A, model_name,paramAdjs);                
    elseif hist10ms
        if varyCtrlMode == 1
%       [K, A, rates, tijs, Nstates, Nparam, param_strings, boundsArray, maxMutationArray,maxMutationCountsArray, minMutationCountsArray,atEdges, boundsArray_noAdj] = model_builder_Polz_v4_10ms(tijs, A, model_name,paramAdjs, c2ctrlAmp_0,c2ctrlTau_0,c4ctrlAmp_0);            
[K, rates, tijs, Nstates, Nparam, param_strings, boundsArray, maxMutationArray,maxMutationCountsArray, minMutationCountsArray, atEdges, boundsArray_noAdj, model_name,dbCell,rates_str_woDB,rates_idx_woDB,rates_idx_DB] = libraryLoader(loadMode, model_set, modelNum , Nstates, tijs, c2ctrlAmp, c2ctrlTau, c4ctrlAmp, K, model_name,dbCell,rates_str_woDB,rates_idx_woDB,rates_idx_DB,Nparam,param_strings, boundsArray, maxMutationArray,maxMutationCountsArray, minMutationCountsArray, atEdges, boundsArray_noAdj,model_lin,model_loop1,model_loopN,model_brch) ;     
        else
      [K, A, rates, tijs, Nstates, Nparam, param_strings, boundsArray, maxMutationArray,maxMutationCountsArray, minMutationCountsArray,atEdges, boundsArray_noAdj] = model_builder_Polz_v3_10ms(tijs, A, model_name,paramAdjs);            
        end
        end
        
    
    %//////////////////////////////////////////////////////////////////////////
    %% PART 2: Load the target data (and plot if opted)
    %\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
    %--------------------------------------------------------------------------
    % Filenames of the data (PART 2: Load the target data)
    %--------------------------------------------------------------------------   
    switch computer_terminal_str
        case 'computer_baimi_mode' %Work Desktop
            genAlgDataLocationPrefix =  ['C:\Users\baimi\Dropbox\MarcusLab\Data\smData_Processed\' constructName '\gp32_0p0uM\ChosenMolecules\genAlgFiles\targetData'];
        case 'computer_bisraels_mode' %Macbook pro
            genAlgDataLocationPrefix =  ['/Users/bisraels/Dropbox/MarcusLab/Data/smData_Processed/' constructName '/gp32_0p0uM/ChosenMolecules/genAlgFiles/targetData'];
        case 'computer_claire_mode'
%             genAlgDataLocationPrefix =  ['/Users/clairealbrecht/Dropbox/MarcusLab/Data/smData_Processed/' constructName '/gp32_0p0uM/ChosenMolecules/genAlgFiles/targetData'];
            genAlgDataLocationPrefix =  ['/Users/clairealbrecht/Dropbox/MATLAB_programs/claire_programs/from_Jack/updatedgenalgcodes/targetData'];
        case 'computer_calbrecht_mode'
            genAlgDataLocationPrefix =  ['/Users/calbrecht/Dropbox/MATLAB_programs/claire_programs/from_Jack/updatedgenalgcodes/targetData'];
            
        case 'computer_JackLaptop_mode' %Jacks Laptop
            if normalizeMode
            genAlgDataLocationPrefix =  ['C:\Users\Ryzen 5\Dropbox\chosenAPD2mat_output\' constructName '\gp32_0p0uM\FPGA\' saltConditions '\genAlgData\targetData'];
            else
                if testTimesMode                    
                    genAlgDataLocationPrefix =  ['C:\Users\Ryzen 5\Dropbox\chosenAPD2mat_output\' constructName '\gp32_0p0uM\FPGA\' saltConditions '\genAlgData_UNorm_wCtrlwWeights_' fileSuffix '\targetData'];                   
                else                    
                    genAlgDataLocationPrefix =  ['C:\Users\Ryzen 5\Dropbox\chosenAPD2mat_output\' constructName '\gp32_0p0uM\FPGA\' saltConditions '\genAlgData_UNorm\targetData'];
                end
%             genAlgDataLocationPrefix =  ['C:\Users\Jack\Dropbox\chosenAPD2mat_output\' constructName '\gp32_0p0uM\FPGA\' saltConditions '\genAlgData_UNorm\targetData'];
            end            
        case 'computer_labTower_mode' %Jacks Laptop
            if normalizeMode
            genAlgDataLocationPrefix =  ['C:\Users\jmaurer3\Documents\KineticModeling\Data\' constructName filesep() saltConditions '\genAlgData\targetData'];
            else
                if testTimesMode                    
                    genAlgDataLocationPrefix =  ['C:\Users\jmaurer3\Documents\KineticModeling\Data\' constructName filesep() saltConditions '\genAlgData_UNorm_testwCtrl_' fileSuffix '\targetData'];                    
                else                    
                    genAlgDataLocationPrefix =  ['C:\Users\jmaurer3\Documents\KineticModeling\Data\' constructName filesep() saltConditions '\genAlgData_UNorm\targetData'];
                end
            end
            
          case 'computer_vostro1_mode'             
             if normalizeMode
                 genAlgDataLocationPrefix =  ['C:\Users\Public\Documents\MATLAB\GenAlg\Data\' constructName '\' saltConditions  '\genAlgData\targetData'];
            else
                if testTimesMode                    
                    genAlgDataLocationPrefix =  ['C:\Users\Public\Documents\MATLAB\GenAlg\Data\' constructName '\' saltConditions  '\genAlgData_UNorm_testwCtrl_' fileSuffix '\targetData'];                   
                else                    
                    genAlgDataLocationPrefix =  ['C:\Users\Public\Documents\MATLAB\GenAlg\Data\' constructName '\' saltConditions  '\genAlgData_UNorm\targetData'];
                end
%             genAlgDataLocationPrefix =  ['C:\Users\Jack\Dropbox\chosenAPD2mat_output\' constructName '\gp32_0p0uM\FPGA\' saltConditions '\genAlgData_UNorm\targetData'];
             end  
            
          case 'computer_labMAC_mode'             
             if normalizeMode
                 genAlgDataLocationPrefix =  ['/Users/ahmarcus/Documents/MATLAB/GenAlg/Data/' constructName '/' saltConditions   '/genAlgData/targetData'];
            else
                if testTimesMode                    
                    genAlgDataLocationPrefix =  ['/Users/ahmarcus/Documents/MATLAB/GenAlg/Data/' constructName '/' saltConditions   '/genAlgData_UNorm_testwCtrl_' fileSuffix '/targetData'];                   
                else                    
                    genAlgDataLocationPrefix =  ['/Users/ahmarcus/Documents/MATLAB/GenAlg/Data/' constructName '/' saltConditions   '/genAlgData_UNorm/targetData'];
                end
%             genAlgDataLocationPrefix =  ['C:\Users\Jack\Dropbox\chosenAPD2mat_output\' constructName '\gp32_0p0uM\FPGA\' saltConditions '\genAlgData_UNorm\targetData'];
            end
            
        case 'computer_talapasDesktop_mode' %Jacks Laptop
            if normalizeMode
            genAlgDataLocationPrefix =  ['C:\Users\Jack\Dropbox\chosenAPD2mat_output\' constructName '\gp32_0p0uM\FPGA\' saltConditions '\genAlgData\targetData'];
            else
                if testTimesMode                    
                    genAlgDataLocationPrefix =  ['/gpfs/home/jmaurer3/Documents/MATLAB/GenAlg_Data/' constructName '/' saltConditions '/genAlgData_UNorm_wCtrlwWeights_' fileSuffix '/targetData'];                   
                else                    
                    genAlgDataLocationPrefix =  ['C:\Users\Jack\Dropbox\chosenAPD2mat_output\' constructName '\gp32_0p0uM\FPGA\' saltConditions '\genAlgData_UNorm\targetData'];
                end
%             genAlgDataLocationPrefix =  ['C:\Users\Jack\Dropbox\chosenAPD2mat_output\' constructName '\gp32_0p0uM\FPGA\' saltConditions '\genAlgData_UNorm\targetData'];
            end    
        otherwise
            error('Cannot detect computer');
    end
    if exist(genAlgDataLocationPrefix,'dir') ~= 7
        error(['Cannot locate data folder ' genAlgDataLocationPrefix]);
    end
    % disp(['Going to ' dataLocationPrefix])
    % cd(dataLocationPrefix);
    if verboseMode == 1
        disp(['Looking for data in ' genAlgDataLocationPrefix]);
    end
    
    
    HistfNameKeyWord = '*HistData.mat';
    if FRETmode==1
    HistfNameKeyWord = '*_0p0uMgp32_001000us_FREThistogram.mat';
    end
    fileNames = dir([genAlgDataLocationPrefix filesep() HistfNameKeyWord]);
    if isempty(fileNames) == 1
        error(['Could not Find any files with the name '  HistfNameKeyWord]);
    end
    histogram_FileName = fileNames(1).name;
    histogram_FilePath = [genAlgDataLocationPrefix filesep() histogram_FileName];
    
    C2fNameKeyWord = '*StitchTCF__wAVG.mat';
    if FRETmode==1
    C2fNameKeyWord = '*_0p0uMgp32_tcfavg_2exp_fitresult.mat';
    end
    fileNames = dir([genAlgDataLocationPrefix filesep() C2fNameKeyWord]);
    if isempty(fileNames) == 1
        error(['Could not Find any files with the name '  C2fNameKeyWord]);
    end
    
    C2_FileName = fileNames(1).name;
    C2_FilePath = [genAlgDataLocationPrefix filesep() C2_FileName];
    
    % fnamekeyword = '000010us_tau2-000000us_Neq035_fourptTCFavg.mat';
    C4fNameKeyWord = 'C4*.mat';
    if FRETmode==1
    C4fNameKeyWord = '*0p0uMgp32_tau2eq000000us_C4.mat';
    end
    fileNames = dir([genAlgDataLocationPrefix filesep() C4fNameKeyWord]);
    if isempty(fileNames) == 1
        error(['Could not Find any files with the name '  C4fNameKeyWord]);
    end
    C4_FileName = fileNames(1).name;
    C4_FilePath = [genAlgDataLocationPrefix filesep() C4_FileName];
    
    %--------------------------------------------------------------------------
    % Set up the figure positions once and done
    %--------------------------------------------------------------------------
    if plotMode == 1
        clear figPosn
        if useFigurePosnMode == 1
%             figPosn = setFigPosn_JM(terminalID);
    fig_1_posn = [ 72 ;  575  ; 560   ; 420];
    fig_2_posn = [675  ; 576 ;  560  ; 420];
    fig_3_posn = [1340 ; 574  ; 560  ; 420];
    fig_4_posn = [16  ; 48 ; 560 ; 420];
    fig_5_posn = [412  ;  52  ; 560 ;  420];
    fig_6_posn = [858  ;  55  ; 560  ; 420];
    fig_7_posn = [1415 ;  56  ;  560 ;  420];
    figPosn = [fig_1_posn, fig_2_posn, fig_3_posn, fig_4_posn, fig_5_posn, fig_6_posn, fig_7_posn];
            if fitHistMode == 1
                figure(1);
                clf;
                set(gcf,'Name','Amp Histogram');
                set(gcf,'Color','w');
                set(gcf,'Position',figPosn(:,1));
            end
            %-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.
            if fitC2Mode == 1
                figure(2);
                clf;
                set(gcf,'Name','Two-Point TCF: C2');
                set(gcf,'Position',figPosn(:,2));
            end
            %-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.
            if fitC4Mode == 1
                figure(3)
                clf;
                set(gcf,'Name','Four-Point TCF: C4');
                set(gcf,'Position',figPosn(:,3));
            end
            %-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.
            if plot_model_GenerationalUpdatesMode == 1
                figure(4);
                clf;
                set(gcf,'Name','Fitting Progress');
                set(gcf,'Position',figPosn(:,4));
            end
        end
        %-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.
        if plot_paramater_histogram_mode == 1
            figure(5);
            clf;
            set(gcf,'Name',' Distribution of initial paramaters : Gen 1');
            set(gcf,'Position',figPosn(:,5));
            set(gcf,'Color','w');
        end
        %-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.
        % Fig 6: Model of the network with rates and Amp values
        if plot_model_GenerationalUpdatesMode == 1
            figure(6);
            clf;
            set(gcf,'Color','w');
            set(gcf,'Name','Model: 3state123 cyclical');
            if useFigurePosnMode == 1
                set(gcf,'Position',figPosn(:,6));
            end
        end
        %-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.
        if plot_GenerationalChiSquaredHistogramMode == 1
            figure(7);
            clf;
            set(gcf,'Name','Chi-squared values guesses');
            set(gcf,'Position',figPosn(:,7));
            set(gcf,'Color','w');
        end
        
        %-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.
    end
    %--------------------------------------------------------------------------
    % (1) Optimization Target #1: 1-D Amp HISTOGRAM (PART 2: Load the target data)
    %--------------------------------------------------------------------------
    % if fitHistMode == 1
    %     load(histogram_FileName,'xData','yData','xFit','yFit');
    if FRETmode==1
    load(histogram_FilePath,'xData','yData');
    else
    load(histogram_FilePath,'edgesAmp','totalAmpHist');
    end
    if fitHistData_mode == 1
        if FRETmode==1
        Amp_bins = xData;
        targetHistogram = yData; 
        edgesAmp=xData;
        else
        Amp_bins = edgesAmp;
        targetHistogram = totalAmpHist;
        end
    elseif fitHistData_mode == 0
        Amp_bins = xFit;
        targetHistogram = yFit;
    end
    
    if strcmp(model_name, 'linear2341_testCutofff_Invert')
        targetHistogram=fliplr(targetHistogram);
    end
    
    targetHistogram = targetHistogram./sum(targetHistogram);
    
%     weight function based on inverse of curve height, suffers due ot
%     edging effects
%     weightHistFunc=1./(targetHistogram);
%     weightHistFunc(isinf(weightHistFunc))=1; 
%     weightHistFunc=weightHistFunc./(sum(weightHistFunc)); 
    
% stepwise weight function (2 windows)
    weightHistFunc=ones(1,length(targetHistogram));
    
     if strcmp(model_name, 'linear2341_testCutofff_Invert')
    weightHistFunc((Amp_bins<=(1-0.21)))=3;   
    else 
    weightHistFunc((Amp_bins>=0.21))=3;
     end   
    weightHistFunc=weightHistFunc./(sum(weightHistFunc));
    
    if plotMode == 1
        figure(1);
        clf
        set(gcf,'Name','Amp Histogram');
        set(gcf,'Color','w');
        
        data_hist_Plot = plot(Amp_bins,targetHistogram);
        data_hist_Plot.LineWidth = 2;
        data_hist_Plot.Color = 'blue';
        data_hist_Plot.DisplayName = 'Data';
        
        xlabel('Amp','FontSize',14);
        ylabel('Frequency','FontSize',14);
        
%         [sample_description, ~] = sample_descriptionGetter();
        sample_description=[constructName ' ' saltConditions];
        title_str = ['Experimental vs Simulated Histograms' ...
            10 sample_description];
        title(title_str,'fontsize',14);
        
        lgd = legend('show');
        lgd.Location = 'northwest';
        lgd.FontSize = 12;
        drawnow();
    end
    % end
    
    %--------------------------------------------------------------------------
    % (2) Optimization Target #2: 2-point TCF (C2)           (20 uSec and on)  (PART 2: Load the target data)
    %--------------------------------------------------------------------------
    
        if FRETmode==1
        load(C2_FilePath,'time','yData','y','yoff');
        else
        load(C2_FilePath,'stitchTaus','stitchTCF');
        end
        %         load(C2_FilePath,'time','yData','y');
        fitC2Data_mode = 1;%If 0 you will fit to the histogram fit
        if fitC2Data_mode == 1
            if FRETmode==1
            C2_exp_x = time;
            C2_exp_y = yData;
            else
%             startIndex=length(find(stitchTaus>=
            startIndList=find((stitchTaus>=startTime));
            endIndList=find((stitchTaus<=cutOffTime));
            C2_exp_x=stitchTaus(startIndList(1):endIndList(end));
%             stitchTaus=stitchTaus((stitchTaus<=cutOffTime));
%             C2_exp_x = stitchTaus;
            C2_exp_y = real(stitchTCF(startIndList(1):endIndList(end)));
            
            C2_exp_x=reshape(C2_exp_x,1,length(C2_exp_x)); 
            end
        elseif fitC2Data_mode == 0
            C2_exp_x = time;
            C2_exp_y = y;
        end
                if normalizeMode == 1
            C2_exp_y = C2_exp_y./C2_exp_y(1);
                end
        c2Control = c2Exp(C2_exp_x,c2AmpFitVal,c2TauFitVal)*C2_exp_y(1);
        controlTimeZero=c2Exp(0,c2AmpFitVal,c2TauFitVal)*C2_exp_y(1);
        
        c2CorrectedSurf = C2_exp_y - c2Control; 
        yoff = abs(mean(c2CorrectedSurf(end-6:end)));%New declaration of y^offset
%         yoff=0;
        
        if  strcmp(model_name ,'linear2314_Ribbon')
        C2_exp_y = [C2_exp_y, ones(1,ribbonPoints)*yoff]; 
        C2_exp_x = [C2_exp_x, linspace(C2_exp_x(end)*1.1,endRibbon,ribbonPoints)];
        end
%          yoff = abs(C2_exp_y(end:end))/2;
%          yoff = 0;%New declaration of y^offset
      if fitC2Mode == 1
%         weightC2func = 1./(sqrt(C2_exp_x));% 2-pt. TCF weighting function
%         changed 12/13/21 to normalize
%         weightC2func=weightC2func./(sum(weightC2func)); 
        if plotMode == 1
            figure(2)
            set(gcf,'Color','w');
            hold on;
            plot(C2_exp_x,C2_exp_y,'b.','MarkerSize',10,'DisplayName','C^{(2)}(\tau) Data');
            
%             [sample_description, ~] = sample_descriptionGetter();
              sample_description=[constructName ' ' saltConditions];

            title_str = ['Experimental vs Simulated C2' ...
                10 sample_description];
            title(title_str,'fontsize',14);
            
            xlabel('Time (sec)','FontSize',14);
            ylabel('C^{(2)}(\tau)','FontSize',14);
            set(gca,'xscale','log');
            drawnow();
        end
      end
      
      
    if addControlMode == 1 && varyCtrlMode == 0
      c2Control = c2Exp(C2_exp_x,c2AmpFitVal,c2TauFitVal);%*C2_exp_y(1);
      c2Control = reshape(c2Control, size(C2_exp_y));
      c2CorrectedSurf = C2_exp_y - c2Control;
      yoff = abs(mean(c2CorrectedSurf(end-6:end)));%New declaration of y^offset
      %     yoff=0;
    else
      c2Control = zeros(size(C2_exp_x));
    end
    
    %---------------------------------------------------------------------------------------------
    % (3) Optimization Target #3: 4-point TCF (C4)    (10 uSec and on) (PART 2: Load the target data)
    %---------------------------------------------------------------------------------------------
    if fitC4Mode == 1
        %         load(C4_FilePath,'FourPtTCF_avg','tau1arrayUsec','tau3arrayUsec','tau2ValUsec');
        %  C4_tau2eq0_exp = FourPtTCF_avg;
        %         C4_tau1range = tau1arrayUsec*1e-6;
        %         C4_tau3range = tau3arrayUsec*1e-6';
        if FRETmode==1
        load(C4_FilePath,'C4','tau1arraySec','tau3arraySec','tau2ValSec');
        C4_tau2eq0_exp = C4;
        C4_tau1range = tau1arraySec;
        C4_tau3range = tau3arraySec';
        tau2 = tau2ValSec;
        tau2usec=tau2*1e6; 
        else
        load(C4_FilePath,'C4finalArray','time','tau2usec');
        time=time(find(time<=cutOffTime));
        C4_tau2eq0_exp = real(C4finalArray(1:length(time),1:length(time)));
%         C4_tau2eq0_exp = real(C4finalArray(3:length(time),3:length(time)));
        C4_tau1range = time;
        C4_tau3range = time';
%         C4_tau1range = time(3:end);
%         C4_tau3range = time(3:end)';
        tau2 = tau2usec;
        end
        
        %
        if  strcmp(model_name ,'linear2314_Ribbon')
        C4_tau2eq0_exp = [C4_tau2eq0_exp, ones(size(C4_tau2eq0_exp,1),ribbonPoints)*zoff]; 
        C4_tau2eq0_exp = [C4_tau2eq0_exp; (ones(size(C4_tau2eq0_exp,2),ribbonPoints)*zoff)'];
        
        C4_tau1range = [C4_tau1range ; linspace(C4_tau1range(end)*1.1,endRibbon,ribbonPoints)'];
        C4_tau3range = [C4_tau3range , linspace(C4_tau3range(end)*1.1,endRibbon,ribbonPoints)];
        end          

        if normalizeMode == 1
            C4_tau2eq0_exp = C4_tau2eq0_exp./C4_tau2eq0_exp(1,1);
        end
        
        c4Control = c4Exp(time,c2AmpFitVal,c2TauFitVal)*C4_tau2eq0_exp(1,1);
        
        c4CorrectedSurf = C4_tau2eq0_exp - c4Control; 
        zoff = abs(mean(mean(c4CorrectedSurf(end-6:end,end-6:end))));
%         zoff = abs(C4finalArray(end:end))/2;
%           zoff = 0;  
        
        % 4-pt. TCF weighting function
%         wC4func = 1./(sqrt(C4_tau1range)).*(1./(sqrt(C4_tau1range)));
%         wC4func=wC4func./(sum(sum(wC4func))); 
        if plotMode == 1
            figure(3)
            clf;
            set(gcf,'Name','Four-Point TCF: C^{(4)}');
            set(gcf,'Color','w');
            
            C4dataPlot = mesh(C4_tau1range,C4_tau3range,C4_tau2eq0_exp);
            C4dataPlot.DisplayName = 'C^{(4)} Data';
            
            xlabel('\tau_1 (sec)','FontSize',14);
            ylabel('\tau_3 (sec)','FontSize',14);
            ylabel('C^{(4)}','FontSize',14);
            %         title('Experimental vs Simulated C4','FontSize',14);
            
%             [sample_description, save_prefix] = sample_descriptionGetter();
            title_str = ['C^{(4)}(\tau_1, \tau_2 = ' tau2usec ', \tau_3)' ...
                10 sample_description];
            title(title_str,'fontsize',14);
            xlabel('\tau_1 (sec)','fontsize',14);
            ylabel('\tau_3 (sec)','fontsize',14);
            
            zlabel('C^{(4)}(\tau_1, \tau_2, \tau_3)','fontsize',18);
            
            %----------| Clean up the plot (non-specific --> specific) |-----------
            grid on;
            set(gca,'FontSize',12);
            axis tight;
            axis square;
            colormap jet;
            colorbar;
            set(gca,'xscale','log');
            set(gca,'yscale','log');
            %         set(gca,'zscale','log');
            %         zlim([1e-4,inf]);
            view(55.3,33.2);
            drawnow();
        end
    end
    
    if addControlMode == 1 && varyCtrlMode == 0
      c4Control = c4Exp(C4_tau1range,c2AmpFitVal,c2TauFitVal); %*C4_tau2eq0_exp(1,1);
      c4CorrectedSurf = C4_tau2eq0_exp - c4Control;
      zoff = abs(mean(mean(c4CorrectedSurf(end-6:end,end-6:end))));
  %     zoff = abs(C4finalArray(end:end))/2;
  %      zoff = 0;
    else
      c4Control = zeros(size(C4_tau2eq0_exp));
    end
    
    if verboseMode == 1
        disp('     Part1: Done Loading and (optionally) plotting the data.');
    end
    % disp('Press enter to continue');
    % pause();
    
    folderInName = [GenAlgFilesLocationPrefix filesep() 'lowestChiSquare'];
    finName =  'BestFitResults.mat';
    finPath = [folderInName filesep() finName];
    if exist(finPath,'file') == 2
        disp(['Detected the file ' finName]);
        load(finPath,'iter');
        fprintf('Num iterations = %d/%d\r',iter,maxIterations);
        if forceMoreIterationsMode == 1
            if iter >= maxIterations
                maxIterations = iter + NumberToExtend;
                fprintf('Extending the maxIterations by %d\r',NumberToExtend);
            end
        end
    else
        iter = 0;
    end
    
%     CALL THE WEIGHTS CALCULATOR HERE
% weightingFactor_Amphist = 5e-1;% Weighting for Amp hist comparison
% weightingFactor_C2 = 0.3;
% weightingFactor_C4_t0 = 9.95e4;

% [wC4func,weightC2func,weightingFactor_C2,weightingFactor_C4_t0,weightingFactor_Amphist] = weightsCalculator_v2(noiseFloorHist,noiseFloorC2,noiseFloorC4,noisePercentDiffC2,noisePercentDiffC4,C2_exp_x, C2_exp_y, C4_tau1range, C4_tau2eq0_exp, targetHistogram,c2Mag,c4Mag,histMag,decadeSegments,c4Slider);
% [wC4func,weightC2func,weightingFactor_C2,weightingFactor_C4_t0,weightingFactor_Amphist] = weightsCalculator_v3_chi2(noiseFloorHist,noiseFloorC2,noiseFloorC4,noisePercentDiffC2,noisePercentDiffC4,C2_exp_x, C2_exp_y, C4_tau1range, C4_tau2eq0_exp, targetHistogram,c2Mag,c4Mag,histMag,decadeSegments,c4Slider);
[wC4func,weightC2func,weightingFactor_C2,weightingFactor_C4_t0,weightingFactor_Amphist] = weightsCalculator_v4_chi2(noiseFloorHist,noiseFloorC2,noiseFloorC4,noisePercentDiffC2,noisePercentDiffC4,C2_exp_x, C2_exp_y, C4_tau1range, C4_tau2eq0_exp, targetHistogram,c2Mag,c4Mag,histMag,decadeSegments,c4Slider,c4WeightTime,Amp_bins)

      %disp(['Nparam = ',num2str(Nparam)])
     % disp(['Nstates = ', num2str(Nstates)])
     % disp(['len(tijs) = ', num2str(length(tijs))])
     % disp(['len(tijs) + 2 * Nstates = ', num2str(length(tijs) + 2*Nstates)])
     % disp(['This leaves ',num2str(Nparam - (length(tijs) + 2*Nstates)), ' control parameters.'])
    
    while iter < maxIterations
        saveMode = saveModeDefault; %Saves the output of the program to a folder
        
        iter = iter + 1;
        disp(['iter = ' num2str(iter) '/' num2str(maxIterations) ]);
        %%
        %//////////////////////////////////////////////////////////////////////////
        % PART 3: Make the initial generation of guesses           (Part3: Gen1)
        %//////////////////////////////////////////////////////////////////////////
        if randomizeModeMode == 1
            percentReproduce = randi([5,50],1,1);%25;%40 works well%5 works with 50% to keep.
            percentToKeep = randi([5,50],1,1);%40;%10 works well. so does 50
            %-------------------DERIVED PARAMATERS------------------
            Nreproduce = round(NmembersInitPop*percentReproduce/100);% Number of individuals who reproduce per generation.
            membersToKeep = NmembersInitPop*percentToKeep/100;
        end
        %Give each member of the population some value between the LB and UB
%         set a while loop to fill Nmembers until full
        membersInList=0;
        population = zeros(NmembersInitPop,Nparam);
        NumLoops=0;
        gen1Guesses= zeros(gen1Nmembers,Nparam);
        C2time= C2_exp_x;
        C4time= C4_tau1range;
        tic
%         while membersInList< NmembersInitPop
%           NumLoops=NumLoops+1;  
%          curMember=zeros(1,Nparam);    
        for param_idx = 1:Nparam
            %To pick a random number in the interval of LB to UB:
            % num = LB + rand*(UB - LB); %If rand = 0 then num = LB. If rand = 1, then num = UB.
%             population(:,param_idx) = boundsArray(param_idx,1) + population(:,param_idx)*(boundsArray(param_idx,2) - boundsArray(param_idx,1));
            [valuesArray,finalLogArray] = logRandUniform(boundsArray(param_idx,1),boundsArray(param_idx,2),pointsPerDecade,gen1Nmembers); 
            gen1Guesses(:,param_idx) = valuesArray; 
        end
        
        pop_chisquared_array_gen1 = zeros(gen1Nmembers,1);

        for k=1:gen1Nmembers
            
            % updated 2022/06/06: vary control amplitude and time values
            if varyCtrlMode == 1
                Nparam_ctrlAdj = Nparam - NctrlParams;
                tijs = gen1Guesses(k,1:(Nparam_ctrlAdj - 2*Nstates));
                A = gen1Guesses(k,(Nparam_ctrlAdj - 2*Nstates+1):Nparam_ctrlAdj-Nstates);
                sigmas = gen1Guesses(k,((Nparam_ctrlAdj-Nstates)+1:Nparam_ctrlAdj));
                
                ctrlParams = gen1Guesses(k,Nparam_ctrlAdj+1:Nparam);
                c2ctrlAmp = ctrlParams(1);
                c2ctrlTau = ctrlParams(2);
                c4ctrlAmp = ctrlParams(3);
                
                % define control surfaces
                c2Control = c2Exp(C2_exp_x, c2ctrlAmp, c2ctrlTau)*C2_exp_y(1);
                c2Control = reshape(c2Control, size(C2_exp_y));
                c2CorrectedSurf = C2_exp_y - c2Control;
                
                yoff = abs(mean(c2CorrectedSurf(end-6:end)));
                c4Control = c4Exp(C4_tau1range, c4ctrlAmp, c2ctrlTau)*C4_tau2eq0_exp(1,1);
                c4CorrectedSurf = C4_tau2eq0_exp - c4Control;
                zoff = abs(mean(mean(c4CorrectedSurf(end-6:end,end-6:end))));
            else
                
                tijs = gen1Guesses(k,1:(Nparam - 2*Nstates));
                A = gen1Guesses(k,(Nparam - 2*Nstates+1):Nparam-Nstates);
                sigmas = gen1Guesses(k,((Nparam-Nstates)+1:Nparam));
                
            end
           
            if hist500ms
            [K, A, rates, tijs, Nstates, Nparam, sigma_A, param_strings, boundsArray, maxMutationArray,maxMutationCountsArray, minMutationCountsArray,atEdges, boundsArray_noAdj] = model_builder_Polz_v3_500ms(tijs, A, model_name,paramAdjs);            
            elseif hist10ms
                if varyCtrlMode
%             [K, A, rates, tijs, Nstates, Nparam, param_strings, boundsArray, maxMutationArray,maxMutationCountsArray, minMutationCountsArray,atEdges, boundsArray_noAdj] = model_builder_Polz_v4_10ms(tijs, A, model_name,paramAdjs,c2ctrlAmp_0,c2ctrlTau_0,c4ctrlAmp_0);            
[K, rates, tijs, Nstates, Nparam, param_strings, boundsArray, maxMutationArray,maxMutationCountsArray, minMutationCountsArray, atEdges, boundsArray_noAdj, model_name,dbCell,rates_str_woDB,rates_idx_woDB,rates_idx_DB] = libraryLoader(loadMode, model_set, modelNum , Nstates, tijs, c2ctrlAmp, c2ctrlTau, c4ctrlAmp, K, model_name,dbCell,rates_str_woDB,rates_idx_woDB,rates_idx_DB,Nparam,param_strings, boundsArray, maxMutationArray,maxMutationCountsArray, minMutationCountsArray, atEdges, boundsArray_noAdj,model_lin,model_loop1,model_loopN,model_brch);           
                else
            [K, A, rates, tijs, Nstates, Nparam, param_strings, boundsArray, maxMutationArray,maxMutationCountsArray, minMutationCountsArray,atEdges, boundsArray_noAdj] = model_builder_Polz_v3_10ms(tijs, A, model_name,paramAdjs);            
                end
            end
            
                if varyCtrlMode
           [chisquared, chisquared_array, chisquared_unweighted, chisquared_unweighted_array,chisquared_Weighted_array] = chiSqCalc_v4_polz(K, A, C2time, C4time,addControlMode,c2Control,c4Control,sigmas,yoff,zoff,c4WeightTime);
                else
           [chisquared, chisquared_array, chisquared_unweighted, chisquared_unweighted_array,chisquared_Weighted_array] = chiSqCalc_v3_polz(K, A, C2time, C4time,addControlMode,c2Control,c4Control,sigmas);
                end
           pop_chisquared_array_gen1(k) = chisquared;
        end 
        
        [pop_chisquare_sorted,index_arr] = sort(pop_chisquared_array_gen1);%First in array have the lowest chisquare value
%         this should take the top 200 guesses from the large gen1 set and
%         propagate them through the rest of the iteration 
        population = gen1Guesses(index_arr(1:NmembersInitPop),:);
        pop_chisquare_sorted = pop_chisquare_sorted(1:NmembersInitPop);
        
         newpopulation = zeros(size(population));
%         pop_chisquared_array = zeros(NmembersInitPop,1);
        guess = zeros(maxGenerations,Nparam+1);
        genNum = 1;
        
        toc       
        
       
        tic

        
        if verboseMode == 1
            disp('     Part2 : Done simulating the first generation. Selecting the best guess to plot.');
        end
        elapsedTime = toc;
        if clockMode == 1
            disp(['     It took ' num2str(elapsedTime) ' seconds to compute the '...
                'RMS of ' num2str(NmembersInitPop) ' guesses.(' num2str(elapsedTime/NmembersInitPop) ' sec/guess)']);
        end
        %------------------------------------------------------------------------------
        % (4) Plot the chisquared value of each member of the population (PART 3: Gen1)
        %-------------------------------------------------------------------------------
        if plotMode == 1
            if plot_Gen1ChiSquaredMode == 1
                figure(4);
                
                set(gcf,'Color','w');
                set(gca,'yscale','log');
                for pop_idx = 1:length(pop_chisquared_array)
                    semilogy(pop_idx,pop_chisquared_array(pop_idx),'--gs',...
                        'LineWidth',2,...
                        'MarkerSize',10,...
                        'MarkerEdgeColor','b',...
                        'MarkerFaceColor',[0.5,0.5,0.5]);
                    hold on;
                end
                hold on;
                axis tight;
                xlim([1,inf]);
                
                xlabel('Member # (Gen1)','FontSize',14);
                ylabel('RMS Fitness (\chi^2)','FontSize',14);
                title('\chi^2 of Gen 1','FontSize',14)
                drawnow();
            end
        end
        
        %-------------------------------------------------------------------------------
        % (7) Make a histogram of the chisquared values of each member   (PART 3: Gen1)
        %--------------------------------------------------------------------------------
        if plotMode == 1
            if plot_GenerationalChiSquaredHistogramMode == 1
                figure(7)
                set(gcf,'Name','Chi-squared values of gen1 guesses');
                
                initialHist = histogram(pop_chisquared_array);
                xlabel('\chi^2 Value');
                ylabel('Frequency');
                title('Distribution of Gen1 guesses chi-squared');
            end
        end
        
        %--------------------------------------------------------------------------
        % Sort Generation 1 as a function of chisquare (PART 3: Gen1)
        %--------------------------------------------------------------------------
        % Sort from Lowest rms to the highest rms (ascending order)
        
%         [pop_chisquare_sorted,index_arr] = sort(pop_chisquared_array);%First in array have the lowest chisquare value
%         
%         population = population(index_arr,:);
        
        guess(genNum,1:Nparam) = population(1,:); % Current best guess
        
        Best_chisquared = pop_chisquare_sorted(1);
        guess(genNum,Nparam+1) = Best_chisquared; % Current best guess' fit value
        
        genNum_array = genNum;
        chisquaredVsGen_array = Best_chisquared;
        
        % Display best guess to screen
%         UPDATE: 03-16-2022 changing the indexing to reflect addition of
%         sigmas 
            if varyCtrlMode == 1
                Nparam_ctrlAdj = Nparam - NctrlParams;
                tijs = guess(genNum,1:(Nparam_ctrlAdj - 2*Nstates));
                A = guess(genNum,(Nparam_ctrlAdj - 2*Nstates+1):Nparam_ctrlAdj-Nstates);
                sigmas = guess(genNum,((Nparam_ctrlAdj-Nstates)+1:Nparam_ctrlAdj));
                
                % define control parameters from guesses
                ctrlParams = guess(genNum,Nparam_ctrlAdj+1:Nparam);
                c2ctrlAmp = ctrlParams(1);
                c2ctrlTau = ctrlParams(2);
                c4ctrlAmp = ctrlParams(3);
                
                % define control surfaces
                c2Control = c2Exp(C2_exp_x, c2ctrlAmp, c2ctrlTau);
                c2Control = reshape(c2Control, size(C2_exp_y))*C2_exp_y(1);
                c2CorrectedSurf = C2_exp_y - c2Control;
                yoff = abs(mean(c2CorrectedSurf(end-6:end)));
                
                c4Control = c4Exp(C4_tau1range, c4ctrlAmp, c2ctrlTau)*C4_tau2eq0_exp(1,1);
                c4CorrectedSurf = C4_tau2eq0_exp - c4Control;
                zoff = abs(mean(mean(c4CorrectedSurf(end-6:end,end-6:end))));
            else
                tijs = guess(genNum, 1:(Nparam - 2*Nstates));
                A = guess(genNum, (Nparam - 2*Nstates +1):(Nparam - Nstates));
                sigma_A = guess(genNum, (Nparam - Nstates)+1:end);
            end
         
         if hist500ms
            [K, A, rates, tijs, Nstates, Nparam, sigma_A, param_strings, boundsArray, maxMutationArray,maxMutationCountsArray, minMutationCountsArray,atEdges, boundsArray_noAdj] = model_builder_Polz_v3_500ms(tijs, A, model_name,paramAdjs);            
         elseif hist10ms
             if varyCtrlMode == 1
%             [K, A, rates, tijs, Nstates, Nparam, param_strings, boundsArray, maxMutationArray,maxMutationCountsArray, minMutationCountsArray,atEdges, boundsArray_noAdj] = model_builder_Polz_v4_10ms(tijs, A, model_name,paramAdjs,c2ctrlAmp_0,c2ctrlTau_0,c4ctrlAmp_0);
[K, rates, tijs, Nstates, Nparam, param_strings, boundsArray, maxMutationArray,maxMutationCountsArray, minMutationCountsArray, atEdges, boundsArray_noAdj, model_name,dbCell,rates_str_woDB,rates_idx_woDB,rates_idx_DB] = libraryLoader(loadMode, model_set, modelNum , Nstates, tijs, c2ctrlAmp, c2ctrlTau, c4ctrlAmp, K, model_name,dbCell,rates_str_woDB,rates_idx_woDB,rates_idx_DB,Nparam,param_strings, boundsArray, maxMutationArray,maxMutationCountsArray, minMutationCountsArray, atEdges, boundsArray_noAdj,model_lin,model_loop1,model_loopN,model_brch) ;        
             else            
            [K, A, rates, tijs, Nstates, Nparam, param_strings, boundsArray, maxMutationArray,maxMutationCountsArray, minMutationCountsArray,atEdges, boundsArray_noAdj] = model_builder_Polz_v3_10ms(tijs, A, model_name,paramAdjs);
             end
             end
%         [K, A, rates, tijs,Nstates, Nparam, sigma_A, param_strings, boundsArray, maxMutationArray,maxMutationCountsArray, minMutationCountsArray,atEdges, boundsArray_noAdj] = model_builder_Polz_v3(tijs, A, model_name,paramAdjs);
        
        [P, V, K, C2time] = K2P(K,C2time);
        
        chisquared = guess(genNum,Nparam+1);
        
        if verboseMode == 1
            fprintf(['Best fit from initial generation was member #%d:\n t12 = %f, t13 = %f, t21 = %f, t23 = %f, t31 = %f, t32 = %f'...
                '\n A1 = %f, A2 = %f, A3 = %f, Best_chisquared = %f\r\n'],...
                index_arr(1),t12,t13,t21,t23,t31,t32,A1,A2,A3,Best_chisquared);
        end
        
        %--------------------------------------------------------------------------
        % (1) Make a histogram which corresponds to the best guess (PART 3: Gen1)
        %--------------------------------------------------------------------------
        if plotMode == 1
            %--------------------------------------------------------------------------
            % (1) Optimization Target #1: 1-D Amp HISTOGRAM  (PART 3: Gen1)
            %--------------------------------------------------------------------------
            if fitHistMode == 1
                figure(1)
                
                hold on;
                
                if useMatrixMethodMode == 1
                    [Peq, denom_hist_sim, hist_sim, Peq_names] = histMaker_Nstate_v3normCorrect(P, A, sigmas, Amp_bins);
                end
                
                if useMatrixMethodMode == 1                   
               
                    figure(1); clf;
                    
                    set(gcf,'Color','w');
                    set(gcf,'Name','Amp Histogram');
                    lineStyle = char('g--','c--', 'm--','y--','k--','r-.','g-.','b-.','c-.','m-.','y-.','k-.'); % Define a list of colors to loop over
                    
                    Amp_bins = linspace(0,1,length(edgesAmp));
                    
                    if exist('data_hist_Plot','var') == 1
                        delete(data_hist_Plot)
                        %                                 disp('data_hist_Plot deleted')
                    end
                    % Replot data Histogram
                    data_hist_Plot = plot(Amp_bins,targetHistogram, 'DisplayName', 'Data');
                    data_hist_Plot.LineWidth = 2;
                    data_hist_Plot.Color = 'blue';
                    %                             x.DisplayName = 'Data';
                    
                    xlabel('Amp','FontSize',14);
                    ylabel('Frequency','FontSize',14);
                    
%                     [sample_description, ~] = sample_descriptionGetter();
                    title_str = ['Experimental vs Simulated Histograms' ...
                        10 sample_description];
                    title(title_str,'fontsize',14);
                    
                    hold on;
                    
                    % Clear plots from previous run
                    if exist('histPlot','var') == 1
                        delete(histPlot)
                        %                         disp('histPlot deleted')
                    end
                    
                    % Plot histogram of each state
                    for i = 1:length(Peq)
%                         histPlot = plot(Amp_bins, Peq(i) * exp(-((Amp_bins - A(i))/sigma_A(i)).^2)./denom_hist_sim,...
                            histPlot = plot(Amp_bins,(Peq(i)/(sqrt(2*pi)*sigmas(i)))*exp(-((Amp_bins - A(i))/sigmas(i)).^2)./denom_hist_sim,...
                            lineStyle(i,:),'LineWidth',1,'DisplayName',[Peq_names(i,:) ' = ' num2str(Peq(i))]);
                        %                                 lgd = legend('show');
                        % lgd.Location = 'northwest';
                        
                        hold on
                    end
                    set(gca, 'FontSize', 14);
                    if exist('histPlotTot','var') == 1
                        delete(histPlotTot)
                    end
                    hold on
                    
                    % Plot overall histogram
                    histPlotTot = plot(Amp_bins, hist_sim,...
                        lineStyle(i),'LineWidth',1,'DisplayName','Total Fit', 'Color','r');
                    lgd_tot = legend('show');
                    lgd_tot.FontSize = 14;
                end
            end
            
            %---------------------------------------------------------------------------------
            % (2) Optimization Target #2: 2-point TCF (C2)  (20 uSec and on)  (PART 3: GEN1)
            %---------------------------------------------------------------------------------
            if fitC2Mode == 1
                figure(2)
                
                if useMatrixMethodMode == 1
                    [C2_sim,C2time] = PA2C2(P,A,C2time,yoff,addControlMode,c2Control);
                end
                
                if normalizeMode == 1
                    C2_sim = C2_sim./C2_sim(1);
                end
                
                hold on;
                if exist('C2_sim_plot','var') == 1
                    delete(C2_sim_plot)
                end
                C2_sim_plot = plot(C2_exp_x,C2_sim);
                set(gca,'xscale','log');
                drawnow();
            end
            
            
            %-----------------------------------------------------------------------------
            % (3) Optimization Target #3: 4-point TCF (C4) (20 uSec and on) (PART 3: GEN1)
            %------------------------------------------------------------------------------
            if fitC4Mode == 1
                
                if useMatrixMethodMode == 1
                    C4time = C4_tau1range;
                    tau2 = 0;
%                     c2 vs c4 time fix
%                     [P, V, K, C4time] = K2P(K,C4time);
                    [C4_sim,C4time] = PAK2C4(P,A,K,C4time,tau2,zoff,addControlMode,c4Control);
                end
                C4_sim = C4_sim; %Add the y-offset^2 to C4_sim
                if normalizeMode == 1
                    C4_sim = C4_sim./C4_sim(1,1);
                end
                
                if exist('C4_plot','var') == 1
                    delete(C4_plot)
                end
                figure(3)
                hold on;
                C4_plot = surf(C4_tau1range,C4_tau1range,real(C4_sim));
                set(gca,'xscale','log');
                set(gca,'yscale','log');
                drawnow();
            end
            
            %-------------------------------------------------------------------------------------------
            % (4) Make a red circle around the best fit and begin optimization from there (PART 3: Gen1)
            %-------------------------------------------------------------------------------------------
            if plot_Gen1ChiSquaredMode == 1
                figure(4);
                plot(index_arr(1),chisquared,'--go','LineWidth',3,...
                    'MarkerSize',15,'MarkerEdgeColor','r','MarkerFaceColor','none')
                y = [chisquared chisquared];
                x = [index_arr(1) 1];
                line1 = line(x,y);
                line1.LineWidth = 2;
                line1.Color = 'r';
                line1.LineStyle = '--';
            end
            
            
            %------------------------------------------------------------------------------
            % (4) Keep track of the chisquared  for each member          (PART 3: Gen1)
            %------------------------------------------------------------------------------
            if plot_GenerationalProgressMode == 1
                figure(4);
                clf;
                plot(genNum,chisquared,'bs',...
                    'MarkerSize',10,...
                    'MarkerEdgeColor','k',...
                    'MarkerFaceColor',[0.5,0.5,0.5]);
                axis tight;
                xlim([1,inf]);
                set(gcf,'Color','w');
                
                xlabel('Generation Number','FontSize',14);
                ylabel('RMS Fitness (\chi^2)','FontSize',14);
                title('Fitness vs Generation Number','FontSize',14)
                drawnow();
            end
            %--------------------------------------------------------------------------
            % (5) Plot histograms of initial guesses for each paramater (PART 3: Gen1)
            %--------------------------------------------------------------------------
            if plot_paramater_histogram_mode == 1
                figure(5)
                clf;
                %                 param_strings = {'t12','t13','t21','t23','t31','A1','A2','A3'};
                maxCountArray = zeros(1,Nparam);
                for param_idx = 1:Nparam
                    numPlots=Nparam;
                    DimRow= floor(sqrt(numPlots));
                    DimCol= floor(sqrt(numPlots));
                    
                    if ((DimRow*DimCol)< numPlots)
                        DimRow=DimRow+1;
                        if(DimRow*DimCol)< numPlots
                            DimCol=DimCol+1;
                        end
                    end
                    subplot(DimRow,DimCol,param_idx)
                    
%                     if ismember(param_idx,1:5)%Rate histograms
%                   UPDATE 3-16-2022 change indexing for sigmas inclusion
                    if ismember(param_idx,1:(Nparam - 2*Nstates))%Rate histograms
                        nbins = 100;
                        
                        % updated 2022-03-02: create histogram from log spaced bins
                        [~, log_bins] = logRandUniform(boundsArray(param_idx,1),boundsArray(param_idx,2),50,1);
                        histogram(population(:,param_idx),'BinEdges',log_bins);
                        [N,edges] = histcounts(population(:,param_idx),'BinEdges',log_bins);
                        set(gca, 'xscale','log')
                        lowExp=floor(log10(boundsArray(param_idx,1)));
                        highExp=ceil(log10(boundsArray(param_idx,2)));
                        lowLim= 1*10^(lowExp); 
                        highLim= 1*10^(highExp); 
                        numDecades=log10(highLim/lowLim);
                        xticks=10.^linspace(lowExp, highExp, numDecades+1);
                        XTickLabels = cellstr(num2str(round(log10(xticks(:))), '10^{%d}'));
                        set(gca,'Xtick',xticks, 'XTickLabels',XTickLabels)

                  
%                         histogram(population(:,param_idx),nbins);
%                         [N,edges] = histcounts(population(:,param_idx),nbins);
                        
%                         [Nwind,edges] = histcounts(population(1:(NmembersInitPop*0.5) ,param_idx),nbins);
                        if windowIter==0
                        histCell{1,param_idx}=N; 
                        else
                        histCell{1,param_idx}=histCell{1,param_idx}+N;    
                        end
                        maxCountArray(param_idx) = max(N);
                        xlabel('Time (Sec)');
%                     elseif ismember(param_idx,[6,7,8]) == 1%FRET Historgrams
                    elseif ismember(param_idx,[(Nparam - 2*Nstates +1):Nparam-Nstates]) == 1%FRET Historgrams
                        nbins = 20;
                        histogram(population(:,param_idx),nbins);
                        [N,edges] = histcounts(population(:,param_idx),nbins);
                        maxCountArray(param_idx) = max(N);
                        xlabel('FRET');
                        xlim([0,1]);
                    end
                    %Plot a blue line for the best guess of that gneration
                    progressLine = line([guess(genNum,param_idx) guess(genNum,param_idx) ],[0  maxCountArray(param_idx)],...
                        'Color','blue','LineStyle','-','LineWidth',4);
                    
                    
                    ylabel('Frequency');
                    title(cellstr(param_strings(param_idx)));
                end
            end
        end%If plot Mode
        
        % error('Quitting before generation 2');
        if pauseBetweenGenerationMode == 1
            disp(['     Gen' num2str(genNum) ' Chi^2 = ' num2str(chisquared) ': Press Enter to proceed to the next generation']);
            pause();
        end
        
        %%
        %///////////////////////////////////////////////////////////////////////////////////////
        % PART 4: Refine the guesses by mutating certain paramaters randomly (PART 4: GenX)
        %///////////////////////////////////////////////////////////////////////////////////////
        %Bookmarking for keeping up with the number of trials
        Ntrials = 0;
        for genNum = 2:maxGenerations
            if clockMode == 1, tic; end
            
            if BrettGenX==1
            % Mix genes of top rated individuals (population is sorted by chisquared: The lowest chisquared are first)
            for n = 1:NmembersInitPop - membersToKeep
%                 newpopulation(n,:) = [population(ceil(Nreproduce*rand(1,1)),1) ,population(ceil(Nreproduce*rand(1,1)),2), population(ceil(Nreproduce*rand(1,1)),3), population(ceil(Nreproduce*rand(1,1)),4),population(ceil(Nreproduce*rand(1,1)),5),population(ceil(Nreproduce*rand(1,1)),6),population(ceil(Nreproduce*rand(1,1)),7),population(ceil(Nreproduce*rand(1,1)),8)];
                newpopulation(n,:) = population(ceil(Nreproduce*rand(1,1)),1:Nparam);
            end

%                 newpopulation(1:membersToKeep,:) = population(1:membersToKeep,1:Nparam);
% 
%                 newpopulation(membersToKeep:NmembersInitPop,:) = population(1:(NmembersInitPop-membersToKeep)+1,1:Nparam);
            % Keep top 5% best parents (These are the largest indices)
            %     newpopulation(NmembersInitPop - 1,:) = population(2,:);
            %     newpopulation(NmembersInitPop,:) = population(1,:);
            pop_idx = 0;
            while pop_idx <= membersToKeep
                %         Let the highest member#s be the best guesses from the first
                %         generation
                newpopulation(NmembersInitPop - pop_idx,:) = population(pop_idx+1,:);
                pop_idx = pop_idx + 1;
            end

            population = newpopulation;
            end 
            
             if JackGenX==1
                 normSelect = ceil(abs(random(makedist('HalfNormal','mu',0,'sigma',Nreproduce/2),(NmembersInitPop - membersToKeep),Nparam)));
                    if normSelect > NmembersInitPop % prevent choosing a value > # of members in pop
                     [row, col] = find(normSelect > NmembersInitPop);
                     normSelect(row,col) = Nreproduce;
                    end
                    %         first retain the top N members of the previous generation
                    %         (sorted by chisquared) -
                    newpopulation(1:membersToKeep,:) = population(1:membersToKeep,1:Nparam);
                    % %       now ???reproduce??? the remaining M members of the
                    %        newpopulation using only random mixtures of the parameters
                    %        from the top Nreproduce members of the previous generation.
                    for j=(membersToKeep+1):NmembersInitPop
                        mixArray=round(rand(1,Nparam)*Nreproduce+0.5);
                        for m=1:Nparam
                            %             newpopulation(j,m)=population(mixArray(m),m);
                            newpopulation(j,m)=population(normSelect(j-(membersToKeep),m),m);
                        end
                    end
                    population=newpopulation;
             end
%                 first retain the top N members of the previous generation
%                 (sorted by chisquared) - 
%               newpopulation(1:membersToKeep,:) = population(1:membersToKeep,1:Nparam);
%             
% % %             now "reproduce" the remaining M members of the
% %               newpopulation using only random mixtures of the parameters
% %               from the top Nreproduce members of the previous generation. 
%                 for j=(membersToKeep+1):NmembersInitPop
%                     mixArray=round(rand(1,Nparam)*Nreproduce+0.5);
%                     for m=1:Nparam
%                     newpopulation(j,m)=population(mixArray(m),m); 
%                     end
%                 end
%             
% 
%                 population=newpopulation; 
%              end
            
            % Mutate some individuals
            for q = 1:Nmutations
                %Let each paramater be mutated by a different amount across entire
                %population
                rndArr = rand(Nparam,1);
                for param_idx = 1:Nparam
                    % Move each value of the population somewhere inbetween (uniformly
                    % distributed) up or down by its  maximum allowable mutatin value
                    if BrettGenX==1
                    popMemberID = ceil((NmembersInitPop - 2)*rndArr(param_idx));
                    end
%                   make it so that only the 3rd through Nth member can be mutated (retaining
%                   the top 2 current guesses as they are) 
                    if JackGenX==1
                    popMemberID = ceil((NmembersInitPop - 2)*rndArr(param_idx)+2);
                    end
%                     UPDATE 3-14-2022 : make the upper bound of the
%                     sampled mutation range 10X maxMut so that the lower
%                     and upper bound are never equal (zero index case)
                    [valuesArray, finalLogArray] = logRandUniform(1e-8,maxMutationArray(param_idx)*10,50,1); 
                    newValue = population(popMemberID,param_idx) + (2*(rand(1,1) - .5))*valuesArray(1);
%                     newValue = population(popMemberID,param_idx) + (2*(rand(1,1) - .5))*maxMutationArray(param_idx);
                    population(popMemberID,param_idx) = newValue;
                    if guessUpdateMode == 1
                        fprintf(' Mutating member #%d''s %s param to %f\r',popMemberID,cellstr(param_strings(param_idx)),newValue)
                    end
                    % If the population has exceeded the maximum allowable amount, set it to
                    % that amount. If it is below the minumum amont, set it to the minimum
                    %************ How often does this happen? This is definitely not what
                    %we want to occur*********
                    if population(popMemberID,param_idx) > boundsArray(param_idx,2)
                        population(popMemberID,param_idx) = boundsArray(param_idx,2);
                        if guessUpdateMode == 1
                            fprintf('     Setting param %d (%s) to the MAX value: %f.\n',param_idx,char(param_strings(param_idx)),boundsArray(param_idx,2));
                        end
                    elseif population(popMemberID,param_idx) < boundsArray(param_idx,1)
                        population(popMemberID,param_idx) = boundsArray(param_idx,1);
                        if guessUpdateMode == 1
                            fprintf('     Setting param %d (%s) to the MIN value: %f.\n',param_idx,char(param_strings(param_idx)),boundsArray(param_idx,1));
                            maxMutationCountsArray(param_idx) = maxMutationCountsArray(param_idx) + 1;
                        end
                    else
                        if guessUpdateMode == 1
                            fprintf('     Setting param %d (%s) to the NEW value: %f.\n',param_idx,char(param_strings(param_idx)),newValue);
                            minMutationCountsArray(param_idx) = minMutationCountsArray(param_idx) + 1;
                        end
                    end
                end
                
            end%End of mutations
            
            %--------------------------------------------------------------------------------
            %After Mutations, Recalculate the chisquared values of each member (PART 4: GenX)
            %--------------------------------------------------------------------------------
            pop_chisquared_array = zeros(NmembersInitPop,1);

            for pop_idx = 1:NmembersInitPop
                
               if varyCtrlMode == 1
                Nparam_ctrlAdj = Nparam - NctrlParams;
                tijs = population(pop_idx,1:(Nparam_ctrlAdj - 2*Nstates));
                A = population(pop_idx,(Nparam_ctrlAdj - 2*Nstates+1):Nparam_ctrlAdj-Nstates);
                sigmas = population(pop_idx,((Nparam_ctrlAdj-Nstates)+1:Nparam_ctrlAdj));
                
                % define control parameters from guesses
                ctrlParams = population(pop_idx,Nparam_ctrlAdj+1:Nparam);
                c2ctrlAmp = ctrlParams(1);
                c2ctrlTau = ctrlParams(2);
                c4ctrlAmp = ctrlParams(3);
                
                % define control surfaces
                c2Control = c2Exp(C2_exp_x, c2ctrlAmp, c2ctrlTau)*C2_exp_y(1);
                c2Control = reshape(c2Control, size(C2_exp_y));
                c2CorrectedSurf = C2_exp_y - c2Control;
                yoff = abs(mean(c2CorrectedSurf(end-6:end)));
                
                c4Control = c4Exp(C4_tau1range, c4ctrlAmp, c2ctrlTau)*C4_tau2eq0_exp(1,1);
                c4CorrectedSurf = C4_tau2eq0_exp - c4Control;
                zoff = abs(mean(mean(c4CorrectedSurf(end-6:end,end-6:end))));
              else
                tijs = population(pop_idx, 1:(Nparam - 2*Nstates));
                A = population(pop_idx, (Nparam - 2*Nstates +1):(Nparam - Nstates));
                sigma_A = population(pop_idx, (Nparam - Nstates)+1:end);
              end
      
                
                if hist500ms
                [K, A, rates, tijs, Nstates, Nparam, sigma_A, param_strings, boundsArray, maxMutationArray,maxMutationCountsArray, minMutationCountsArray,atEdges, boundsArray_noAdj] = model_builder_Polz_v3_500ms(tijs, A, model_name,paramAdjs);            
                elseif hist10ms
                    if varyCtrlMode == 1
%                 [K, A, rates, tijs, Nstates, Nparam, param_strings, boundsArray, maxMutationArray,maxMutationCountsArray, minMutationCountsArray,atEdges, boundsArray_noAdj] = model_builder_Polz_v4_10ms(tijs, A, model_name,paramAdjs,c2ctrlAmp_0,c2ctrlTau_0,c4ctrlAmp_0);            
[K, rates, tijs, Nstates, Nparam, param_strings, boundsArray, maxMutationArray,maxMutationCountsArray, minMutationCountsArray, atEdges, boundsArray_noAdj, model_name,dbCell,rates_str_woDB,rates_idx_woDB,rates_idx_DB] = libraryLoader(loadMode, model_set, modelNum , Nstates, tijs, c2ctrlAmp, c2ctrlTau, c4ctrlAmp, K, model_name,dbCell,rates_str_woDB,rates_idx_woDB,rates_idx_DB,Nparam,param_strings, boundsArray, maxMutationArray,maxMutationCountsArray, minMutationCountsArray, atEdges, boundsArray_noAdj,model_lin,model_loop1,model_loopN,model_brch) ;          
                    else
                [K, A, rates, tijs, Nstates, Nparam, param_strings, boundsArray, maxMutationArray,maxMutationCountsArray, minMutationCountsArray,atEdges, boundsArray_noAdj] = model_builder_Polz_v3_10ms(tijs, A, model_name,paramAdjs);            
                    end
                    end
%                 [K, A, rates, tijs,Nstates, Nparam, sigma_A, param_strings, boundsArray, maxMutationArray,maxMutationCountsArray, minMutationCountsArray,atEdges, boundsArray_noAdj] = model_builder_Polz_v3(tijs, A, model_name,paramAdjs);                
                
                if useMatrixMethodMode == 1
                    C2time = C2_exp_x;
                    C4time = C4_tau1range;
                    Amp_bins = linspace(0,1,length(edgesAmp));
                    if varyCtrlMode == 1
                   [chisquared, chisquared_array, chisquared_unweighted, chisquared_unweighted_array,chisquared_Weighted_array] = chiSqCalc_v4_polz(K, A, C2time, C4time,addControlMode,c2Control,c4Control,sigmas,yoff,zoff,c4WeightTime);  
                    else
                    [chisquared, chisquared_array, chisquared_unweighted, chisquared_unweighted_array,chisquared_Weighted_array] = chiSqCalc_v3_polz(K, A, C2time, C4time,addControlMode,c2Control,c4Control,sigmas);
                    end
                    end
                %Add the chisquared value to its corresponding slot in the
                %population
                pop_chisquared_array(pop_idx) = chisquared;
            end
            
            % Evaluate Generation
            % Sort from highest rated to lowest rated individuals
            [pop_chisquare_sorted,index_arr] = sort(pop_chisquared_array);
            population = population(index_arr,:);
            
            %Add to the best guess array the best guess of the population
%             guess(genNum,1:8) = population(1,:); % Current best guess
%             guess(genNum,9) = pop_chisquare_sorted(1); % Current best guess' fit value
            guess(genNum,1:Nparam) = population(1,:); % Current best guess
            guess(genNum,Nparam+1) = pop_chisquare_sorted(1); % Current best guess' fit value
            
            %population(1,:) % Display best guess to screen
            %     popfit(1);
            if guessUpdateMode == 1
                disp(['     Gen ' num2str(genNum) ': Best fitness guess so far: ' num2str(pop_chisquare_sorted(1))]);
            end
            
            %----------------------------------------------------------------------
            % (7) Plot the distribution of guesses (PART 4: genX)
            %----------------------------------------------------------------------
            if plot_GenerationalChiSquaredHistogramMode == 1
                figure(7)
                set(gcf,'Name','Chi-squared values of initial guesses');
                set(gcf,'Color','w');
                %         switch computer_terminal_str
                %             case 'computer_baimi_mode'
                %                 set(gcf,'Position',[2535 573 560 420]);%Work computer
                %             case 'computer_bisraels_mode'
                %                 %                 set(gcf,'Position',[1121 41 560 420]);%Macbook computer
                %                 disp('Set an option for placement of this figure');
                %         end
                histogram(pop_chisquared_array);
                xlabel('\chi^2 Value');
                ylabel('Frequency');
                title(['Gen' num2str(genNum) ' Distribution of guesses chi-squared']);
            end
            
            
            %--------------------------------------------------------------------------
            % Update the plots with the current best guesses (PART 4: genX)
            %--------------------------------------------------------------------------
            if plot_GenerationalProgressMode == 1
                
%                 now that the population is chi-square sorted, take the
%                 first entry
               if varyCtrlMode == 1
                Nparam_ctrlAdj = Nparam - NctrlParams;
                tijs = population(1,1:(Nparam_ctrlAdj - 2*Nstates));
                A = population(1,(Nparam_ctrlAdj - 2*Nstates+1):Nparam_ctrlAdj-Nstates);
                sigmas = population(1,((Nparam_ctrlAdj-Nstates)+1:Nparam_ctrlAdj));

                % define control parameters from guesses
                ctrlParams = population(1,Nparam_ctrlAdj+1:Nparam);
                c2ctrlAmp = ctrlParams(1);
                c2ctrlTau = ctrlParams(2);
                c4ctrlAmp = ctrlParams(3);

                % define control surfaces
                c2Control = c2Exp(C2_exp_x, c2ctrlAmp, c2ctrlTau);
                c2Control = reshape(c2Control, size(C2_exp_y))*C2_exp_y(1);
                c2CorrectedSurf = C2_exp_y - c2Control;
                yoff = abs(mean(c2CorrectedSurf(end-6:end)));

                c4Control = c4Exp(C4_tau1range, c4ctrlAmp, c2ctrlTau)*C4_tau2eq0_exp(1,1);
                c4CorrectedSurf = C4_tau2eq0_exp - c4Control;
                zoff = abs(mean(mean(c4CorrectedSurf(end-6:end,end-6:end))));
              else
                tijs = population(1,1:(Nparam - 2*Nstates));
                A = population(1,(Nparam - 2*Nstates+1):Nparam-Nstates);
                sigma_A = population(1, (Nparam-Nstates+1):end);
              end
                  
                  if hist500ms
                 [K, A, rates, tijs, Nstates, Nparam, sigma_A, param_strings, boundsArray, maxMutationArray,maxMutationCountsArray, minMutationCountsArray,atEdges, boundsArray_noAdj] = model_builder_Polz_v3_500ms(tijs, A, model_name,paramAdjs);            
                  elseif hist10ms
                      if varyCtrlMode == 1
%                  [K, A, rates, tijs, Nstates, Nparam, param_strings, boundsArray, maxMutationArray,maxMutationCountsArray, minMutationCountsArray,atEdges, boundsArray_noAdj] = model_builder_Polz_v4_10ms(tijs, A, model_name,paramAdjs,c2ctrlAmp_0,c2ctrlTau_0,c4ctrlAmp_0);              
[K, rates, tijs, Nstates, Nparam, param_strings, boundsArray, maxMutationArray,maxMutationCountsArray, minMutationCountsArray, atEdges, boundsArray_noAdj, model_name,dbCell,rates_str_woDB,rates_idx_woDB,rates_idx_DB] = libraryLoader(loadMode, model_set, modelNum , Nstates, tijs, c2ctrlAmp, c2ctrlTau, c4ctrlAmp, K, model_name,dbCell,rates_str_woDB,rates_idx_woDB,rates_idx_DB,Nparam,param_strings, boundsArray, maxMutationArray,maxMutationCountsArray, minMutationCountsArray, atEdges, boundsArray_noAdj,model_lin,model_loop1,model_loopN,model_brch) ;           
                      else
                 [K, A, rates, tijs, Nstates, Nparam, param_strings, boundsArray, maxMutationArray,maxMutationCountsArray, minMutationCountsArray,atEdges, boundsArray_noAdj] = model_builder_Polz_v3_10ms(tijs, A, model_name,paramAdjs);            
                      end
                      end
%                 [K, A, rates, tijs, Nstates, Nparam, sigma_A, param_strings, boundsArray, maxMutationArray,maxMutationCountsArray, minMutationCountsArray,atEdges, boundsArray_noAdj] = model_builder_Polz_v3(tijs, A, model_name,paramAdjs);
                
%                 time = C2_exp_x;
                [P, V, K, C2time] = K2P(K,C2time);
                
                
                chisquared = guess(genNum,Nparam+1);
                
                genNum_array = [genNum_array; genNum];
                chisquaredVsGen_array = [chisquaredVsGen_array; chisquared];
                
                %--------------------------------------------------------------------------
                % (1) Optimization Target #1: 1-D Amp HISTOGRAM (PART 4: genX)
                %--------------------------------------------------------------------------
                if fitHistMode == 1                
                    
                   if useMatrixMethodMode == 1
                        [Peq, denom_hist_sim, hist_sim, Peq_names] = histMaker_Nstate_v3normCorrect(P, A, sigmas, Amp_bins);
                    end
                    
                    if useMatrixMethodMode == 1
                        
                        figure(1); clf;
                        set(gcf,'Color','w');
                        set(gcf,'Name','Amp Histogram');
                        lineStyle = char('g--','c--', 'm--','y--','k--','r-.','g-.','b-.','c-.','m-.','y-.','k-.'); % Define a list of colors to loop over
                        Amp_bins = linspace(0,1,length(edgesAmp));
                        
                        if exist('data_hist_Plot','var') == 1
                            delete(data_hist_Plot)
                            %                                 disp('data_hist_Plot deleted')
                        end
                        % Replot data Histogram
                        data_hist_Plot = plot(Amp_bins,targetHistogram, 'DisplayName', 'Data');
                        data_hist_Plot.LineWidth = 2;
                        data_hist_Plot.Color = 'blue';
                        %                             x.DisplayName = 'Data';
                        xlabel('Amp','FontSize',14);
                        ylabel('Frequency','FontSize',14);
                        
%                         [sample_description, ~] = sample_descriptionGetter();
                        title_str = ['Experimental vs Simulated Histograms' ...
                            10 sample_description];
                        title(title_str,'fontsize',14);
                        hold on;
                        % Clear plots from previous run
                        if exist('histPlot','var') == 1
                            delete(histPlot)
                            %                             disp('histPlot deleted')
                        end
                        % Plot histogram of each state
                        for i = 1:length(Peq)
%                             histPlot = plot(Amp_bins, Peq(i) * exp(-((Amp_bins - A(i))/sigma_A(i)).^2)./denom_hist_sim,...
                                histPlot = plot(Amp_bins,(Peq(i)/(sqrt(2*pi)*sigmas(i)))*exp(-((Amp_bins - A(i))/sigmas(i)).^2)./denom_hist_sim,...                                
                                lineStyle(i,:),'LineWidth',1,'DisplayName',[Peq_names(i,:) ' = ' num2str(Peq(i))]);
                            %                                 lgd = legend('show');
                            % lgd.Location = 'northwest';
                            hold on
                        end
                        set(gca, 'FontSize', 14);
                        if exist('histPlotTot','var') == 1
                            delete(histPlotTot)
                        end
                        hold on
                        % Plot overall histogram
                        histPlotTot = plot(Amp_bins, hist_sim,...
                            lineStyle(i),'LineWidth',1,'DisplayName','Total Fit', 'Color','r');
                        lgd_tot = legend('show');
                        lgd_tot.FontSize = 14;
                    end
                end
                %--------------------------------------------------------------------------
                % (2) Optimization Target #2: 2-point TCF (C2) (20 uSec and on) (PART 4: genX)
                %--------------------------------------------------------------------------
                if fitC2Mode == 1
                    
                    if useMatrixMethodMode == 1
                        C2time = C2_exp_x;
                        [C2_sim,C2time] = PA2C2(P,A,C2time,yoff,addControlMode,c2Control);
                    end
                    
                    if normalizeMode == 1
                        C2_sim = C2_sim./C2_sim(1);
                    end
                    %Plot the correlation function
                    figure(2)
                    if exist('C2_sim_plot','var') == 1
                        delete(C2_sim_plot)
                    end
                    C2_sim_plot = plot(C2_exp_x,C2_sim);
                    set(gca,'xscale','log');
                    drawnow();
                end
                
                %--------------------------------------------------------------------------
                % (3) Optimization Target #3: 4-point TCF (C4) (10 uSec and on) (PART 4: genX)
                %--------------------------------------------------------------------------
                if fitC4Mode == 1
                    %Plot the correlation function
                    figure(3)
                    if exist('C4_plot','var') == 1
                        delete(C4_plot)
                    end
                    
                    if useMatrixMethodMode == 1
%                         C4time = C4_tau1range;
                        tau2 = 0;
%                         time fix
%                         [P, V, K, C4time] = K2P(K,C4time);
                        [C4_sim,C4time] = PAK2C4(P,A,K,C4time,tau2,zoff,addControlMode,c4Control);
                    end
                    C4_sim = C4_sim; %Add the y-offset^2 to C4_sim
                    if normalizeMode == 1
                        C4_sim = C4_sim./C4_sim(1,1);
                    end
                    
                    if exist('C4_plot','var') == 1
                        delete(C4_plot)
                    end
                    figure(3)
                    hold on;
                    C4_plot = surf(C4_tau1range,C4_tau1range,real(C4_sim));
                    set(gca,'xscale','log');
                    set(gca,'yscale','log');
                    drawnow();
                end
                
                %--------------------------------------------------------------------------
                % (4) Plotting the progress as a function of Generation
                %--------------------------------------------------------------------------
                figure(4);
                
                hold on;
                plot(genNum_array,chisquaredVsGen_array,'o',...
                    'MarkerSize',10,...
                    'MarkerEdgeColor','b',...
                    'MarkerFaceColor',[0.5,0.5,0.5]);
                drawnow();
                %------------------------------------------------------------------------------------------
                % (5) Update paramter-histograms as a function of the current best guesses  (PART 4: genX)
                %------------------------------------------------------------------------------------------
                if plot_paramater_histogram_mode == 1
                    figure(5)
                    set(gcf,'Name',['Distribution of paramaters: Gen ' num2str(genNum)]);
                    set(gcf,'Color','w');
                    
%                     param_strings = {'t12','t13','t21','t23','t31','A1','A2','A3'};
                    maxCountArray = zeros(1,Nparam);
                    for param_idx = 1:Nparam
                        subplot(DimRow,DimCol,param_idx)
                        
                        if ismember(param_idx,1:Nparam-2*Nstates)%Rate histograms
%                             nbins = 100;
%                             histogram(population(:,param_idx),nbins);
%                             [N,edges] = histcounts(population(:,param_idx),nbins);
%                             [Nwind,edges] = histcounts(population(1:(NmembersInitPop*0.2) ,param_idx),nbins);

                            % updated 2022-03-02: create histogram from log spaced bins
                            [~, log_bins] = logRandUniform(boundsArray(param_idx,1),boundsArray(param_idx,2),50,1);
                            histogram(population(:,param_idx),'BinEdges',log_bins);
                            [Nwind,edges] = histcounts(population(:,param_idx),'BinEdges',log_bins);
                            set(gca, 'xscale','log')
                            lowExp=floor(log10(boundsArray(param_idx,1)));
                            highExp=ceil(log10(boundsArray(param_idx,2)));
                            lowLim= 1*10^(lowExp); 
                            highLim= 1*10^(highExp); 
                            numDecades=log10(highLim/lowLim);
                            xticks=10.^linspace(lowExp, highExp, numDecades+1);
                            XTickLabels = cellstr(num2str(round(log10(xticks(:))), '10^{%d}'));
                            set(gca,'Xtick',xticks, 'XTickLabels',XTickLabels)

                            
                            histCell{1,param_idx}=histCell{1,param_idx}+Nwind;
                            maxCountArray(param_idx) = max(N);
                            xlabel('Time (Sec)');
                        elseif ismember(param_idx,(Nparam-2*Nstates+1):Nparam-Nstates) == 1%Amp Historgrams
                            nbins = 20;
                            histogram(population(:,param_idx),nbins);
                            [N,edges] = histcounts(population(:,param_idx),nbins);
                            maxCountArray(param_idx) = max(N);
                            xlabel('Amp');
                            xlim([0,1]);
                        end
                        if exist('progressLine','var') == 1
                            delete(progressLine)
                        end
                        %Plot a blue line for the best guess of that gneration
                        progressLine = line([guess(genNum,param_idx) guess(genNum,param_idx) ],[0  maxCountArray(param_idx)],...
                            'Color','blue','LineStyle','-','LineWidth',4);
                        
                        ylabel('Frequency');
                        title(cellstr(param_strings(param_idx)));
                    end
                    histChart=gcf;
           
%%                   plot the upper 20% of guesses 
          figure(8)
%           clf;
          set(gcf,'Name',['Upper 20% Distribution of paramaters: Iter ' num2str(iter) ', windowIter ' num2str(windowIter) ', Gen ' num2str(genNum) ]);
          set(gcf,'Color','w');
%           maxCountArray = zeros(1,Nparam);

        % updated 2022-03-02: only create tij plots for top 20%
          for param_idx = 1:Nparam-2*Nstates
              
              numPlots_top20=Nparam-2*Nstates;
              DimRow_top20= floor(sqrt(numPlots_top20));
              DimCol_top20= floor(sqrt(numPlots_top20));
              
              if ((DimRow_top20*DimCol_top20)< numPlots_top20)
                  DimRow_top20=DimRow_top20+1;
                  if(DimRow_top20*DimCol_top20)< numPlots_top20
                      DimCol_top20=DimCol_top20+1;
                  end
              end
              subplot(DimRow_top20,DimCol_top20,param_idx)
              
              %             subplot(DimRow,DimCol,param_idx)
%               if ismember(param_idx,1:Nparam-Nstates)%Rate histograms
                  %               [~, HistLogArray]=logRandUniform(boundsArray(param_idx,1),boundsArray(param_idx,2),pointsPerDecade,100);
                  nbins = 100;
                  %               h = histogram(histCell{param_idx},nbins);
                  %               nbins = 50;
                  %               histcell_edges = linspace(boundsArray(param_idx,1),boundsArray(param_idx,2),nbins+1);
                  %               h = histogram('BinCounts',histCell{param_idx},'BinEdges',histcell_edges);%(param_idx,:));
                  %               [~, finalLogArray_hist]=logRandUniform(boundsArray(param_idx,1),boundsArray(param_idx,2),pointsPerDecade,1);%length(histCell{param_idx})+1);
                  
                  % updated 2022-03-02: plot histcell array using the logRandUniform bins
                  lowerBound_hist = boundsArray(param_idx,1);
                  upperBound_hist = boundsArray(param_idx,2);
                  [~, bins_hist] = logRandUniform(lowerBound_hist, upperBound_hist,50,1);
                  h = histogram('BinCounts',histCell{param_idx},'BinEdges', bins_hist);
                  set(gca, 'xscale','log')
                  lowExp=floor(log10(boundsArray(param_idx,1)));
                  highExp=ceil(log10(boundsArray(param_idx,2))); 
                  lowLim= 1*10^(lowExp); 
                    highLim= 1*10^(highExp); 
                    numDecades=log10(highLim/lowLim);
                  xticks=10.^linspace(lowExp, highExp, numDecades+1);
                  XTickLabels = cellstr(num2str(round(log10(xticks(:))), '10^{%d}'));
                  set(gca,'Xtick',xticks, 'XTickLabels',XTickLabels)

                  
                  %               histogram(population(1:(NmembersInitPop*0.2),param_idx),HistLogArray);
                  %               h = histogram(histCell{param_idx},'BinEdges',HistLogArray);%(param_idx,:));
                  %               [N,edges] = histcounts(population(:,param_idx),nbins);
                  %               maxCountArray(param_idx) = max(N);
                  xlabel('Time (Sec)');
                  
                  maxCountArray20(param_idx) = max(h.Values);% max(histCell{param_idx});
                  
                  hold on
                  if exist('progressLine','var') == 1
                      delete(progressLine)
                  end
                  line([guess(genNum,param_idx) guess(genNum,param_idx) ],[0 maxCountArray20(param_idx)],...
                      'Color','red','LineStyle',':','LineWidth',2);
                  set(gca, 'xscale','log')
                  hold off
%               end
              
              ylabel('Frequency');
              title(cellstr(param_strings(param_idx)));
              hold off
          end
                end
                %--------------------------------------------------------------------------
                % (6) Display a graphic of the network                  (PART 4: genX)
                %--------------------------------------------------------------------------
%                 if plot_model_GenerationalUpdatesMode == 1
%                     figure(6);
%                     clf
%                     
%                     set(gcf,'Name',['Gen' num2str(genNum) ': Model Results Gen Alg 3state123 cyclical']);
%                     
%                     xlim([0 1.1]);
%                     ylim([0 1]);
%                     
%                     %Designate spots for the states
%                     state1_loc = [0.5 1];
%                     state2_loc = [1 0];
%                     state3_loc = [0 0];
%                     
%                     hold on;
%                     %Plot the gen number
%                     text(state1_loc(1) + (state2_loc(1) - state1_loc(1))/2,state1_loc(2),['Gen# = ' num2str(genNum)],'FontSize',14);
%                     
%                     %Plot the state symbols
%                     text(state1_loc(1)-0.04,state1_loc(2)+0.05,'1','FontSize',24)
%                     text(state2_loc(1),state2_loc(2),'2','FontSize',24);
%                     text(state3_loc(1)-0.05,state3_loc(2),'3','FontSize',24);
%                     
%                     %Plot the Amp values
%                     text(state1_loc(1),state1_loc(2) + 0.05,['=' num2str(A1,'%.2f')],'FontSize',16)
%                     text(state2_loc(1)+ 0.05,state2_loc(2),['=' num2str(A2,'%.2f')],'FontSize',16);
%                     text(state3_loc(1),state3_loc(2),['=' num2str(A3,'%.2f')],'FontSize',16);
%                     
%                     %Plot Lines Between states
%                     line([state1_loc(1) state2_loc(1)],[state1_loc(2) state2_loc(2)],'Color','k','LineStyle','-');
%                     line([state2_loc(1) state3_loc(1)],[state2_loc(2) state3_loc(2)],'Color','k','LineStyle','-');
%                     line([state3_loc(1) state1_loc(1)],[state3_loc(2) state1_loc(2)],'Color','k','LineStyle','-');
%                     
%                     %Plot the inverse of the rates
%                     t12_loc = [state1_loc(1)+(state2_loc(1)- state1_loc(1))/2,state1_loc(2)+(state2_loc(2)- state1_loc(2))/2];
%                     t12_msg = ['t_{12} = ' num2str(t12) 10 't_{21} = ' num2str(t21)];
%                     text(t12_loc(1),t12_loc(2),t12_msg,'FontSize',14);
%                     
%                     t23_loc = [state3_loc(1)+(state2_loc(1)- state3_loc(1))/2,state3_loc(2)+(state2_loc(2)- state3_loc(2))/2];
%                     t23_msg = ['t_{23} = ' num2str(t23) 10 't_{32} = ' num2str(t32)];
%                     text(t23_loc(1),t23_loc(2),t23_msg,'FontSize',14);
%                     
%                     t31_loc = [state1_loc(1)+(state3_loc(1)- state1_loc(1))/2,state1_loc(2)+(state3_loc(2)- state1_loc(2))/2];
%                     t31_msg = ['t_{31} = ' num2str(t31) 10 't_{13} = ' num2str(t13)];
%                     text(t31_loc(1),t31_loc(2),t31_msg,'FontSize',14);
%                     
%                     %Plot the chi-squared value
%                     text(state3_loc(1),state1_loc(2),['\chi^2 = ' num2str(chisquared)],'FontSize',14);
%                     
%                     %Plot a title with the information in it
% %                     [sample_description, save_prefix] = sample_descriptionGetter();
%                     title_str = ['3state Cyclical: ' sample_description];
%                     text(state3_loc(1),state1_loc(2)+0.15,title_str,'fontsize',16);
%                     ylim([0,1.1]);
%                     
%                     axis off;
%                 end
%             end
            %----------------------------------------------------------------------
            % Determine if the guesses have converged (PART 4: genX)
            %----------------------------------------------------------------------
            if (guess(genNum-1,Nparam+1) - guess(genNum,Nparam+1))/guess(genNum,Nparam+1) <= threshold
                Ntrials = Ntrials + 1;
            else
                Ntrials = 0;
            end
            if Ntrials == maxRepeats
                guess = guess(1:genNum,:);
                if verboseMode == 1
                    disp(['The fitness has no longer improved after ' num2str(maxRepeats) ' repeats. Ending search.']);
                end
                break;
            end
            if clockMode == 1
                elapsedTime = toc;
                disp(['The amount of time per generation is: ' num2str(elapsedTime) ' seconds']);
            end
            if pauseBetweenGenerationMode == 1
                disp(['     Gen' num2str(genNum) ' Chi^2 = ' num2str(chisquared) ': Press Enter to proceed to the next generation']);
                pause();
            end
            end
        end
        %%
        %//////////////////////////////////////////////////////////////////////////
        %--------------------------------------------------------------------------
        % PART 5: Display the Final values to the user
        %--------------------------------------------------------------------------
        %//////////////////////////////////////////////////////////////////////////
        
% %         tijs = population(pop_idx,1:(Nparam - Nstates));
% %         A = population(pop_idx,(Nparam - Nstates+1):Nparam);
%             tijs = population(1,1:(Nparam - 2*Nstates));
%             A = population(1,(Nparam - 2*Nstates+1):Nparam-Nstates);
%             sigmas = population(1,((Nparam-Nstates)+1:Nparam));
% %         tijs = population(1,1:(Nparam - Nstates));
% %         A = population(1,(Nparam - Nstates+1):Nparam);
%          if hist500ms
%             [K, A, rates, tijs, Nstates, Nparam, sigma_A, param_strings, boundsArray, maxMutationArray,maxMutationCountsArray, minMutationCountsArray,atEdges, boundsArray_noAdj] = model_builder_Polz_v3_500ms(tijs, A, model_name,paramAdjs);            
%           elseif hist10ms
%             [K, A, rates, tijs, Nstates, Nparam, param_strings, boundsArray, maxMutationArray,maxMutationCountsArray, minMutationCountsArray,atEdges, boundsArray_noAdj] = model_builder_Polz_v3_10ms(tijs, A, model_name,paramAdjs);            
%           end
% %         [K, A, rates, tijs,Nstates, Nparam, sigma_A, param_strings, boundsArray, maxMutationArray,maxMutationCountsArray, minMutationCountsArray,atEdges, boundsArray_noAdj] = model_builder_Polz_v3(tijs, A, model_name,paramAdjs);
%         
%         [P, V, K, C2time] = K2P(K,C2time);
%         
%         if useMatrixMethodMode == 1
%             [Peq, denom_hist_sim, hist_sim, Peq_names] = histMaker_Nstate_v3normCorrect(P, A, sigmas, Amp_bins);
%             p1_eq = Peq(1);
%             p2_eq = Peq(2);
%             p3_eq = Peq(3);
%         end
%         
%         chisquared = guess(genNum,Nparam+1);
        
        if verboseMode == 1
            fprintf(['     Final values at the end of the GenAlg:\n t12 = %f, t13 = %f, t21 = %f, t23 = %f, t31 = %f, t32 = %f\n'...
                'A1 = %f, A2 = %f, A3 = %f, chisquared = %f\n'...
                'p1_eq = %f, p2_eq = %f, p3_eq = %f\r\n'],...
                t12,t13,t21,t23,t31,t32,A1,A2,A3,chisquared,p1_eq,p2_eq,p3_eq);
        end
        %%
        %--------------------------------------------------------------------------
        % Save the data (PART 5: Final Values)
        %--------------------------------------------------------------------------
        
        %--------------------------------------------------------------------------
        % 5.1: Figure out where to save the data
        %--------------------------------------------------------------------------
        if sum(abs(imag(Peq))>0)~=0 || sum(Peq<0)~=0
        
        else
        if saveMode == 1
            clear('outputFolderName')
            %Make output folder if it doesnt exist
            lowestChiSquareFolderName = [GenAlgFilesLocationPrefix filesep() 'lowestChiSquare'];
            secondLowestChiSquareFolderName = [GenAlgFilesLocationPrefix filesep() 'secondLowestChiSquare'];
            thirdLowestChiSquareFolderName = [GenAlgFilesLocationPrefix filesep() 'thirdLowestChiSquare'];
            if exist(lowestChiSquareFolderName,'dir') ~= 7
                mkdir(lowestChiSquareFolderName);
                disp(['Making a folder called "' lowestChiSquareFolderName '" to hold the output']);
            end
            
            %If the folder exist and there are files in it, check the new fit
            %against them.
            
            bestfitFileName = [lowestChiSquareFolderName filesep() 'BestFitResults' '.mat'];
            if exist(bestfitFileName,'file') == 2
                
                % % %                 iter = iter + 1;
                save(bestfitFileName,'iter','-append');%***one place we update
                disp([save_prefix ': Updated the best fit with another iteration: #' num2str(iter)]);
                
                %Update the fitresults with the iteration number no matter what
                %Save the data in a higher chisquare folder
                outputFolderName = [GenAlgFilesLocationPrefix filesep() 'higherChiSquares'];
                if exist(outputFolderName,'dir') ~= 7
                    mkdir(outputFolderName);
                    disp(['Making a folder called "' outputFolderName '" to hold the output']);
                end
                tempFitfoutName = ['fitResults_iter' num2str(iter,'%04d') '.mat'];
                disp(['Saving the results as ' tempFitfoutName '. ChiSqaure = ' num2str(chisquared)]);
                if varyCtrlMode == 1
                save([outputFolderName filesep() tempFitfoutName],...
                    'model_name','boundsArray','param_strings','Nparam','Nstates',...
                    'maxMutationArray','maxMutationCountsArray','minMutationCountsArray',...
                    'A','K','P','Peq','tijs','sigmas','rates',...
                    'chisquared','guess','genNum','iter',...
                    'constructName','sample_description','save_prefix','yoff','zoff','c2ctrlAmp_0','c2ctrlTau_0', ...
                    'c2ctrlAmp_0', 'c2ctrlAmp','c2ctrlTau','c4ctrlAmp','chisquared_Weighted_array');
                else
                  save([outputFolderName filesep() tempFitfoutName],...
                    'model_name','boundsArray','param_strings','Nparam','Nstates',...
                    'maxMutationArray','maxMutationCountsArray','minMutationCountsArray',...
                    'A','K','P','Peq','tijs','sigmas','rates',...
                    'chisquared','guess','genNum','iter',...
                    'constructName','sample_description','save_prefix','yoff','zoff','chisquared_Weighted_array');  
                end
                
                current_LowestChisquared = chisquared;
                load(bestfitFileName,'chisquared','iter');
                prev_LowestChisquared = chisquared;%Save the previous X2 as prev_LowestChisquared
                chisquared = current_LowestChisquared;%reset chisquared to its current value
                if current_LowestChisquared < prev_LowestChisquared
                    %Just found a new best fit. Need to save it.
                    
                    fprintf('Congratulations! %f < %f (#1 fit). Will update best fit .\r\n',current_LowestChisquared,prev_LowestChisquared)
                    
                    %Moving the previous 2nd best fit to the number 3 slot
                    
                    if exist(secondLowestChiSquareFolderName,'dir') == 7
                        copyfile(secondLowestChiSquareFolderName,thirdLowestChiSquareFolderName);
                        disp('     Copied the previous second best fit into the third place slot');
                    end
                    
                    %Moving the previous best fit to the number 2 slot
                    copyfile(lowestChiSquareFolderName,secondLowestChiSquareFolderName);
                    disp('     Copied the previous best fit into the second place slot');
                    
                    outputFolderName = lowestChiSquareFolderName; %Assign the save folder to the lowest chi square value
                    
                elseif current_LowestChisquared == prev_LowestChisquared
                    
                    fprintf('Nothing changed.%f = %f. \r\n',current_LowestChisquared,prev_LowestChisquared);
                    saveMode = 0;
                elseif current_LowestChisquared > prev_LowestChisquared
                    fprintf('Oh no! %f > %f (#1 fit). Will check to see if it is lower than #2 best fit .\n',current_LowestChisquared,prev_LowestChisquared)
                    % Two options:
                    % 1) The current x2 is higher than the number 2 and 3. Update iteration number on the best fit
                    % 2) current x2 is higher than #1 but less than #3. Save it as #2.
                    
                    %------------------------------------------------------------------
                    % Check if  the new chisquared deserves spot number 2
                    %------------------------------------------------------------------
                    if exist(secondLowestChiSquareFolderName,'dir') == 7
                        load([secondLowestChiSquareFolderName filesep() 'BestFitResults' '.mat'],'chisquared');
                        save([secondLowestChiSquareFolderName filesep() 'BestFitResults' '.mat'],'iter','-append');
                        prev_2ndBestChisquared = chisquared;
                        chisquared = current_LowestChisquared;%reset chisquared to its current value
                        if current_LowestChisquared < prev_2ndBestChisquared
                            fprintf('Congratulations! %f < %f (#2 fit). You found a new #2 fit.  Updating 2nd best fit.\r\n',current_LowestChisquared,prev_2ndBestChisquared)
                            %Move the number 2 fit to the number 3 slot
                            copyfile(secondLowestChiSquareFolderName,thirdLowestChiSquareFolderName);
                            disp('     Copied the previous second best fit into the third place slot');
                            outputFolderName = secondLowestChiSquareFolderName;
                        elseif current_LowestChisquared == prev_2ndBestChisquared
                            fprintf('Nothing changed.%f = %f (#2 fit).  \r\n',current_LowestChisquared,prev_2ndBestChisquared);
                            saveMode = 0;
                        elseif current_LowestChisquared > prev_2ndBestChisquared
                            fprintf('Oh No! %f > %f (#2 fit). Checking to see if fit is better than #3:\n',current_LowestChisquared,prev_2ndBestChisquared);
                            %------------------------------------------------------------------
                            % Check if  the new chisquared deserves spot number 3
                            %------------------------------------------------------------------
                            if exist(thirdLowestChiSquareFolderName,'dir') == 7
                                current_iter = iter;
                                load([thirdLowestChiSquareFolderName filesep() 'BestFitResults' '.mat'],'chisquared','iter');
                                prev_iter = iter;
                                iter = current_iter;%Reset the iteration to the current one
                                save([thirdLowestChiSquareFolderName filesep() 'BestFitResults' '.mat'],'iter','-append');
                                prev_3rdBestChisquared = chisquared;
                                chisquared = current_LowestChisquared;%reset chisquared to its current value;
                                if current_LowestChisquared < prev_3rdBestChisquared
                                    fprintf('Congratulations! %f < %f (#3fit). You found a new #3 fit.  Updating 3rd best fit.\r\n',current_LowestChisquared,prev_3rdBestChisquared)
                                    outputFolderName = thirdLowestChiSquareFolderName;
                                    movefile([thirdLowestChiSquareFolderName filesep() 'BestFitResults' '.mat'],[GenAlgFilesLocationPrefix filesep() 'higherChiSquares'  filesep() 'fitResults_iter' num2str(prev_iter,'%04d') '.mat']);
                                    
                                elseif current_LowestChisquared == prev_3rdBestChisquared
                                    fprintf('Nothing changed.%f = %f (#3 fit). \r\n',current_LowestChisquared,prev_3rdBestChisquared);
                                    saveMode = 0;
                                elseif current_LowestChisquared > prev_3rdBestChisquared
                                    fprintf('Oh No! %f > %f (#3 fit).\n',current_LowestChisquared,prev_3rdBestChisquared);
                                    
                                    saveMode = 0;
                                    
                                end
                            elseif exist(thirdLowestChiSquareFolderName,'dir') ~= 7
                                outputFolderName = thirdLowestChiSquareFolderName;
                                
                                mkdir(outputFolderName);
                                disp(['Making a folder called "' outputFolderName '" to hold the output']);
                                
                            end
                        end
                    elseif exist(secondLowestChiSquareFolderName,'dir') ~= 7
                        
                        outputFolderName = secondLowestChiSquareFolderName;
                        
                        mkdir(outputFolderName);
                        disp(['Making a folder called "' outputFolderName '" to hold the output']);
                        
                    end
                end
            elseif exist(bestfitFileName,'file') ~= 2
                iter = 1;
                outputFolderName = lowestChiSquareFolderName;
                
            end
            
        end
       end
        
        %------------------------------------------------------------------
        % Update the folders with the correct information
        %------------------------------------------------------------------
      if sum(abs(imag(Peq))>0)~=0 || sum(Peq<0)~=0
          
      else
        if saveMode == 1
            
            disp(['Will save the data in ' outputFolderName]);
            
            foutName = [outputFolderName filesep() 'fitInputData.mat'];
            if fitHistMode == 1
                save(foutName,'histogram_FileName','Amp_bins','targetHistogram')
            end
            if fitC2Mode == 1
                if exist(foutName,'file') == 2
                    save(foutName,'C2_FileName','C2_exp_x','C2_exp_y','yoff','weightC2func','-append');
                else
                    save(foutName,'C2_FileName','C2_exp_x','C2_exp_y','yoff','weightC2func');
                end
                
            end
            if fitC4Mode == 1
                if exist(foutName,'file') == 2
                    save(foutName,'C4_FileName','C4_tau2eq0_exp','C4_tau1range','C4_tau3range','wC4func','zoff','-append');
                else
                    save(foutName,'C4_FileName','C4_tau2eq0_exp','C4_tau1range','C4_tau3range','wC4func','zoff');
                end
            end
            
            
            if useMatrixMethodMode == 1
                [Peq, denom_hist_sim, hist_sim, Peq_names] = histMaker_Nstate_v3normCorrect(P, A, sigmas, Amp_bins);
                p1_eq = Peq(1);
                p2_eq = Peq(2);
                p3_eq = Peq(3);
            end      
                
                if varyCtrlMode == 1
                save([outputFolderName filesep() 'BestFitResults' '.mat'],...
                    'model_name','boundsArray','param_strings','Nparam','Nstates',...
                    'maxMutationArray','maxMutationCountsArray','minMutationCountsArray',...
                    'A','K','P','Peq','tijs','sigmas','rates',...
                    'chisquared','guess','genNum','iter',...
                    'constructName','sample_description','save_prefix','yoff','zoff','c2ctrlAmp_0','c2ctrlTau_0', ...
                    'c2ctrlAmp_0', 'c2ctrlAmp','c2ctrlTau','c4ctrlAmp','chisquared_Weighted_array');
                else
                  save([outputFolderName filesep() 'BestFitResults' '.mat'],...
                    'model_name','boundsArray','param_strings','Nparam','Nstates',...
                    'maxMutationArray','maxMutationCountsArray','minMutationCountsArray',...
                    'A','K','P','Peq','tijs','sigmas','rates',...
                    'chisquared','guess','genNum','iter',...
                    'constructName','sample_description','save_prefix','yoff','zoff','chisquared_Weighted_array');  
                end
            
            %Save the paramaters used for the genetic algorithm
            save([outputFolderName filesep()  'genAlgParamaters.mat'],'programName','threshold',...
                'NmembersInitPop','maxGenerations','maxRepeats','Nreproduce','maxIterations',...
                'Nmutations','boundsArray','maxMutationArray',...
                'percentReproduce','percentToKeep',...
                'fitC2Mode','fitC4Mode','fitHistMode','fitHistData_mode',...
                'constructName','sample_description', 'save_prefix');
            
            %Save paramatres used to plot
            save([outputFolderName filesep() 'plottingParamaters.mat'],'Amp_bins','sigmas',...
                'C2_exp_x','C2_exp_y','constructName','sample_description', 'save_prefix',...
                'C2_sim','C4_sim','Peq','sigma_A','c2Control','c4Control');
            
        end
      end
        
        
        %--------------------------------------------------------------------------
        % Save the data PLOTS     (PART 5: Final Values)
        %--------------------------------------------------------------------------
        if sum(abs(imag(Peq))>0)~=0 || sum(Peq<0)~=0
           
        else
        if plotMode == 1
            if showBestFitMode == 1
                disp('     Loading the best fit from all runs.');
                load([lowestChiSquareFolderName filesep() 'BestFitResults.mat'],...                    
                    'A', 'sigmas','tijs','rates','Peq',...
                    'chisquared','guess',...
                    'iter','param_strings','Nparam','Nstates');
                saveMode = 1;
                %         disp('Setting saveMode to 0 because showBestFitMode = 1');
                outputFolderName = lowestChiSquareFolderName;
            end
            %//////////////////////////////////////////////////////////////////////////
            % Plot final results of the algorithm   (Part 5: Final Values)
            %//////////////////////////////////////////////////////////////////////////
            %     if fitHistMode == 1
            figure(1)
            hold on;
            
            if useMatrixMethodMode == 1
                [Peq, denom_hist_sim, hist_sim, Peq_names] = histMaker_Nstate_v3normCorrect(P, A, sigmas, Amp_bins);
                p1_eq = Peq(1);
                p2_eq = Peq(2);
                p3_eq = Peq(3);
            end
            
            if useMatrixMethodMode == 1
                figure(1); clf;
                set(gcf,'Color','w');
                set(gcf,'Name','Amp Histogram');
                lineStyle = char('g--','c--', 'm--','y--','k--','r-.','g-.','b-.','c-.','m-.','y-.','k-.'); % Define a list of colors to loop over
                Amp_bins = linspace(0,1,length(edgesAmp));
                
                if exist('data_hist_Plot','var') == 1
                    delete(data_hist_Plot)
                    %                                 disp('data_hist_Plot deleted')
                end
                % Replot data Histogram
                data_hist_Plot = plot(Amp_bins,targetHistogram, 'DisplayName', 'Data');
                data_hist_Plot.LineWidth = 2;
                data_hist_Plot.Color = 'blue';
                %                             x.DisplayName = 'Data';
                xlabel('Amp','FontSize',14);
                ylabel('Frequency','FontSize',14);
                
%                 [sample_description, ~] = sample_descriptionGetter();
                title_str = ['Experimental vs Simulated Histograms' ...
                    10 sample_description];
                title(title_str,'fontsize',14);
                hold on;
                % Clear plots from previous run
                if exist('histPlot','var') == 1
                    delete(histPlot)
                    %                             disp('histPlot deleted')
                end
                % Plot histogram of each state
                for i = 1:length(Peq)
%                     histPlot = plot(Amp_bins, Peq(i) * exp(-((Amp_bins - A(i))/sigma_A(i)).^2)./denom_hist_sim,...
                        histPlot = plot(Amp_bins,(Peq(i)/(sqrt(2*pi)*sigmas(i)))*exp(-((Amp_bins - A(i))/sigmas(i)).^2)./denom_hist_sim,...
                        lineStyle(i,:),'LineWidth',1,'DisplayName',[Peq_names(i,:) ' = ' num2str(Peq(i))]);
                    %                                 lgd = legend('show');
                    % lgd.Location = 'northwest';
                    hold on
                end
                set(gca, 'FontSize', 14);
                if exist('histPlotTot','var') == 1
                    delete(histPlotTot)
                end
                hold on
                % Plot overall histogram
                histPlotTot = plot(Amp_bins, hist_sim,...
                    lineStyle(i),'LineWidth',1,'DisplayName','Total Fit', 'Color','r');
                lgd_tot = legend('show');
                lgd_tot.FontSize = 14;
            end
            
            
            
            if saveMode == 1
                saveas(gcf,[outputFolderName filesep() 'hist_bestFit.fig'])
            end
            %     end
            
            %---------------------------------------------------------------------------------------
            % (2) Optimization Target #2: 2-point TCF (C2)  (20 uSec and on)  (Part 5: Final Values)
            %----------------------------------------------------------------------------------------
            %     if fitC2Mode == 1
            figure(2)
            
            if useMatrixMethodMode == 1
                [C2_sim,C2time] = PA2C2(P,A,C2time,yoff,addControlMode,c2Control);
            end
            
            
            if normalizeMode == 1
                C2_sim = C2_sim./C2_sim(1);
            end
            
            hold on;
            if exist('C2_sim_plot','var') == 1
                delete(C2_sim_plot)
            end
            C2_sim_plot = plot(C2time,C2_sim,'r-',...
                'LineWidth',2,'DisplayName','Final Fit');
            set(gca,'xscale','log');
            drawnow();
            
            if saveMode == 1
                saveas(gcf,[outputFolderName filesep() 'C2_bestFit.fig'])
            end
            %     end
            
            %---------------------------------------------------------------------------------------
            % (3) Optimization Target #3: 4-point TCF (C4) (20 uSec and on)  (Part 5: Final Values)
            %---------------------------------------------------------------------------------------
                if fitC4Mode == 1
            
            if useMatrixMethodMode == 1
                C4time = C4_tau1range;
                tau2 = 0;
                [C4_sim,C4time] = PAK2C4(P,A,K,C4time,tau2,zoff,addControlMode,c4Control);
            end
            C4_sim = C4_sim; %Add the y-offset^2 to C4_sim
            if normalizeMode == 1
                C4_sim = C4_sim./C4_sim(1,1);
            end
            
            if exist('C4_plot','var') == 1
                delete(C4_plot)
            end
            figure(3)
            hold on;
            C4_plot = surf(C4_tau1range,C4_tau1range,real(C4_sim));
            set(gca,'xscale','log');
            set(gca,'yscale','log');
            drawnow();
            
            if saveMode == 1
                saveas(gcf,[outputFolderName filesep() 'C4_bestFit.fig'])
            end
            %     end
            
                end
                
            if plot_paramater_histogram_mode
               if saveMode ==1
                  saveas(histChart,[outputFolderName filesep() 'ParamHist_bestIter.fig']) 
               end
            end
            %--------------------------------------------------------------------------
            % (6) Display a graphic of the network        (Part 5: Final Values)
            %--------------------------------------------------------------------------
         if plot_model_GenerationalUpdatesMode
            figure(6);
            clf;
            set(gcf,'Name','Model: 3state123 cyclical');
            if useFigurePosnMode == 1
                switch computer_terminal_str
                    case 'coputer_baimi_mode' %Work Desktop
                        set(gcf,'Position',fig_6_posn);
                    case 'computer_bisraels_mode' %Macbook pro
                        set(gcf,'Position',fig_6_posn);
                end
                
                
            end
            set(gcf,'Color','w');
            
            set(gcf,'Name','Model final Results:  GenAlg');
            set(gcf,'Color',[1 1 1]);
            hold on;
            xlim([0 1.1]);
            ylim([0 1]);
            
            %Designate spots for the states
            state1_loc = [0.5 1];
            state2_loc = [1 0];
            state3_loc = [0 0];
            
            %Plot the state symbols
            text(state1_loc(1)-0.04,state1_loc(2)+0.05,'1','FontSize',24)
            text(state2_loc(1),state2_loc(2),'2','FontSize',24);
            text(state3_loc(1)-0.05,state3_loc(2),'3','FontSize',24);
            
            %Plot the Amp values
            text(state1_loc(1),state1_loc(2) + 0.05,['=' num2str(A1,'%.3f')],'FontSize',16)
            text(state2_loc(1)+ 0.05,state2_loc(2),['=' num2str(A2,'%.3f')],'FontSize',16);
            text(state3_loc(1),state3_loc(2),['=' num2str(A3,'%.3f')],'FontSize',16);
            
            %Plot Lines Between states
            line([state1_loc(1) state2_loc(1)],[state1_loc(2) state2_loc(2)],'Color','k','LineStyle','-');
            line([state2_loc(1) state3_loc(1)],[state2_loc(2) state3_loc(2)],'Color','k','LineStyle','-');
            line([state3_loc(1) state1_loc(1)],[state3_loc(2) state1_loc(2)],'Color','k','LineStyle','-');
            
            %Plot the inverse of the rates
            t12_loc = [state1_loc(1)+(state2_loc(1)- state1_loc(1))/2,state1_loc(2)+(state2_loc(2)- state1_loc(2))/2];
            t12_msg = ['t_{12} = ' num2str(t12) 10 't_{21} = ' num2str(t21)];
            text(t12_loc(1),t12_loc(2),t12_msg,'FontSize',14);
            
            t23_loc = [state3_loc(1)+(state2_loc(1)- state3_loc(1))/2,state3_loc(2)+(state2_loc(2)- state3_loc(2))/2];
            t23_msg = ['t_{23} = ' num2str(t23) 10 't_{32} = ' num2str(t32)];
            text(t23_loc(1),t23_loc(2),t23_msg,'FontSize',14);
            
            t31_loc = [state1_loc(1)+(state3_loc(1)- state1_loc(1))/2,state1_loc(2)+(state3_loc(2)- state1_loc(2))/2];
            t31_msg = ['t_{31} = ' num2str(t31) 10 't_{13} = ' num2str(t13)];
            text(t31_loc(1),t31_loc(2),t31_msg,'FontSize',14);
            
            %Plot the chi-squared value
            text(state3_loc(1),state1_loc(2),['\chi^2 = ' num2str(chisquared)],'FontSize',14);
            
            %Plot a title with the information in it
%             [sample_description, save_prefix] = sample_descriptionGetter();
            title_str = ['3state Cyclical: ' sample_description];
            text(state3_loc(1),state1_loc(2)+0.15,title_str,'fontsize',16);
            ylim([0,1.1]);
            
            
            axis off;
            if saveMode == 1
                saveas(gcf,[outputFolderName filesep() 'ModelResultsFigure.fig'])
            end
        end
        
        end
       end
        
        if windowIter==maxWindIters
        
%         set the boolean array up based onthe cuulative histograms - then
%         reset them both
        for j=1:length(histCell)
%         [v,I]=maxk(histCell{1,j},2);
%             if I(1)==1 && v(2)<(v(1)/2) && atEdges(j)==0
%            paramAdjs(j)=paramAdjs(j)-1;    
%             elseif I(1)==numel(histCell{1,j}) && v(2)<(v(1)/2) && atEdges(j)==0
%             paramAdjs(j)=paramAdjs(j)+1;
%             end
%             UPDATE: 3/10/2022 - make the shift criteria that the lower/upper 1/4 of the
%             LogHistArray contains 75% or more of the occurences
% UPDATE : 04 -11 - 22 turned off the paramAdjs sliding windows to allow
% wide bounds to be set in modelBuilder rather than slide narrow bounds
%             curHist=histCell{1,j};
%             if sum(curHist(1:end/4)) >= (0.75)*(sum(curHist)) && atEdges(j)==0
%                paramAdjs(j)=paramAdjs(j)-1; 
%             elseif sum(curHist(3*(end/4):end)) >= (0.75)*(sum(curHist)) && atEdges(j)==0
%                paramAdjs(j)=paramAdjs(j)+1;  
%             end
        end
        
          outputFolderName = [GenAlgFilesLocationPrefix filesep() 'ParameterAdjustments'];
%           if exist(outputFolderName,'dir') ~= 7
%                     mkdir(outputFolderName);
% %                     disp(['Making a folder called "' outputFolderName '" to hold the output']);
%           end
          boundsArrayCurrent=boundsArray; 
           save(outputFolderName,'paramAdjs','histCell','boundsArray_noAdj','boundsArrayCurrent'); 
                
        histCell=cell(1,Nparam-2*Nstates);
        windowIter=0;
    else
        windowIter=windowIter+1;  
        end
    
        
    end
    
    
    iter = iter + 1;
    %--------------------------------------------------------------------------
    %The end of the program
    %--------------------------------------------------------------------------
    disp(['iter = ' num2str(iter) ': ' programName ' complete.']);
end

disp([programName ' complete. All ' num2str(NconstructFolderNames) ' completed.']);

end
% end
