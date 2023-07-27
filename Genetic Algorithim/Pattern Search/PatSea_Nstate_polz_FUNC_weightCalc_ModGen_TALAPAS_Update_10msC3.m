%__________________________________________________________________________
% AUTHOR:Jack Maurer +  Claire Albrecht
%
% NAME: PatSea_Nstate_CA.m
%
% FUNCTION: % Fits (1) the FRET Histogram, (2) the 2pt TCF and (3) the 4-point TCF

%--------------------------------------------------------------------------
% PROCEDURE:
%--------------------------------------------------------------------------
% (1) Load the experimental histogram
% (2) Load the experimental 2pt time correlation functions to fit
% (3) Load the experimental 4pt TCFs for various tau2 values
% (4) Make an initial guess of the model paramaters:
%       x = [t12,t13,t21,t23,t31,t24,t42,t25,t52,t56,t63,t36,y37,t73,t68,t86,A1,A2,A3,A4,A5,A6,A7,A8];
%   for 9 state: x = [t12,t13,t21,t23,t31,t24,t42,t25,t52,t56,t63,t36,y37,t73,t68,t86,t89,t98,t69,A1,A2,A3,A4,A5,A6,A7,A8,A9];
% (5) Make an array of lower bounds
% (6) Make an ar1ray of upper bounds
% (7) feed the guess into the optimization algorithm: x = patternsearch(fun,x0,A,b,Aeq,beq,lb,ub);

%--------------------------------------------------------------------------
% INPUT: ############## EDIT FROM 3 state #############
%--------------------------------------------------------------------------
% (1) FRET Histogram :
% (2) Two-Point TCF :
% (3) Four-Point TCF :
% (4) Solutions to Diff Eqns:  [DONT NEED THIS ANYMORE?]

% OUTPUT:
% (1) BestFitResults.mat         %All the fitting paramaters
% (2) BestFitRestults_hist.fig   %
% (3) fitInputData.mat
% (4) genAlgParamaters.mat
% (5) ModelResultsFigure.fig
% (6) plottingParamaters.mat

%--------------------------------------------------------------------------
% EXTERNAL PROGRAMS CALLS: (needs these codes in MATLAB PATH to function)
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
% MODIFICATION LOG:
%--------------------------------------------------------------------------
% close all
% clear all
%%
function [] = PatSea_Nstate_polz_wGActrl_FUNC_multiOut_v2_weightCalc_ModGen(model_idx,setIdx)
%__________________________________________________________________________
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
programName = 'PatSea_Nstate_polz';
global NbinsHist histResUsec
global normalizeMode testTimesMode Nparam

NbinsHist = 25;
pointsPerDecade=25;
histResUsec = 1000;
% disp(['Now Running ' programName '.m']);

% addpath statements for talapas use case - 06012022
  addpath(genpath('/gpfs/home/jmaurer3/Documents/MATLAB/MatlabCodes/KineticModeling/KineticNetwork'));


%% PART 1: Set the program up with various options
constructFolderNames = {'E3E4'};
%constructFolderNames = {'E1E2'};
% constructFolderNames = {'E5E6'};
%constructFolderNames = {'(+15)Dimer'};

saltConditions='100mM_NaCl';
%saltConditions='20mM_NaCl';
%saltConditions='300mM_NaCl';
%saltConditions='20mM_NaCl_0mm_MgCl2';
% saltConditions='1uMgp4462_1p5uMgp45'; 

altWeightScheme=1;
Nstates = 3;
% PatSea output to load and pick up on for optimization 
loadPatSeaMode=0;

testTimesMode=1; %will adjust the C2 and C4 time range so it can be matched
if testTimesMode
    time_LB_usec=10000;
%     time_LB_usec=250;
    time_UB_sec= 2.5; 
end

%label on patSea folder containing output##'s'
patSeaWeightNum=6;

%weight num for the scheme being used at top level Data folder
%*****CHANGE THIS FOR DATA OTHER THAN -2 100mM******
weightVersionNum=2;

pauseBetweenConstructsMode = 0;

% model_set = {'model_loopN'};
model_list = { 'model_lin', 'model_loop1','model_loopN','model_brch'}; 
model_set = model_list(setIdx); 
modelNum=model_idx; 
% loadMode=1;


% noise increase on the c2 and c4 at the latest tau compared to th first
% tau (across decades) - in terms of % noise change, i.e. last decade 15%
% noisier than first decade . Base on 1/sqrt(NtotalPairs) of TCF wAVG

%******************NEW WEIGHT SCHEME FOR V2 CALC*****************

noisePercentDiffC2=0.1;
noisePercentDiffC4=0.25;
decadeSegments=10;

% c4 slider less than 1 is flatter, greate than 1 is steeper (at a fixed
% noise floor and ramp) - lower limit of ~0.15 to avoid going below
% secondary maxima on c4 weight func surface.

%3 state data
c4Slider=1.5;

%+15 Dimer
% c4Slider=2.0;


% set the base line noise on the c2/c4 as well as hist. less baseline noise
% will steepen the weight function for the c2/c4
noiseFloorHist=0.02;
noiseFloorC2=0.01;
noiseFloorC4=0.05;

% set the relative contributions of the 3 surfaces 
%c2Mag=2;
%c4Mag=25;
%histMag=10;

%new intial scalars for the GenAlg with WeightCalc_v3_chi2
%c2Mag=4000;
%c4Mag=1e6;
%histMag=1e5;

c2Mag=3000;
c4Mag=4e3;
histMag=1e5;
%****************************************************************

%based on 10-13-2022 reuslts 
% set the relative contributions of the 3 surfaces given the results from each data set

if strcmp(constructFolderNames,'E5E6')

    if strcmp(saltConditions,'100mM_NaCl')
    
    %new scalars for the true chi2 calc 1/2/2023
    c2Mag=c2Mag;
    c4Mag=c4Mag;
    histMag=histMag;
    outputToLoad=1;

    elseif strcmp(saltConditions,'20mM_NaCl')
    c2Mag=c2Mag;
    c4Mag=c4Mag;
    histMag=histMag;
        outputToLoad=1;

     elseif strcmp(saltConditions,'20mM_NaCl_0mm_MgCl2')
    
    c2Mag=c2Mag;
    c4Mag=c4Mag;
    histMag=histMag;
        outputToLoad=1;

    elseif strcmp(saltConditions,'300mM_NaCl')
    c2Mag=c2Mag;
    c4Mag=c4Mag;
    histMag=histMag;
        outputToLoad=1;

    end

elseif strcmp(constructFolderNames,'E1E2')

    c2Mag=c2Mag;
    c4Mag=c4Mag;
    histMag=histMag;
        outputToLoad=1;


elseif strcmp(constructFolderNames,'E3E4')
    
    if strcmp(saltConditions,'100mM_NaCl')
        c2Mag=c2Mag;
        c4Mag=c4Mag;
        histMag=histMag;
            outputToLoad=1;

      elseif strcmp(saltConditions,'20mM_NaCl')
        c2Mag=c2Mag;
        c4Mag=c4Mag;
        histMag=histMag;
            outputToLoad=1;

    elseif strcmp(saltConditions,'20mM_NaCl_0mm_MgCl2')
        c2Mag=c2Mag;
        c4Mag=c4Mag;
        histMag=histMag;
            outputToLoad=1;

     elseif strcmp(saltConditions,'300mM_NaCl')
        c2Mag=c2Mag;
        c4Mag=c4Mag;
        histMag=histMag;
            outputToLoad=1;

      end

elseif strcmp(constructFolderNames,'(+15)Dimer')

    c2Mag=c2Mag;
    c4Mag=c4Mag;
    histMag=histMag;
        outputToLoad=1;

    
end


% ******************** SELECT MODEL YOU WANT TO USE ******************** %
% get model name from genAlgParameters.mat based on what was used in
% GenAlg_Nstate
% ********************************************************************************
% global model_name
%%%%% Which GenAlg_Nstate output do you want? %%%%%
% Nstates = 3;
% GenAlg_progName = 'GenAlg_Nstate_test3state_output';


%SET OF GLOBAL PARAMETERS SO THAT MODEL GENERAOTR AND LIBRARY LOADER ARE
%COMPATIBLE WITH THE FORM OF MULTIGOAL_NSTATE FUNCTION - IN PRINCIPLE
%EVERYTHIGN HERE IS "STATIC" AFTER FIRST CALL, SO SHOULD NOT BE PROBLEMATIC
%FOR GLOBAL DECLARATION
global model_set modelNum model_name dbCell rates_str_woDB rates_idx_woDB rates_idx_DB 
global param_strings  boundsArray maxMutationArray maxMutationCountsArray minMutationCountsArray atEdges
global boundsArray_noAdj model_lin model_loop1 model_loopN model_brch


GenAlg_progName = ['GenAlg_Nstate_updated2'];
folder2load = 'lowestChiSquare';
MaxIterations = 10000;

% add the extra resolution suffix to the file path for variable C2/C4 resol
addC2resSuffix=1;

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
    fileSuffix=[fileSuffix '_c2_' num2str(time_LB_usec) 'usec_v' num2str(weightVersionNum)];
end
%--------------------------------------------------------------------------
% Declare global variables
%--------------------------------------------------------------------------
global normalizeMode verboseMode guessUpdateMode diagnoseMode
global fitHistMode fitC2Mode fitC4Mode
global targetHistogram weightingFactor_FREThist
global  FRET_bins  %sigma_A  Histogram optimization
% global sigma_A1 sigma_A2 sigma_A3 sigma_A4 sigma_A5 sigma_A6 sigma_A7 sigma_A8 sigma_A9
global C2_exp_x C2_exp_y weightingFactor_C2  weightC2func 
global C4_tau1range C4_tau2eq0_exp weightingFactor_C4_t0 wC4func 
global showProgressOnFit_mode overwriteWeightsMode
global Amp_bins weightingFactor_Amphist sample_description ctrlScaling ctrlTime

%--------------------------------------------------------------------------
% User Options                                                  (Part 1)
%--------------------------------------------------------------------------
%%%%% What protein condition are you dealing with? %%%%%
global protein_str Nstates polzMode
protein_str = 'gp32_0p0uM';
% protein_str = 'gp32_0p5uM';
polzMode = 1;

% *******parameters and functions for addition of control timescales*******
global addControlMode
addControlMode=1;

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
global varyCtrlMode NctrlParams   % added to vary the control timescale
global c2ctrlAmp_GenAlg c4ctrlAmp_GenAlg c2ctrlTau_GenAlg hardSetCtrlMode

varyCtrlMode=1;
hardSetCtrlMode=1;
NctrlParams=3;

if time_LB_usec==1000
c2AmpFitVal = 0.509;
c2TauFitVal = 3; 

elseif time_LB_usec==250 || time_LB_usec==10000
% % based on 4 ExP fit 
% c2AmpFitVal = 0.30837;
% c2TauFitVal = 2.5; 

% % based on 5exp fit
% c2AmpFitVal = 0.329;
% c2TauFitVal = 3.5; 

% arbitrary but closer? +1 Dimer
%c2AmpFitVal = 0.3; 
%c2TauFitVal = 3; 

% ***** (-2) Dimer Values (100mM) ****** 
c2AmpFitVal = 0.1;
c2TauFitVal = 5; 

% ***** (-2) Dimer Values (20mM) + 300mM ****** 
%c2AmpFitVal = 0.125;
%c2TauFitVal = 4; 

% ***** (+1) Dimer 20mM NaCl ******
%c2AmpFitVal = 0.2;
%c2TauFitVal = 3; 

% ***** (+1) Dimer 300mM NaCl ******
%c2AmpFitVal = 0.1;
%c2TauFitVal = 3.5; 

% ***** (-1) Dimer Values ******
%c2AmpFitVal = 0.23;
%c2TauFitVal = 6; 

% +15 values - also run by accident with +1 minus23
%c2AmpFitVal = 0.155;
%c2TauFitVal = 4; 

% Pat values for P1% 
% c2AmpFitVal = 0.19;
% c2TauFitVal = 3.6; 

% ***** (-2) Dimer Values 20mm/0mm NaMg ******
%c2AmpFitVal = 0.1;
%c2TauFitVal = 3; 

% ***** (+1) Dimer Values 20mm/0mm NaMg ******
%c2AmpFitVal = 0.1;
%c2TauFitVal = 4; 

end
% ***********************************************************************
% first range is for the amplitude of the exp, second is for the coefficent
% on the 4pt outer product
% NctrlParams = 3; 
% ctrlScaling=[0.75*c2AmpFitVal 1.25*c2AmpFitVal ;...
%             0.25*c2AmpFitVal 1.75*c2AmpFitVal];
% ctrlTime=[0.5*c2TauFitVal 1.5*c2TauFitVal];

% Nstates = 5;
% GenAlg_progName = 'GenAlg_Nstate_Brett_varySig_5stateNoLoops';


% outputProgramName = [GenAlg_progName(1:end-6) programName];
outputProgramName = [GenAlg_progName '_' programName];
% tag the output with the GenAlg model and the patternsearch version

normalizeMode = 0;%Normalizes the 2pt and the 4pt correlation functions
verboseMode = 1;
guessUpdateMode = 0; %Very detailed: shows changes to any paramaters
diagnoseMode = 1; %Shows relative contribution of the 2pt TCF to the histogram at the end of each chi square calculation
showProgressOnFit_mode = 1;%Updates fit during each iteration of optimization

clockMode = 0; %Times various features of the algorithm
saveModeDefault = 1;
saveMode = saveModeDefault;
datestr_mode = 0;
plotMode = 1;                               %Makes plots.
useFigurePosnMode = 1;

fitHistMode = 1; fitHistData_mode = 1;%If 0 you will fit to the histogram fit
fitC2Mode = 1;
fitC4Mode = 1;
overwriteWeightsMode=1; 

showBestFitMode = 0;
% updateExcelModeDefault = 0;
% updateExcelMode = updateExcelModeDefault;
if sum([fitHistMode,fitC2Mode,fitC4Mode]) ==0
    error('No surfaces to optimize to! Pick atleast one.');
end

%----------------------------------------------------------------------
%% Start in the single molecule folder (smData_Processed): comp specific
%----------------------------------------------------------------------

%Start in the single molecule folder (smData_Processed): comp specific
% [computer_terminal_str, terminalID, DropboxLocationPrefix] = computerMode();
global computer_terminal_str terminalID
% [computer_terminal_str, terminalID] = computerMode_JM();
computer_terminal_str ='computer_JackLaptop_mode';
terminalID = 'C:\Users\Ryzen 5\';
%% -----------------------------------------------------------------------
%  Set up the figure positions once and done
% -----------------------------------------------------------------------
if plotMode == 1
    clear figPosn
    if useFigurePosnMode == 1
%         figPosn = setFigPosn_JM(terminalID);
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
            set(gcf,'Name','FRET Histogram');
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
    end
    
    
    %-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.
end
%----------------------------------------------------------------------
%% Loop Over the different constructs
%----------------------------------------------------------------------
NconstructFolderNames = numel(constructFolderNames);

for construct_idx = 1:NconstructFolderNames
    %         for construct_idx = NconstructFolderNames:-1:1
    constructName = char(constructFolderNames(construct_idx));
    disp(['       Construct ' num2str(construct_idx) '/' num2str(NconstructFolderNames) ': ' constructName]);
    
    if exist('constructName','var') ~= 1
        error('Not sure what constrcut to analyze');
    end
     
    
    %--------------------------------------------------------------------------
    %% Navigate to the correct folder & load data file paths
    %--------------------------------------------------------------------------
    
     %********************START LIBRARY LOADER SECTION ******************
 
%     intialize and load the relevant parameters for the current model -
%     then load best fit results so the tijs, K , A, sigmas are all correct
%     from the GenAlg best fit. 

    [model_lin,model_loop1,model_loopN,model_brch] = model_generator_v5(Nstates,model_set);
   
   % initialize dummy variables
    c2ctrlAmp_0= c2AmpFitVal;
    c4ctrlAmp_0= c2AmpFitVal;
    c2ctrlTau_0= c2TauFitVal; 
    c2ctrlAmp= c2AmpFitVal;
    c4ctrlAmp= c2AmpFitVal;
    c2ctrlTau= c2TauFitVal; 
   
    loadMode = 1;
    A = zeros([Nstates, 1]);
    tijs = ones([5*Nstates, 1]);
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
    
 %***************************************************************************
    
    
    % Find GenAlg files and load target data:
%     protein_str = 'gp32_0p5uM';
%     modelName = 'GenAlg8or9_9state_1000usec25binHist_output'; %%%% Figure out a better way to select an file to optimize than putting in by hand?
%     modelName = 'GenAlg_Nstate_test3state_output';    

%     [GenAlgFilesLocationPrefix, genAlgDataLocationPrefix] = fmincon_fileLocator(terminalID, constructName, GenAlg_progName, protein_str);
    
    normalizeMode = 0 ;
%     [GenAlgFilesLocationPrefix, genAlgDataLocationPrefix] = polz_fileLocator(terminalID, constructName, model_name, GenAlg_progName, protein_str);
    
    [GenAlgFilesLocationPrefix, genAlgDataLocationPrefix, histogram_FilePath, C2_FilePath,C4_FilePath] = polz_fileLocator_wCtrl_ModGen(terminalID, constructName, model_name, GenAlg_progName, protein_str, saltConditions, fileSuffix);

    cd(genAlgDataLocationPrefix);
    cd ..   %Go to the GenAlgFiles Folder
    wd = pwd;
    if verboseMode == 1
        disp(['Part1 complete: You are now in the folder: ' wd]);
    end
    
    %//////////////////////////////////////////////////////////////////////////
    %% PART 2: Load the target data (and plot if opted)
    %\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

    %--------------------------------------------------------------------------
    % (2.1) Optimization Target #1: 1-D FRET HISTOGRAM (PART 2: Load the target data)
    %--------------------------------------------------------------------------
    % if fitHistMode == 1
    %     load(histogram_FileName,'xData','yData','xFit','yFit');
    load(histogram_FilePath,'edgesAmp','totalAmpHist', 'numMolecules');
   
%     FRET_bins = xData;
     Amp_bins = edgesAmp;
    targetHistogram = totalAmpHist; %yData;
    
    targetHistogram = targetHistogram./sum(targetHistogram);
    
    if plotMode == 1
        figure(1);
        clf;
        set(gcf,'Name','Visibility Histogram');
        set(gcf,'Color','w');
        
        data_hist_Plot = plot(Amp_bins,targetHistogram);
        data_hist_Plot.LineWidth = 2;
        data_hist_Plot.Color = 'blue';
        data_hist_Plot.DisplayName = 'Data';
        
        xlabel('visibility','FontSize',14);
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
        hold on;
    end
    % end
    
    %--------------------------------------------------------------------------
    % (2) Optimization Target #2: 2-point TCF (C2)           (20 uSec and on)  (PART 2: Load the target data)
    %--------------------------------------------------------------------------
    if fitC2Mode == 1

        if strcmp(protein_str, 'gp32_0p5uM') == 1
            load(C2_FilePath,'C2avg','tauArraySec','res','Nfiles')
            C2_exp_x = tauArraySec;
            C2_exp_y = C2avg;
        elseif polzMode == 0
            load(C2_FilePath,'time','yData','y','yoff');
            C2_exp_x = time;
            C2_exp_y = yData;
        elseif polzMode == 1
            load(C2_FilePath, 'stitchTaus','stitchTCF')
        end
        C2_exp_x = stitchTaus;
        C2_exp_y = stitchTCF;
        
        tauArraySec = stitchTaus;%tauArraySec;

        %--------------------------------------------------------------
        %Take out the lower bound time limit and cut off at the upper bound
        time_usec = tauArraySec*1e6;
        startIndList=find((time_usec>=time_LB_usec));
        endIndList=find((time_usec<=time_UB_sec*1e6));
        C2_exp_x=stitchTaus(startIndList(1):endIndList(end));
        %             stitchTaus=stitchTaus((stitchTaus<=cutOffTime));
        %             C2_exp_x = stitchTaus;
        C2_exp_y = real(stitchTCF(startIndList(1):endIndList(end)));
        
        C2_exp_x=reshape(C2_exp_x,1,length(C2_exp_x)); 


        
        c2Control = c2Exp(C2_exp_x,c2AmpFitVal,c2TauFitVal)*C2_exp_y(1);
        
        c2CorrectedSurf = C2_exp_y - c2Control; 
        yoff = abs(mean(c2CorrectedSurf(end-6:end)));%New declaration of y^offset
%             yoff=0;
        
%         load(C2_FilePath,'tauArraySec','C2avg');
%         yoff = abs(mean(C2avg(end-3:end)));
        fprintf('Setting yoff = abs(mean(C2_exp_y(end-3:end))) = %f\r\n',yoff);
        
      
        
        
        if normalizeMode == 1
            C2_exp_y = C2_exp_y./C2_exp_y(1);
        end



        if plotMode == 1
            figure(2); clf;
            set(gcf,'Color','w');
            hold on;
            plot(C2_exp_x,C2_exp_y,'b.','MarkerSize',10,'DisplayName','C^{(2)}(\tau) Data');
            
           sample_description=[constructName ' ' saltConditions];
            title_str = ['Experimental vs Simulated C2' ...
                10 sample_description];
            title(title_str,'fontsize',14);
            
            xlabel('Time (sec)','FontSize',14);
            ylabel('C^{(2)}(\tau)','FontSize',14);
            set(gca,'xscale','log');
            drawnow();
            hold on;
        end
    end
    
    %---------------------------------------------------------------------------------------------
    % (3) Optimization Target #3: 4-point TCF (C4)    (10 uSec and on) (PART 2: Load the target data)
    %---------------------------------------------------------------------------------------------
    if fitC4Mode == 1
        
%         load(C4_FilePath,'C4','tau1arraySec','tau3arraySec','tau2ValSec');
        load(C4_FilePath,'C4finalArray','tauArrayUsecTotal','time','tau2usec');
        
        load(C4_FilePath,'C4finalArray','time','tau2usec');
        startIndList=find((time>=time_LB_usec*1e-6));
        endIndList=find((time<=time_UB_sec));
        time=time(startIndList(1):endIndList(end));
        C4_tau2eq0_exp = real(C4finalArray(startIndList(1):endIndList(end),startIndList(1):endIndList(end)));        
        

        C4_tau1range = time;
        C4_tau3range = time';
     
%         time = C2_exp_x;
        tau2ValSec = tau2usec;
        tau2 = tau2usec;
        
        %         load(C2_FilePath,'time','yData','y','yoff');
        %         zoff = yoff*yoff;
%         zoff = abs(mean(mean(C4_tau2eq0_exp(end-3:end-3))));
%         fprintf('setting zoff = abs(mean(mean(C4_tau2eq0_exp(end-5:end-5)))) = %f\r\n',zoff);
        if normalizeMode == 1
            C4_tau2eq0_exp = C4_tau2eq0_exp./C4_tau2eq0_exp(1,1);
        end
        
        % 4-pt. TCF weighting function
%         wC4func = 1./(C4_tau1range).*(1./((C4_tau1range)));
        % new weight function based on both lack of outer product and less steep
        % slope at early times - 06-22-22
%         if altWeightScheme
%          wC4func =((1./nthroot(C4_tau1range,1.25)).*(1./nthroot(C4_tau1range,1.25)'));  
%         else
%           wC4func =((1./nthroot(C4_tau1range,2)).*(1./nthroot(C4_tau1range,2)')); 
%         end

          c4Control = c4Exp(C4_tau1range,c2AmpFitVal,c2TauFitVal)*C4_tau2eq0_exp(1,1);
        
        c4CorrectedSurf = C4_tau2eq0_exp - c4Control; 
        zoff = abs(mean(mean(c4CorrectedSurf(end-6:end,end-6:end))));
        
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
            
            sample_description=[constructName ' ' saltConditions];
            title_str = ['C^{(4)}(\tau_1, \tau_2 = ' num2str(tau2ValSec) '\musec, \tau_3)' ...
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
            hold on;
        end
    end
    
    
    if verboseMode == 1
        disp('     Part1: Done Loading and (optionally) plotting the data.');
    end
    
    %% PART 3: LOAD THE BEST FIT (Determine starting paramaters)
% Go to the folder with the best GenAlg outputs and load lowest chi squared
% best fit results
%     cd([GenAlgFilesLocationPrefix filesep() '3state123_cyclical' filesep() ); 
%     GenAlg_best_output_folder = ['3state123_cyclical' filesep() GenAlg_prog_Name filesep() 'lowestChiSquare'];
    
%     GenAlg_best_output_folder = [GenAlgFilesLocationPrefix filesep() GenAlg_progName filesep() 'lowestChiSquare'];
%     GenAlg_best_output_folder = ['/Users/clairealbrecht/Dropbox/MarcusLab/Data/smData_Processed/S1S2/gp32_0p5uM/ChosenMolecules/genAlgFiles/outputs' filesep() GenAlg_progName filesep() 'lowestchisquare'];
%     GenAlg_best_output_folder =[ '/Users/clairealbrecht/Dropbox/MarcusLab/Data/smData_Processed/S1S2/gp32_0p5uM/ChosenMolecules/genAlgFiles/outputs/GenAlg_Nstate_Brett_varySig_5stateNoLoops/lowestchisquare'];
    
    GenAlg_best_output_folder = [GenAlgFilesLocationPrefix filesep() folder2load];
    PatSea_best_output_folder = [GenAlgFilesLocationPrefix filesep() programName '_outputs_v' num2str(patSeaWeightNum) ...
                                 filesep() 'output' num2str(outputToLoad)];
    
    if loadPatSeaMode
        cd(PatSea_best_output_folder)
        if exist(['BestFitResults_' programName '.mat']);
            load(['BestFitResults_' programName '.mat'])
            sigma_A = sigmas;
        else
            error(["Can't find the BestFitResults file for this PatSea Output Number"])
        end
        
    else
    
    cd(GenAlg_best_output_folder)
    if exist('BestFitResults.mat')
        load('BestFitResults.mat')
        sigma_A = sigmas;
    else
        error(["Can't find the BestFitResults.mat file in: " GenAlg_best_output_folder])
    end
    
    end
%     RIGHT HERE: if the values within Peq are negative or imag then auto
%     generate a guess for the intial points in PatSea
    paramAdjs = zeros(length(tijs));
%     [wC4func,weightC2func,weightingFactor_C2,weightingFactor_C4_t0,weightingFactor_Amphist] = weightsCalculator_v1(noiseFloorHist,noiseFloorC2,noiseFloorC4,noisePercentDiffC2,noisePercentDiffC4,C2_exp_x, C2_exp_y, C4_tau1range, C4_tau2eq0_exp, targetHistogram,c2Mag,c4Mag,histMag);
% [wC4func,weightC2func,weightingFactor_C2,weightingFactor_C4_t0,weightingFactor_Amphist] = weightsCalculator_v2(noiseFloorHist,noiseFloorC2,noiseFloorC4,noisePercentDiffC2,noisePercentDiffC4,C2_exp_x, C2_exp_y, C4_tau1range, C4_tau2eq0_exp, targetHistogram,c2Mag,c4Mag,histMag,decadeSegments,c4Slider);    
[wC4func,weightC2func,weightingFactor_C2,weightingFactor_C4_t0,weightingFactor_Amphist] = weightsCalculator_v3_chi2(noiseFloorHist,noiseFloorC2,noiseFloorC4,noisePercentDiffC2,noisePercentDiffC4,C2_exp_x, C2_exp_y, C4_tau1range, C4_tau2eq0_exp, targetHistogram,c2Mag,c4Mag,histMag,decadeSegments,c4Slider);

weightHistFunc = ones(size(Amp_bins));


if sum(abs(imag(Peq))>0)~=0 || sum(Peq<0)~=0
%     first define the control params for the model builder call
   if hardSetCtrlMode==1
    c2ctrlAmp_GenAlg = c2AmpFitVal;
    c4ctrlAmp_GenAlg = c2AmpFitVal;
    c2ctrlTau_GenAlg = c2TauFitVal; 
    c2ctrlAmp_x0 = c2ctrlAmp;
    c4ctrlAmp_x0 = c4ctrlAmp;
    c2ctrlTau_x0 = c2ctrlTau;
    else
    c2ctrlAmp_GenAlg = c2ctrlAmp;
    c4ctrlAmp_GenAlg = c4ctrlAmp;
    c2ctrlTau_GenAlg = c2ctrlTau;
    c2ctrlAmp_x0 = c2ctrlAmp;
    c4ctrlAmp_x0 = c4ctrlAmp;
    c2ctrlTau_x0 = c2ctrlTau;
    end 
    
[K, rates, tijs, Nstates, Nparam, param_strings, boundsArray, maxMutationArray,maxMutationCountsArray, minMutationCountsArray, atEdges, boundsArray_noAdj, model_name,dbCell,rates_str_woDB,rates_idx_woDB,rates_idx_DB] = libraryLoader(loadMode, model_set, modelNum , Nstates, tijs, c2ctrlAmp, c2ctrlTau, c4ctrlAmp, K, model_name,dbCell,rates_str_woDB,rates_idx_woDB,rates_idx_DB,Nparam,param_strings, boundsArray, maxMutationArray,maxMutationCountsArray, minMutationCountsArray, atEdges, boundsArray_noAdj,model_lin,model_loop1,model_loopN,model_brch) ;   
randGuess=zeros(1,Nparam); 
 % now with the proper boundsArray fill a guess and use as the initial
% points for PatSea 
    for param_idx = 1:Nparam
            %To pick a random number in the interval of LB to UB:
            % num = LB + rand*(UB - LB); %If rand = 0 then num = LB. If rand = 1, then num = UB.
%             population(:,param_idx) = boundsArray(param_idx,1) + population(:,param_idx)*(boundsArray(param_idx,2) - boundsArray(param_idx,1));
            [valuesArray,finalLogArray] = logRandUniform(boundsArray(param_idx,1),boundsArray(param_idx,2),pointsPerDecade,1); 
            randGuess(1,param_idx) = valuesArray; 
    end
    Nparam_ctrlAdj = Nparam - NctrlParams;
    tijs = randGuess(1,1:(Nparam_ctrlAdj - 2*Nstates));
    A = randGuess(1,(Nparam_ctrlAdj - 2*Nstates+1):Nparam_ctrlAdj-Nstates);
    sigmas = randGuess(1,((Nparam_ctrlAdj-Nstates)+1:Nparam_ctrlAdj));
    
    tijs_GenAlg = tijs;
    A_GenAlg = A;
    sigmaA_GenAlg = sigma_A;
    chisquared_GenAlg = chisquared; 

else

    tijs_GenAlg = tijs;
    A_GenAlg = A;
    sigmaA_GenAlg = sigma_A;
    chisquared_GenAlg = chisquared;
    if hardSetCtrlMode==1
    c2ctrlAmp_GenAlg = c2AmpFitVal;
    c4ctrlAmp_GenAlg = c2AmpFitVal;
    c2ctrlTau_GenAlg = c2TauFitVal; 
    c2ctrlAmp_x0 = c2ctrlAmp;
    c4ctrlAmp_x0 = c4ctrlAmp;
    c2ctrlTau_x0 = c2ctrlTau;
    else
    c2ctrlAmp_GenAlg = c2ctrlAmp;
    c4ctrlAmp_GenAlg = c4ctrlAmp;
    c2ctrlTau_GenAlg = c2ctrlTau;
    c2ctrlAmp_x0 = c2ctrlAmp;
    c4ctrlAmp_x0 = c4ctrlAmp;
    c2ctrlTau_x0 = c2ctrlTau;
    end
    
end
%     output_loc = [genAlgDataLocationPrefix(1:end-10) GenAlg_best_output_folder];
if altWeightScheme
     output_loc = [GenAlgFilesLocationPrefix filesep() programName '_outputs_v' num2str(patSeaWeightNum)];   
else
    output_loc = [GenAlgFilesLocationPrefix filesep() programName '_outputs_v2'];
end

    if exist(output_loc,'dir')
        cd(output_loc)
    else
        mkdir(output_loc)
        cd(output_loc)
    end
    
    
    % See if there are pattern search results to start from
%     outputFolderName_test = [outputProgramName '_output'];
%     if datestr_mode == 1
%         outputFolderName_test = [outputFolderName_test '_' date];
%     end
%     

    if overwriteWeightsMode
    weightingFactor_AmphistNew = weightingFactor_Amphist;
    weightingFactor_C2New = weightingFactor_C2;
    weightingFactor_C4_t0New = weightingFactor_C4_t0;
    end

    PatSea_filename = ['BestFitResults_' programName '.mat'];
    if exist(PatSea_filename,'file') == 2
        load(PatSea_filename)
        disp(['... Loading best fit results from ' programName '_output folder'])
    end       
    
    if overwriteWeightsMode
    weightingFactor_Amphist = weightingFactor_AmphistNew;
    weightingFactor_C2 = weightingFactor_C2New;
    weightingFactor_C4_t0 = weightingFactor_C4_t0New;
    end  
    %//////////////////////////////////////////////////////////////////////////
    %--------------------------------------------------------------------------
    %% Determine the bounds of the paramaters of the model
    %--------------------------------------------------------------------------
    %//////////////////////////////////////////////////////////////////////////

%     modelName = 'ssDNA_3state_cyclical';
%     model_name = '5stateNoLoops';

%     [K, A, rates, tijs, Nstates, Nparam, sigma_A, param_strings, boundsArray, maxMutationArray,maxMutationCountsArray, minMutationCountsArray] = model_builder(tijs, A, modelName);
%     SelectedTimeRange = 'full';
    %     varySigmaMode = 1;

    %     paramAdjs_path='/Users/clairealbrecht/Dropbox/MarcusLab/Data/smData_Processed/S1S2/gp32_0p5uM/ChosenMolecules/genAlgFiles/outputs/GenAlg_Nstate_Brett_varySig_5stateNoLoops/ParameterAdjustments.mat';
    %     paramAdjs_path = '/Users/clairealbrecht/Dropbox/MarcusLab/Data/smData_Processed/S1S2/gp32_0p5uM/ChosenMolecules/genAlgFiles/outputs/GenAlg_Nstate_Brett_varySig_5stateNoLoops';
    %     paramAdjs_path = '/Users/clairealbrecht/Dropbox/MATLAB_programs/claire_programs/from_Jack/updatedgenalgcodes/outputs/GenAlg_Nstate_updated_linear3124';
%     paramAdjs_path = [GenAlgFilesLocationPrefix filesep() 'ParameterAdjustments.mat'];
%     if exist(paramAdjs_path,'file')
%         paramAdj_vars = load(paramAdjs_path, 'paramAdjs', 'histCell');
%         paramAdjs = paramAdj_vars.paramAdjs;
%         histCell = paramAdj_vars.histCell;
%     else
%     end

    %     [K, A, rates, tijs, Nstates, Nparam, sigma_A, param_strings, boundsArray, maxMutationArray,maxMutationCountsArray, minMutationCountsArray, atEdges, boundsArray_noAdj] = model_builderBI_fullRange(tijs, A, sigma_A, modelName,SelectedTimeRange, paramAdjs);
    loadMode=1;
    [K, rates, tijs, Nstates, Nparam, param_strings, boundsArray, maxMutationArray,maxMutationCountsArray, minMutationCountsArray, atEdges, boundsArray_noAdj, model_name,dbCell,rates_str_woDB,rates_idx_woDB,rates_idx_DB] = libraryLoader(loadMode, model_set, modelNum , Nstates, tijs, c2ctrlAmp, c2ctrlTau, c4ctrlAmp, K, model_name,dbCell,rates_str_woDB,rates_idx_woDB,rates_idx_DB,Nparam,param_strings, boundsArray, maxMutationArray,maxMutationCountsArray, minMutationCountsArray, atEdges, boundsArray_noAdj,model_lin,model_loop1,model_loopN,model_brch) ;   
    loadMode=0;
    % define separate variables for the GenAlg outputs to compare later on

%     boundsArray = [boundsArray ; ctrlScaling ; ctrlTime]
    LBarray = boundsArray(:,1);
    UBarray = boundsArray(:,2);
    
    %//////////////////////////////////////////////////////////////////////////
    %--------------------------------------------------------------------------
    %% PART 4: Use the optimization algorithm patternsearch
    %--------------------------------------------------------------------------
    %//////////////////////////////////////////////////////////////////////////

    fun = @multigoaltcf_Nstate_polz_v4_ctrl_ModGen;
    tijs = reshape(tijs,[1, length(tijs)]);
    A = reshape(A, [1, length(A)]);
    sigma_A = reshape(sigma_A, [1,length(sigma_A)]);
    x0 = [tijs, A, sigma_A, ...
        c2ctrlAmp_x0, c2ctrlTau_x0, c4ctrlAmp_x0];
    disp(['x0: ' num2str(x0)]);

    Nparam = length(x0);
    global SelectedTimeRange  paramAdjs varySigmaMode param_strings weightHistFunc
    varySigmaMode =1;
%     [chisquared, chisquared_array, chisquared_unweighted, chisquared_unweighted_array] = multigoaltcf_Nstate(x0);
    [chisquared, chisquared_array, chisquared_unweighted, chisquared_unweighted_array] = multigoaltcf_Nstate_polz_v4_ctrl_ModGen(x0);
 
    LinearInequalityConstraintsMatrix = [];
    LinearInequalityConstraintsRealVector = [];

    a = [];
    b = [];

    aeq = [];
    beq = [];
    lb = LBarray;
    ub = UBarray;
    
    hold off;
    

    if showProgressOnFit_mode == 1
        nonlcon = [];
        figure(4);
        options = optimoptions('patternsearch','Display','iter','PlotFcn',{@psplotbestf, @psplotfuncount},'MaxIterations',MaxIterations,'CompletePoll','on','PollMethod','GPSPositiveBasis2N','MaxFunEvals',25); %10000);
        %         [x,fval] = patternsearch(fun,x0,A,b,Aeq,beq,lb,ub,nonlcon,options);
        [x,fval,exitflag,output] = patternsearch(fun,x0,a,b,aeq,beq,lb,ub,nonlcon,options);
    else
        [x,fval] = patternsearch(fun,x0,a,b,aeq,beq,lb,ub);
    end
    
    %//////////////////////////////////////////////////////////////////////////
    %--------------------------------------------------------------------------
    %% PART 5: Display the Final values to the user
    %--------------------------------------------------------------------------
    %//////////////////////////////////////////////////////////////////////////

%         tijs = x(1:(Nparam - Nstates));
%         A = x((Nparam - Nstates +1):Nparam);

        [chisquared, chisquared_array, chisquared_unweighted, chisquared_unweighted_array] = multigoaltcf_Nstate_polz_v4_ctrl_ModGen(x);
        
%         disp(['param_strings size: ',num2str(size(param_strings))])
        disp(['    Final values at the end of the GenAlg: (before patternsearch opt)'])
        params_GA = [tijs_GenAlg, A_GenAlg, sigmaA_GenAlg, chisquared_GenAlg,c2ctrlAmp_GenAlg,c4ctrlAmp_GenAlg,c2ctrlTau_GenAlg];
        param_strings2 = [param_strings,{'chi_sq'}];
        for i = 1:length(params_GA)
            fprintf(['    ',param_strings2{i},' = %f \n'], params_GA(i));
        end
        
                Nparam_ctrlAdj = Nparam - NctrlParams;
                tijs = x(1:(Nparam_ctrlAdj - 2*Nstates));
                A = x((Nparam_ctrlAdj - 2*Nstates+1):Nparam_ctrlAdj-Nstates);
                sigmas = x(((Nparam_ctrlAdj-Nstates)+1:Nparam_ctrlAdj));

                % define control parameters from guesses
                ctrlParams = x(Nparam_ctrlAdj+1:Nparam);
                c2CtrlAmp = ctrlParams(1);
                C2C4_ctrlTime = ctrlParams(2); 
                c4CtrlAmp = ctrlParams(3);
                
                c2Control = c2Exp(C2_exp_x,c2CtrlAmp,C2C4_ctrlTime)*C2_exp_y(1);

                c2CorrectedSurf = C2_exp_y - c2Control;
                yoff = abs(mean(c2CorrectedSurf(end-6:end)));
                
                c4Control = c4Exp(C4_tau1range,c4CtrlAmp,C2C4_ctrlTime)*C4_tau2eq0_exp(1,1);
                c4CorrectedSurf = C4_tau2eq0_exp - c4Control;
                zoff = abs(mean(mean(c4CorrectedSurf(end-6:end,end-6:end))));
%         tijs = x(1:(Nparam - 2*Nstates));
%         A = x((Nparam - 2*Nstates +1):Nparam-Nstates);
%         sigma_A = x(Nparam-Nstates+1:Nparam); 
%         c2CtrlAmp=x(end-2);
%         c4CtrlAmp=x(end-1);
%         C2C4_ctrlTime=x(end);

        %         [K, A, rates, tijs, ~, ~, ~, ~, ~, ~,~, ~] = model_builder(tijs, A, modelName);
%         [K, A, rates, tijs, Nstates, Nparam, sigma_A, ~, ~, ~,~, ~, ~, ~] = model_builderBI_fullRange(tijs, A, sigma_A, modelName,SelectedTimeRange, paramAdjs);
% [pass the GenAlg values for the amps and times into model_builder() so the bounds for patSea are constant from the start        
[K, rates, tijs, Nstates, Nparam, param_strings, boundsArray, maxMutationArray,maxMutationCountsArray, minMutationCountsArray, atEdges, boundsArray_noAdj, model_name,dbCell,rates_str_woDB,rates_idx_woDB,rates_idx_DB] = libraryLoader(loadMode, model_set, modelNum , Nstates, tijs, c2ctrlAmp, c2ctrlTau, c4ctrlAmp, K, model_name,dbCell,rates_str_woDB,rates_idx_woDB,rates_idx_DB,Nparam,param_strings, boundsArray, maxMutationArray,maxMutationCountsArray, minMutationCountsArray, atEdges, boundsArray_noAdj,model_lin,model_loop1,model_loopN,model_brch) ;   

    disp(['Chisquared after patternsearch: X^2(x) = ' num2str(chisquared)]);
    
    %Calculate the rawchisquared value
    rawchisquared = chisquared_array(1)*fitHistMode/weightingFactor_Amphist + ...
        chisquared_array(2)*fitC2Mode/weightingFactor_C2+...
        chisquared_array(3)*fitC4Mode/weightingFactor_C4_t0;
    disp(['The unweighted chisquared value is ' num2str(rawchisquared)]);
    
    chisquared_Weighted_array = zeros(size(chisquared_array));
        chisquared_Weighted_array(1) = chisquared_array(1)*weightingFactor_Amphist;
        chisquared_Weighted_array(2) = chisquared_array(2)*weightingFactor_C2;
        chisquared_Weighted_array(3) = chisquared_array(3)*weightingFactor_C4_t0;
    
    [P, V, K, C2_exp_x] = K2P(K,C2_exp_x);
    
    [Peq, denom_hist_sim, hist_sim, Peq_names] = histMaker_Nstate_v3normCorrect(P, A, sigma_A, Amp_bins);
   
    
        disp(['    Final values at the end of the patternsearch:'])
        
        params = [tijs, A, sigma_A, chisquared];
        param_strings2 = [param_strings,{'chi_sq'}];
        for i = 1:length(params) %Nparam+1
            fprintf(['    ',param_strings2{i},' = %f \n'], params(i));
        end

    %--------------------------------------------------------------------------
    %% Save the data (PART 5: Final Values)
    %--------------------------------------------------------------------------
    
    %--------------------------------------------------------------------------
    % 5.1: Figure out where to save the data
    %--------------------------------------------------------------------------
    if saveMode == 1
        clear('outputFolderName')
        foutName = ['BestFitResults_' programName '.mat'];
        
%         outputFolderName = [outputProgramName '_output'];
        if altWeightScheme
        outputFolderName = [programName '_outputs_v' num2str(patSeaWeightNum)];   
        else
        outputFolderName = [programName '_outputs_v2'];
        end
        if datestr_mode == 1
            outputFolderName = [outputFolderName '_' date];
        end
        disp(['Checking for existance of ' outputFolderName]);
        if exist(outputFolderName,'dir') ~= 7
            mkdir(outputFolderName);
            disp(['Making a folder called "' outputFolderName '" to hold the output']);
        else
            startFold=pwd;
%             cd(outputFolderName);
            listOutputs = dir('output*');
            numOutputsCur=length(listOutputs);
            cd(startFold); 
        end      
       
        
        curOutputDir=['output' num2str(numOutputsCur+1)];
        if exist(curOutputDir,'dir') ~= 7
            mkdir(curOutputDir);
        end
        bestfitFilePath = ['output' num2str(numOutputsCur+1) filesep() foutName];
        if exist(bestfitFilePath,'file') ~= 2 %If the best fit DNE
            disp(['Saving the results for the first time as "'  'output' num2str(numOutputsCur+1) filesep() foutName '". ChiSqaure = ' num2str(chisquared)]);
            saveFitMode = 1;
        elseif exist(['output' num2str(numOutputsCur+1) filesep() foutName],'file') == 2
            current_LowestUWChisquared = chisquared_unweighted;
            load(['output' num2str(numOutputsCur+1) filesep() foutName],'chisquared_unweighted');
            
            prev_LowestUWChisquared = chisquared_unweighted;%Save the previous X2 as prev_LowestChisquared
            chisquared_unweighted = current_LowestUWChisquared;%reset chisquared to its current value
            
            if current_LowestUWChisquared < prev_LowestUWChisquared
                %Just found a new best fit. Need to save it.\
                fprintf('Congratulations! %f < %f (#1 fit). Will update best fit for %s .\r\n',current_LowestUWChisquared,prev_LowestUWChisquared,constructName)
                saveFitMode = 1;
            elseif current_LowestUWChisquared == prev_LowestUWChisquared
                fprintf('Nothing changed.%f = %f. \r\n',current_LowestUWChisquared,prev_LowestUWChisquared);
                saveMode = 0;
                saveFitMode = 0;
                updateExcelMode = 0;
                
            elseif current_LowestUWChisquared > prev_LowestUWChisquared
                fprintf('Oh no! %f > %f (#1 fit).\n',current_LowestUWChisquared,prev_LowestUWChisquared)
                saveMode = 0;
                saveFitMode = 0;
                updateExcelMode = 0;
            end
        end
        
        if saveFitMode == 1
%             if exist(outputFolderName, 'dir') ~= 7
%                 disp(['Making folder ' outputFolderName ' to hold output']);
%                 mkdir(outputFolderName);
%             end
            
%             if exist([outputFolderName filesep() foutName],'file') == 1
%              if exist( foutName,'file') == 1
%                 save(foutName, 'x','model_name','protein_str', 'Nstates',...
%                     'rates','tijs','A','K','P','sigma_A','Amp_bins','yoff','zoff',...
%                     'chisquared','chisquared_unweighted','weightingFactor_Amphist',...
%                     'weightingFactor_C2','weightingFactor_C4_t0','-append');
%             else
                save(bestfitFilePath, 'x','model_name','protein_str','Nstates',...
                    'rates','tijs','A','K','P','sigma_A','Amp_bins','yoff','zoff',...
                    'chisquared','chisquared_unweighted','weightingFactor_Amphist',...
                    'weightingFactor_C2','weightingFactor_C4_t0','c2CtrlAmp','c4CtrlAmp','C2C4_ctrlTime', ...
                    'wC4func','weightC2func','chisquared_Weighted_array');
%             end
        end
        
        if saveMode == 1
            disp(['Will save the data in ' outputFolderName filesep() 'output' num2str(numOutputsCur+1) filesep() foutName]);
            
%             foutName = [outputFolderName filesep() 'fitInputData.mat'];
            foutName = ['fitInputData.mat'];
            fitInputFilePath = ['output' num2str(numOutputsCur+1) filesep() foutName];

            if fitHistMode == 1
%                 save(foutName,'histogram_FileName','FRET_bins','targetHistogram')
                save(fitInputFilePath,'histogram_FilePath','Amp_bins','targetHistogram')
            end
            if fitC2Mode == 1
                if exist(fitInputFilePath,'file') == 2
%                     save(foutName,'C2_FileName','C2_exp_x','C2_exp_y','yoff','weightC2func','-append');
                    save(fitInputFilePath,'C2_FilePath','C2_exp_x','C2_exp_y','yoff','weightC2func','-append');
                
                else
%                     save(foutName,'C2_FileName','C2_exp_x','C2_exp_y','yoff','weightC2func');
                    save(fitInputFilePath,'C2_FilePath','C2_exp_x','C2_exp_y','yoff','weightC2func');
                end
                
            end
            if fitC4Mode == 1
                if exist(fitInputFilePath,'file') == 2
%                     save(foutName,'C4_FileName','C4_tau2eq0_exp','C4_tau1range','C4_tau3range','wC4func','zoff','-append');
                    save(fitInputFilePath,'C4_FilePath','C4_tau2eq0_exp','C4_tau1range','C4_tau3range','wC4func','zoff','-append');
                else
%                     save(foutName,'C4_FileName','C4_tau2eq0_exp','C4_tau1range','C4_tau3range','wC4func','zoff');
                    save(fitInputFilePath,'C4_FilePath','C4_tau2eq0_exp','C4_tau1range','C4_tau3range','wC4func','zoff');
                end
            end
        
        end

        
%         
%         if showBestFitMode == 1
%             disp('     Loading the best fit from all runs.');
%             if model == 8
%                 load([outputFolderName filesep() foutName],...
%                     'A1','A2','A3','A4','A5','A6','A7','A8',...
%                     't12','t13','t21','t23','t31','t32',...
%                     't24','t42','t25','t52','t56','t65','t36','t63','t37','t73','t68','t86',...
%                     'chisquared','chisquared_array',... %'rawchisquared',...
%                     'k12','k13','k21','k23','k31','k32',...
%                     'k24','k42','k25','k52','k56','k65','k36','k63','k37','k73','k68','k86');
%             elseif model == 9
%                 load([outputFolderName filesep() foutName],...
%                     'A1','A2','A3','A4','A5','A6','A7','A8','A9',...
%                     't12','t13','t21','t23','t31','t32',...
%                     't24','t42','t25','t52','t56','t65','t36','t63','t37','t73','t68','t86',...
%                     't69','t96','t89','t98',...
%                     'chisquared','chisquared_array',... %'rawchisquared',...
%                     'k12','k13','k21','k23','k31','k32',...
%                     'k24','k42','k25','k52','k56','k65','k36','k63','k37','k73','k68','k86',...
%                     'k69','k96','k89','k98');
%             end
%             saveMode = saveModeDefault;
%             updateExcelMode = updateExcelModeDefault;
%             disp('showBestFitMode: Turning save mode on');
%         end
    end
    
    
    %--------------------------------------------------------------------------
    %% Write Data to excel file
    %--------------------------------------------------------------------------
%     if saveMode == 1
%         if updateExcelMode == 1
%             excelOutFolder = [DropboxLocationPrefix filesep() 'BrettsProjects' filesep() 'SSB_gp32'];
%             excelFileName = 'SSBgp32_numericalResults.xlsx';
%             %         excelFileName = 'Results.xlsx';
%             excelFilePath = [excelOutFolder filesep() excelFileName];
%             if model == 8
%                 sheetOUT = 'TimeResults_8state';
%             elseif model == 9
%                 sheetOUT = 'TimeResults_9state';
%             end
%             disp(['Will update the excel file: ' excelFilePath ...
%                 ', and excel sheet: ' sheetOUT]);
%             
%             if exist(excelFilePath,'file') == 2
%                 if strcmp(constructName,'S1S2') == 1
%                     if model == 8
%                         rangeOUT = 'D13:AE13';
%                     elseif model == 9
%                         rangeOUT = 'D13:AJ13';
%                     end
%                 elseif strcmp(constructName,'S18S2') == 1
%                     if model == 8
%                         rangeOUT = 'D14:AE14';
%                     elseif model == 9
%                         rangeOUT = 'D14:AJ14';
%                     end
%                 elseif strcmp(constructName,'S4S5') == 1
%                     if model == 8
%                         rangeOUT = 'D15:AE15';
%                     elseif model == 9
%                         rangeOUT = 'D15:AJ15';
%                     end
%                 elseif strcmp(constructName,'S19S5') == 1
%                     if model == 8
%                         rangeOUT = 'D16:AE16';
%                     elseif model == 9
%                         rangeOUT = 'D16:AJ16';
%                     end
%                 end
%                 
%                 if model == 8
%                     outCell = {t12,t13,t21,t23,t31,t32,...
%                         t24,t42,t25,t52,t56,t65,t36,t63,t37,t73,t68,t86,...
%                         A1,A2,A3,A4,A5,A6,A7,A8,chisquared,rawchisquared};
%                 elseif model == 9
%                     outCell = {t12,t13,t21,t23,t31,t32,...
%                         t24,t42,t25,t52,t56,t65,t36,t63,t37,t73,t68,t86,...
%                         t69,t96,t89,t98,...
%                         A1,A2,A3,A4,A5,A6,A7,A8,A9,chisquared,rawchisquared};
%                 end
%                 writecell(outCell,excelFilePath,'Sheet',sheetOUT,'Range',rangeOUT);
%                 disp(['Updated excel file ' excelFileName ' with the fits']);
%             else
%                 disp('Coud not find excel file to update');
%             end
%         end
%     end
    
    %--------------------------------------------------------------------------
    %% Save the data PLOTS     (PART 5: Final Values)
    %--------------------------------------------------------------------------
    if plotMode == 1
        
        %//////////////////////////////////////////////////////////////////////////
        % Plot final results of the algorithm   (Part 5: Final Values)
        %//////////////////////////////////////////////////////////////////////////
        %     if fitHistMode == 1
        %--------------------------------------------------------------------------
        % (1) Final Histogram Result
        %--------------------------------------------------------------------------
        
        figure(1);
        clf;
                                   

        lineStyle = char('g--','c--', 'm--','y--','k--','r-.','g-.','b-.','c-.','m-.','y-.','k-.'); % Define a list of colors to loop over
                            
        set(gcf,'Color','w');
        set(gcf,'Name','Amp Histogram');
        %         lineStyle = char('r--','g--','b--','c--', 'm--','y--','k--','r-.','g-.','b-.','c-.','m-.','y-.','k-.'); % Define a list of colors to loop over
%         lineColor = char('g','#D95319','c','m','#EDB120','k','#7E2F8E','#A2142F','#4DBEEE');
        lineColor = char('g','b','c','m','y','k','orange','A2142F','4DBEEE');

        if exist('data_hist_Plot','var') == 1
            delete(data_hist_Plot)
        end
        % Replot data Histogram
        data_hist_Plot = plot(Amp_bins,targetHistogram);
        data_hist_Plot.LineWidth = 2;
        data_hist_Plot.Color = 'blue';
        %                             x.DisplayName = 'Data';
        xlabel('Visibility','FontSize',14);
        ylabel('Frequency','FontSize',14);
        %         [sample_description, ~] = sample_descriptionGetter();
        sample_description=[constructName ' ' saltConditions];
        title_str = ['Experimental vs Simulated Histograms' ...
            10 sample_description];
        title(title_str,'fontsize',14);
        hold on;
        % Clear plots from previous run
        if exist('histPlot','var') == 1
            delete(histPlot)
            %                                 disp('histPlot deleted')
        end
        
        [Peq, denom_hist_sim, hist_sim, Peq_names] = histMaker_Nstate_v3normCorrect(P, A, sigma_A, Amp_bins);

        for i = 1:length(Peq)
            %                                 Peq(i)/(sqrt(2*pi)*sigma_A(i))
            %                                 histPlot = plot(Amp_bins, Peq(i)/(sqrt(2*pi)*sigma_A(i))*exp(-((Amp_bins - A(i))/sigma_A(i)).^2)./denom_hist_sim,...
            %
            histPlot = plot(Amp_bins,(Peq(i)/(sqrt(2*pi)*sigma_A(i)))*exp(-((Amp_bins - A(i))/(sqrt(2)*sigma_A(i))).^2)./denom_hist_sim,...
                lineStyle(i,:),'LineWidth',1,'DisplayName',[Peq_names(i,:) ' = ' num2str(Peq(i))]);
            %                                 lgd = legend('show');
            % lgd.Location = 'northwest';
            
            hold on
        end
        
        %         N = length(K);
%         % Plot histogram of each state
%         for state_idx = 1:N
%             histPlot = plot(Amp_bins, Peq(state_idx)/(sqrt(2*pi)*sigma_A(state_idx))* exp(-((Amp_bins - A(state_idx))/sigma_A(state_idx)).^2)./denom_hist_sim,...
%                 'color',deblank(lineColor(state_idx,:)),'LineStyle','--','LineWidth',1,'DisplayName',Peq_names(state_idx,:));
%             lgd = legend('show');
%             lgd.Location = 'northwest';
%             lgd.FontSize = 14;
%             hold on
%         end
%         set(gca, 'FontSize', 14);
%         if exist('histPlotTot','var') == 1
%             delete(histPlotTot)
%         end
%         hold on
%         
        % Plot overall histogram
        histPlotTot = plot(Amp_bins, hist_sim,...
            'color','r','LineWidth',2,'DisplayName','Total Fit');
        lgd_tot = legend('show');
        
        if saveMode == 1
           % saveas(gcf,['output' num2str(numOutputsCur+1) filesep() 'hist_bestFit.fig'])
        end
        
        %% --------------------------------------------------------------------------
        % (1) Final C2 Result
        %--------------------------------------------------------------------------
        if fitC2Mode == 1
            [C2,C2_exp_x] = PA2C2(P,A,C2_exp_x,yoff,addControlMode,c2Control);
            
            if plotMode == 1
                figure(2)
                clf;
                hold on;
                set(gcf,'Color','w');
                
                if exist('C2_plot','var')  == 1
                    delete(C2_plot)
                end
                % Plot data C2
                plot(C2_exp_x,C2_exp_y,'b.','MarkerSize',10,'DisplayName','C^{(2)}(\tau) Data');
                
                %                 [sample_description, ~] = sample_descriptionGetter();
                sample_description=[constructName ' ' saltConditions];
                title_str = ['Experimental vs (Best fit) Simulated C2' ...
                    10 sample_description];
                title(title_str,'fontsize',14);
                
                drawnow();
                hold on;
                C2_plot = plot(C2_exp_x,C2,'color','r','LineWidth',2.5);
                
                xlabel('\tau (sec)','fontsize',16);
                ylabel('C^{(2)}(\tau)','fontsize',16);
                set(gca,'yscale','linear');
                set(gca,'xscale','log');
                set(gca,'FontSize',14);
                grid on
                axis tight;
                
                drawnow();
            end
            
            if saveMode == 1
               % saveas(gcf,['output' num2str(numOutputsCur+1) filesep() 'C2_bestFit.fig'])
            end
            
        end
        
        %% --------------------------------------------------------------------------
        % (1) Final C4 Result
        %--------------------------------------------------------------------------
        if fitC4Mode == 1
            tau2 = 0;
            [C4_sim,time] = PAK2C4(P,A,K,time,tau2,zoff,addControlMode,c4Control);
            
            if plotMode == 1
                figure(3);
                clf;
                hold on;
                if exist('C4_plot','var')  == 1
                    delete(C4_plot)
                end
                
                C4dataPlot = mesh(C4_tau1range,C4_tau3range,C4_tau2eq0_exp);
                C4dataPlot.DisplayName = 'C^{(4)} Data';
                set(gcf,'Color','w');
                hold on;
                C4_plot = surf(time,time,C4_sim);
                title_str = ['Four point time correlation function'];
                title(title_str,'FontSize',18);
                xlabel('\tau_1 (sec)','fontsize',16);
                ylabel('\tau_3 (sec)','fontsize',16);
                zlabel('C^{(4)}(\tau)','fontsize',16);
                grid on;
                set(gca,'FontSize',12);
                axis tight;
                axis square;
                colormap jet;
                colorbar;
                set(gca,'xscale','log');
                set(gca,'yscale','log');         
                view(55.3,33.2);  
                
                drawnow();
                
            end
            
            if saveMode == 1
               % saveas(gcf,['output' num2str(numOutputsCur+1) filesep() 'C4_bestFit.fig'])
            end
        end
        
        %% --------------------------------------------------------------------------
        % (4) Display a graphic of the network        (Part 5: Final Values)
        %--------------------------------------------------------------------------
%         figure(4);
%         if model == 8
%             networkPlotter_8stateModel1(sample_description,chisquared,...
%                 A1,A2,A3,A4,A5,A6,A7,A8,...
%                 t12,t13,t21,t23,t24,t25,t31,t32,t36,t37,t42,t52,t56,t63,t65,t68,t73,t86)
%         elseif model == 9
%             networkPlotter_9stateModel(sample_description,chisquared,...
%                 A1,A2,A3,A4,A5,A6,A7,A8,A9,...
%                 t12,t13,t21,t23,t24,t25,t31,t32,t36,t37,t73,t42,t52,t56,t63,t65,t68,t86,t69,t96,t89,t98)
%         end
%         
        if pauseBetweenConstructsMode == 1
            disp('Press any button to continue');
            pause();
        end
        
        disp([programName ' complete. All ' num2str(NconstructFolderNames) ' constructs completed.']);
        
    end
end%construct_idx loop

%% *********************** Internal Function Calls ***********************
%% Calculate the rate matrix
% function K = rates2K_5stateModel(k12,k21,k23,k24,k32,k42,k45,k54);
% K = [-(k12),    k21,             0,            0,         0;...
%     k12    ,  -(k21+k23+ k24), k32,          k42,         0;....
%     0      ,    k23,        -(k32),            0,         0;...
%     0      ,    k24,             0, -(k42 + k54),       k45;...
%     0      ,      0,             0,          k54,    -(k45);];
% end
% 
% % ********* 8 state *********
% function K = rates2K_8stateModel1(k12,k13,k21,k23,k24,k25,k31,k32,k36,k37,k42,k52,k56,k63,k65,k68,k73,k86)
% K = [-(k12+k13),  k21, k31, 0, 0, 0, 0, 0;...
%     k12, -(k21+k23+k24+k25), k32, k42, k52, 0, 0, 0;...
%     k13, k23, -(k31+k32+k36+k37), 0, 0, k63, k73, 0;...
%     0, k24, 0, -k42, 0, 0, 0, 0;...
%     0, k25, 0, 0, -(k52+k56), k65, 0, 0;...
%     0, 0, k36, 0, k56, -(k65+k63+k68), 0, k86;...
%     0, 0, k37, 0, 0, 0, -k73, 0;...
%     0, 0, 0, 0, 0, k68, 0, -k86;];
% end

% ********* 9 state *********
% function K = rates2K_9stateModel(k12,k13,k21,k23,k24,k25,k31,k32,k36,k37,k42,k52,k56,k63,k65,k68,k73,k86,k89,k69,k96,k98)
% K = [-(k12+k13),  k21, k31, 0, 0, 0, 0, 0, 0;...
%     k12, -(k21+k23+k24+k25), k32, k42, k52, 0, 0, 0, 0;...
%     k13, k23, -(k31+k32+k36+k37), 0, 0, k63, k73, 0, 0;...
%     0, k24, 0, -k42, 0, 0, 0, 0, 0;...
%     0, k25, 0, 0, -(k52+k56), k65, 0, 0, 0;...
%     0, 0, k36, 0, k56, -(k65+k63+k68+k69), 0, k86, k96;...
%     0, 0, k37, 0, 0, 0, -k73, 0, 0;...
%     0, 0, 0, 0, 0, k68, 0, -(k86+k89),k98;...
%     0, 0, 0, 0, 0, k69, 0, k89, -(k96+k98)];
% end

%% Convert guess to paramaters
% ********* 8 state *********
% function   [t12,t13,t21,t23,t31,t24,t42,t25,t52,t56,t65,t36,t37,t73,t68,t86,...
%     A1,A2,A3,A4,A5,A6,A7,A8,chisquared,chisquared_unweighted] = guess2param_8state(guessIN)
% 
% % invRates = [t12,t13,t21,t23,t31,t24,t42,t25,t52,t56,t65,t36,t37,t73,t68,t86];
% % paramNum = [  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15, 16];
% % derivedParam = [k63,t63,k32,t32];
% 
% % Display best guess to screen
% t12 = guessIN(1);
% t13 = guessIN(2);
% t21 = guessIN(3);
% t23 = guessIN(4);
% t31 = guessIN(5);
% % * t32 will be determined by others in loop
% t24 = guessIN(6);
% t42 = guessIN(7);
% t25 = guessIN(8);
% t52 = guessIN(9);
% t56 = guessIN(10);
% t65 = guessIN(11);
% % * t63 will be determined by others in loop
% t36 = guessIN(12);
% t37 = guessIN(13);
% t73 = guessIN(14);
% t68 = guessIN(15);
% t86 = guessIN(16);
% 
% A1 = guessIN(17);        %21);
% A2 = guessIN(18);        %22);
% A3 = guessIN(19);        %23);
% A4 = guessIN(20);        %24);
% A5 = guessIN(21);        %25);
% A6 = guessIN(22);        %26);
% A7 = guessIN(23);        %27);
% A8 = guessIN(24);
% 
% chisquared = guessIN(25);
% chisquared_unweighted = guessIN(26);
% end%end of guess2param
% 
% % ********* 9 state *********
% function   [t12,t13,t21,t23,t31,t24,t42,t25,t52,t56,t65,t36,t37,t73,t68,t86,t89, t98, t69,...
%     A1,A2,A3,A4,A5,A6,A7,A8,A9,chisquared,chisquared_unweighted] = guess2param_9state(guessIN)
% 
% % invRates = [t12,t13,t21,t23,t31,t24,t42,t25,t52,t56,t65,t36,t37,t73,t68,t86];
% % paramNum = [  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15, 16];
% % derivedParam = [k63,t63,k32,t32];
% 
% % Display best guess to screen
% t12 = guessIN(1);
% t13 = guessIN(2);
% t21 = guessIN(3);
% t23 = guessIN(4);
% t31 = guessIN(5);
% % * t32 will be determined by others in loop
% t24 = guessIN(6);
% t42 = guessIN(7);
% t25 = guessIN(8);
% t52 = guessIN(9);
% t56 = guessIN(10);
% t65 = guessIN(11);
% % * t63 will be determined by others in loop
% t36 = guessIN(12);
% t37 = guessIN(13);
% t73 = guessIN(14);
% t68 = guessIN(15);
% t86 = guessIN(16);
% t89 = guessIN(17);
% t98 = guessIN(18);
% % * t96 will be determined by others in loop
% t69 = guessIN(19);
% 
% A1 = guessIN(20);        %21);
% A2 = guessIN(21);        %22);
% A3 = guessIN(22);        %23);
% A4 = guessIN(23);        %24);
% A5 = guessIN(24);        %25);
% A6 = guessIN(25);        %26);
% A7 = guessIN(26);        %27);
% A8 = guessIN(27);
% A9 = guessIN(28);
% 
% chisquared = guessIN(29);
% chisquared_unweighted = guessIN(30);
% end%end of guess2param_9state
% 
% %%  Convert times to rates
% % ********* 8 state *********
% function [k12,k21,k23,k24,k32,k42,k45,k54]...
%     = times2rates_5stateModel(t12,t21,t23,t24,t32,t42,t45,t54)
% 
% k12 = 1/t12;
% k21 = 1/t21;
% k23 = 1/t23;
% k24 = 1/t24;
% k32 = 1/t32;
% k42 = 1/t42;
% k45 = 1/t45;
% k54 = 1/t54;
% 
% end%end of times2rates_5stateModel
% 
% % ********* 8 state *********
% function [t32,t63,k12,k13,k21,k23,k24,k25,k31,k32,k36,k37,k42,k52,k56,k63,k65,k68,k73,k86]...
%     = times2rates_8stateModel1(t12,t13,t21,t23,t24,t25,t31,t36,t37,t42,t52,t56,t65,t68,t73,t86)
% 
% k12 = 1/t12;
% k13 = 1/t13;
% k21 = 1/t21;
% k23 = 1/t23;
% k31 = 1/t31;
% k24 = 1/t24;
% k42 = 1/t42;
% k25 = 1/t25;
% k52 = 1/t52;
% k56 = 1/t56;
% k65 = 1/t65;
% k36 = 1/t36;
% k37 = 1/t37;
% k73 = 1/t73;
% k68 = 1/t68;
% k86 = 1/t86;
% 
% % Loop conditions
% k32 = (k12 * k23 * k31)/(k13 * k21);
% k63 = (k23 * k36 * k65* k52)/(k32 * k25 * k56);
% t32 = 1/k32;
% t63 = 1/k63;
% end%end of times2rates_8stateModel1
% 
% % ********* 9 state *********
% function [t32,t63,t96,k12,k13,k21,k23,k24,k25,k31,k32,k36,k37,k42,k52,k56,k63,k65,k68,k73,k86,k89,k98,k96,k69] = ...
%     times2rates_9stateModel(t12,t13,t21,t23,t24,t25,t31,t36,t37,t42,t52,t56,t65,t68,t73,t86,t89,t98,t69)
% 
% k12 = 1/t12;
% k13 = 1/t13;
% k21 = 1/t21;
% k23 = 1/t23;
% k31 = 1/t31;
% k24 = 1/t24;
% k42 = 1/t42;
% k25 = 1/t25;
% k52 = 1/t52;
% k56 = 1/t56;
% k65 = 1/t65;
% k36 = 1/t36;
% k37 = 1/t37;
% k73 = 1/t73;
% k68 = 1/t68;
% k86 = 1/t86;
% 
% k69 = 1/t69;
% k89 = 1/t89;
% k98 = 1/t98;
% 
% % Loop conditions
% k32 = (k12 * k23 * k31)/(k13 * k21);
% k63 = (k23 * k36 * k65* k52)/(k32 * k25 * k56);
% k96 = (k69 * k98 * k86)/(k68 * k89);
% t32 = 1/k32;
% t63 = 1/k63;
% t96 = 1/k96;
% end%end of times2rates_9stateModel1
% 
% 
% %% Plot a model of the 8 state network
% % ********* 8 state *********
% function networkPlotter_8stateModel1(sample_description,chisquared,...
%     A1,A2,A3,A4,A5,A6,A7,A8,...
%     t12,t13,t21,t23,t24,t25,t31,t32,t36,t37,t42,t52,t56,t63,t65,t68,t73,t86)
% 
% clf;
% set(gcf,'Name','Model: 8 state');
% set(gcf,'Color','white');
% xlim([0 4.1]);
% ylim([0 4.1]);
% 
% %Designate spots for the states
% state1_loc = [0.5 2];
% state2_loc = [1.5 2.5];
% state3_loc = [1.5 1.5];
% state4_loc = [2.5 3.5];
% state5_loc = [2.5 2.5];
% state6_loc = [2.5 1.5];
% state7_loc = [2.5 0.5];
% state8_loc = [3.5 1.5];
% % or for 9 state
% %              state8_loc = [4 2.5];
% %              state9_loc = [4 1.5];
% hold on;
% % text(0.1,4,['Gen# = ' num2str(genNum)],'FontSize',12);
% 
% % Plot the state symbols
% text(state1_loc(1)-0.04,state1_loc(2)+0.05,'1','FontSize',24)
% text(state2_loc(1),state2_loc(2),'2','FontSize',24);
% text(state3_loc(1)-0.05,state3_loc(2),'3','FontSize',24);
% text(state4_loc(1)-0.04,state4_loc(2)+0.05,'4','FontSize',24)
% text(state5_loc(1),state5_loc(2),'5','FontSize',24);
% text(state6_loc(1)-0.05,state6_loc(2),'6','FontSize',24);
% text(state7_loc(1)-0.04,state7_loc(2)+0.05,'7','FontSize',24)
% text(state8_loc(1),state8_loc(2),'8','FontSize',24);
% %              text(state9_loc(1),state9_loc(2),'9','FontSize',24);
% 
% 
% %Plot the FRET values
% text(state1_loc(1),state1_loc(2)-0.2,['=' num2str(A1,'%.2f')],'FontSize',16, 'Color','b');
% text(state2_loc(1),state2_loc(2)-0.2,['=' num2str(A2,'%.2f')],'FontSize',16, 'Color','b');
% text(state3_loc(1),state3_loc(2)-0.2,['=' num2str(A3,'%.2f')],'FontSize',16, 'Color','b');
% text(state4_loc(1),state4_loc(2)-0.2,['=' num2str(A4,'%.2f')],'FontSize',16, 'Color','b');
% text(state5_loc(1),state5_loc(2)-0.2,['=' num2str(A5,'%.2f')],'FontSize',16, 'Color','b');
% text(state6_loc(1),state6_loc(2)-0.2,['=' num2str(A6,'%.2f')],'FontSize',16, 'Color','b');
% text(state7_loc(1),state7_loc(2)-0.2,['=' num2str(A7,'%.2f')],'FontSize',16, 'Color','b');
% text(state8_loc(1),state8_loc(2)-0.2,['=' num2str(A8,'%.2f')],'FontSize',16, 'Color','b');
% %              text(state9_loc(1),state8_loc(2)-0.2,['=' num2str(A9,'%.2f')],'FontSize',16);
% 
% 
% %Plot Lines Between states
% line([state1_loc(1) state2_loc(1)],[state1_loc(2) state2_loc(2)],'Color','k','LineStyle','-');
% line([state2_loc(1) state3_loc(1)],[state2_loc(2) state3_loc(2)],'Color','k','LineStyle','-');
% line([state3_loc(1) state1_loc(1)],[state3_loc(2) state1_loc(2)],'Color','k','LineStyle','-');
% line([state2_loc(1) state4_loc(1)],[state2_loc(2) state4_loc(2)],'Color','k','LineStyle','-');
% line([state2_loc(1) state5_loc(1)],[state2_loc(2) state5_loc(2)],'Color','k','LineStyle','-');
% line([state5_loc(1) state6_loc(1)],[state5_loc(2) state6_loc(2)],'Color','k','LineStyle','-');
% line([state3_loc(1) state6_loc(1)],[state3_loc(2) state6_loc(2)],'Color','k','LineStyle','-');
% line([state3_loc(1) state7_loc(1)],[state3_loc(2) state7_loc(2)],'Color','k','LineStyle','-');
% line([state6_loc(1) state8_loc(1)],[state6_loc(2) state8_loc(2)],'Color','k','LineStyle','-');
% 
% %              % For 9 states
% %              line([state8_loc(1) state9_loc(1)],[state8_loc(2) state9_loc(2)],'Color','k','LineStyle','-');
% %              line([state6_loc(1) state9_loc(1)],[state6_loc(2) state9_loc(2)],'Color','k','LineStyle','-');
% %
% %Plot the inverse of the rates
% t12_loc = [state1_loc(1)+(state2_loc(1)- state1_loc(1))/2,state1_loc(2)+(state2_loc(2)- state1_loc(2))/2];
% t12_msg = ['t_{12} = ' num2str(round(t12,6)) 10 't_{21} = ' num2str(round(t21,6))];
% text(t12_loc(1),t12_loc(2),t12_msg,'FontSize',10);
% 
% t23_loc = [state3_loc(1)+(state2_loc(1)- state3_loc(1))/2,state3_loc(2)+(state2_loc(2)- state3_loc(2))/2];
% t23_msg = ['t_{23} = ' num2str(round(t23,6)) 10 't_{32} = ' num2str(round(t32,6))];
% text(t23_loc(1),t23_loc(2),t23_msg,'FontSize',10);
% 
% t31_loc = [state1_loc(1)+(state3_loc(1)- state1_loc(1))/2,state1_loc(2)+(state3_loc(2)- state1_loc(2))/2];
% t31_msg = ['t_{31} = ' num2str(round(t31,6)) 10 't_{13} = ' num2str(round(t13,6))];
% text(t31_loc(1),t31_loc(2),t31_msg,'FontSize',10);
% 
% t24_loc = [state2_loc(1)+abs((state2_loc(1)- state4_loc(1))/2),state2_loc(2)+abs((state2_loc(2)- state4_loc(2))/2)];
% t24_msg = ['t_{24} = ' num2str(round(t24,6)) 10 't_{42} = ' num2str(round(t42,6))];
% text(t24_loc(1),t24_loc(2),t24_msg,'FontSize',10);
% 
% t37_loc = [state7_loc(1) - abs((state3_loc(1) - state7_loc(1))/2),state7_loc(2)+ abs((state3_loc(2)- state7_loc(2))/2)];
% t37_msg = ['t_{37} = ' num2str(round(t37,6)) 10 't_{73} = ' num2str(round(t73,6))];
% text(t37_loc(1),t37_loc(2),t37_msg,'FontSize',10);
% 
% t25_loc = [state2_loc(1)+abs((state2_loc(1)- state5_loc(1))/2),state2_loc(2)+abs((state2_loc(2)- state5_loc(2))/2)];
% t25_msg = ['t_{25} = ' num2str(round(t25,6)) 10 't_{52} = ' num2str(round(t52,6))];
% text(t25_loc(1),t25_loc(2),t25_msg,'FontSize',10);
% 
% t36_loc = [state3_loc(1)+abs((state3_loc(1)- state6_loc(1))/2),state3_loc(2)+abs((state3_loc(2)- state6_loc(2))/2)];
% t36_msg = ['t_{36} = ' num2str(round(t36,6)) 10 't_{63} = ' num2str(round(t63,6))];
% text(t36_loc(1),t36_loc(2),t36_msg,'FontSize',10);
% 
% t56_loc = [state5_loc(1)+abs((state5_loc(1)- state6_loc(1))/2),state5_loc(2)-abs((state5_loc(2)- state6_loc(2))/2)];
% t56_msg = ['t_{56} = ' num2str(round(t56,6)) 10 't_{65} = ' num2str(round(t65,6))];
% text(t56_loc(1),t56_loc(2),t56_msg,'FontSize',10);
% 
% t68_loc = [state6_loc(1)+abs((state8_loc(1)- state6_loc(1))/2),state6_loc(2)+abs((state8_loc(2)- state6_loc(2))/2)];
% t68_msg = ['t_{68} = ' num2str(round(t68,6)) 10 't_{86} = ' num2str(round(t86,6))];
% text(t68_loc(1),t68_loc(2),t68_msg,'FontSize',10);
% 
% %Plot the chi-squared value
% text(0.5,3.7,['\chi^2 = ' num2str(chisquared)],'FontSize',12);
% 
% %Plot a title with the information in it
% 
% title_str = ['8 state: ' sample_description];
% text(2.5,4,title_str,'fontsize',16);
% 
% axis off;
% end% end of networkPlotter_8stateModel1
% 
% %% Plot a model of the 9 state network
% % ********* 9 state *********
% function networkPlotter_9stateModel(sample_description,chisquared,...
%     A1,A2,A3,A4,A5,A6,A7,A8,A9,...
%     t12,t13,t21,t23,t24,t25,t31,t32,t36,t37,t73,t42,t52,t56,t63,t65,t68,t86,t69,t96,t89,t98)
% 
% clf;
% set(gcf,'Name','Model: 8 state');
% set(gcf,'Color','white');
% xlim([0 4.1]);
% ylim([0 4.1]);
% 
% %Designate spots for the states
% state1_loc = [0.5 2];
% state2_loc = [1.5 2.5];
% state3_loc = [1.5 1.5];
% state4_loc = [2.5 3.5];
% state5_loc = [2.5 2.5];
% state6_loc = [2.5 1.5];
% state7_loc = [2.5 0.5];
% % state8_loc = [3.5 1.5];
% % or for 9 state
% state8_loc = [3.5 2];
% state9_loc = [3.5 1];
% hold on;
% % text(0.1,4,['Gen# = ' num2str(genNum)],'FontSize',12);
% 
% % Plot the state symbols
% text(state1_loc(1)-0.04,state1_loc(2)+0.05,'1','FontSize',24)
% text(state2_loc(1),state2_loc(2),'2','FontSize',24);
% text(state3_loc(1)-0.05,state3_loc(2),'3','FontSize',24);
% text(state4_loc(1)-0.04,state4_loc(2)+0.05,'4','FontSize',24)
% text(state5_loc(1),state5_loc(2),'5','FontSize',24);
% text(state6_loc(1)-0.05,state6_loc(2),'6','FontSize',24);
% text(state7_loc(1)-0.04,state7_loc(2)+0.05,'7','FontSize',24)
% text(state8_loc(1),state8_loc(2),'8','FontSize',24);
% text(state9_loc(1),state9_loc(2),'9','FontSize',24);
% 
% 
% %Plot the FRET values
% text(state1_loc(1),state1_loc(2)-0.2,['=' num2str(A1,'%.2f')],'FontSize',16, 'Color','b');
% text(state2_loc(1),state2_loc(2)-0.2,['=' num2str(A2,'%.2f')],'FontSize',16, 'Color','b');
% text(state3_loc(1),state3_loc(2)-0.2,['=' num2str(A3,'%.2f')],'FontSize',16, 'Color','b');
% text(state4_loc(1),state4_loc(2)-0.2,['=' num2str(A4,'%.2f')],'FontSize',16, 'Color','b');
% text(state5_loc(1),state5_loc(2)-0.2,['=' num2str(A5,'%.2f')],'FontSize',16, 'Color','b');
% text(state6_loc(1),state6_loc(2)-0.2,['=' num2str(A6,'%.2f')],'FontSize',16, 'Color','b');
% text(state7_loc(1),state7_loc(2)-0.2,['=' num2str(A7,'%.2f')],'FontSize',16, 'Color','b');
% text(state8_loc(1),state8_loc(2)-0.2,['=' num2str(A8,'%.2f')],'FontSize',16, 'Color','b');
% text(state9_loc(1),state9_loc(2)-0.2,['=' num2str(A9,'%.2f')],'FontSize',16, 'Color','b');
% 
% 
% %Plot Lines Between states
% line([state1_loc(1) state2_loc(1)],[state1_loc(2) state2_loc(2)],'Color','k','LineStyle','-');
% line([state2_loc(1) state3_loc(1)],[state2_loc(2) state3_loc(2)],'Color','k','LineStyle','-');
% line([state3_loc(1) state1_loc(1)],[state3_loc(2) state1_loc(2)],'Color','k','LineStyle','-');
% line([state2_loc(1) state4_loc(1)],[state2_loc(2) state4_loc(2)],'Color','k','LineStyle','-');
% line([state2_loc(1) state5_loc(1)],[state2_loc(2) state5_loc(2)],'Color','k','LineStyle','-');
% line([state5_loc(1) state6_loc(1)],[state5_loc(2) state6_loc(2)],'Color','k','LineStyle','-');
% line([state3_loc(1) state6_loc(1)],[state3_loc(2) state6_loc(2)],'Color','k','LineStyle','-');
% line([state3_loc(1) state7_loc(1)],[state3_loc(2) state7_loc(2)],'Color','k','LineStyle','-');
% % line([state6_loc(1) state8_loc(1)],[state6_loc(2) state8_loc(2)],'Color','k','LineStyle','-');
% % For 9 states
% line([state8_loc(1) state9_loc(1)],[state8_loc(2) state9_loc(2)],'Color','k','LineStyle','-');
% line([state6_loc(1) state9_loc(1)],[state6_loc(2) state9_loc(2)],'Color','k','LineStyle','-');
% line([state6_loc(1) state8_loc(1)],[state6_loc(2) state8_loc(2)],'Color','k','LineStyle','-');
% 
% 
% %Plot the inverse of the rates
% t12_loc = [state1_loc(1)+(state2_loc(1)- state1_loc(1))/2,state1_loc(2)+(state2_loc(2)- state1_loc(2))/2];
% t12_msg = ['t_{12} = ' num2str(round(t12,6)) 10 't_{21} = ' num2str(round(t21,6))];
% text(t12_loc(1),t12_loc(2),t12_msg,'FontSize',10);
% 
% t23_loc = [state3_loc(1)+(state2_loc(1)- state3_loc(1))/2,state3_loc(2)+(state2_loc(2)- state3_loc(2))/2];
% t23_msg = ['t_{23} = ' num2str(round(t23,6)) 10 't_{32} = ' num2str(round(t32,6))];
% text(t23_loc(1),t23_loc(2),t23_msg,'FontSize',10);
% 
% t31_loc = [state1_loc(1)+(state3_loc(1)- state1_loc(1))/2,state1_loc(2)+(state3_loc(2)- state1_loc(2))/2];
% t31_msg = ['t_{31} = ' num2str(round(t31,6)) 10 't_{13} = ' num2str(round(t13,6))];
% text(t31_loc(1),t31_loc(2),t31_msg,'FontSize',10);
% 
% t24_loc = [state2_loc(1)+abs((state2_loc(1)- state4_loc(1))/2),state2_loc(2)+abs((state2_loc(2)- state4_loc(2))/2)];
% t24_msg = ['t_{24} = ' num2str(round(t24,6)) 10 't_{42} = ' num2str(round(t42,6))];
% text(t24_loc(1),t24_loc(2),t24_msg,'FontSize',10);
% 
% t37_loc = [state7_loc(1) - abs((state3_loc(1) - state7_loc(1))/2),state7_loc(2)+ abs((state3_loc(2)- state7_loc(2))/2)];
% t37_msg = ['t_{37} = ' num2str(round(t37,6)) 10 't_{73} = ' num2str(round(t73,6))];
% text(t37_loc(1),t37_loc(2),t37_msg,'FontSize',10);
% 
% t25_loc = [state2_loc(1)+abs((state2_loc(1)- state5_loc(1))/2),state2_loc(2)+abs((state2_loc(2)- state5_loc(2))/2)];
% t25_msg = ['t_{25} = ' num2str(round(t25,6)) 10 't_{52} = ' num2str(round(t52,6))];
% text(t25_loc(1),t25_loc(2),t25_msg,'FontSize',10);
% 
% t36_loc = [state3_loc(1)+abs((state3_loc(1)- state6_loc(1))/2),state3_loc(2)+abs((state3_loc(2)- state6_loc(2))/2)];
% t36_msg = ['t_{36} = ' num2str(round(t36,6)) 10 't_{63} = ' num2str(round(t63,6))];
% text(t36_loc(1),t36_loc(2),t36_msg,'FontSize',10);
% 
% t56_loc = [state5_loc(1)+abs((state5_loc(1)- state6_loc(1))/2),state5_loc(2)-abs((state5_loc(2)- state6_loc(2))/2)];
% t56_msg = ['t_{56} = ' num2str(round(t56,6)) 10 't_{65} = ' num2str(round(t65,6))];
% text(t56_loc(1),t56_loc(2),t56_msg,'FontSize',10);
% 
% t68_loc = [state6_loc(1)+abs((state8_loc(1)- state6_loc(1))/2),state6_loc(2)+abs((state8_loc(2)- state6_loc(2))/2)];
% t68_msg = ['t_{68} = ' num2str(round(t68,6)) 10 't_{86} = ' num2str(round(t86,6))];
% text(t68_loc(1),t68_loc(2),t68_msg,'FontSize',10);
% 
% t69_loc = [state9_loc(1)-abs((state6_loc(1)- state9_loc(1))/2),state9_loc(2)+abs((state6_loc(2)- state9_loc(2))/2)];
% t69_msg = ['t_{69} = ' num2str(round(t69,6)) 10 't_{96} = ' num2str(round(t96,6))];
% text(t69_loc(1),t69_loc(2),t69_msg,'FontSize',10);
% 
% t89_loc = [state9_loc(1)-abs((state8_loc(1)- state9_loc(1))/2),state9_loc(2)+abs((state8_loc(2)- state9_loc(2))/2)];
% t89_msg = ['t_{89} = ' num2str(round(t89,6)) 10 't_{98} = ' num2str(round(t98,6))];
% text(t89_loc(1),t89_loc(2),t89_msg,'FontSize',10);
% 
% 
% %Plot the chi-squared value
% text(0.5,3.7,['\chi^2 = ' num2str(chisquared)],'FontSize',12);
% 
% %Plot a title with the information in it
% 
% title_str = ['8 state: ' sample_description];
% text(2.5,4,title_str,'fontsize',16);
% 
% axis off;
% end% end of networkPlotter_8stateModel1
% 
% 
% %%
% %//////////////////////////////////////////////////////////////////////////
% % OPTIMIZATION FUNCTION: MINIMIZE chisquared
% %//////////////////////////////////////////////////////////////////////////
% % FUNCTION: multigoaltcf
% % PURPOSE: returns the rms (root mean square) deviation from the model
% INPUT: (1) Set of state-to-state times ("rates" {t_ij} where k_ij = 1/t_ij}
%        (2) FRET Values for each state
% ********* 5 state *********
% function [chisquared,chisquared_array, chisquared_unweighted, chisquared_unweighted_array] = multigoaltcf_5state(x)
% global normalizeMode diagnoseMode
% global targetHistogram weightingFactor_FREThist
% global C2_exp_x C2_exp_y weightingFactor_C2 weightC2func
% global C4_tau1range C4_tau2eq0_exp weightingFactor_C4_t0 wC4func
% global fitHistMode fitC2Mode fitC4Mode
% global yoff zoff FRET_bins
% global sigma_A
% 
% t12 = x(1);
% t21 = x(2);
% t23 = x(3);
% t24 = x(4);
% t32 = x(5);
% t42 = x(6);
% t45 = x(7);
% t54 = x(8);
% A1 = x(9);
% A2 = x(10);
% A3 = x(11);
% A5 = x(12);
% 
% A4 = A3;    %A3 and A4 should have the same FRET value because they are essentially the same state
% 
% A = [A1, A2, A3, A4, A5];
% [k12,k21,k23,k24,k32,k42,k45,k54] = times2rates_5stateModel(t12,t21,t23,t24,t32,t42,t45,t54);
% K = rates2K_5stateModel(k12,k21,k23,k24,k32,k42,k45,k54);
% 
% C2_time = C2_exp_x;
% P = K2P(K, C2_time);
% C4_time = C4_tau1range;
% 
% [chisquared,chisquared_array, chisquared_unweighted, chisquared_unweighted_array] = chiSqCalc(K, A, P, sigma_A, FRET_bins, C2_time, C4_time);
% 
% end

% % ********* 8 state *********
% function [chisquared,chisquared_array, chisquared_unweighted, chisquared_unweighted_array] = multigoaltcf_8state(x0)
% global normalizeMode diagnoseMode
% global targetHistogram weightingFactor_FREThist
% global C2_exp_x C2_exp_y weightingFactor_C2 weightC2func
% global C4_tau1range C4_tau2eq0_exp weightingFactor_C4_t0 wC4func
% global fitHistMode fitC2Mode fitC4Mode
% global yoff zoff FRET_bins
% global sigma_A
% 
% t12 = x0(1);
% t13 = x0(2);
% t21 = x0(3);
% t23 = x0(4);
% t31 = x0(5);
% 
% t24 = x0(6);
% t42 = x0(7);
% t25 = x0(8);
% t52 = x0(9);
% t56 = x0(10);
% t65 = x0(11);
% t36 = x0(12);
% t37 = x0(13);
% t73 = x0(14);
% t68 = x0(15);
% t86 = x0(16);
% 
% A1 = x0(17);
% A2 = x0(18);
% A3 = x0(19);
% A4 = x0(20);
% A5 = x0(21);
% A6 = x0(22);
% A7 = x0(23);
% A8 = x0(24);
% 
% [t32,t63,k12,k13,k21,k23,k24,k25,k31,k32,k36,k37,k42,k52,k56,k63,k65,k68,k73,k86]...
%     = times2rates_8stateModel1(t12,t13,t21,t23,t24,t25,t31,t36,t37,t42,t52,t56,t65,t68,t73,t86);    % Function at end of code
% 
% %------------------------------------------------------------------
% % Calculate chisquared array and value
% %------------------------------------------------------------------
% % Define K matrix
% % 8 state model:
% K = rates2K_8stateModel1(k12,k13,k21,k23,k24,k25,k31,k32,k36,k37,k42,k52,k56,k63,k65,k68,k73,k86);  % Function at end of code
% 
% A = [A1,A2,A3,A4,A5,A6,A7,A8];
% 
% C2_time = C2_exp_x;
% P = K2P(K, C2_time);
% C4_time = C4_tau1range;
% 
% [chisquared,chisquared_array, chisquared_unweighted, chisquared_unweighted_array] = chiSqCalc(K, A, P, sigma_A, FRET_bins, C2_time, C4_time);
% end
% 
% % ********* 9 state *********
% function [chisquared,chisquared_array, chisquared_unweighted, chisquared_unweighted_array] = multigoaltcf_9state(x0)
% global normalizeMode diagnoseMode
% global targetHistogram weightingFactor_FREThist
% global C2_exp_x C2_exp_y weightingFactor_C2 weightC2func
% global C4_tau1range C4_tau2eq0_exp weightingFactor_C4_t0 wC4func
% global fitHistMode fitC2Mode fitC4Mode
% global yoff zoff FRET_bins
% global sigma_A
% 
% t12 = x0(1);
% t13 = x0(2);
% t21 = x0(3);
% t23 = x0(4);
% t31 = x0(5);
% 
% t24 = x0(6);
% t42 = x0(7);
% t25 = x0(8);
% t52 = x0(9);
% t56 = x0(10);
% t65 = x0(11);
% t36 = x0(12);
% t37 = x0(13);
% t73 = x0(14);
% t68 = x0(15);
% t86 = x0(16);
% t89 = x0(17);
% t98 = x0(18);
% t69 = x0(19);
% 
% 
% A1 = x0(20);
% A2 = x0(21);
% A3 = x0(22);
% A4 = x0(23);
% A5 = x0(24);
% A6 = x0(25);
% A7 = x0(26);
% A8 = x0(27);
% A9 = x0(28);
% 
% [t32,t63,t96,k12,k13,k21,k23,k24,k25,k31,k32,k36,k37,k42,k52,k56,k63,k65,k68,k73,k86,k89,k98,k96,k69] = ...
%     times2rates_9stateModel(t12,t13,t21,t23,t24,t25,t31,t36,t37,t42,t52,t56,t65,t68,t73,t86,t89,t98,t69);
% 
% %------------------------------------------------------------------
% % Calculate chisquared array and value
% %------------------------------------------------------------------
% % Calculate K matrix
% K = rates2K_9stateModel(k12,k13,k21,k23,k24,k25,k31,k32,k36,k37,k42,k52,k56,k63,k65,k68,k73,k86,k89,k69,k96,k98);
% 
% A = [A1,A2,A3,A4,A5,A6,A7,A8,A9];
% 
% C2_time = C2_exp_x;
% P = K2P(K, C2_time);
% C4_time = C4_tau1range;
% 
% [chisquared,chisquared_array, chisquared_unweighted, chisquared_unweighted_array] = chiSqCalc(K, A, P, sigma_A, FRET_bins, C2_time, C4_time);
% end%end of multigoal9state

%% Save the information
% ********* 8 state *********
% function saveInfo_5state(outputFolderName,foutName,...
%     k12,k21,k23,k24,k32,k42,k45,k54,...
%     t12,t21,t23,t24,t32,t42,t45,t54,...
%     A1,A2,A3,A4,A5,time,Peq,K,A,P,sigma_A,FRET_bins,C2_time,C4_time,yoff,zoff,...
%     guess,genNum_array,chisquared_vs_Gen_array,chisquaredUnweighted_vs_Gen_array,iter,chisquared,...
%     chisquared_unweighted,weightingFactor_FREThist,weightingFactor_C2,weightingFactor_C4_t0)
% 
% if exist(outputFolderName,'dir')~= 7
%     disp(['Making folder ' outputFolderName ' to hold output']);
%     mkdir(outputFolderName);
% end
% 
% if exist([outputFolderName filesep() foutName],'file') == 1
%     save([outputFolderName filesep() foutName],...
%         'k12','k21','k23','k24','k32','k42','k45','k54',...
%         't12','t21','t23','t24','t32','t42','t45','t54',...
%         'A1','A2','A3','A4','A5','time','Peq',...
%         'K','A','P','sigma_A','FRET_bins','C2_time','C4_time','yoff','zoff',...
%         'guess','genNum_array','chisquared_vs_Gen_array','chisquaredUnweighted_vs_Gen_array',...
%         'iter','chisquared','chisquared_unweighted','weightingFactor_FREThist',...
%         'weightingFactor_C2','weightingFactor_C4_t0','-append');
% else
%     save([outputFolderName filesep() foutName],...
%         'k12','k21','k23','k24','k32','k42','k45','k54',...
%         't12','t21','t23','t24','t32','t42','t45','t54',...
%         'A1','A2','A3','A4','A5','time','Peq',...
%         'K','A','P','sigma_A','FRET_bins','C2_time','C4_time','yoff','zoff',...
%         'guess','genNum_array','chisquared_vs_Gen_array','chisquaredUnweighted_vs_Gen_array',...
%         'iter','chisquared','chisquared_unweighted','weightingFactor_FREThist',...
%         'weightingFactor_C2','weightingFactor_C4_t0');
% end
% end%end of saveInfo_5state function
% 
% % ********* 8 state *********
% function saveInfo_8state(outputFolderName,foutName,k12,k13,k21,k23,k24,k25,k31,k32,k36,k37,k42,k52,k56,k63,k65,k68,k73,k86,...
%     t12,t13,t21,t23,t24,t25,t31,t32,t36,t37,t42,t52,t56,t63,t65,t68,t73,t86,...
%     A1,A2,A3,A4,A5,A6,A7,A8,time,Peq,K,A,P,sigma_A,FRET_bins,C2_time,C4_time,yoff,zoff,...
%     guess,genNum_array,chisquared_vs_Gen_array,chisquaredUnweighted_vs_Gen_array,iter,chisquared,...
%     chisquared_unweighted,weightingFactor_FREThist,weightingFactor_C2,weightingFactor_C4_t0)
% 
% if exist(outputFolderName,'dir')~= 7
%     disp(['Making folder ' outputFolderName ' to hold output']);
%     mkdir(outputFolderName);
% end
% 
% if exist([outputFolderName filesep() foutName],'file') == 1
%     save([outputFolderName filesep() foutName],...
%         'k12','k13','k21','k23','k24','k25','k31','k32','k36','k37','k42','k52','k56','k63','k65','k68','k73','k86',...
%         't12','t13','t21','t23','t24','t25','t31','t32','t36','t37','t42','t52','t56','t63','t65','t68','t73','t86',...
%         'A1','A2','A3','A4','A5','A6','A7','A8','time','Peq',...
%         'K','A','P','sigma_A','FRET_bins','C2_time','C4_time','yoff','zoff',...
%         'guess','genNum_array','chisquared_vs_Gen_array','chisquaredUnweighted_vs_Gen_array',...
%         'iter','chisquared','chisquared_unweighted','weightingFactor_FREThist',...
%         'weightingFactor_C2','weightingFactor_C4_t0','-append');
% else
%     save([outputFolderName filesep() foutName],...
%         'k12','k13','k21','k23','k24','k25','k31','k32','k36','k37','k42','k52','k56','k63','k65','k68','k73','k86',...
%         't12','t13','t21','t23','t24','t25','t31','t32','t36','t37','t42','t52','t56','t63','t65','t68','t73','t86',...
%         'A1','A2','A3','A4','A5','A6','A7','A8','time','Peq',...
%         'K','A','P','sigma_A','FRET_bins','C2_time','C4_time','yoff','zoff',...
%         'guess','genNum_array','chisquared_vs_Gen_array','chisquaredUnweighted_vs_Gen_array',...
%         'iter','chisquared','chisquared_unweighted','weightingFactor_FREThist',...
%         'weightingFactor_C2','weightingFactor_C4_t0');
% end
% end%end of saveInfo_8state function
% 
% % ********* 9 state *********
% function saveInfo_9state(outputFolderName,foutName,k12,k13,k21,k23,k24,k25,k31,k32,k36,k37,k42,k52,k56,k63,k65,k68,k73,k86,...
%     k69,k96,k89,k98,...
%     t12,t13,t21,t23,t24,t25,t31,t32,t36,t37,t42,t52,t56,t63,t65,t68,t73,t86,...
%     t69,t96,t89,t98,...
%     A1,A2,A3,A4,A5,A6,A7,A8,A9,time,Peq,K,A,P,sigma_A,FRET_bins,C2_time,C4_time,yoff,zoff,...
%     guess,genNum_array,chisquared_vs_Gen_array,chisquaredUnweighted_vs_Gen_array,iter,chisquared,...
%     chisquared_unweighted,weightingFactor_FREThist,weightingFactor_C2,weightingFactor_C4_t0)
% 
% if exist(outputFolderName,'dir')~= 7
%     disp(['Making folder ' outputFolderName ' to hold output']);
%     mkdir(outputFolderName);
% end
% 
% if exist([outputFolderName filesep() foutName],'file') == 1
%     save([outputFolderName filesep() foutName],...
%         'k12','k13','k21','k23','k24','k25','k31','k32','k36','k37','k42','k52','k56','k63','k65','k68','k73','k86',...
%         't12','t13','t21','t23','t24','t25','t31','t32','t36','t37','t42','t52','t56','t63','t65','t68','t73','t86',...
%         'k69','k96','k89','k98',...
%         't69','t96','t89','t98',...
%         'A1','A2','A3','A4','A5','A6','A7','A8','A9','time','Peq',...
%         'K','A','P','sigma_A','FRET_bins','C2_time','C4_time','yoff','zoff',...
%         'guess','genNum_array','chisquared_vs_Gen_array','chisquaredUnweighted_vs_Gen_array',...
%         'iter','chisquared','chisquared_unweighted','weightingFactor_FREThist',...
%         'weightingFactor_C2','weightingFactor_C4_t0','-append');
% else
%     save([outputFolderName filesep() foutName],...
%         'k12','k13','k21','k23','k24','k25','k31','k32','k36','k37','k42','k52','k56','k63','k65','k68','k73','k86',...
%         't12','t13','t21','t23','t24','t25','t31','t32','t36','t37','t42','t52','t56','t63','t65','t68','t73','t86',...
%         'k69','k96','k89','k98',...
%         't69','t96','t89','t98',...
%         'A1','A2','A3','A4','A5','A6','A7','A8','A9','time','Peq',...
%         'K','A','P','sigma_A','FRET_bins','C2_time','C4_time','yoff','zoff',...
%         'guess','genNum_array','chisquared_vs_Gen_array','chisquaredUnweighted_vs_Gen_array',...
%         'iter','chisquared','chisquared_unweighted','weightingFactor_FREThist',...
%         'weightingFactor_C2','weightingFactor_C4_t0');
% end
% end%end of saveInfo_9state function
% 
end