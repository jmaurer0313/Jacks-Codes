%__________________________________________________________________________
% AUTHOR: Brett Israels
%
% NAME: fmincon_general_Corrected
%
% FUNCTION: % Fits (1) the FRET Histogram, (2) the 2pt TCF and (3) the 4-point TCF

%--------------------------------------------------------------------------
% PROCEDURE:
%--------------------------------------------------------------------------
% (1) Load the experimental histogram
% (2) Load the experimental 2pt time correlation functions to fit
% (3) Load the experimental 4pt TCFs for   various tau2 values
% (4) Make an initial guess of the model paramaters:
% (5) Make an array of lower bounds
% (6) Make an array of upper bounds
% (7) feed the guess into the optimization algorithm: x = ga(fun,x0,A,b,Aeq,beq,lb,ub);

%--------------------------------------------------------------------------
% INPUT:
%--------------------------------------------------------------------------
% (1) FRET Histogram :
% (2) Two-Point TCF :
% (3) Four-Point TCF :

%--------------------------------------------------------------------------
% OUTOUT:
%--------------------------------------------------------------------------
% (1) Histogram Fit
% (2) C2 Fit
% (3) C4 Fit
% (4) Final fit values
% (5) Model with FRET state values & Rates


%--------------------------------------------------------------------------
% EXTERNAL PROGRAMS CALLS: (needs these codes in MATLAB PATH to function)
%--------------------------------------------------------------------------
% writecell.m (only in matlab 2020 not 2018)

%--------------------------------------------------------------------------
% MODIFICATION LOG:
%--------------------------------------------------------------------------
% 4-20-20 BI: Now works on models 5,6,6.1 (loop)
% Nov 14 uses new optimization
% Global optimization

programName = 'fmincon_general_Corrected';
disp(['Beggining ' programName ' at ' char(datetime) ]);
% disp(['Now Running ' programName '.m']);

%% Declare global variables
global normalizeMode verboseMode  diagnoseMode
global fitHistMode fitC2Mode fitC4Mode
global targetHistogram weightingFactor_FREThist
global sigma_A FRET_bins  % Histogram optimization
global C2_exp_x C2_exp_y weightingFactor_C2  weightC2func yoff
global C4_tau1range C4_tau2eq0_exp weightingFactor_C4_t0 wC4func zoff
global showProgressOnFit_mode
global NbinsHist histResUsec

%% PART 1: Set the program up with various options

constructName = getConstructFolderName();

%% Load the best fit from the genetic algorithm
fileBeginName = 'BestFitResults_ga_general_';
fileEndName = 'state.mat';
fileKeyWord = [fileBeginName '*' fileEndName];
fileNames = dir(['*' fileKeyWord]);
if length(fileNames) < 1
    error(['Cannot find any files ending with the name ' fileKeyWord ' : Exiting program.']);
else
    ga_bestFitFileName = fileNames(1).name;
    disp(['Loading ' ga_bestFitFileName ' to get the values']);
    load(ga_bestFitFileName);
end


%% Determine the model the user is using
strStart = strfind(ga_bestFitFileName,fileBeginName)+length(fileBeginName);
modelstr = ga_bestFitFileName(strStart:length(ga_bestFitFileName) - length(fileEndName));
modelNum = str2double(modelstr);
if isnan(modelNum)
    error('Cannot detect state name');
end
fprintf('Detected Model = %d\r',modelNum);

switch modelNum
    case 5
        disp('Optimizing the 5 state model');
    case 6
        disp('Optimizing the 6 state model');
    case 6.1
        disp('Optimizing the 6state Loop model');
end

NbinsHist = 25;
histResUsec = 1000;

%% User Options                                                  (Part 1)
plotDataOnlyMode = 0;
plotDataFirstMode = 1;

normalizeMode = 1; %Normalizes the 2pt and the 4pt correlation functions
verboseMode = 1;
diagnoseMode = 0 ; %Shows relative contribution of the 2pt TCF to the histogram at the end of each chi square calculation
showProgressOnFit_mode = 0;%Updates fit during each iteration of optimization
showFminconUpdates = 0;

clockMode = 0; %Times various features of the algorithm

saveModeDefault = 0;
updateGlobalFitMode_default = 1;

datestr_mode = 1;

plotMode = 1;                               %Makes plots.
useFigurePosnMode = 0;

%Choose which surfaces to fit
fitHistMode = 1;
fitC2Mode = 1;
fitC4Mode = 1;

%Choose Weighting Amount
weightingFactor_FREThist = 100;
weightingFactor_C2 = 1;
weightingFactor_C4_t0 = 1;

% time_LB = 20e-6;
time_LB_usec = 20;

if sum([fitHistMode,fitC2Mode,fitC4Mode]) == 0
    error('No surfaces to optimize to! Pick atleast one.');
end

%% Start in the single molecule folder (smData_Processed): comp specific

%Start in the single molecule folder (smData_Processed): comp specific
% [computer_terminal_str, terminalID, DropboxLocationPrefix] = computerMode();
[computer_terminal_str, terminalID] = computerMode();
%%  Set up the figure positions once and done
if plotMode == 1
    clear figPosn
    if useFigurePosnMode == 1
        figPosn = setFigPosn(terminalID);
        if fitHistMode == 1
            figure(11);
            clf;
            set(gcf,'Name','FRET Histogram');
            set(gcf,'Color','w');
            set(gcf,'Position',figPosn(:,1));
        end
        %-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.
        if fitC2Mode == 1
            figure(22);
            clf;
            set(gcf,'Name','Two-Point TCF: C2');
            set(gcf,'Position',figPosn(:,2));
        end
        %-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.
        if fitC4Mode == 1
            figure(33)
            clf;
            set(gcf,'Name','Four-Point TCF: C4');
            set(gcf,'Position',figPosn(:,3));
        end
    end
    %-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.
end

%% SELECT MODEL YOU WANT TO USE

outputProgramName = [programName '_' modelstr 'state']; % rename program output depending on which model you run

% Reset some paramaters
updateGlobalFitMode = updateGlobalFitMode_default;
saveMode = saveModeDefault;
%% Navigate to the correct folder

% Find GenAlg files and load target data:
protein_str = 'gp32_0p5uM';
[GenAlgFilesLocationPrefix, genAlgDataLocationPrefix] = ga_fileLocator(terminalID, constructName, protein_str);
%         cd(GenAlgFilesLocationPrefix);

[sample_description, save_prefix] = sample_descriptionGetter();
wd = pwd;
if verboseMode == 1
    disp(['Part1 complete: You are now in the folder: ' wd]);
end

%% PART 2: Load the target data (and plot if opted)
[histogram_FilePath, C2_FilePath, C4_FilePath] = loadFilePath_dataPlots(genAlgDataLocationPrefix, constructName);

%% Optimization Target #1: 1-D FRET HISTOGRAM
load(histogram_FilePath,'xData','yData');
FRET_bins = xData;
targetHistogram = yData;

% 'Normalize' the histogram so that the sum of frequency = 1;
targetHistogram = targetHistogram./sum(targetHistogram);

if plotMode == 1 && plotDataFirstMode == 1
    figure(11);
    clf;
    FRETDataHistPlot(FRET_bins,targetHistogram);
end

%% (2) Optimization Target #2: 2-point TCF (C2)           (20 uSec and on)  (PART 2: Load the target data)
if fitC2Mode == 1
    %         load(C2_FilePath,'time','yData','y','yoff');
    load(C2_FilePath,'tauArraySec','C2avg');
    yoff = abs(mean(C2avg(end-3:end)));
    if verboseMode == 1
        fprintf('Setting C2 yoff = abs(mean(C2_exp_y(end-3:end))) = %f\r\n',yoff);
    end
    
    C2_exp_x = tauArraySec;
    C2_exp_y = C2avg';
    
    %Take out the lower bound time limit
    time_usec = tauArraySec*1e6;
    timeStartIndex = find(time_usec >= time_LB_usec,1,'first');
    C2_exp_y = C2_exp_y(timeStartIndex:end);
    C2_exp_x = C2_exp_x(timeStartIndex:end);
    
    time = C2_exp_x;
    
    if normalizeMode == 1
        C2_exp_y = C2_exp_y./C2_exp_y(1);
    end
    
    weightC2func = 1./power(C2_exp_x,0.5);% 2-pt. TCF weighting function
    
    if plotMode == 1 && plotDataFirstMode == 1
        figure(22);
        clf;
        C2dataPlotter(C2_exp_x,C2_exp_y);
    end
end


%% (3) Optimization Target #3: 4-point TCF (C4)
if fitC4Mode == 1
    % load(C4_FilePath,'C4','tau1arraySec','tau3arraySec','tau2ValSec');
    % C4_tau2eq0_exp = C4;
    
    load(C4_FilePath,'C4_tau2eq0_avg','tau1arraySec','tau3arraySec','tau2ValSec');
    C4_tau2eq0_exp = C4_tau2eq0_avg;
    
    zoff = abs(mean(mean(C4_tau2eq0_exp(end-3:end-3))));
    if verboseMode == 1
        fprintf('setting zoff = abs(mean(mean(C4_tau2eq0_exp(end-3:end-3)))) = %f\r\n',zoff);
    end
    
    C4_tau1range = tau1arraySec;
    C4_tau3range = tau3arraySec';
    
    %Take out the lower bound time limit
    time_usec = tau1arraySec*1e6;
    timeStartIndex = find(time_usec >= time_LB_usec,1,'first');
    C4_tau1range = C4_tau1range(timeStartIndex:end);
    C4_tau3range = C4_tau3range(timeStartIndex:end);
    C4_tau2eq0_exp = C4_tau2eq0_exp(timeStartIndex:end,timeStartIndex:end);
    
    tau2 = tau2ValSec;
    
    if normalizeMode == 1
        C4_tau2eq0_exp = C4_tau2eq0_exp./C4_tau2eq0_exp(1,1);
    end
    
    
    % 4-pt. TCF weighting function
    %         wC4func = 1./sqrt(C4_tau1range).*(1./(sqrt(C4_tau1range)))';
    wC4func = 1./power(C4_tau1range,.5).*(1./(power(C4_tau1range,.5)))';
    
    if plotMode == 1 && plotDataFirstMode == 1
        figure(33)
        clf;
        C4dataPlotter(C4_tau1range,C4_tau3range,C4_tau2eq0_exp,tau2ValSec);
    end
end

if verboseMode == 1
    disp('     Part1: Done Loading and (optionally) plotting the data.');
end
if plotDataOnlyMode == 1
    disp('plotDataOnlyMode == 1, so quitting');
    return;
end
%% PART 3: Display THE BEST FIT (Determine starting paramaters)
if verboseMode == 1
    switch modelNum
        case 5
            fprintf(['Starting points: t12 = %f, t21 = %f, t23 = %f, t24 = %f, t32 = %f, t42 = %f'...
                '\n t45 = %f, t54 = %f, A1 = %f, A2 = %f, A3 = %f, A4 = %f, A5 = %f\r\n'],...
                t12,t21,t23,t24,t32,t42,t45,t54,A1,A2,A3,A4,A5);
        case 6
            fprintf(['Starting points: t12 = %f, t21 = %f, t23 = %f, t24 = %f, t32 = %f, t42 = %f'...
                '\n t45 = %f, t54 = %f, t56 = %f, t65 = %f\r\n'...
                'A1 = %f, A2 = %f, A3 = %f, A4 = %f, A5 = %f, A6 = %f\r\n'],...
                t12,t21,t23,t24,t32,t42,t45,t54,t56,t65,A1,A2,A3,A4,A5,A6);
        case 6.1
            fprintf(['Starting points: t12 = %f, t21 = %f, t23 = %f, t24 = %f, t32 = %f, t42 = %f'...
                '\n t45 = %f, t46 = %f, t54 = %f, t56 = %f, t64 = %f'...
                'tA1 = %f, A2 = %f, A3 = %f, A4 = %f, A5 = %f, A6 = %f\r\n'],...
                t12,t21,t23,t24,t32,t42,t45,t46,t54,t56,t64,A1,A2,A3,A4,A5,A6);
    end
end

%% Set the conformational width of each state

sigma_A1 = 0.15;
sigma_A2 = 0.15;
sigma_A3 = 0.1;
sigma_A4 = 0.1;
sigma_A5 = 0.1;
sigma_A6 = 0.1;

switch modelNum
    case 5
        sigma_A = [sigma_A1; sigma_A2; sigma_A3; sigma_A4; sigma_A5];
    case {6,6.1}
        sigma_A = [sigma_A1; sigma_A2; sigma_A3; sigma_A4; sigma_A5; sigma_A6];
end

%% Determine the bounds of the paramaters of the model
switch modelNum
    case {5,6,6.1}
        %5statebounds_
        percent_window = 0.25;
        
        t12_bounds = [t12*(1-percent_window)*1,t12*(1+percent_window)*1];
        t21_bounds = [t21*(1-percent_window)*1,t21*(1+percent_window)*1];
        
        %         percent_window = 0.1;
        t23_bounds = [t23*(1-percent_window)*1,t23*(1+percent_window)*1];
        t24_bounds = [t24*(1-percent_window)*1,t24*(1+percent_window)*1];
        t32_bounds = [t32*(1-percent_window)*1,t32*(1+percent_window)*1];
        t42_bounds = [t42*(1-percent_window)*1,t42*(1+percent_window)*1];
        t45_bounds = [t45*(1-percent_window)*1,t45*(1+percent_window)*1];
        t54_bounds = [t54*(1-percent_window)*1,t54*(1+percent_window)*1];
        
        %         percent_window = 0.1;
        A1_bounds = [A1*(1-percent_window),A1*(1+percent_window)];%
        A2_bounds = [A2*(1-percent_window),A2*(1+percent_window)];%
        A3_bounds = [A3*(1-percent_window),A3*(1+percent_window)];%
        A4_bounds = [A4*(1-percent_window),A4*(1+percent_window)];%
        A5_bounds = [A5*(1-percent_window),A5*(1+percent_window)];%
        
        
        if modelNum == 6
            t56_bounds = [t56*(1-percent_window)*1,t56*(1+percent_window)*1];
            t65_bounds = [t65*(1-percent_window)*1,t65*(1+percent_window)*1];
            A6_bounds = [A6*(1-percent_window),A6*(1+percent_window)];%
        end
        
        if modelNum == 6.1
            t46_bounds = [t46*(1-percent_window)*1,t46*(1+percent_window)*1];
            t56_bounds = [t56*(1-percent_window)*1,t56*(1+percent_window)*1];
            t64_bounds = [t64*(1-percent_window)*1,t64*(1+percent_window)*1];
            A6_bounds = [A6*(1-percent_window)*1,A6*(1+percent_window)*1];
        end
        
end

switch modelNum
    case  5
        boundsArray = [t12_bounds; t21_bounds; t23_bounds; t24_bounds; t32_bounds;...
            t42_bounds; t45_bounds; t54_bounds;...
            A1_bounds; A2_bounds;A3_bounds; A4_bounds; A5_bounds];
    case  6
        boundsArray = [t12_bounds; t21_bounds; t23_bounds; t24_bounds; t32_bounds;...
            t42_bounds; t45_bounds; t54_bounds; t56_bounds; t65_bounds;...
            A1_bounds; A2_bounds;A3_bounds; A4_bounds; A5_bounds; A6_bounds];
    case 6.1
        boundsArray = [t12_bounds; t21_bounds; t23_bounds; t24_bounds; t32_bounds;...
            t42_bounds; t45_bounds; t46_bounds; t54_bounds; t56_bounds; t64_bounds;...
            A1_bounds; A2_bounds;A3_bounds; A4_bounds; A5_bounds; A6_bounds];
end

lb = boundsArray(:,1);
ub = boundsArray(:,2);

%% PART 4: Use the optimization algorithms
switch  modelNum
    case 5
        fun = @multigoaltcf_5state;
        invrates = [t12,t21,t23,t24,t32,t42,t45,t54];%8 rates
        FRETstates = [A1,A2,A3,A4,A5];%5 FRET states
        x0 = [invrates,FRETstates];% Start point (row vector)
    case 6
        fun = @multigoaltcf_6state;
        invrates = [t12,t21,t23,t24,t32,t42,t45,t54,t56,t65];%10 rates
        FRETstates = [A1,A2,A3,A4,A5,A6];%6 fret states
        x0 = [invrates,FRETstates];% Start point (row vector)
        [chisquared,~] = multigoaltcf_6state(x0);
        
    case 6.1
        fun = @multigoaltcf_6stateLoop;
        invrates = [t12,t21,t23,t24,t32,t42,t45,t46,t54,t56,t64];%11 rates
        FRETstates = [A1,A2,A3,A4,A5,A6];%6 fret states
        x0 = [invrates,FRETstates];% Start point (row vector)
        [chisquared,~] = multigoaltcf_6stateLoop(x0);
        
end

disp(['chisquared after ga: X^2(x0) = ' num2str(chisquared)]);
disp(['chisquared_unweighted after ga: X^2(x0) = ' num2str(chisquared_unweighted)]);

if showFminconUpdates == 1
    nonlcon = [];
    options = optimoptions('fmincon','Display','iter');
    [x,fval] = fmincon(fun,x0,[],[],[],[],lb,ub,[],options);
else
    [x,fval] = fmincon(fun,x0,[],[],[],[],lb,ub);
end

%% Whats the point of the next section?
% switch modelNum
%     case 5
%         [chisquared, chisquared_array, chisquared_unweighted, chisquared_unweighted_array] = multigoaltcf_5state(x);
%     case 6
%         [chisquared, chisquared_array, chisquared_unweighted, chisquared_unweighted_array] = multigoaltcf_6state(x);
%     case 6.1
%         [chisquared, chisquared_array, chisquared_unweighted, chisquared_unweighted_array]= multigoaltcf_6stateLoop(x);
% end

%% PART 5: Display the Final values to the user
switch modelNum
    case 5
        
        t12 = x(1);
        t21 = x(2);
        t23 = x(3);
        t24 = x(4);
        t32 = x(5);
        t42 = x(6);
        t45 = x(7);
        t54 = x(8);
        A1 = x(9);
        A2 = x(10);
        A3 = x(11);
        A4 = x(12);
        A5 = x(13);
        
        A = [A1, A2, A3, A4, A5];
        [k12,k21,k23,k24,k32,k42,k45,k54] = times2rates_5stateModel(t12,t21,t23,t24,t32,t42,t45,t54);
        K = rates2K_5stateModel(k12,k21,k23,k24,k32,k42,k45,k54);
        
        [chisquared,chisquared_array, chisquared_unweighted, chisquared_unweighted_array]= multigoaltcf_5state(x);
    case 6
        
        t12 = x(1);
        t21 = x(2);
        t23 = x(3);
        t24 = x(4);
        t32 = x(5);
        t42 = x(6);
        t45 = x(7);
        t54 = x(8);
        t56 = x(9);
        t65 = x(10);
        
        A1 = x(11);
        A2 = x(12);
        A3 = x(13);
        A4 = x(14);
        A5 = x(15);
        A6 = x(16);
        
        A = [A1, A2, A3, A4, A5, A6];
        [k12,k21,k23,k24,k32,k42,k45,k54,k56,k65] = times2rates_6stateModel(t12,t21,t23,t24,t32,t42,t45,t54,t56,t65);
        K = rates2K_6stateModel(k12,k21,k23,k24,k32,k42,k45,k54,k56,k65);
        
        [chisquared, chisquared_array, chisquared_unweighted, chisquared_unweighted_array]= multigoaltcf_6state(x);
    case 6.1
        t12 = x(1);
        t21 = x(2);
        t23 = x(3);
        t24 = x(4);
        t32 = x(5);
        t42 = x(6);
        t45 = x(7);
        t46 = x(8);
        t54 = x(9);
        t56 = x(10);
        t64 = x(11);
        A1 = x(12);
        A2 = x(13);
        A3 = x(14);
        A4 = x(15);
        A5 = x(16);
        A6 = x(17);
        A = [A1, A2, A3, A4, A5, A6];
        [t65,k12,k21,k23,k24,k32,k42,k45,k46,k54,k56,k64,k65] = times2rates_6stateLoopModel(t12,t21,t23,t24,t32,t42,t45,t46,t54,t56,t64);
        K = rates2K_6stateLoopModel(k12,k21,k23,k24,k32,k42,k45,k46,k54,k56,k64,k65);
        
        [chisquared, chisquared_array, chisquared_unweighted, chisquared_unweighted_array]= multigoaltcf_6stateLoop(x);
        
        
end

disp(['chisquared after fmincon: X^2(x) = ' num2str(chisquared)]);
disp(['chisquared_unweighted fmincon: X^2(x) = ' num2str(chisquared_unweighted)]);

[P, V, K, time] = K2P(K,time);

% [Peq, denom_hist_sim, hist_sim, Peq_names] = histMaker_Nstate(P, A, sigma_A, FRET_bins);
[Peq, denom_hist_sim, hist_sim, Peq_names] = histMaker_Nstate_v3normCorrect(P, A, sigma_A, FRET_bins);

if verboseMode == 1
    switch modelNum
        case 5
            fprintf(['     Final values at the end of fmincon:\n t12 = %f, t21 = %f, t23 = %f, t24 = %f, t32 = %f, t42 = %f\n '...
                't45 = %f, t54 = %f \n A1 = %f, A2 = %f, A3 = %f, A4 = %f, A5 = %f\n'...
                'chisquared = %f, chisquared_unweighted = %f, p1_eq = %f,p2_eq = %f,p3_eq = %f, p4_eq = %f, p5_eq = %f\r\n'],...
                t12,t21,t23,t24,t32,t42,t45,t54,...
                A1,A2,A3,A4,A5,chisquared,chisquared_unweighted,...
                Peq(1),Peq(2),Peq(3),Peq(4),Peq(5));
        case 6
            fprintf(['     Final values at the end of fmincon:\n t12 = %f, t21 = %f, t23 = %f, t24 = %f, t32 = %f, t42 = %f\n '...
                't45 = %f, t54 = %f, t56 = %f, t65 = %f \n'...
                'A1 = %f, A2 = %f, A3 = %f, A4 = %f, A5 = %f, A6 = %f\n'...
                'chisquared = %f, chisquared_unweighted = %f\r'...
                'p1_eq = %f,p2_eq = %f,p3_eq = %f, p4_eq = %f, p5_eq = %f, p6_eq = %f\r\n'],...
                t12,t21,t23,t24,t32,t42,t45,t54,t56,t65,...
                A1,A2,A3,A4,A5,A6,chisquared,chisquared_unweighted,...
                Peq(1),Peq(2),Peq(3),Peq(4),Peq(5),Peq(6));
        case 6.1
            fprintf(['     Final values at the end of fmincon:\n t12 = %f, t21 = %f, t23 = %f, t24 = %f, t32 = %f, t42 = %f\n '...
                't45 = %f, t46 = %f, t54 = %f, t56 = %f, t64 = %f, t65 = %f \n'...
                'A1 = %f, A2 = %f, A3 = %f, A4 = %f, A5 = %f, A6 = %f\n'...
                'chisquared = %f, chisquared_unweighted = %f\r'...
                'p1_eq = %f,p2_eq = %f,p3_eq = %f, p4_eq = %f, p5_eq = %f, p6_eq = %f\r\n'],...
                t12,t21,t23,t24,t32,t42,t45,t46,t54,t56,t64,t65,...
                A1,A2,A3,A4,A5,A6,chisquared,chisquared_unweighted,...
                Peq(1),Peq(2),Peq(3),Peq(4),Peq(5),Peq(6));
            
    end
end

%--------------------------------------------------------------------------
%% Save the data (PART 5: Final Values)
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% 5.1: Figure out where to save the data
%--------------------------------------------------------------------------
if saveMode == 1
    clear('LocalOutputFolderName')
    foutName = [save_prefix 'BestFitResults_' outputProgramName '.mat'];
    
    LocalOutputFolderName = [outputProgramName '_output'];
    
    disp(['Checking for existance of ' LocalOutputFolderName]);
    if exist(LocalOutputFolderName,'dir') ~= 7
        mkdir(LocalOutputFolderName);
        disp(['     Making a folder called "' LocalOutputFolderName '" to hold the output']);
    end
    
    bestfitFilePath = [LocalOutputFolderName filesep() foutName];
    if exist(LocalOutputFolderName,'file') ~= 2 %If the best fit DNE
        disp(['     Will save the local fmincon results as "' LocalOutputFolderName filesep() foutName '"']);
        saveFitMode = 1;
    elseif exist([LocalOutputFolderName filesep() foutName],'file') == 2
        current_LowestUWChisquared = chisquared_unweighted;
        load([LocalOutputFolderName filesep() foutName],'chisquared_unweighted');
        
        prev_LowestUWChisquared = chisquared_unweighted;%Save the previous X2 as prev_LowestChisquared
        chisquared_unweighted = current_LowestUWChisquared;%reset chisquared to its current value
        
        if current_LowestUWChisquared < prev_LowestUWChisquared
            %Just found a new best fit. Need to sPave it.\
            fprintf('Congratulations! %f < %f (#1 fit). Will update best fit for %s .\r\n',current_LowestUWChisquared,prev_LowestUWChisquared,constructName)
            saveFitMode = 1;
            saveMode = 1;
        elseif current_LowestUWChisquared == prev_LowestUWChisquared
            fprintf('Nothing changed.%f = %f. \r\n',current_LowestUWChisquared,prev_LowestUWChisquared);
            saveMode = 1;
            saveFitMode = 1;
            
        elseif current_LowestUWChisquared > prev_LowestUWChisquared
            fprintf('Oh no! %f > %f (#1 fit).\n',current_LowestUWChisquared,prev_LowestUWChisquared)
            saveMode = 0;
            saveFitMode = 0;
        end
    end
    
    
    if saveFitMode == 1
        switch modelNum
            case 5
                saveInfo_5state(LocalOutputFolderName,foutName,x,...
                    k12,k21,k23,k24,k32,k42,k45,k54,t12,t21,t23,t24,t32,t42,t45,t54,...
                    A1,A2,A3,A4,A5,time,Peq,K,A,P,sigma_A,FRET_bins,yoff,zoff,...
                    chisquared,chisquared_unweighted,weightingFactor_FREThist,weightingFactor_C2,weightingFactor_C4_t0);
            case 6
                saveInfo_6state(LocalOutputFolderName,foutName,x,...
                    k12,k21,k23,k24,k32,k42,k45,k54,k56,k65,...
                    t12,t21,t23,t24,t32,t42,t45,t54,t56,t65,...
                    A1,A2,A3,A4,A5,A6,time,Peq,K,A,P,sigma_A,FRET_bins,yoff,zoff,...
                    chisquared,chisquared_unweighted,weightingFactor_FREThist,weightingFactor_C2,weightingFactor_C4_t0);
            case 6.1
                saveInfo_6stateLoop(LocalOutputFolderName,foutName,x,...
                    k12,k21,k23,k24,k32,k42,k45,k46,k54,k56,k64,k65,...
                    t12,t21,t23,t24,t32,t42,t45,t46,t54,t56,t64,t65,...
                    A1,A2,A3,A4,A5,A6,time,Peq,K,A,P,sigma_A,FRET_bins,yoff,zoff,...
                    chisquared,chisquared_unweighted,weightingFactor_FREThist,weightingFactor_C2,weightingFactor_C4_t0);
        end
    end
end




%% Make a folder to hold the Global best fit

%Retain the current chisquared unweighted value
Current_chisquared_unweighted = chisquared_unweighted;%hold the value
if saveMode == 1
    BestFitFolder = ['ga_general' filesep() 'BestFits' filesep() 'BestFit_' modelstr 'state'];
    BestFitFolder_Path = [GenAlgFilesLocationPrefix filesep() BestFitFolder];
    if exist(BestFitFolder_Path,'dir') ~= 7
        disp('Making folder to store the best fit');
        mkdir(BestFitFolder_Path)
        GlobalLowestUWChiSq = Current_chisquared_unweighted;
    else %If the folder does exist
        finName = [save_prefix 'GlobalBestFitResults_fmincon_general_Corrected_' modelstr 'state.mat'];
        load([BestFitFolder_Path filesep() finName],...
            'chisquared_unweighted');
        GlobalLowestUWChiSq = chisquared_unweighted;
        chisquared_unweighted = Current_chisquared_unweighted;%reassign
        disp(['The GlobalLowestUWChiSqGlobalLowestUWChiSq value to beat for the ' modelstr ...
            '-state model is ' num2str(GlobalLowestUWChiSq)]);
        disp(['The current chisquared_unweighted is ' num2str(chisquared_unweighted)]);
    end
    
    %% Check to see if new fit has a lower chi squared than the global lowest
    
    if Current_chisquared_unweighted < GlobalLowestUWChiSq
        fprintf(['Congratulations!!!! %f < %f. You found a new global minimum '...
            'for the %d state model. Will update excel and best fit folder \r\n'],...
            Current_chisquared_unweighted,GlobalLowestUWChiSq,modelNum);
        updateGlobalFitMode = 1;
    elseif Current_chisquared_unweighted == GlobalLowestUWChiSq
        disp(['chisquared_unweighted =  current_GlobalLowestUWChiSq.'...
            ' Will update the excel sheet and Best fit folder.']);
        updateGlobalFitMode = 1;
    elseif Current_chisquared_unweighted > GlobalLowestUWChiSq
        fprintf(['Oh NO! %f > %f. The fit was not better. Will not update excel' ...
            ' with any new fits\r\n'],Current_chisquared_unweighted,GlobalLowestUWChiSq);
        updateGlobalFitMode = 0;
    end
end


%% Write Data to excel file
if saveMode == 1
    if updateGlobalFitMode == 1
        [computer_terminal_str, terminalID, DropboxLocationPrefix] = computerMode();
        %         excelOutFolder = [DropboxLocationPrefix filesep() 'BrettsProjects' filesep() 'SSB_gp32'];
        excelOutFolder = [DropboxLocationPrefix filesep() 'BrettsProjects' ...
            filesep() 'Writing' filesep() 'gp32_Paper'];
        
        %         excelFileName = 'SSBgp32_numericalResults.xlsx';
%         excelFileName = 'gp32_numericalResults.xlsx';
        excelFileName = 'gp32_numericalResults_Corrected.xlsx';
        
        excelFilePath = [excelOutFolder filesep() excelFileName];
        switch modelNum
            case 5
                sheetOUT = 'TimeResults_5state';
            case 6
                sheetOUT = 'TimeResults_6state';
            case 6.1
                sheetOUT = 'TimeResults_6stateLoop';
        end
        disp(['Will update the excel file: ' excelFilePath ...
            ', and excel sheet: ' sheetOUT]);
        
        if exist(excelFilePath,'file') == 2
            
            switch constructName
                case 'S1S2'
                    ConstructRow = 2;
                case 'S18S2'
                    ConstructRow = 3;
                case 'S4S5'
                    ConstructRow = 4;
                case 'S19S5'
                    ConstructRow = 5;
            end
            ConstructRow = ConstructRow + 10;
            switch modelNum
                case 5
                    startExcelCol = 'D';
                    finalExcelCol = 'S';
                case 6
                    startExcelCol = 'D';
                    finalExcelCol = 'V';
                case 6.1
                    startExcelCol = 'D';
                    finalExcelCol = 'X';
            end
            rangeOUT = [startExcelCol num2str(ConstructRow) ':' finalExcelCol num2str(ConstructRow)];
            
            switch modelNum
                case 5
                    outCell = {t12,t21,t23,t24,t32,t42,t45,t54,A1,A2,A3,A4,A5,...
                        chisquared,chisquared_unweighted,date};
                case 6
                    outCell = {t12,t21,t23,t24,t32,t42,t45,t54,t56,t65,...
                        A1,A2,A3,A4,A5,A6,chisquared,chisquared_unweighted,date};
                case 6.1
                    outCell = {t12,t21,t23,t24,t32,t42,t45,t46,t54,t56,t64,t65,...
                        A1,A2,A3,A4,A5,A6,chisquared,chisquared_unweighted,date};
            end
            writecell(outCell,excelFilePath,'Sheet',sheetOUT,'Range',rangeOUT);
            disp(['Updated excel file ' excelFileName ' with the fits']);
        else
            disp('Coud not find excel file to update');
        end
    end %end of update excel mode
    
end %End of save mode

%% Update the global best fit
if updateGlobalFitMode == 1
    foutName = [save_prefix 'GlobalBestFitResults_' outputProgramName '.mat'];
    switch modelNum
        case 5
            saveInfo_5state(BestFitFolder_Path,foutName,x,...
                k12,k21,k23,k24,k32,k42,k45,k54,t12,t21,t23,t24,t32,t42,t45,t54,...
                A1,A2,A3,A4,A5,time,Peq,K,A,P,sigma_A,FRET_bins,yoff,zoff,...
                chisquared,chisquared_unweighted,weightingFactor_FREThist,weightingFactor_C2,weightingFactor_C4_t0);
            
        case 6
            saveInfo_6state(BestFitFolder_Path,foutName,x,...
                k12,k21,k23,k24,k32,k42,k45,k54,k56,k65,...
                t12,t21,t23,t24,t32,t42,t45,t54,t56,t65,...
                A1,A2,A3,A4,A5,A6,time,Peq,K,A,P,sigma_A,FRET_bins,yoff,zoff,...
                chisquared,chisquared_unweighted,weightingFactor_FREThist,weightingFactor_C2,weightingFactor_C4_t0);
        case 6.1
            saveInfo_6stateLoop(BestFitFolder_Path,foutName,x,...
                k12,k21,k23,k24,k32,k42,k45,k46,k54,k56,k64,k65,...
                t12,t21,t23,t24,t32,t42,t45,t46,t54,t56,t64,t65,...
                A1,A2,A3,A4,A5,A6,time,Peq,K,A,P,sigma_A,FRET_bins,yoff,zoff,...
                chisquared,chisquared_unweighted,weightingFactor_FREThist,weightingFactor_C2,weightingFactor_C4_t0);
    end
end


%%  Plot final results of the algorithm   (Part 5: Final Values)
if plotMode == 1
    %% Final Histogram Result
    if fitHistMode == 1
        figure(11);
        clf;
        
        FRETDataHistPlot(FRET_bins,targetHistogram);
        hold on;
        
        % Plot FRET Fit Histogram
        NumStates = length(Peq);
        % Plot histogram of each state
        %         lineColor = char('g','#D95319','c','m','#EDB120','k','#7E2F8E','#A2142F','#4DBEEE');
        lineColor = char('g','r','c','m','b','k');
        for state_idx = 1:NumStates
            hist_sim_num = Peq(state_idx)/(sqrt(2*pi)*sigma_A(state_idx))*exp(-((FRET_bins-A(state_idx))/(sqrt(2)*sigma_A(state_idx))).^2);
            
            histPlot = plot(FRET_bins, hist_sim_num/denom_hist_sim,...
                'color',deblank(lineColor(state_idx,:)),'LineStyle','--','LineWidth',1,'DisplayName',Peq_names(state_idx,:));
            
            hold on
        end
        lgd = legend('show');
        lgd.Location = 'eastoutside';
        lgd.FontSize = 12;
        set(gca, 'FontSize', 14);
        hold on
        
        % Plot overall histogram
        histPlotTot = plot(FRET_bins, hist_sim,...
            'color','r','LineWidth',2,'DisplayName','Total Fit');
        lgd_tot = legend('show');
        
        if saveMode == 1
            foutName = [save_prefix 'Histfit'];
           saveas(gcf,[LocalOutputFolderName filesep() foutName]);
        end
        
        if updateGlobalFitMode == 1
            vectorFigSaver(foutName,BestFitFolder_Path);
            figSaver(foutName,BestFitFolder_Path); 
        end
    end
    %% (1) Final C2 Result
    if fitC2Mode == 1
        figure(22);
        clf;
        P = K2P(K,C2_exp_x);
        [C2_sim,C2_exp_x] = PA2C2(P,A,C2_exp_x);
        C2_sim = C2_sim + yoff;
        if normalizeMode == 1
            C2_sim = C2_sim./C2_sim(1);
        end
        % Add the yoffset back in (March/April 2021)
        C2_sim = C2_sim + yoff;
        C2dataPlotter(C2_exp_x,C2_exp_y);
        C2_plot = plot(C2_exp_x,C2_sim,'r-','LineWidth',3);
        drawnow();
        axis tight;
        axis square;
        if saveMode == 1
            foutName = [save_prefix 'C2fit'];
            saveas(gcf,[LocalOutputFolderName filesep() foutName]);
        end
        
        if updateGlobalFitMode == 1
            vectorFigSaver(foutName,BestFitFolder_Path);
            figSaver(foutName,BestFitFolder_Path);
            
        end
    end
    
    
    
    %% Final C4 Result
    if fitC4Mode == 1
        tau2 = 0;
        P = K2P(K,C4_tau1range);
        
        [C4_sim,C4_tau1range] = PAK2C4(P,A,K,C4_tau1range,tau2,zoff);
        if normalizeMode == 1
            C4_sim = C4_sim./C4_sim(1,1);
        end
        figure(33);
        clf;
        C4dataPlotter(C4_tau1range,C4_tau3range,C4_tau2eq0_exp,tau2ValSec);
        
        hold on;
        C4_plot = surf(C4_tau1range,C4_tau1range,C4_sim);
        drawnow();
        axis tight;
        axis square;
        
        if saveMode == 1
            foutName = [save_prefix 'C4fit'];
            saveas(gcf,[LocalOutputFolderName filesep() foutName]);
        end
        
        if updateGlobalFitMode == 1
            vectorFigSaver(foutName,BestFitFolder_Path);
            figSaver(foutName,BestFitFolder_Path);
            
        end
        
    end
    
    
    
    %%  Final: Display a graphic of the network        (Part 5: Final Values)
    figure(5);
    clf;
    switch modelNum
        case 5
            networkPlotter_5stateModel(x,chisquared_unweighted,Peq)
        case 6
            networkPlotter_6stateModel(x,chisquared_unweighted,Peq)
        case 6.1
            networkPlotter_6stateLoopModel(x,chisquared_unweighted,Peq)
    end
    if saveMode == 1
        foutName = [save_prefix 'Network'];
        saveas(gcf,[LocalOutputFolderName filesep() foutName]);
    end
    
    if updateGlobalFitMode == 1
            vectorFigSaver(foutName,BestFitFolder_Path);
            figSaver(foutName,BestFitFolder_Path); 
    end
    
end %End of plotMode

%% End of the program
disp(['Ending ' outputProgramName ' at ' char(datetime) ]);
