% C2 fitter for stretched exponential functions 


%__________________________________________________________________________
% AUTHOR: Jack 
%
% NAME: C2fitter_StretchExp_JM.m
%
% FUNCTION: % Fits the 2-point time correlation functions of simulated data.
%
% MODIFICATION LOG:
% BI 20181128 Creation.
% BI 20181206 Modified program to produce datestr yyyymmdd
% BI 20181227 Adds labels to the plot for transition rates
% BI 20190304 Removes the "step labels", adds option for docking
% BI 20190329 Removing the starting points to get better fits
% BI 20190406 Fixed labeling error and added a red font to text
% BI 20190415 Uses an foutName instead of fileKeywordope to fit
% BI 20190607 Making a 'findBestFit_mode' routine
% BI Taking over tcfFitter_v5.m
% CA 20200305 changed upper and lower bounds to fit within 20us-1ms (only 1-4 exps - skipped 5 exp for now)
%             changed raw data to plot from the first x value to 1ms
% JM 20210216 edited code to fit stretched exponentials rather than pure
% exponentials
%__________________________________________________________________________
% close all
clf

% tempMsg = 'Going to the folder 1ms res 3`15mer';
% disp(tempMsg);
% cd('/Users/bisraels/Dropbox/MarcusLab/Data/smData_Processed/S1S2/gp32_0p5uM/ChosenMolecules/001000us/trace2bf2ptTCF_fluc_output/tcfmat2tcfavg_output_Neq038');
% cd('/Users/bisraels/Dropbox/MarcusLab/Data/smData_Processed/S4S5/gp32_0p0uM/ChosenMolecules/000010us/trace2bf2ptTCF_fluc_output/');

%--------------------------------------------------------------------------
% Introduce the program
%--------------------------------------------------------------------------
functionName = 'C2fitS_';
disp(['Now executing the function: ' functionName]);

%--------------------------------------------------------------------------
% User Options (Plotting)
%--------------------------------------------------------------------------
docked_mode = 0;

labelPlotWithFit_mode = 1;
verboseMode = 0;

add_date_str_mode = 0;
special_folder_ending = '';


plot_raw_mode = 1;
findBestFit_mode = 1;

plot_current_fit_mode = 0;
showProgressMode = 0;

%--------------------------------------------------------------------------
% User Options (Save Preferences)
%--------------------------------------------------------------------------
saveMode = 1;
saveFigMode = 1;

%--------------------------------------------------------------------------
% User Options (Fitting)
%--------------------------------------------------------------------------
fit_PowerLaw_mode = 0;
fit_1_exp_mode = 0;
fit_2_exp_mode = 1;
fit_3_exp_mode = 0;
fit_4_exp_mode = 0;
fit_5_exp_mode = 0;

offsetMode=1; 
% time in seconds for the powerLaw fit
cutOffTime=0.4;

% Set time resolution
usres = 1000;
msres = usres*10^-3;
res = usres*1e-6;

% Set information for this construct
sample_description= ['E3cE4'];
save_prefix= ['E3cE4_' num2str(usres) 'usec'];

% STITCHED TCF PLOT OR NOT
stitchFit=0;  
simFit=1;

 %-------- Global Optimization -----------|
maxTrials = 200;
%Make accomodations based on choices
if findBestFit_mode == 0, maxTrials = 1; end

%--------------------------------------------------------------------------
% Determine Resolution based on chosen molecules folder (June 2019)
%--------------------------------------------------------------------------
% wd = pwd;
% usres_str = wd(strfind(wd,['ChosenMolecules' filesep()])+ length(['ChosenMolecules' filesep()]):strfind(wd,[filesep() 'trace2bf2ptTCF'])-3);
% usres = str2double(usres_str);
% msres = usres*10^-3;
% res = usres*1e-6;
disp(['     Setting the resolution equal to ' num2str(usres) ' microseconds (' num2str(msres) ' msec) (' num2str(res) ' sec)']);

%--------------------------------------------------------------------------
% Fit Paramters
%--------------------------------------------------------------------------
if usres == 0.1
    tcf_min_time_to_fit_sec = 1.5e-6;%1.1e-6;
    tcf_max_time_to_fit_sec = 20e-6;
elseif usres == 10
    tcf_min_time_to_fit_sec = res; % start at 20us (instead of res)
%     tcf_max_time_to_fit_sec = .001;
    tcf_max_time_to_fit_sec = 5;  % Only fit the data up to 1ms
elseif usres == 1
    tcf_min_time_to_fit_sec = 2e-6;
    tcf_max_time_to_fit_sec = 5;
elseif usres == 1000
    tcf_min_time_to_fit_sec = res;
    tcf_max_time_to_fit_sec = 5;%3;
    
    elseif usres == 100
        tcf_min_time_to_fit_sec = res;
    tcf_max_time_to_fit_sec = 5;%3;
else
    error('Set paramaters for resolution'); 
end

if tcf_min_time_to_fit_sec >= tcf_max_time_to_fit_sec
    error('Min time to fit (%f) is less than max time (%f). Remedy this')
else
    fprintf('Fitting from %f seconds to %f seconds.\r',tcf_min_time_to_fit_sec,tcf_max_time_to_fit_sec);
end
%--------------------------------------------------------------------------
% Prepare for plotting
%--------------------------------------------------------------------------
if docked_mode
    set(0,'DefaultFigureWindowStyle','docked')
else
    set(0,'DefaultFigureWindowStyle','normal')
end

%--------------------------------------------------------------------------
% Make a Sample Desctiption based on the path (11/27/17)
%--------------------------------------------------------------------------
% [sample_description,save_prefix] = sample_descriptionGetter();

%--------------------------------------------------------------------------
% Find all files matching fileEndName.
%--------------------------------------------------------------------------
fileKeyWord = 'wAVG*';
fileEndName = [fileKeyWord '.mat'];
fileNames = dir(['*' fileEndName]);
if length(fileNames) < 1
    disp(['Cannot find any files ending with the name ' fileEndName ' : Exiting program.']);
    return
else
    tcfavg_filename = fileNames(1).name;
end
parentFileName = tcfavg_filename(1:end-length('.mat'));



%--------------------------------------------------------------------------
% Make an output filename
%--------------------------------------------------------------------------
newStr = erase(fileKeyWord,'*');
foutName = newStr;


%--------------------------------------------------------------------------
% Find number of files based on the filename
%--------------------------------------------------------------------------
% Neq_loc = strfind(tcfavg_filename,'Neq');
% Neq_str = tcfavg_filename(Neq_loc+3:Neq_loc+5);
% Nfiles = str2double(Neq_str);

%--------------------------------------------------------------------------
% Create a folder to store the output
%--------------------------------------------------------------------------
% date_str = dateGetter(); %Code is below for other users.
if add_date_str_mode
    date_str = ['_' datestr(now,'yyyymmdd')];
else
    date_str = '';
end

special_folder_ending = ['_' special_folder_ending];
if length(special_folder_ending) < 2
    special_folder_ending = '';
end

if saveMode
    %     outputFolderName = [functionName '_output' '_Neq' num2str(Nfiles,'%03d') date_str special_folder_ending];
    outputFolderName = [functionName '_out' date_str special_folder_ending];
    %If the folder doesn't exist, create it
%     exist(outputFolderName,'dir') 
%     if 7 ~= exist(outputFolderName,'dir') 
        fprintf('     ***Making a folder called %s to store output.\r\n',outputFolderName);
        mkdir(outputFolderName);
%     end
end

%--------------------------------------------------------------------------
% Create a log file to hold information about the fitting procedure
%--------------------------------------------------------------------------
if saveMode
    
    fid = fopen([outputFolderName filesep() functionName '.log'],'w');
    fprintf(fid,'%s\r\n',datestr(now));
    fprintf(fid,'Fitting the file: %s\r\n',tcfavg_filename);
    
    fprintf(fid,'Xmin = %f, Xmax = %f\r\n',tcf_min_time_to_fit_sec,tcf_max_time_to_fit_sec);
end

%--------------------------------------------------------------------------
% Read the tcf file and load in the data
%--------------------------------------------------------------------------
disp(['     Opening file ' tcfavg_filename ' to retrieve raw data.']);
%tcfavg_filename is a matfile holding several variables
if stitchFit~=1
    if simFit ~=1
load(tcfavg_filename,'tauArraySec','UrightTCF','numMolecules');
[r,c] = size(tauArraySec); if c > r, tauArraySec = tauArraySec'; end %Make Array a Column Vector
timeFull = tauArraySec;
    elseif simFit==1
      load(tcfavg_filename,'tauArraySec','pairsTcfFinal','totalInList');
      numMolecules=totalInList; 
    [r,c] = size(tauArraySec); if c > r, tauArraySec = tauArraySec'; end %Make Array a Column Vector
    timeFull = tauArraySec;  
        
    end
if offsetMode
    if simFit==1
tcfavgFull = pairsTcfFinal + abs(min(pairsTcfFinal));
tcfavgFull=tcfavgFull./((tcfavgFull(2)));
tcfavgFull=tcfavgFull(1:end-3);
timeFull=timeFull(1:end-3);  
    else
tcfavgFull = UrightTCF + abs(min(UrightTCF));
tcfavgFull=tcfavgFull./((tcfavgFull(2)));
tcfavgFull=tcfavgFull(1:end-3);
timeFull=timeFull(1:end-3); 
    end
else
tcfavgFull = UrightTCF;
end
Nfiles= numMolecules; 

elseif stitchFit==1
load(tcfavg_filename);
tauArraySec=stitchTaus; 
[r,c] = size(tauArraySec); if c > r, tauArraySec = tauArraySec'; end %Make Array a Column Vector
timeFull = stitchTaus;
tcfavgFull = stitchTCF;
Nfiles= totalMatches; 
end

%--------------------------------------------------------------------------
% Normalization
%--------------------------------------------------------------------------
% tcfavg = tcfavg./tcfavg(1);

%%
%**************************************************************************
%**************************************************************************
%************************** Prepare Fit Input *****************************
%**************************************************************************
%**************************************************************************

%--------------------------------------------------------------------------
% Curtail the TCF to tcf_max_time_to_fit_sec
%--------------------------------------------------------------------------
point_to_start_tcfavg = find(tauArraySec >= tcf_min_time_to_fit_sec,1,'first');
time_to_start_fitting = tauArraySec(point_to_start_tcfavg);

point_to_end_tcfavg = find(timeFull <= tcf_max_time_to_fit_sec,1,'last');
time_to_end_fitting = timeFull(point_to_end_tcfavg);

fprintf('     The 2point TCF is being fit from t = %f sec to t = %f seconds.\r\n',time_to_start_fitting,time_to_end_fitting);
if saveMode
    fprintf(fid,'     The 2point TCF is being fit from t = %f sec to t = %f seconds.\r\n',time_to_start_fitting,time_to_end_fitting);
end
time = tauArraySec(point_to_start_tcfavg:point_to_end_tcfavg);
tcfavg = tcfavgFull(point_to_start_tcfavg:point_to_end_tcfavg);


%--------------------------------------------------------------------------
%  Save the TCF as a two column file: time & C^(2)(\tau)
%--------------------------------------------------------------------------
[xData, yData] = prepareCurveData( time,tcfavg );
if saveMode
    foutName_prefix =  [save_prefix foutName '_FitInput'];
    fprintf(fid,['\t Saving the tcfavg as %s.dat in %s\r\n'],foutName_prefix,outputFolderName);
    save([outputFolderName filesep() foutName_prefix '.mat'],...
        'time','tcfavg','timeFull','tcfavgFull','xData','yData','Nfiles');
end
%%
%--------------------------------------------------------------------------
%  Plot the Raw Data
%--------------------------------------------------------------------------
if plot_raw_mode
    if ~docked_mode
        %         fig  = figure;
        fig.Position = [4 442 579 550];%Left side of screen
    end
%     set(gcf,'Color',[1 1 1]);
    set(gcf,'Name','Step 1: time_tcfavg to be fit: Curtailed Raw Data');
    
    
    scatter(xData,yData,'b.')
    xlim([min(xData) max(xData)])
    title_str = [sample_description '{\color{Black}  N = ' num2str(Nfiles) '  }{\color{red}Res = ' num2str(num2str(usres)) '\musec}'];
    title(title_str,'FontSize',18);
    xlabel('\tau (sec)','fontsize',16);
    ylabel('C^{(2)}(\tau)','fontsize',16);
    set(gca,'yscale','log');
    set(gca,'xscale','log');
    set(gca,'FontSize',14);
    grid on
    axis tight;
    xlim([min(xData) max(xData)]);
    drawnow();
    %-------------------------------------------------------------------------
    % Save the Raw TCF Figure
    %-------------------------------------------------------------------------
    if saveFigMode == 1
        fig = gcf;
        fig.PaperPositionMode = 'auto';
        extension = 'png';
        foutName_prefix =  [foutName '_FitInput'];
        disp(['     Saving figure as ' foutName_prefix '.' extension ' in ' outputFolderName]);
        print(fig,[outputFolderName filesep() foutName_prefix],['-d' extension],'-r0');
        saveas(gcf,[outputFolderName filesep() foutName_prefix],['fig']);
    end
end


%**************************************************************************
%****************** Fit to a power law decay *********************
%**************************************************************************
%**************************************************************************
if fit_PowerLaw_mode
    % Set up fittype and options.
    [Ival, ~]= find(xData>cutOffTime); 
    lastPt=Ival(1); 
    xData=xData(1:lastPt);
    yData=yData(1:lastPt);
    [xData, yData] = prepareCurveData( time(1:lastPt),tcfavg(1:lastPt) );
    fitFunction = 'A*(x.^(-k)) + yoff';
    if saveMode
        
        fprintf(fid,'Fitting a function to the raw data: y(x) = %s\r\n',fitFunction);
    end
    fprintf('Fitting a function to the raw data: y(x) = %s\r\n',fitFunction);
    
    ft = fittype( fitFunction, 'independent', 'x', 'dependent', 'y' );
    
    % Choose Fitting options
    chosen_method = 'NonlinearLeastSquares';
    opts = fitoptions( 'Method', chosen_method );
    opts.Display = 'off';
    if saveMode
        fprintf(fid,'\t Fit options: %s\r\n',chosen_method);
    end
    fprintf('\t Fit options: %s\n',chosen_method);
    
    
    A_LB = 0;
    A_UB = 2*max(yData);
    A_SP = (A_LB + A_UB)/2;
    
    k_LB = 0;
    k_UB = 5;
    k_SP = (k_LB + k_UB)/2;
    
    yoff_LB = 0;
    yoff_UB = max(yData);
    yoff_SP = yData(end);
    
    
    %*** THESE MUT BE IN ALPHABETICAL ORDER ***
    opts.Lower =      [A_LB  k_LB  yoff_LB];
    opts.Upper =      [A_UB  k_UB  yoff_UB];
    
    % Fit model to data.
    
    randmult = 5;
    for i = 1:maxTrials
        A_SPm = A_SP + A_SP*randn*randmult;
        k_SPm = k_SP + k_SP*randn*randmult;      
        yoff_SPm = yoff_SP + yoff_SP*randn*randmult;
        
        opts.StartPoint = [A_SPm  k_SPm  yoff_SPm];
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Do the Fit
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %         opts.StartPoint = [A1_SP t1_SP  yoff_SP];
        [fitresult, gof] = fit( xData, yData, ft, opts );
        rsquare = gof.rsquare;
        rsquareNew = rsquare;
        if verboseMode
            disp(['i = ' num2str(i) '/' num2str(maxTrials)]);%
            disp(['    rsquare = ' num2str(rsquare)]);
            disp(fitresult);
        end
        
        if showProgressMode
            figure(10)
            plot(i,rsquare,'--gs',...
                'LineWidth',2,...
                'MarkerSize',10,...
                'MarkerEdgeColor','b',...
                'MarkerFaceColor',[0.5,0.5,0.5])
            ylim([.9,1]);
            hold on;
        end
        %-------------------------------------------------------------------------
        % Save the fit values
        %-------------------------------------------------------------------------
        
        A = fitresult.A;
        k= fitresult.k;
        yoff = fitresult.yoff;
        
        y = A*(time.^(-k)) + yoff;
        
        if saveMode
            fitFunctionName = 'PowerLaw_fitresult';
            foutName_prefix = [save_prefix foutName '_' fitFunctionName];
            filepath_fit = [outputFolderName filesep() foutName_prefix '.mat'];
            if findBestFit_mode
                if exist(filepath_fit,'file') ~= 2
                    disp('Did not find the file. Saving the initial fit');
                    iterNumber = 1;
                    save(filepath_fit,'fitresult','gof','A','k','yoff','time','y','rsquare','iterNumber');
                else
                    load(filepath_fit,'rsquare','iterNumber');
                    iterNumber = iterNumber + 1 ;
                    %                     disp(['Already found a fit. IterNumber = ' num2str(iterNumber)]);
                    if rsquareNew > rsquare
                        disp(['*** Found a better fit: ' num2str(rsquareNew) ' > ' num2str(rsquare) '. Saving.']);
                        rsquare = rsquareNew;
                        save(filepath_fit,'fitresult','gof','A','k','yoff','yoff','time','y','rsquare','iterNumber');
                    else
                        %                         disp(['     appending the iteration number to' filepath]);
                        save(filepath_fit,'iterNumber','-append')
                    end
                end
            else
                disp(['     Saving the fit results as ' foutName_prefix '.mat in ' outputFolderName]);
                save(filepath_fit,'fitresult','gof','A','k','yoff','time','y','rsquare','iterNumber');
            end
        end
        fprintf('%d/%d fitresult #%d: A = %f, k=%f, yoff = %f\n',i,maxTrials,iterNumber,A,k,yoff);
    end
    
     if ~docked_mode
%         fig.Position = [4 442 579 550];%Left side of screen
    end
    fig = figure(1);
    fig.Name = 'Power Law Fit with y-intercept';
    
    plot(xData,yData,'b.','MarkerSize',10,'color','blue','DisplayName','C^{(2)}(\tau)')
    hold on;
    if plot_current_fit_mode
        plot(time,y,'color','green','LineWidth',2,'DisplayName',['Current Fit R^2 = ' num2str(gof.rsquare,'%.3f')]);
        hold on;
        if labelPlotWithFit_mode
            
        text((time(end)/50),0.5,['{\bf\leftarrow}' 'k=' num2str(k,'%2.1f') ' , Amp=' num2str(A)],'Color','black','FontSize',14);      
                       
        end
    end
    load(filepath_fit,'fitresult','gof','A','k','yoff','time','y','rsquare','iterNumber');
    plot(time,y,'color','red','LineWidth',2,'DisplayName',['Best Fit R^2 = ' num2str(gof.rsquare,'%.3f')]);
    hold on;
    
    legend('show');
    lgd = legend;
    lgd.Location = 'NorthEast';
    lgd.FontSize = 12;
    
    
    title_str = [sample_description '{\color{Black}  N = ' num2str(Nfiles) '  }{\color{red}Res = ' num2str(num2str(usres)) '\musec}'...
        10 'y = ' fitFunction];
    title(title_str,'FontSize',16);
    xlabel('\tau (sec)','fontsize',16);
    ylabel('C^{(2)}(\tau)','fontsize',16);
    
    
    % axis tight;
%     whitebg;
    grid on;
    set(gca,'yscale','log');
    set(gca,'xscale','log');
    set(gca,'FontSize',14);
    axis tight;
    axis square
    
    if labelPlotWithFit_mode
                       
        text((time(end)/50),0.5,['{\bf\leftarrow}' 'k=' num2str(k,'%2.1f') ' , Amp=' num2str(A)],'Color','black','FontSize',14);      
        
        drawnow();
        
        %-------------------------------------------------------------------------
        % Save the Average TCF overlayed with the Exp Fit
        %-------------------------------------------------------------------------
        if saveFigMode
            fitFunctionName = 'PowerLaw_fitresult';
            foutName_prefix = [save_prefix foutName '_' fitFunctionName];
            fig = gcf;
            fig.PaperPositionMode = 'auto';
            extension = 'png';
            disp(['     Saving figure as ' foutName_prefix '.' extension ' in ' outputFolderName]);
            print(fig,[outputFolderName filesep() foutName_prefix],['-d' extension],'-r0');
            saveas(gcf,[outputFolderName filesep() foutName_prefix],['fig']);
            
        end
        
        plot_error_mode = 0;
        if plot_error_mode
            rectangle('Position',[x_begin,y_begin,x_width,y_height],'FaceColor',[0 .5 .5])
        end
    end
end 
%%
%**************************************************************************
%**************************************************************************
%****************** Fit to a single stretched exponential decay *********************
%**************************************************************************
%**************************************************************************
if fit_1_exp_mode
    %----------------------------------------------------------------------
    %  Fit the Time Correlation Function to a biexponential
    %----------------------------------------------------------------------
    
    % Set up fittype and options.
    [xData, yData] = prepareCurveData( time,tcfavg );
    fitFunction = 'A1*exp(-(x/t1).^b1) + yoff';
    if saveMode
        
        fprintf(fid,'Fitting a function to the raw data: y(x) = %s\r\n',fitFunction);
    end
    fprintf('Fitting a function to the raw data: y(x) = %s\r\n',fitFunction);
    
    ft = fittype( fitFunction, 'independent', 'x', 'dependent', 'y' );
    
    % Choose Fitting options
    chosen_method = 'NonlinearLeastSquares';
    opts = fitoptions( 'Method', chosen_method );
    opts.Display = 'off';
    if saveMode
        fprintf(fid,'\t Fit options: %s\r\n',chosen_method);
    end
    fprintf('\t Fit options: %s\n',chosen_method);
    
    %--------------------------------------------------------------------------
    % Assign the Fit Paramaters: fit options
    % LowerBound (LB), Upper Bound (UB) and StartingPoint (SP)
    % fitFunction = 'A1*exp(-(x/t1).^b1)+e';
    %--------------------------------------------------------------------------
    
    %A1 is the amplitude of the first decay, t1 is the time constant, b1 is
    %the stretching exponent
    A1_LB = 0;
    A1_UB = 2*max(yData);%2*max(yData);
    A1_SP = max(yData);%(A1_LB + A1_UB)/2;
    
    t1_LB = tcf_min_time_to_fit_sec;
    t1_UB = tcf_max_time_to_fit_sec;
    t1_SP = (t1_LB + t1_UB)/2;
    
    b1_LB=0.1;
    b1_UB=0.99;
    b1_SP=(b1_LB + b1_UB)/2;
    
    yoff_LB = 0;
    yoff_UB = max(yData);
    yoff_SP = yData(end);
    
    %*** THESE MUT BE IN ALPHABETICAL ORDER ***
    opts.Lower =      [A1_LB  b1_LB  t1_LB  yoff_LB];
    opts.Upper =      [A1_UB  b1_UB  t1_UB  yoff_UB];
    
    % Fit model to data.
    
    randmult = 1;
    for i = 1:maxTrials
        A1_SPm = A1_SP + A1_SP*randn*randmult;
        t1_SPm = t1_SP + t1_SP*randn*randmult;
        b1_SPm = b1_SP + b1_SP*randn*randmult;
        yoff_SPm = yoff_SP + yoff_SP*randn*randmult;
        
        opts.StartPoint = [A1_SPm  b1_SPm  t1_SPm  yoff_SPm];
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Do the Fit
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %         opts.StartPoint = [A1_SP t1_SP  yoff_SP];
        [fitresult, gof] = fit( xData, yData, ft, opts );
        rsquare = gof.rsquare;
        rsquareNew = rsquare;
        if verboseMode
            disp(['i = ' num2str(i) '/' num2str(maxTrials)]);%
            disp(['    rsquare = ' num2str(rsquare)]);
            disp(fitresult);
        end
        
        if showProgressMode
            figure(10)
            plot(i,rsquare,'--gs',...
                'LineWidth',2,...
                'MarkerSize',10,...
                'MarkerEdgeColor','b',...
                'MarkerFaceColor',[0.5,0.5,0.5])
            ylim([.9,1]);
            hold on;
        end
        %-------------------------------------------------------------------------
        % Save the fit values
        %-------------------------------------------------------------------------
        
        A1 = fitresult.A1;
        b1= fitresult.b1;
        t1 = fitresult.t1;
        yoff = fitresult.yoff;
        
        y = A1*exp(-(time/t1).^b1) + yoff;
        
        if saveMode
            fitFunctionName = '1exp_fitresult';
            foutName_prefix = [save_prefix foutName '_' fitFunctionName];
            filepath_fit = [outputFolderName filesep() foutName_prefix '.mat'];
            if findBestFit_mode
                if exist(filepath_fit,'file') ~= 2
                    disp('Did not find the file. Saving the initial fit');
                    iterNumber = 1;
                    save(filepath_fit,'fitresult','gof','A1','b1','t1','yoff','time','y','rsquare','iterNumber');
                else
                    load(filepath_fit,'rsquare','iterNumber');
                    iterNumber = iterNumber + 1 ;
                    %                     disp(['Already found a fit. IterNumber = ' num2str(iterNumber)]);
                    if rsquareNew > rsquare
                        disp(['*** Found a better fit: ' num2str(rsquareNew) ' > ' num2str(rsquare) '. Saving.']);
                        rsquare = rsquareNew;
                        save(filepath_fit,'fitresult','gof','A1','b1','t1','yoff','time','y','rsquare','iterNumber');
                    else
                        %                         disp(['     appending the iteration number to' filepath]);
                        save(filepath_fit,'iterNumber','-append')
                    end
                end
            else
                disp(['     Saving the fit results as ' foutName_prefix '.mat in ' outputFolderName]);
                save(filepath_fit,'fitresult','gof','A1','b1','t1','yoff','time','y','rsquare','iterNumber');
            end
        end
        fprintf('%d/%d fitresult #%d: A1 = %f, b1=%f, t1 = %f, yoff = %f\n',i,maxTrials,iterNumber,A1,b1,t1,yoff);
    end
%     save(filepath_fit,'xData','yData','opts','-append');
    %-------------------------------------------------------------------------
    % Plot the Fit with Data
    %-------------------------------------------------------------------------
    
    if ~docked_mode
%         fig.Position = [4 442 579 550];%Left side of screen
    end
    fig = figure(1);
    fig.Name = '1exp Fit with y-intercept';
    
    plot(xData,yData,'b.','MarkerSize',10,'color','blue','DisplayName','C^{(2)}(\tau)')
    hold on;
    if plot_current_fit_mode
        plot(time,y,'color','green','LineWidth',2,'DisplayName',['Current Fit R^2 = ' num2str(gof.rsquare,'%.3f')]);
        hold on;
        if labelPlotWithFit_mode
            if t1 < 1e-3
                text(t1,fitresult(t1),['{\bf\leftarrow}' num2str(t1 *1e6,'%2.1f') ' \musec'],'Color','black','FontSize',14);
            else
                text(t1,fitresult(t1),['{\bf\leftarrow}' num2str(t1*1e3,'%2.1f') ' msec'],'Color','black','FontSize',14);
            end
            
        end
    end
    load(filepath_fit,'fitresult','gof','A1','b1','t1','yoff','time','y','rsquare','iterNumber');
    plot(time,y,'color','red','LineWidth',2,'DisplayName',['Best Fit R^2 = ' num2str(gof.rsquare,'%.3f')]);
    hold on;
    
    legend('show');
    lgd = legend;
    lgd.Location = 'NorthEast';
    lgd.FontSize = 12;
    
    
    title_str = [sample_description '{\color{Black}  N = ' num2str(Nfiles) '  }{\color{red}Res = ' num2str(num2str(usres)) '\musec}'...
        10 'y = ' fitFunction];
    title(title_str,'FontSize',16);
    xlabel('\tau (sec)','fontsize',16);
    ylabel('C^{(2)}(\tau)','fontsize',16);
    
    
    % axis tight;
%     whitebg;
    grid on;
    set(gca,'yscale','log');
    set(gca,'xscale','log');
    set(gca,'FontSize',14);
    axis tight;
    axis square
    
    if labelPlotWithFit_mode
        if t1 < 1e-3
            text(t1,fitresult(t1),['{\bf\leftarrow}' num2str(t1 *1e6,'%2.1f') ' \musec, Amp=' num2str(A1)],'Color','black','FontSize',14);
        else
            text(t1,fitresult(t1),['{\bf\leftarrow}' num2str(t1*1e3,'%2.1f') ' msec, Amp=' num2str(A1)],'Color','black','FontSize',14);
        end
        
        drawnow();
        
        %-------------------------------------------------------------------------
        % Save the Average TCF overlayed with the Exp Fit
        %-------------------------------------------------------------------------
        if saveFigMode
            fitFunctionName = '1exp_fitresult';
            foutName_prefix = [save_prefix foutName '_' fitFunctionName];
            fig = gcf;
            fig.PaperPositionMode = 'auto';
            extension = 'png';
            disp(['     Saving figure as ' foutName_prefix '.' extension ' in ' outputFolderName]);
            print(fig,[outputFolderName filesep() foutName_prefix],['-d' extension],'-r0');
            saveas(gcf,[outputFolderName filesep() foutName_prefix],['fig']);
            
        end
        
        plot_error_mode = 0;
        if plot_error_mode
            rectangle('Position',[x_begin,y_begin,x_width,y_height],'FaceColor',[0 .5 .5])
        end
    end
end
%%
%**************************************************************************
%**************************************************************************
%******************** Fit to a biexponential decay ************************
%**************************************************************************
%**************************************************************************
if fit_2_exp_mode
    fitFunction = 'A1*exp(-(x/t1).^(b1)) + A2*exp(-(x/t2).^(b2)) + yoff';
    %----------------------------------------------------------------------
    % Fit the Time Correlation Function to a biexponential + yoffset
    %----------------------------------------------------------------------
    ft = fittype( fitFunction, 'independent', 'x', 'dependent', 'y' );
    chosen_method = 'NonlinearLeastSquares';
    opts = fitoptions( 'Method', chosen_method );
    opts.Display = 'off';
    
    if saveMode
        fprintf(fid,'Fitting data to the function: y(x) = %s\r\n',fitFunction);
        fprintf(fid,'\t Fit options: %s\r\n',chosen_method);
    end
    fprintf('Will fit raw data to the function: y(x) = %s\r\n',fitFunction);
    fprintf('\t Fit options: %s\n',chosen_method);
    
    %--------------------------------------------------------------------------
    % Assign the Fit Paramaters: fit options
    % LowerBound (LB), Upper Bound (UB) and StartingPoint (SP)
    % fitFunction = 'a*exp(-x/b)+c*exp(-x/d)+e';
    %--------------------------------------------------------------------------
    
    %A1 is the amplitude of the first decay, t1 is the time constant
    A1_LB = 0;
    A1_UB = max(yData);
    A1_SP = (A1_LB + A1_UB)/2;
    
    t1_LB = tcf_min_time_to_fit_sec;
    t1_UB = 100*tcf_min_time_to_fit_sec;  %t1_UB = abs(t1_UB + .2*randn);
    t1_SP = (t1_LB + t1_UB)/2;
    
    b1_LB=0.1;
    b1_UB=1.0;
    b1_SP=(b1_LB + b1_UB)/2;
    
    %A2 is the amplitude of the second decay, t2 is its time constant
    A2_LB = 0;
    A2_UB = max(yData);
    A2_SP = (A2_LB + A2_UB)/2;
    
    t2_LB = t1_UB;
    t2_UB = tcf_max_time_to_fit_sec;
    t2_SP = (t2_LB + t2_UB)/2;
    
    b2_LB=0.1;
    b2_UB= 1.0;
    b2_SP=(b2_LB + b2_UB)/2;
    
    yoff_LB = 0;
    yoff_UB = 0.05*(max(yData));
    yoff_SP = yData(end);
    %*** THESE MUT BE IN ALPHABETICAL ORDER ***
    opts.Lower =      [A1_LB A2_LB b1_LB b2_LB t1_LB t2_LB yoff_LB];
    %        opts.StartPoint = [A1_SP A2_SP t1_SP t2_SP yoff_SP];
    opts.Upper =      [A1_UB A2_UB b1_UB b2_UB t1_UB t2_UB yoff_UB];
    
    randmult = 5;
    for i = 1:maxTrials
        A1_SPm = A1_SP + A1_SP*randn*randmult;
        t1_SPm = t1_SP + t1_SP*randn*randmult;
        b1_SPm = b1_SP + b1_SP*randn*randmult;
        A2_SPm = A2_SP + A2_SP*randn*randmult;
        t2_SPm = t2_SP + t2_SP*randn*randmult;
        b2_SPm = b2_SP + b2_SP*randn*randmult;
        yoff_SPm = yoff_SP + yoff_SP*randn*randmult;
        
        opts.StartPoint = [A1_SPm b1_SPm t1_SPm  A2_SPm b2_SPm t2_SPm yoff_SPm];
        
        %*** THESE MUT BE IN ALPHABETICAL ORDER ***
        %         opts.StartPoint = [A1_SP A2_SP t1_SP t2_SP yoff_SP]*randn;
        [fitresult, gof] = fit( xData, yData, ft, opts );
        rsquare = gof.rsquare;
        rsquareNew = rsquare;
        if verboseMode
            disp(['i = ' num2str(i) '/' num2str(maxTrials)]);%
            %             disp(['    rsquare = ' num2str(rsquare)]);
            disp(fitresult);
        end
        
        if showProgressMode
            figure(10)
            plot(i,rsquare,'--gs','LineWidth',2,'MarkerSize',10,'MarkerEdgeColor','b','MarkerFaceColor',[0.5,0.5,0.5]);
            ylim([.99,1]);
            disp(['    rsquare = ' num2str(rsquare)]);
            hold on;
        end
        %-------------------------------------------------------------------------
        % Save the fitresult and fit-values
        %-------------------------------------------------------------------------
        A1 = fitresult.A1;
        A2 = fitresult.A2;
        b1 = fitresult.b1;
        b2 = fitresult.b2;
        t1 = fitresult.t1;
        t2 = fitresult.t2;
        yoff = fitresult.yoff;
        
        y = A1*exp(-(time/t1).^b1) + A2*exp(-(time/t2).^b2) + yoff;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % NEW JUNE 9 2019 randomized global optimization
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if saveMode
            fitFunctionName = '2exp_fitresult';
            foutName_prefix = [save_prefix foutName '_' fitFunctionName];
            filepath_fit = [outputFolderName filesep() foutName_prefix '.mat'];
            if findBestFit_mode
                if exist(filepath_fit,'file') ~= 2
                    disp('Did not find the file. Saving the initial fit');
                    iterNumber = 1;
                    save(filepath_fit,'fitresult','gof','A1','t1','b1','A2','b2','t2','yoff','time','y','rsquare','iterNumber');
                else
                    load(filepath_fit,'rsquare','iterNumber');
                    iterNumber = iterNumber + 1 ;
                    %                     disp(['Already found a fit. IterNumber = ' num2str(iterNumber)]);
                    if rsquareNew > rsquare
                        disp(['***Found a better fit: ' num2str(rsquareNew) ' > ' num2str(rsquare) '. Saving.']);
                        rsquare = rsquareNew;
                        save(filepath_fit,'fitresult','gof','A1','t1','b1','A2','b2','t2','yoff','time','y','rsquare','iterNumber');
                    else
                        save(filepath_fit,'iterNumber','-append')
                    end
                end
                fprintf('%d/%d fitresult #%d: A1 = %f, b1 = %f, t1 = %f,  A2 = %f, b2 = %f, t2 = %f, yoff = %f\n',i,maxTrials,iterNumber,A1,b1,t1,A2,b2,t2,yoff);
                
            else
                disp(['     Saving the fit results as ' foutName_prefix '.mat in ' outputFolderName]);
                save(filepath_fit,'fitresult','gof','A1','t1','b1','A2','b2','t2','yoff','time','y','rsquare','iterNumber');
                fprintf('\t fitresults: A1 = %f, b1 = %f, t1 = %f,  A2 = %f, b2 = %f, t2 = %f, yoff = %f\n',A1,b1,t1,A2,b2,t2,yoff);
            end
        end
    end
    if saveMode
    save(filepath_fit,'xData','yData','opts','-append');
    end
    %-------------------------------------------------------------------------
    % Plot the Fit with Data
    %-------------------------------------------------------------------------
    if ~docked_mode
%         fig.Position = [4 442 579 550];%Left side of screen
    end
    fig = figure(2);
    fig.Name = '2exp Fit with y-intercept';
    
    plot(xData,yData,'b.','MarkerSize',10,'color','blue','DisplayName','C^{(2)}(\tau)')
    hold on;
    if plot_current_fit_mode == 1
        plot(time,y,'color','green','LineWidth',2,'DisplayName',['Current Fit R^2 = ' num2str(gof.rsquare,'%.3f')]);
        hold on;
        if labelPlotWithFit_mode
            if t1 < 1e-3
                text(t1,fitresult(t1),['{\bf\leftarrow}' num2str(t1 *1e6,'%2.1f') ' \musec, Amp=' num2str(A1) ' b=' num2str(b1)],'Color','black','FontSize',14);
            else
                text(t1,fitresult(t1),['{\bf\leftarrow}' num2str(t1*1e3,'%2.1f') ' msec, Amp=' num2str(A1) ' b=' num2str(b1)],'Color','black','FontSize',14);
            end
            if t2 < 1e-3
                text(t2,fitresult(t2),['{\bf\leftarrow}' num2str(t2*1e6,'%2.1f') ' \musec, Amp=' num2str(A2) ' b=' num2str(b2)],'Color','black','FontSize',14);
            else
                text(t2,fitresult(t2),['{\bf\leftarrow}' num2str(t2*1e3,'%2.1f') ' msec, Amp=' num2str(A2) ' b=' num2str(b2)],'Color','black','FontSize',14);
            end
        end
    end
    
     fitFunctionName = '2exp_fitresult';
            foutName_prefix = [save_prefix foutName '_' fitFunctionName];
outputFolderName = [functionName '_out' date_str special_folder_ending];
            filepath_fit = [outputFolderName filesep() foutName_prefix '.mat'];
            
    load(filepath_fit,'fitresult','gof','A1','b1','t1','A2','b2','t2','yoff','time','y','rsquare','iterNumber');
    
    plot(time,y,'color','red','LineWidth',2,'DisplayName',['Best Fit R^2 = ' num2str(gof.rsquare,'%.3f')]);
    hold on;
    
    legend('show');
    lgd = legend;
    lgd.Location = 'NorthEast';
    lgd.FontSize = 12;
    
    
    title_str = [sample_description '{\color{Black}  N = ' num2str(Nfiles) '  }{\color{red}Res = ' num2str(num2str(usres)) '\musec}'...
        10 'y = ' fitFunction];
    title(title_str,'FontSize',16);
    xlabel('\tau (sec)','fontsize',16);
    ylabel('C^{(2)}(\tau)','fontsize',16);
    grid on;
    
    % axis tight;
%     whitebg;
    grid on;
    set(gca,'yscale','log');
    set(gca,'xscale','log');
    set(gca,'FontSize',14);
    axis tight;
    axis square
    
    if labelPlotWithFit_mode
        if t1 < 1e-3
            text(t1,fitresult(t1),['{\bf\leftarrow}' num2str(t1 *1e6,'%2.1f') ' \musec, A=' num2str(A1) ' b=' num2str(b1)],'Color','black','FontSize',14);
        else
            text(t1,fitresult(t1),['{\bf\leftarrow}' num2str(t1*1e3,'%2.1f') ' msec, A=' num2str(A1) ' b=' num2str(b1)],'Color','black','FontSize',14);
        end
        if t2 < 1e-3
            text(t2,fitresult(t2),['{\bf\leftarrow}' num2str(t2*1e6,'%2.1f') ' \musec, A=' num2str(A2) ' b=' num2str(b2)],'Color','black','FontSize',14);
        else
            text(t2,fitresult(t2),['{\bf\leftarrow}' num2str(t2*1e3,'%2.1f') ' msec, A=' num2str(A2) ' b=' num2str(b2)],'Color','black','FontSize',14);
        end
    end
    
    drawnow();
    
    %-------------------------------------------------------------------------
    % Save the Average TCF overlayed with the BiExpYoffset Fit
    %-------------------------------------------------------------------------
    if saveFigMode
        foutName_prefix = [save_prefix foutName '_2exp_Fit'];
        fig = gcf;
        fig.PaperPositionMode = 'auto';
        extension = 'png';
        disp(['     Saving figure as ' foutName_prefix '.' extension ' in ' outputFolderName]);
        print(fig,[outputFolderName filesep() foutName_prefix],['-d' extension],'-r0');
        saveas(gcf,[outputFolderName filesep() foutName_prefix],['fig']);
    end
end
%%
%**************************************************************************
%**************************************************************************
%******************* Fit to a Tri- exp with a y offset ********************
%**************************************************************************
%**************************************************************************
if fit_3_exp_mode == 1
    fitFunction = 'A1*exp(-(x/t1).^(b1)) + A2*exp(-(x/t2).^(b2)) + A3*exp(-(x/t3).^(b3)) + yoff';
    %----------------------------------------------------------------------
    %  Fit the Time Correlation Function to a triiexponential + yoffset
    %----------------------------------------------------------------------
    ft = fittype( fitFunction, 'independent', 'x', 'dependent', 'y' );
    chosen_method = 'NonlinearLeastSquares';
    opts = fitoptions( 'Method', chosen_method );
    opts.Display = 'Off';
    
    if saveMode
        fprintf(fid,'Fitting data to the function: y(x) = %s\r\n',fitFunction);
        fprintf(fid,'\t Fit options: %s\r\n',chosen_method);
    end
    fprintf('Will fit raw data to the function: y(x) = %s\r\n',fitFunction);
    fprintf('\t Fit options: %s\n',chosen_method);
    
    %--------------------------------------------------------------------------
    % Assign the Fit Paramaters: fit options
    % LowerBound (LB), Upper Bound (UB) and StartingPoint (SP)
    % fitFunction = 'a*exp(-x/b)+c*exp(-x/d)+e*exp(-x/f)+g';
    %--------------------------------------------------------------------------
    
    %A1 is the amplitude of the first decay, t1 is the time constant
    A1_LB = 0;
    A1_UB = max(yData);
    A1_SP = (A1_LB + A1_UB)/2;
    
    b1_LB=0.1;
    b1_UB=1.0;
    b1_SP=(b1_LB + b1_UB)/2;
    
    t1_LB = tcf_min_time_to_fit_sec;
    t1_UB = 4*tcf_min_time_to_fit_sec;
    t1_SP = (t1_LB + t1_UB)/2;
    
    %A2 is the amplitude of the second decay, t2 is its time constant
    A2_LB = 0;
    A2_UB = max(yData);
    A2_SP = (A2_LB + A2_UB)/2;
    
    b2_LB=0.1;
    b2_UB=1.0;
    b2_SP=(b2_LB + b2_UB)/2;
    
    t2_LB = t1_UB;
    t2_UB = 100*tcf_min_time_to_fit_sec;
    t2_SP = (t2_LB + t2_UB)/2;
    
    %A3 is the amplitude of the second decay, t3 is its time constant
    A3_LB = 0;
    A3_UB = max(yData);
    A3_SP = (A3_LB + A3_UB)/2;
    
    b3_LB=0.1;
    b3_UB=1.0;
    b3_SP=(b3_LB + b3_UB)/2;
    
    t3_LB = t2_UB;
    t3_UB = tcf_max_time_to_fit_sec;
    t3_SP = (t3_LB + t3_UB)/2;
    
    yoff_LB = 0;
    yoff_UB = max(yData);
    yoff_SP = yData(end);
    
    %*** THESE MUST BE IN ALPHABETICAL ORDER ***
    opts.Lower =      [A1_LB A2_LB A3_LB b1_LB b2_LB b3_LB t1_LB t2_LB t3_LB yoff_LB];
    opts.Upper =       [A1_UB A2_UB A3_UB b1_UB b2_UB b3_UB t1_UB t2_UB t3_UB yoff_UB];
    
    randmult = 5;
    for i = 1:maxTrials
        A1_SPm = A1_SP + A1_SP*randn*randmult;
        t1_SPm = t1_SP + t1_SP*randn*randmult;
        b1_SPm = b1_SP + b1_SP*randn*randmult;
        A2_SPm = A2_SP + A2_SP*randn*randmult;
        t2_SPm = t2_SP + t2_SP*randn*randmult;
        b2_SPm = b2_SP + b2_SP*randn*randmult;
        A3_SPm = A3_SP + A3_SP*randn*randmult;
        t3_SPm = t3_SP + t3_SP*randn*randmult;
        b3_SPm = b3_SP + b3_SP*randn*randmult;
        yoff_SPm = yoff_SP + yoff_SP*randn*randmult;
        
        opts.StartPoint = [A1_SPm b1_SPm t1_SPm A2_SPm b2_SPm t2_SPm A3_SPm b3_SPm t3_SPm yoff_SPm];
        %*** THESE MUT BE IN ALPHABETICAL ORDER ***
        %         opts.StartPoint =  [A1_SP A2_SP A3_SP t1_SP t2_SP t3_SP yoff_SP];
        [fitresult, gof] = fit( xData, yData, ft, opts );
        rsquare = gof.rsquare;
        rsquareNew = rsquare;
        if verboseMode
            disp(['i = ' num2str(i) '/' num2str(maxTrials)]);%
            %             disp(['    rsquare = ' num2str(rsquare)]);
            disp(fitresult);
        end
        
        if showProgressMode
            figure(10);
            plot(i,rsquare,'--gs','LineWidth',2,'MarkerSize',10,'MarkerEdgeColor','b','MarkerFaceColor',[0.5,0.5,0.5]);
            ylim([.99,1]);
            disp(['    rsquare = ' num2str(rsquare)]);
            hold on;
        end
        %-------------------------------------------------------------------------
        % Save the fitresult and fit-values
        %-------------------------------------------------------------------------
        A1 = fitresult.A1;
        A2 = fitresult.A2;
        A3 = fitresult.A3;
        t1 = fitresult.t1;
        t2 = fitresult.t2;
        t3 = fitresult.t3;
        b1 = fitresult.b1;
        b2 = fitresult.b2;
        b3 = fitresult.b3;
        yoff = fitresult.yoff;
        
        y = A1*exp(-(time/t1).^b1) + A2*exp(-(time/t2).^b2) + A3*exp(-(time/t3).^b3)+yoff;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % NEW JUNE 9 2019 randomized global optimization
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if saveMode
            fitFunctionName = '3exp_fitresult';
            foutName_prefix = [save_prefix foutName '_' fitFunctionName];
            filepath_fit = [outputFolderName filesep() foutName_prefix '.mat'];
            if findBestFit_mode
                if exist(filepath_fit,'file') ~= 2
                    disp('Did not find the file. Saving the initial fit');
                    iterNumber = 1;
                    save(filepath_fit,'fitresult','gof','A1','t1','b1','A2','t2','b2','A3','t3','b3','yoff','time','y','rsquare','iterNumber');
                else
                    load(filepath_fit,'rsquare','iterNumber');
                    iterNumber = iterNumber + 1 ;
                    %                     disp(['Already found a fit. IterNumber = ' num2str(iterNumber)]);
                    if rsquareNew > rsquare
                        disp(['***Found a better fit: ' num2str(rsquareNew) ' > ' num2str(rsquare) '. Saving.']);
                        rsquare = rsquareNew;
                        save(filepath_fit,'fitresult','gof','A1','t1','b1','A2','t2','b2','A3','t3','b3','yoff','time','y','rsquare','iterNumber');
                    else
                        save(filepath_fit,'iterNumber','-append')
                    end
                end
                fprintf('%d/%d fitresult #%d: A1 = %f, b1 = %f, t1 = %f, A2 = %f, b2 = %f, t2 = %f, A3 = %f, b3 = %f, t3 = %f,yoff = %f\n',i,maxTrials,iterNumber,A1,b1,t1,A2,b2,t2,A3,b3,t3,yoff);
            else
                disp(['     Saving the fit results as ' foutName_prefix '.mat in ' outputFolderName]);
                fprintf('\tfitresults: A1 = %f, b1 = %f, t1 = %f, A2 = %f, b2 = %f, t2 = %f, A3 = %f, b3 = %f, t3 = %f, yoff = %f\r',A1,b1,t1,A2,b2,t2,A3,b3,t3,yoff);
                save(filepath_fit,'fitresult','gof','A1','t1','b1','A2','t2','b2','A3','t3','b3','yoff','time','y','rsquare','iterNumber');
            end
        end
    end
    save(filepath_fit,'xData','yData','opts','-append');
    %-------------------------------------------------------------------------
    % Plot the Fit with Data
    %-------------------------------------------------------------------------
    if ~docked_mode
%         fig.Position = [4 442 579 550];%Left side of screen
    end
    fig = figure(3);
    set(gcf,'Name','3exp Fit with y-intercept');
    
    plot(xData,yData,'b.','MarkerSize',10,'color','blue','DisplayName','C^{(2)}(\tau)')
    hold on;
    if plot_current_fit_mode == 1
        currFitPlot = plot(time,y,'color','green','LineWidth',2,'DisplayName',['Current Fit R^2 = ' num2str(gof.rsquare,'%.3f')]);
        hold on;
        if labelPlotWithFit_mode
            if t1 < 1e-3
                text(t1,fitresult(t1),['{\bf\leftarrow}' num2str(t1 *1e6,'%2.1f') ' \musec'],'Color','black','FontSize',14);
            else
                text(t1,fitresult(t1),['{\bf\leftarrow}' num2str(t1*1e3,'%2.1f') ' msec'],'Color','black','FontSize',14);
            end
            if t2 < 1e-3
                text(t2,fitresult(t2),['{\bf\leftarrow}' num2str(t2*1e6,'%2.1f') ' \musec'],'Color','black','FontSize',14);
            else
                text(t2,fitresult(t2),['{\bf\leftarrow}' num2str(t2*1e3,'%2.1f') ' msec'],'Color','black','FontSize',14);
            end
        end
    end
    
    load(filepath_fit,'fitresult','gof','A1','t1','b1','A2','t2','b2','A3','t3','b3','yoff','time','y','rsquare','iterNumber');
    
    y = A1*exp(-(time/t1).^b1) + A2*exp(-(time/t2).^b2) + A3*exp(-(time/t3).^b3)+yoff;
    plot(time,y,'color','red','LineWidth',2,'DisplayName',['Best Fit R^2 =' num2str(gof.rsquare,'%.3f')]);
    hold on;
    
    legend('show');
    lgd = legend;
    lgd.Location = 'SouthWest';
    lgd.FontSize = 12;
    
    title_str = [sample_description '{\color{Black}  N = ' num2str(Nfiles) '  }{\color{red}Res = ' num2str(num2str(usres)) '\musec}'...
        10 'y = ' fitFunction];
    title(title_str,'FontSize',16);
    xlabel('\tau (sec)','fontsize',16);
    ylabel('C^{(2)}(\tau)','fontsize',16);
    grid on;
    
%     whitebg;
    grid on;
    set(gca,'yscale','log');
    set(gca,'xscale','log');
    set(gca,'FontSize',15);
    axis tight;
    axis square
    
    if labelPlotWithFit_mode
        if t1 < 1e-3
            text(t1,fitresult(t1),['{\bf\leftarrow}' num2str(t1 *1e6,'%2.1f') ' \musec, A=' num2str(A1) ' b=' num2str(b1)],'Color','black','FontSize',14);
        else
            text(t1,fitresult(t1),['{\bf\leftarrow}' num2str(t1*1e3,'%2.1f') ' msec, A=' num2str(A1) ' b=' num2str(b1)],'Color','black','FontSize',14);
        end
        
        if t2 < 1e-3
            text(t2,fitresult(t2),['{\bf\leftarrow}' num2str(t2*1e6,'%2.1f') ' \musec, A=' num2str(A2) ' b=' num2str(b2)],'Color','black','FontSize',14);
        else
            text(t2,fitresult(t2),['{\bf\leftarrow}' num2str(t2*1e3,'%2.1f') ' msec, A=' num2str(A2) ' b=' num2str(b2)],'Color','black','FontSize',14);
        end
        
        if t3 < 1e-3
            text(t3,fitresult(t3),['{\bf\leftarrow}' num2str(t3*1e6,'%2.1f') ' \musec, A=' num2str(A3) ' b=' num2str(b3)],'Color','black','FontSize',14);
        else
            text(t3,fitresult(t3),['{\bf\leftarrow}' num2str(t3*1e3,'%2.1f') ' msec, A=' num2str(A3) ' b=' num2str(b3)],'Color','black','FontSize',14);
        end
        
    end
    
    drawnow();
    %-------------------------------------------------------------------------
    % Save the Average TCF overlayed with the Triexponential Fit
    %-------------------------------------------------------------------------
    if saveFigMode == 1
        foutName_prefix = [save_prefix foutName '_3exp_Fit'];
        fig = gcf;
        fig.PaperPositionMode = 'auto';
        extension = 'png';
        disp(['     Saving figure as ' foutName_prefix '.' extension ' in ' outputFolderName]);
        print(fig,[outputFolderName filesep() foutName_prefix],['-d' extension],'-r0');
        saveas(gcf,[outputFolderName filesep() foutName_prefix],['fig']);
    end
end
% disp(fitresult);
%%
%**************************************************************************
%**************************************************************************
%**************** Fit to Four exponentials witt a y offset ****************
%**************************************************************************
%**************************************************************************
if fit_4_exp_mode
    %----------------------------------------------------------------------
    %  Fit the Time Correlation Function to a triiexponential + yoffset
    %----------------------------------------------------------------------
    
    % Set up fittype and options.
    [xData, yData] = prepareCurveData( time, tcfavg );
    fitFunction = 'A1*exp(-(x/t1).^(b1)) + A2*exp(-(x/t2).^(b2)) + A3*exp(-(x/t3).^(b3)) + A4*exp(-(x/t4).^b4) + yoff';
    if saveMode
        fprintf(fid,'Will fit raw data to the function: y(x) = %s\r\n',fitFunction);
    end
    fprintf('Will fit raw data to the function: y(x) = %s\r\n',fitFunction);
    
    ft = fittype( fitFunction, 'independent', 'x', 'dependent', 'y' );
    
    % Choose Fitting options
    chosen_method = 'NonlinearLeastSquares';
    opts = fitoptions( 'Method', chosen_method );
    opts.Display = 'Off';
    if saveMode
        fprintf(fid,'\t Fit options: %s\r\n',chosen_method);
    end
    fprintf('\t Fit options: %s\r\n',chosen_method);
    
    %--------------------------------------------------------------------------
    % Assign the Fit Paramaters: fit options
    % LowerBound (LB), Upper Bound (UB) and StartingPoint (SP)
    % fitFunction = 'a*exp(-x/b)+c*exp(-x/d)+e*exp(-x/f)+g';
    %--------------------------------------------------------------------------
    
    %A1 is the amplitude of the first decay, t1 is the time constant
    A1_LB = 0;       A1_UB = max(yData);     A1_SP = (A1_LB + A1_UB)/4;
    t1_LB = tcf_min_time_to_fit_sec;       t1_UB = 5*tcf_min_time_to_fit_sec;     t1_SP = (t1_LB + t1_UB)/2;
    b1_LB=0.1;    b1_UB=1.0;    b1_SP=(b1_LB + b1_UB)/2;  
    %t1_UB = 10*tcf_min_time_to_fit_sec;
    
    %A2 is the amplitude of the second decay, t2 is its time constant
    A2_LB = 0;       A2_UB = max(yData);    A2_SP = (A2_LB + A2_UB)/4;
    t2_LB = t1_UB;      t2_UB = 30*tcf_min_time_to_fit_sec;   t2_SP = (t2_LB + t2_UB)/2;
    b2_LB=0.1;    b2_UB=1.0;    b2_SP=(b2_LB + b2_UB)/2;
% t2_LB = 10*tcf_min_time_to_fit_sec; 
    %     t2_UB = 100*tcf_min_time_to_fit_sec;

    %A3 is the amplitude of the second decay, t3 is its time constant
    A3_LB = 0;       A3_UB = max(yData);       A3_SP = (A3_LB + A3_UB)/4;
    t3_LB = t2_UB;  t3_UB = 70000*tcf_min_time_to_fit_sec;   t3_SP = (t3_LB + t3_UB)/2;
    b3_LB=0.1;    b3_UB=1.0;    b3_SP=(b3_LB + b3_UB)/2;
% t3_LB = 100*tcf_min_time_to_fit_sec;
    %     t3_UB = 1000*tcf_min_time_to_fit_sec;

    %A4 is the amplitude of the second decay, t4 is its time constant
    A4_LB = 0;       A4_UB = max(yData);       A4_SP = (A4_LB + A4_UB)/4;
    t4_LB = t3_UB;  t4_UB = tcf_max_time_to_fit_sec;    t4_SP = (t4_LB + t4_UB)/2;
    b4_LB=0.1;    b4_UB=1.0;    b4_SP=(b4_LB + b4_UB)/2;
% t4_LB = 1000*tcf_min_time_to_fit_sec;
    
    %yoff is the y-offset
    %     yoff_LB = -1*max(yData);       yoff_UB = max(yData);   yoff_SP = (yoff_LB + yoff_UB)/2;
    yoff_LB = 0;   yoff_UB = max(yData);  yoff_SP =yData(end);
    
    %*** THESE MUT BE IN ALPHABETICAL ORDER ***
    opts.Lower =      [A1_LB  A2_LB A3_LB A4_LB b1_LB b2_LB b3_LB b4_LB t1_LB t2_LB t3_LB t4_LB yoff_LB];
    opts.Upper =       [A1_UB A2_UB A3_UB A4_UB b1_UB b2_UB b3_UB b4_UB t1_UB t2_UB t3_UB t4_UB yoff_UB];
    
    randmult = 5;
    for i = 1:maxTrials % Loop to do many fits
        A1_SPm = A1_SP + A1_SP*randn*randmult;
        t1_SPm = t1_SP + t1_SP*randn*randmult;
        b1_SPm = b1_SP + b1_SP*randn*randmult;
        A2_SPm = A2_SP + A2_SP*randn*randmult;
        t2_SPm = t2_SP + t2_SP*randn*randmult;
        b2_SPm = b2_SP + b2_SP*randn*randmult;
        A3_SPm = A3_SP + A3_SP*randn*randmult;
        t3_SPm = t3_SP + t3_SP*randn*randmult;
        b3_SPm = b3_SP + b3_SP*randn*randmult;
        A4_SPm = A4_SP + A4_SP*randn*randmult;
        t4_SPm = t4_SP + t4_SP*randn*randmult;
        b4_SPm = b4_SP + b4_SP*randn*randmult;
        yoff_SPm = yoff_SP + yoff_SP*randn*randmult;
        
        opts.StartPoint =  [A1_SPm b1_SPm t1_SPm A2_SPm b2_SPm t2_SPm A3_SPm b3_SPm t3_SPm A4_SPm b4_SPm t4_SPm yoff_SPm];
       
        disp(['i = ' num2str(i) '/' num2str(maxTrials)]);
        [fitresult, gof] = fit( xData, yData, ft, opts );
        rsquare = gof.rsquare;
        rsquareNew = rsquare;
        disp(['    rsquare = ' num2str(rsquare)]);
        disp(fitresult);
        
        %-------------------------------------------------------------------------
        % Save the fitresult and fit-values
        %-------------------------------------------------------------------------
        A1 = fitresult.A1;
        A2 = fitresult.A2;
        A3 = fitresult.A3;
        A4 = fitresult.A4;
        t1 = fitresult.t1;
        t2 = fitresult.t2;
        t3 = fitresult.t3;
        t4 = fitresult.t4;
        b1 = fitresult.b1;
        b2 = fitresult.b2;
        b3 = fitresult.b3;
        b4 = fitresult.b4;
        yoff = fitresult.yoff;
        
        y = A1*exp(-(time/t1).^b1) + A2*exp(-(time/t2).^b2) + A3*exp(-(time/t3).^b3) + A4*exp(-(time/t4).^b4) + yoff;
        
        fprintf('\tfitresults: A1 = %f, b1 = %f, t1 = %f, A2 = %f, b2 = %f, t2 = %f, A3 = %f, b3 = %f, t3 = %f, A4 = %f, b4 = %f, t4 = %f, yoff = %f\n',A1,b1,t1,A2,b2,t2,A3,b3,t3,A4,b4,t4,yoff);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % NEW JUNE 9 2019 randomized global optimization
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if saveMode
            fitFunctionName = '4exp_fitresult';
            foutName_prefix = [save_prefix foutName '_' fitFunctionName];
            filepath_fit = [outputFolderName filesep() foutName_prefix '.mat'];
            if findBestFit_mode
                if exist(filepath_fit,'file') ~= 2
                    disp('Did not find the file. Saving the initial fit');
                    iterNumber = 1;
                    save(filepath_fit,'fitresult','gof','A1','b1','t1','A2','b2','t2','A3','b3','t3','A4','b4','t4','yoff','time','y','rsquare','iterNumber');
                else
                    load(filepath_fit,'rsquare','iterNumber');
                    iterNumber = iterNumber + 1 ;
                    %                     disp(['Already found a fit. IterNumber = ' num2str(iterNumber)]);
                    if rsquareNew > rsquare
                        disp(['***Found a better fit: ' num2str(rsquareNew) ' > ' num2str(rsquare) '. Saving.']);
                        rsquare = rsquareNew;
                        save(filepath_fit,'fitresult','gof','A1','b1','t1','A2','b2','t2','A3','b3','t3','A4','b4','t4','yoff','time','y','rsquare','iterNumber');
                    else
                        save(filepath_fit,'iterNumber','-append')
                    end
                end
            else
                disp(['     Saving the fit results as ' foutName_prefix '.mat in ' outputFolderName]);
                save(filepath_fit,'fitresult','gof','A1','b1','t1','A2','b2','t2','A3','b3','t3','A4','b4','t4','yoff','time','y','rsquare','iterNumber');
            end
        end
    end
    
    %-------------------------------------------------------------------------
    % Plot the Fit with Data
    %-------------------------------------------------------------------------
    if ~docked_mode
        fig = figure;
        fig.Name = '4exp Fit with y-intercept';
        %     fig.Position = [540 405 579 550];%middle of screen
        %         fig.Position = [1102 405 579 550];%Right of screen
        %          fig.Position = [1164 442 579 550];%Right of screen
%         fig.Position = [279 74 579 550];%Bottom middle-left of screen
    end
    % h = plot( fitresult, xData, yData );
    % legend( h, 'C^{(2)}(\tau) vs. time', [fitFunction 10 'Fit R^2=' num2str(gof.rsquare)], 'Location', 'NorthEast' );
    
    
    plot(xData,yData,'b.','MarkerSize',10,'color','blue','DisplayName','C^{(2)}(\tau)')
    hold on;
    plot(time,y,'color','red','LineWidth',2,'DisplayName',['Fit R^2=' num2str(gof.rsquare,'%.3f')]);
    hold on;
    
    legend('show');
    lgd = legend;
    lgd.Location = 'NorthEast';
    lgd.FontSize = 12;
    
    title_str = [sample_description '{\color{Black}  N = ' num2str(Nfiles) '  }{\color{red}Res = ' num2str(num2str(usres)) '\musec}'...
        10 'y = ' fitFunction];
    title(title_str,'FontSize',16);
    xlabel('\tau (sec)','fontsize',16);
    ylabel('C^{(2)}(\tau)','fontsize',16);
    grid on;
%     whitebg;
    grid on;
    set(gca,'yscale','log');
    set(gca,'xscale','log');
    axis tight;
    axis square
    
    if labelPlotWithFit_mode
        if t1 < 1e-3
            text(t1,fitresult(t1),['{\bf\leftarrow}' num2str(t1 *1e6,'%2.1f') ' \musec, Amp=' num2str(A1) ' b=' num2str(b1)],'Color','black','FontSize',14);
        else
            text(t1,fitresult(t1),['{\bf\leftarrow}' num2str(t1*1e3,'%2.1f') ' msec, Amp=' num2str(A1) ' b=' num2str(b1)],'Color','black','FontSize',14);
        end
        
        if t2 < 1e-3
            text(t2,fitresult(t2),['{\bf\leftarrow}' num2str(t2*1e6,'%2.1f') ' \musec, Amp=' num2str(A2) ' b=' num2str(b2)],'Color','black','FontSize',14);
        else
            text(t2,fitresult(t2),['{\bf\leftarrow}' num2str(t2*1e3,'%2.1f') ' msec, Amp=' num2str(A2) ' b=' num2str(b2)],'Color','black','FontSize',14);
        end
        
        if t3 < 1e-3
            text(t3,fitresult(t3),['{\bf\leftarrow}' num2str(t3*1e6,'%2.1f') ' \musec, Amp=' num2str(A3) ' b=' num2str(b3)],'Color','black','FontSize',14);
        else
            text(t3,fitresult(t3),['{\bf\leftarrow}' num2str(t3*1e3,'%2.1f') ' msec, Amp=' num2str(A3) ' b=' num2str(b3)],'Color','black','FontSize',14);
        end
        
        if t4 < 1e-3
            text(t4,fitresult(t4),['{\bf\leftarrow}' num2str(t4*1e6,'%2.1f') ' \musec, Amp=' num2str(A4) ' b=' num2str(b4)],'Color','black','FontSize',14);
        elseif t4 < 1
            text(t4,fitresult(t4),['{\bf\leftarrow}' num2str(t4*1e3,'%2.1f') ' msec, Amp=' num2str(A4) ' b=' num2str(b4)],'Color','black','FontSize',14);
        else
            text(t4,fitresult(t4),['{\bf\leftarrow}' num2str(t4,'%2.1f') ' sec, Amp=' num2str(A4) ' b=' num2str(b4)],'Color','black','FontSize',14);
        end
    end
    
    drawnow();
    %-------------------------------------------------------------------------
    % Plot lines which connect to the rate in question
    %-------------------------------------------------------------------------
    %     line([t3 t3],[0 y_four(find(time <= t3,1,'last'))],'Color','red','LineStyle','-','LineWidth',2)
    % line([t2 t2],[0 y_four(find(time <= t2,1,'last'))],'Color','red','LineStyle','-','LineWidth',2)
    % line([t1 t1],[0 y_four(find(time <= t1,1,'last'))],'Color','red','LineStyle','-','LineWidth',2)
    % line([t4 t4],[0 y_four(find(time <= t4,1,'last'))],'Color','red','LineStyle','-','LineWidth',2)
    % hold on;
    % xarr = ones(1,100)*t1;
    % yarr = linspace(min(yData),y_four(find(time <= t1,1,'last')),100);
    % plot(xarr,yarr,'r-');
    %-------------------------------------------------------------------------
    % Save the Average TCF overlayed with the fourexponential Fit
    %-------------------------------------------------------------------------
    if saveFigMode == 1
        foutName_prefix = [save_prefix foutName '_4exp_Fit'];
        fig = gcf;
        fig.PaperPositionMode = 'auto';
        extension = 'png';
        disp(['     Saving figure as ' foutName_prefix '.' extension ' in ' outputFolderName]);
        print(fig,[outputFolderName filesep() foutName_prefix],['-d' extension],'-r0');
        saveas(gcf,[outputFolderName filesep() foutName_prefix],['fig']);
        
%         vectorFigSaver(foutName_prefix,outputFolderName)
        
    end
end
%%
%**************************************************************************
%**************************************************************************
%**************** Fit to Five exponentials witt a y offset ****************
%**************************************************************************
%**************************************************************************
if fit_5_exp_mode
    %----------------------------------------------------------------------
    %  Fit the Time Correlation Function to a 5-expoentnial + yoffset
    %----------------------------------------------------------------------
    
    % Set up fittype and options.
    [xData, yData] = prepareCurveData( time, tcfavg );
    fitFunction = 'A1*exp(-x/t1) + A2*exp(-x/t2) + A3*exp(-x/t3) + A4*exp(-x/t4) + A5*exp(-x/t5) + yoff';
    if saveMode
        fprintf(fid,'Will fit raw data to the function: y(x) = %s\r\n',fitFunction);
    end
    fprintf('Will fit raw data to the function: y(x) = %s\r\n',fitFunction);
    
    ft = fittype( fitFunction, 'independent', 'x', 'dependent', 'y' );
    
    % Choose Fitting options
    chosen_method = 'NonlinearLeastSquares';
    opts = fitoptions( 'Method', chosen_method );
    opts.Display = 'Off';
    if saveMode
        fprintf(fid,'\t Fit options: %s\r\n',chosen_method);
    end
    fprintf('\t Fit options: %s\r\n',chosen_method);
    
    %--------------------------------------------------------------------------
    % Assign the Fit Paramaters: fit options
    % LowerBound (LB), Upper Bound (UB) and StartingPoint (SP)
    % fitFunction = 'a*exp(-x/b)+c*exp(-x/d)+e*exp(-x/f)+g';
    %--------------------------------------------------------------------------
    
    %A1 is the amplitude of the first decay, t1 is the time constant
    A1_LB = 0;       A1_UB = max(yData);     A1_SP = (A1_LB + A1_UB)/5;
    t1_LB = tcf_min_time_to_fit_sec;       t1_UB = 2*tcf_min_time_to_fit_sec; t1_SP = (t1_LB + t1_UB)/2;
    
    %A2 is the amplitude of the second decay, t2 is its time constant
    A2_LB = 0;       A2_UB = max(yData);    A2_SP = (A2_LB + A2_UB)/5;
    t2_LB = t1_UB;       t2_UB = 30*tcf_min_time_to_fit_sec;   t2_SP = (t2_LB + t2_UB)/2;
    
    %A3 is the amplitude of the third decay, t3 is its time constant
    A3_LB = 0;       A3_UB = max(yData);       A3_SP = (A3_LB + A3_UB)/5;
    t3_LB = t2_UB;  t3_UB = 500*tcf_min_time_to_fit_sec;    t3_SP = (t3_LB + t3_UB)/2;
    
    %A4 is the amplitude of the fourth decay, t4 is its time constant
    A4_LB = 0;       A4_UB = max(yData);       A4_SP = (A4_LB + A4_UB)/5;
    t4_LB = 500*tcf_min_time_to_fit_sec;  t4_UB = 11000*tcf_min_time_to_fit_sec;    t4_SP = (t4_LB + t4_UB)/2;
    
    %A5 is the amplitude of the fifth decay, t4 is its time constant
    A5_LB = 0;       A5_UB = max(yData);       A5_SP = (A4_LB + A4_UB)/5;
    t5_LB = t4_UB;  t5_UB = tcf_max_time_to_fit_sec;    t5_SP = (t4_LB + t4_UB)/2;
    %t5_UB = tcf_max_time_to_fit_sec;
    
    %yoff is the y-offset
    yoff_LB = -1*max(yData);       yoff_UB = max(yData);       yoff_SP = (yoff_LB + yoff_UB)/2;
    
    %*** THESE MUT BE IN ALPHABETICAL ORDER ***
    opts.Lower =      [A1_LB  A2_LB A3_LB A4_LB A5_LB t1_LB t2_LB t3_LB t4_LB t5_LB yoff_LB];
    %     opts.StartPoint =  [A1_SP A2_SP A3_SP A4_SP A5_SP t1_SP t2_SP t3_SP t4_SP t5_SP yoff_SP];
    opts.Upper =       [A1_UB A2_UB A3_UB A4_UB A5_UB t1_UB t2_UB t3_UB t4_UB t5_UB yoff_UB];
    %%
    randmult=5; 
    for i = 1:maxTrials % Loop to do many fits
        opts.StartPoint =  [A1_SP A2_SP A3_SP A4_SP A5_SP t1_SP t2_SP t3_SP t4_SP t5_SP yoff_SP]*randn*randmult;
        disp(['i = ' num2str(i) '/' num2str(maxTrials)]);
        [fitresult, gof] = fit( xData, yData, ft, opts );
        rsquare = gof.rsquare;
        rsquareNew = rsquare;
        disp(['    rsquare = ' num2str(rsquare)]);
        disp(fitresult);
        
        %-------------------------------------------------------------------------
        % Save the fitresult and fit-values
        %-------------------------------------------------------------------------
        A1 = fitresult.A1;
        A2 = fitresult.A2;
        A3 = fitresult.A3;
        A4 = fitresult.A4;
        A5 = fitresult.A5;
        t1 = fitresult.t1;
        t2 = fitresult.t2;
        t3 = fitresult.t3;
        t4 = fitresult.t4;
        t5 = fitresult.t5;
        yoff = fitresult.yoff;
        
        y = A1*exp(-time/t1) + A2*exp(-time/t2) + A3*exp(-time/t3) + A4*exp(-time/t4) + A5*exp(-time/t5) + yoff;
        
        fprintf('\tfitresults: A1 = %f, t1 = %f, A2 = %f, t2= %f, A3 = %f, t3 = %f, A4 = %f, t4 = %f, A5 = %f, t5 = %f, yoff = %f\r\n',A1,t1,A2,t2,A3,t3,A4,t4,A5,t5,yoff);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % NEW JUNE 9 2019 randomized global optimization
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if saveMode
            fitFunctionName = '5exp_fitresult';
            foutName_prefix = [save_prefix foutName '_' fitFunctionName];
            filepath_fit = [outputFolderName filesep() foutName_prefix '.mat'];
            if findBestFit_mode
                if exist(filepath_fit,'file') ~= 2
                    disp('Did not find the file. Saving the initial fit');
                    iterNumber = 1;
                    save(filepath_fit,'fitresult','gof','A1','t1','A2','t2','A3','t3','A4','t4','A5','t5','yoff','time','y','rsquare','iterNumber');
                else
                    load(filepath_fit,'rsquare','iterNumber');
                    iterNumber = iterNumber + 1 ;
                    %                     disp(['Already found a fit. IterNumber = ' num2str(iterNumber)]);
                    if rsquareNew > rsquare
                        disp(['***Found a better fit: ' num2str(rsquareNew) ' > ' num2str(rsquare) '. Saving.']);
                        rsquare = rsquareNew;
                        save(filepath_fit,'fitresult','gof','A1','t1','A2','t2','A3','t3','A4','t4','A5','t5','yoff','time','y','rsquare','iterNumber');
                    else
                        save(filepath_fit,'iterNumber','-append')
                    end
                end
            else
                disp(['     Saving the fit results as ' foutName_prefix '.mat in ' outputFolderName]);
                save(filepath_fit,'fitresult','gof','A1','t1','A2','t2','A3','t3','A4','t4','A5','t5','yoff','time','y','rsquare','iterNumber');
            end
        end
    end
    %-------------------------------------------------------------------------
    % Plot the Fit with Data
    %-------------------------------------------------------------------------
    if ~docked_mode
        fig = figure;
        fig.Name = '5exp Fit with y-intercept';
        %     fig.Position = [540 405 579 550];%middle of screen
        %         fig.Position = [1102 405 579 550];%Right of screen
        %          fig.Position = [1200 408 721 696];%Right of screen
%         fig.Position = [859 74 578 550];%Bottom middle right of screen
    end
    % h = plot( fitresult, xData, yData );
    % legend( h, 'C^{(2)}(\tau) vs. time', [fitFunction 10 'Fit R^2=' num2str(gof.rsquare)], 'Location', 'NorthEast' );
    
    
    plot(xData,yData,'b.','MarkerSize',10,'color','blue','DisplayName','C^{(2)}(\tau)');
    hold on;
    plot(time,y,'color','red','LineWidth',2,'DisplayName',['Fit R^2=' num2str(gof.rsquare,'%.3f')]);
    hold on;
    
    legend('show');
    lgd = legend;
    lgd.Location = 'SouthWest';
    lgd.FontSize = 12;
    
    title_str = [sample_description '{\color{Black}  N = ' num2str(Nfiles) '  }{\color{red}Res = ' num2str(num2str(usres)) '\musec}'...
        10 'y = {\fontsize{12}' fitFunction '}'];
    title(title_str,'FontSize',16);
    xlabel('\tau (sec)','fontsize',16);
    ylabel('C^{(2)}(\tau)','fontsize',16);
    grid on;
%     whitebg;
    grid on;
    set(gca,'yscale','log');
    set(gca,'xscale','log');
    axis tight;
    axis square
    
    if labelPlotWithFit_mode
        if t1 < 1e-3
            text(t1,fitresult(t1),['{\bf\leftarrow}' num2str(t1 *1e6,'%2.1f') ' \musec, Amp=' num2str(A1)]);
        else
            text(t1,fitresult(t1),['{\bf\leftarrow}' num2str(t1*1e3,'%2.1f') ' msec, Amp=' num2str(A1)]);
        end
        
        if t2 < 1e-3
            text(t2,fitresult(t2),['{\bf\leftarrow}' num2str(t2*1e6,'%2.1f') ' \musec, Amp=' num2str(A2)]);
        else
            text(t2,fitresult(t2),['{\bf\leftarrow}' num2str(t2*1e3,'%2.1f') ' msec, Amp=' num2str(A2)]);
        end
        
        if t3 < 1e-3
            text(t3,fitresult(t3),['{\bf\leftarrow}' num2str(t3*1e6,'%2.1f') ' \musec, Amp=' num2str(A3)]);
        else
            text(t3,fitresult(t3),['{\bf\leftarrow}' num2str(t3*1e3,'%2.1f') ' msec, Amp=' num2str(A3)]);
        end
        
        if t4 < 1e-3
            text(t4,fitresult(t4),['{\bf\leftarrow}' num2str(t4*1e6,'%2.1f') ' \musec, Amp=' num2str(A4)]);
        else
            text(t4,fitresult(t4),['{\bf\leftarrow}' num2str(t4*1e3,'%2.1f') ' msec, Amp=' num2str(A4)]);
        end
        
        if t5 < 1e-3
            text(t5,fitresult(t5),['{\bf\leftarrow}' num2str(t5*1e6,'%2.1f') ' \musec, Amp=' num2str(A5)]);
        else
            text(t5,fitresult(t5),['{\bf\leftarrow}' num2str(t5*1e3,'%2.1f') ' msec, Amp=' num2str(A5)]);
        end
    end
    
    drawnow();
    %-------------------------------------------------------------------------
    % Save the Average TCF overlayed with the 5-exponential Fit
    %-------------------------------------------------------------------------
    if saveFigMode == 1
        foutName_prefix = [save_prefix foutName '_fiveExpYoffset_Fit'];
        fig = gcf;
        fig.PaperPositionMode = 'auto';
        extension = 'png';
        disp(['     Saving figure as ' foutName_prefix '.' extension ' in ' outputFolderName]);
        print(fig,[outputFolderName filesep() foutName_prefix],['-d' extension],'-r0');
        saveas(gcf,[outputFolderName filesep() foutName_prefix],['fig']);
        
        
    end
end
