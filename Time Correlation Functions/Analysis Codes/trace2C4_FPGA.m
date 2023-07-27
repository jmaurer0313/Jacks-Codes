% function trace2C4_FPGA()
%__________________________________________________________________________
% Author: Jack Maurer + Brett Israels 
%
% FUNCTION: trace2C4() Makes the 4pt TCF from raw data (no 2pt Subtraction)
%
% INPUT: 1) Be in folder with trace.mat files
%
% OUTPUT: 1) Folder named after the function with a 4-pt TCF for each trace
% file. 2) Figure files
%
% DEPENDENCIES:
%
% MODIFICATION LOG:
% BI/JM 20190620 Speeding up by checking for NaN's before product
% JM 4-14-2021: updated and overhauled code to handle FPGA data
%__________________________________________________________________________

% close all;
%--------------------------------------------------------------------------
% Program information
%--------------------------------------------------------------------------
functionName = 'trace2C4';
% disp(['        >> Entering ' functionName])
% 
% %--------------------------------------------------------------------------
% % Extra information
% %--------------------------------------------------------------------------
% extraDetail = '';

%--------------------------------------------------------------------------
% Set resolution and Construct Name
%--------------------------------------------------------------------------
constructName='(+15)Dimer'; 
usres = 100;
msres = usres*10^-3;%('i.e. 1000 usec = 1 msec')
res = usres*10^-6;%(i.e. 1 million useconds turn into 1 second)
disp(['Setting the resolution equal to ' num2str(usres) ' microseconds ('...
    num2str(msres) ' msec) (' num2str(res) ' sec)']);

%--------------------------------------------------------------------------
% Look for all files ending with fileEndName
%--------------------------------------------------------------------------
original_fileEndName = 'usec.mat';
fileNames = dir(['*' original_fileEndName]);
Nfiles = numel(fileNames);
if Nfiles == 0
    disp(['No files matching *' original_fileEndName ' found.']);
    return
else
    disp(['Found ' num2str(Nfiles) ' ' original_fileEndName ' files']);
end

%--------------------------------------------------------------------------
% User options
%--------------------------------------------------------------------------
plot_opt = 0;

save_opt = 1;
skipMode = 1;
figSaveMode = 0;
verbose_mode = 0;
calcAverage_mode = 0;
testMode = 0;

%--------------------------------------------------------------------------
% Determine size of 4point output
%--------------------------------------------------------------------------
% if usres < 1000
%     botRowMode = 1;
%     extraDetail = '_BotRow';
%     disp('     Calculating the bottom row of the 4point TCF');
% else
%     botRowMode = 0;
%     disp('     Calculating the entire 4point TCF');
% end

% *****BLOCK TO RUN FOR CHECKS oN TAU ARRAY SPACING AND SAMPLING********
%%
logarithmicTauSamplingMode = 1;
if logarithmicTauSamplingMode
    %--------------------------------------------------------------------------
    % Choose values for the tau arrays
    %--------------------------------------------------------------------------
%     tauArrayUsec_Combo = [];
    
%     these parameters here acheive a similar Npts in the final tauArray as
%     Bretts original spacing
    NptsTCFDesired = 325;
    NptsToAdd = NptsTCFDesired - 10;
    tauArrayUsec = [1:9, round(logspace(1, 6,NptsToAdd))]';
    tauArrayUsec = tauArrayUsec(tauArrayUsec>=res*1e6);
    tauArray = unique(floor(tauArrayUsec/(res*1e6)));

    % The tauArray will be in sec, running from 0 to 3 seconds.
    tauArray = [0; tauArray];
    tauArraySec= tauArray*res;
    
    %--------------------------------------------------------------------------
    % Paramaters across all resolutions
    %--------------------------------------------------------------------------
%     NptsTCFDesired = 500;%110;%500  will result in 207 points%50 as 28
%     maxtimeSec = 3;
%     maxtimeUsec = maxtimeSec*1e6;
%     NptsToAdd = NptsTCFDesired;
%     tauArrayUsec_all = [unique(round(logspace(log10(1), log10(maxtimeUsec),NptsToAdd)))];
    
    %--------------------------------------------------------------------------
    % Choose points to add when the resolution is 1000 usec
    %--------------------------------------------------------------------------
%     usresTemp = 1000;
%     tauArrayUsec = tauArrayUsec_all(tauArrayUsec_all >= usresTemp);
% %     convert the times in usec to number of indices at the resolution
%     tau1ArrayUnits = unique(round(tauArrayUsec/usresTemp));
%     tauArrayUsec = tau1ArrayUnits*usresTemp;
%     tauArrayUsec_Combo = unique([tauArrayUsec_Combo tauArrayUsec]);
%     maxTauUsec = tauArrayUsec(1);
    
    %--------------------------------------------------------------------------
    % Reduce the tauArrayUsec_all array to the altered tauArrayUsec_Combo
%     %--------------------------------------------------------------------------
%     last_idx = find(tauArrayUsec_all <= min(tauArrayUsec_Combo),1,'last');
%     tauArrayUsec_all = unique([tauArrayUsec_all(1:last_idx) tauArrayUsec_Combo]);
%     
%     %--------------------------------------------------------------------------
%     % Choose points to add when the resolution is 10 usec
%     %--------------------------------------------------------------------------
%     usresTemp = 10;
%     tauArrayUsec = tauArrayUsec_all(tauArrayUsec_all >= usresTemp);
%     tau1ArrayUnits = unique(round(tauArrayUsec/usresTemp));
%     tauArrayUsec = tau1ArrayUnits*usresTemp;
%     tauArrayUsec_Combo = unique([tauArrayUsec_Combo tauArrayUsec]);
%     tauArrayUsec = tauArrayUsec_Combo(tauArrayUsec_Combo <= maxTauUsec);
%     tauArrayUsec = tauArrayUsec(tauArrayUsec >= usresTemp);
%     maxTauUsec = tauArrayUsec(1);
    
    %--------------------------------------------------------------------------
    % Reduce the tauArrayUsec_all array to the altered tauArrayUsec_Combo
    %--------------------------------------------------------------------------
%     last_idx = find(tauArrayUsec_all <= min(tauArrayUsec_Combo),1,'last');
%     tauArrayUsec_all = unique([tauArrayUsec_all(1:last_idx) tauArrayUsec_Combo]);
    
    %--------------------------------------------------------------------------
    % Choose points to add when the resolution is 1 usec
    %--------------------------------------------------------------------------
    if usres == 1
        usresTemp = 1;
        tauArrayUsec_Combo = tauArrayUsec_all;
        
        tauArrayUsec = tauArrayUsec_Combo(tauArrayUsec_Combo <= maxTauUsec);
        tauArrayUsec = tauArrayUsec(tauArrayUsec >= usresTemp);
    end
else
    tauArrayUsec_Combo = [1000:1000:250000];
end
% %--------------------------------------------------------------------------
% %--------------------------------------------------------------------------
% %--------------------------------------------------------------------------
% %--------------------------------------------------------------------------
%
% %--------------------------------------------------------------------------
% % Limit the array to only points >= the resolution
% %--------------------------------------------------------------------------
% if usres > 1
%     tauArrayUsec_Combo = tauArrayUsec_Combo(tauArrayUsec_Combo >= usres);
% elseif usres == 1
%     tauArrayUsec_Combo = tauArrayUsec_Combo(tauArrayUsec_Combo >= usres);
% end

% tauArraySec_Combo = tauArrayUsec_Combo*1e-6;

NptsTCF = numel(tauArray);
 tau1ArrayUnits = tauArray;
 tau3ArrayUnits = tauArray;

% if botRowMode == 0
%     tau1ArrayUnits = tauArrayUsec_Combo/usres;
%     tau3ArrayUnits = tauArrayUsec_Combo/usres;
% elseif  botRowMode
%     tau1ArrayUnits = tauArrayUsec/usres;
%     tau3ArrayUnits = tauArrayUsec_Combo/usres;
% end

fprintf('     tau1: %d points [%d - %d] usec\r',numel(tau1ArrayUnits),tau1ArrayUnits(1)*usres,tau1ArrayUnits(end)*usres);
fprintf('     tau3: %d points [%d - %d] usec\r',numel(tau3ArrayUnits),tau3ArrayUnits(1)*usres,tau3ArrayUnits(end)*usres);


%%
%--------------------------------------------------------------------------
% Extra information
%--------------------------------------------------------------------------
% extraDetail = [extraDetail '_Npts' num2str(NptsTCF)];
if testMode
    NfilesToIterate = 1;
    %     disp('Testmode: Will only iterate over 1 molecule')
else
    NfilesToIterate = Nfiles;
end

%--------------------------------------------------------------------------
% Make the sample_description based on the information found
%--------------------------------------------------------------------------
% [sample_description,save_prefix,sample_description_plain,conditions,DNA_save_prefix] = sample_descriptionGetter()
% [sample_description,save_prefix] = sample_descriptionGetter();


%--------------------------------------------------------------------------
% Create a log file to hold information about the 2ptTCFs.
%--------------------------------------------------------------------------
%open a file for appending to
% log_exist = exist([folderOutName filesep() functionName,'.log']);
% fid = fopen([folderOutName filesep() functionName,'.log'],'a');
% if log_exist == 0
%     fprintf(fid,'ScanID      Npoints   C(1usec)\r\n');
% end

%--------------------------------------------------------------------------
% For each file, read the donor and acceptor, compute the FRET, then use
% the bf2ptTCF_v3 MATLAB function to calculate the 2pt-TCF.
%--------------------------------------------------------------------------
% if plot_opt
%     fig = figure;
%     fig.Color = 'w';
%     fig.Position = [185 75 1095 921];
%     %     fig.Position = [28.3333 105 1.4767e+03 1250];%Home computer
%     fig.Name = [functionName extraDetail];
% end

%---------------------
% -----------------------------------------------------
%*********** Choose Values for tau2 ******************
%--------------------------------------------------------------------------
% tau2ArrayUnits_us = [0 1 10 100 1e3 1e4 1e5 1e6];
% tau2ArrayUnits_us = [1e6 1e5 1e4 1e3 1e2 1 0];
tau2ArrayUnits_us = [0];

% tau2ArrayUnits_sec = [0 1 5 10 15 20 30 50 75 100]*1e-3;
% tau2ArrayUnits_sec = [50 75 100]*1e-3;
% tau2ArrayUnits_us = tau2ArrayUnits_sec*1e6;
% tau2ArrayUnits_us = [0 1 10 100 1000 3e3 5e3 10e3 50e3 10e4 1e6 3e6 5e6 10e6];
inlcude_tau2eq0_mode = 1;
if inlcude_tau2eq0_mode
    tau2ArrayUnits_us = [0 tau2ArrayUnits_us(tau2ArrayUnits_us >= usres)];
else
    tau2ArrayUnits_us = [tau2ArrayUnits_us(tau2ArrayUnits_us >= usres)];
end
tau2ArrayUnits = tau2ArrayUnits_us./usres;


%--------------------------------------------------------------------------
% Set up the time arrays: tau1 tau2 and tau3
%--------------------------------------------------------------------------
tau1arraySec = tau1ArrayUnits*res;
tau2arraySec = tau2ArrayUnits*res;
tau3arraySec = tau3ArrayUnits*res;

tau1arrayUsec = tau1ArrayUnits*usres;
tau2arrayUsec = tau2ArrayUnits*usres;
tau3arrayUsec = tau3ArrayUnits*usres;

fprintf(['The dimensions of the 4pt will be tau1 = %d:%d usec (%d pts)' ...
    'and tau3 = %d:%d usec (%d pts)\r\n'],tau1arrayUsec(1),tau1arrayUsec(end),length(tau1arrayUsec),tau3arrayUsec(1),tau3arrayUsec(end),length(tau3arrayUsec));


%--------------------------------------------------------------------------
% Vary tau2
%--------------------------------------------------------------------------

for tau2_idx = 1:numel(tau2ArrayUnits)
    tau2_units = tau2ArrayUnits(tau2_idx);
    tau2ValUsec = tau2ArrayUnits(tau2_idx)*usres;
    if verbose_mode
        disp(['Tau2 = ' num2str(tau2ArrayUnits(tau2_idx)*res) ' seconds']);
    end
    
    %--------------------------------------------------------------------------
    % Create a folder to store the output
    %----------------------------------------------------------------------
    tau2_subFolderName = ['tau2-' num2str(tau2ValUsec,'%06d') 'us'];
    folderOutName = ['4ptTCFs_' constructName  filesep() tau2_subFolderName];
    
    if save_opt
        if exist(folderOutName,'dir') ~= 7
            mkdir(folderOutName);
            disp(['Making a folder "' folderOutName '"  to hold the four point TCFs']);
        end
    end
    newFileEndName = 'fourptTCF.mat';
    Nfiles_soFar = 0;
    %         for i = NfilesToIterate:-1:1
    
    for i = 1:NfilesToIterate
        tic
        steps = numel(tau1ArrayUnits)*numel(tau3ArrayUnits);
        step = 0;
        
        Nfiles_soFar = Nfiles_soFar + 1;
        fileName = fileNames(i).name;
%         scanIDFormat = 'YYYYMMDD-###';
        scanID = fileName([1:15]);
        
        foutName = [scanID '_' num2str(usres,'%06d') 'us_tau2eq' num2str(tau2ValUsec,'%06d') 'us_fourptTCF.mat'];
        %         foutName = [fileName([1:end-length(original_fileEndName)]) newFileEndName]
        file_path = [folderOutName  filesep() foutName];
        
% %         perform a check here to see if the file is already in the
%           output folder, if skipMode is on, then skip it
        if  exist(file_path,'file') == 2 && skipMode == 1
            fprintf('%d/%d: Already detected %s in %s. Skipping.\n',i,Nfiles,foutName,folderOutName);
            
            load(file_path,'fourptTCF','NpairsTau_fourptTCF','tau1arrayUsec','tau2ValUsec','tau3arrayUsec','scanID','res');
            tau1arraySec = tau1arrayUsec*1e-6;
            tau3arraySec = tau3arrayUsec*1e-6;
            
        else
            fprintf('%d/%d: Now computing the 4pt TCF of %s. Saving as %s in %s.\n',i,Nfiles,scanID,foutName,folderOutName);
            
            
           
            load(fileName);
            f = P1PhaseFact;
            fAvg = nanmean(f);
           
            
            %--------------------------------------------------------------------------
            % Calculate the 4-point TCF
            %--------------------------------------------------------------------------
            
            
            
            %/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\
            % bf4ptTCF.m (BEGIN)
            %/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\
           
            
            NptsRaw = numel(f);
            MaxTrajLengthSec = NptsRaw*res;
            disp(['     Trajectory is ' num2str(MaxTrajLengthSec) ' seconds long']);
            
            %--------------------------------------------------------------------------
            % Prepare the TCF arrays.
            %--------------------------------------------------------------------------
            fourptTCF = zeros(numel(tau1ArrayUnits),numel(tau3ArrayUnits));
            NpairsTau_fourptTCF = zeros(size(fourptTCF));
            
            %--------------------------------------------------------------------------
            % Vary Tau 1
            %--------------------------------------------------------------------------
            for  tau1_idx = 1:numel(tau1ArrayUnits)
                tau1_units = tau1ArrayUnits(tau1_idx);
                tau1ValUsec = tau1ArrayUnits(tau1_idx)*res*1e6;
                if verbose_mode
                    disp(['          Tau1 = ' num2str(tau1ValUsec*1e-6) ' seconds']);
                end
                
                %--------------------------------------------------------------------------
                % Vary Tau 3
                %--------------------------------------------------------------------------
                for  tau3_idx = 1:numel(tau3ArrayUnits)
                    tau3_units = tau3ArrayUnits(tau3_idx);
                    step = step + 1;
                    tau3ValUsec = tau3ArrayUnits(tau3_idx)*res*1e6;
                    %                    
                    %--------------------------------------------------------------------------
                    % Step through the trajectories one element at a time.
                    %--------------------------------------------------------------------------
                    %                     
                    for t = 1:NptsRaw - tau1_units - tau2_units - tau3_units
                        %--------------------------------------------------
                        % Speed up the calculation by checking for NaN's
                        % before making a product
                        %--------------------------------------------------
         if isnan(f(t)) 
            continue;
          elseif isnan(f(t+tau1_units))
            continue
          elseif isnan(f(t+tau1_units + tau2_units))
            continue
          elseif isnan(f(t+tau1_units+ + tau2_units + tau3_units))
            continue
         end
         
          %------------------------------------------------------------------
          % If fourptTCF_temp ~= NaN, add it to the fourptTC
          %------------------------------------------------------------------
         fourptTCF_temp = (f(t)-fAvg)*(conj(f(t+tau1_units)-fAvg))*(f(t + tau1_units + tau2_units)-fAvg)*(conj(f(t + tau1_units + tau2_units + tau3_units)-fAvg));
         fourptTCF(tau1_idx,tau3_idx) = fourptTCF(tau1_idx,tau3_idx) + fourptTCF_temp;
         NpairsTau_fourptTCF(tau1_idx,tau3_idx) =  NpairsTau_fourptTCF(tau1_idx,tau3_idx) + 1;               
                       
%                         fourptTCF_temp = f(t)*f(t+tau1_units)*f(t + tau1_units + tau2_units)*f(t + tau1_units + tau2_units + tau3_units);
%                         if isnan(fourptTCF_temp) ~= 1
%                             % Subtract off the product of the 2ptTCF
%                             %                             fourptTCF_temp = fourptTCF_temp - f(t)*f(t+tau1_units)*f(t)*f(t + tau3_units);
%                             %                             fprintf('%d*%d*%d*%d =  %d\r\n',f(t),f(t+tau1_units),f(t + tau1_units + tau2_units),f(t + tau1_units + tau2_units + tau3_units),fourptTCF_temp);
%                             fourptTCF(tau1_idx,tau3_idx) = fourptTCF(tau1_idx,tau3_idx) + fourptTCF_temp;
%                             NpairsTau_fourptTCF(tau1_idx,tau3_idx) =  NpairsTau_fourptTCF(tau1_idx,tau3_idx) + 1;
%                         end
                        
                        
                        
                    end
                    
                    if verbose_mode
                        %                         disp(['          Tau3 = ' num2str(tau3ValUsec*1e-6) ' seconds']);
                        fprintf('          4pt(%d,%d,%d) = %d (Npairs = %d) [step%d/%d]\r',...
                            tau1_units,tau2_units,tau3_units,fourptTCF(tau1_idx,tau3_idx),NpairsTau_fourptTCF(tau1_idx,tau3_idx),step,steps);
                    end
                end
            end
            
            %--------------------------------------------------------------------------
            % Normalize the Arrays
            %---------------------------------------
            fourptTCF = fourptTCF./NpairsTau_fourptTCF;
            
            %/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\
            % bf4ptTCF.m (END)
            %/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\
            if save_opt
%                 foutName = [scanID '_' num2str(usres,'%06d') 'us_tau2eq' num2str(tau2ValUsec,'%06d') 'us_fourptTCF.mat'];
                save([folderOutName filesep() foutName],'fourptTCF','NpairsTau_fourptTCF','tau1arrayUsec','tau2ValUsec','tau3arrayUsec','scanID','res');               
                disp(['     Saving variables in ' folderOutName ' as ' foutName]);
            end
            disp('Done Calculating the 4 pt TCF');
        end
        
        if plot_opt
            figure(2);
            whitebg;           
           surf(tau1arraySec,tau3arraySec,real(fourptTCF));
%             else
%                 surf(tau1arraySec,tau3arraySec,fourptTCF');
%                 hold on;
%                 surf(tau3arraySec',tau1arraySec',fourptTCF);
%                 hold off;
%             end
            
            %             contour3(tau1arraySec,tau3arraySec,fourptTCF,100,'LineWidth',1)
            %             contourf(tau1arraySec,tau3arraySec,fourptTCF,100,'LineWidth',1,'LineStyle','none')
%             [sample_description] = sample_descriptionGetter();
            
            title_str = [constructName '  "' scanID '"' 10 ' (\tau_2 = ' num2str(tau2ValUsec*1e-6) 'sec)    {\color{red}Res = ' num2str(usres) '\musec}'];
            
            title(title_str,'fontsize',20);
            xlabel('\tau_1 (sec)','fontsize',18);
            ylabel('\tau_3 (sec)','fontsize',18);
            zlabel('C^{(4)}(\tau)','fontsize',18);
            grid on;
            axis tight;
            axis square;
            view(0,90)
            colormap jet
            colorbar;
            %             caxis([min(min(fourptTCF)),fourptTCF(2,2)]);
            set(gca,'yscale','log');
            set(gca,'xscale','log');
            drawnow;
            
            %         legend('off')
            %         lgd = legend('show');
            %
            %         lgd.FontSize = 12;
            %         lgd.Location = 'eastoutside';
            
            
%             if figSaveMode
% %                 foutName = [scanID '_' num2str(usres,'%06d') 'us_tau2eq' num2str(tau2ValUsec,'%06d') 'us_fourptTCF'];
%                 save([folderOutName filesep() foutName],'fourptTCF',...
%                     'NpairsTau_fourptTCF','tau1arrayUsec','tau2ValUsec','tau3arrayUsec','scanID','res');
%                 disp(['Saving figure in ' folderOutName ' as ' foutName]);
%                 
%                 extension = 'png';
%                 pdfSaveModeIN= 0;
%                 figSaveModeIN = 0;
%                 figSaver(foutName,folderOutName,extension,pdfSaveModeIN,figSaveModeIN)
%             end
        end
        toc
        dt = datestr(now,'mmmm dd, yyyy HH:MM:SS.FFF AM');
        
    end
    
    
    %--------------------------------------------------------------------------
    % Calculate the average
    %--------------------------------------------------------------------------
    if calcAverage_mode
        avg_filePath = [folderOutName filesep() 'fourPttcfAverager_output' ...
            filesep() num2str(usres,'%06d') 'us_tau2eq' num2str(tau2ValUsec,'%06d') 'us_fourptTCFavg.mat'];
        if exist(avg_filePath,'file') == 2
            disp('Already detected the average');
        else
            disp('Calculating the average correlation function');
            wd = pwd;
            cd(folderOutName);
            fourPttcfAverager_v2();
            cd(wd);
        end
    end
end