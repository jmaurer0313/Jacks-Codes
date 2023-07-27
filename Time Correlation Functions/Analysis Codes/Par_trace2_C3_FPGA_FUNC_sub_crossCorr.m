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
function [totalRunTime,avgFileRunTime,tau1L,tau3L,NfilesToIterate,filesList] = Par_trace2_C3_FPGA_FUNC_sub_crossCorr(usres,constructName,tau2ArrayUnits_us,subsetOpt,subsetLoc,subsetName,quantityName,calcZeroTau,crossCorr) 
% disp(['        >> Entering ' functionName])
% 
% %--------------------------------------------------------------------------
% % Extra information
% %--------------------------------------------------------------------------
% extraDetail = '';
tic
%--------------------------------------------------------------------------
% Set resolution and Construct Name
%--------------------------------------------------------------------------
% constructName='(+15)Dimer'; 
% usres = 100;
msres = usres*10^-3;%('i.e. 1000 usec = 1 msec')
res = usres*10^-6;%(i.e. 1 million useconds turn into 1 second)
disp(['Setting the resolution equal to ' num2str(usres) ' microseconds ('...
    num2str(msres) ' msec) (' num2str(res) ' sec)']);

%--------------------------------------------------------------------------
% Look for all files ending with fileEndName
%--------------------------------------------------------------------------
if strcmp(quantityName,'FRET')
original_fileEndName = 'trace.mat';    
else
original_fileEndName = 'usec.mat';
end
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
    checkTaus=0;    
    verbose_mode = 0;
    testMode = 0;
    save_opt = 1;
    skipMode = 1;
    
    
%     load in the proper subset scanID list if the the option is selected
    if subsetOpt
        load(subsetLoc);
        
        if contains(subsetName,'Total') || contains(subsetName,'total')
            filesList=scanIDarray; 
            
        elseif contains(subsetName,'Uright') || contains(subsetName,'uright')
            filesList=Uright; 
       
        elseif contains(subsetName,'Boxset') || contains(subsetName,'boxset') || contains(subsetName,'boxSet') 
            filesList=boxSet; 
            
        end
        
    else
        filesList=[];
    end
            
   

% Set a value for tau2
% tau2ArrayUnits_us=0; 

 if tau2ArrayUnits_us < usres && tau2ArrayUnits_us~=0
        disp(['Tau2 value less than the resolution']);
        return;
  end 

% *****BLOCK TO RUN FOR CHECKS oN TAU ARRAY SPACING AND SAMPLING********
%%
    %--------------------------------------------------------------------------
    % Choose values for the tau arrays
    %--------------------------------------------------------------------------
% build a "total" Tau array which contains the times for all resolutions to
% be handled, such that the sticthing routine will alwsys match at the same
% tau values. 
% ******This construction is based on 10usec to 100usec to 1ms stiching plan. *****
    NptsTCFDesired = 325;
    NptsToAdd = NptsTCFDesired - 200;
    tauArray_msec = ([1:9, round(logspace(1, 3,NptsToAdd))]')*1000;
    tauArray_msec = tauArray_msec(tauArray_msec>=1000);
    tauArray_msec = (unique(tauArray_msec/1000))*1000;
    tauArrayUsecTotal = [linspace(10,100,10), linspace(100,1000,10),(tauArray_msec)']';
%     all based on starting at 10usec resolution for the stitching...
    tauArrayTotal = unique(floor(tauArrayUsecTotal/(10)));
    tauArrayUsecTotal=tauArrayTotal*10;
    tauArraySecTotal= tauArrayUsecTotal*1e-5;
    
   if usres==10
       tau1S=tauArrayUsecTotal(tauArrayUsecTotal<=100);
       tau3S=tauArrayUsecTotal(tauArrayUsecTotal<=100);
       tau1L=tauArrayUsecTotal;
       tau3L=tauArrayUsecTotal;

       
   elseif usres==100
       
       tau1S=tauArrayUsecTotal(tauArrayUsecTotal>=100);
       tau1S=tau1S(tau1S<=1000);
       tau3S=tauArrayUsecTotal(tauArrayUsecTotal>=100);
       tau3S=tau3S(tau3S<=1000);
       tau1L=tauArrayUsecTotal(tauArrayUsecTotal>=100);
       tau3L=tauArrayUsecTotal(tauArrayUsecTotal>=100);
       
     
       
       
   elseif usres==1000
       tau1S=[];
       tau3S=[];
       tau1L=tauArrayUsecTotal(tauArrayUsecTotal>=1000);
       tau3L=tauArrayUsecTotal(tauArrayUsecTotal>=1000);
       
   elseif usres==10000
       
    tauArray_msec = ([1:9, round(logspace(1, 3,NptsToAdd))]')*10000;
    tauArray_msec = tauArray_msec(tauArray_msec>=10000);
    tauArray_msec = (unique(tauArray_msec/10000))*10000;
    tauArrayUsecTotal = [linspace(10,100,10), linspace(100,1000,10),(tauArray_msec)']';
       tau1S=[];
       tau3S=[];
       tau1L=tauArrayUsecTotal(tauArrayUsecTotal>=10000);
       tau3L=tauArrayUsecTotal(tauArrayUsecTotal>=10000);
       
   end
   
   if calcZeroTau && ~isempty(tau1S) 
       tauArrayUsecTotal=[0;tauArrayUsecTotal];
       tau1S=[0;tau1S];
       tau1L=[0;tau1L];
       tau3S=[0;tau3S];
       tau3L=[0;tau3L];
   elseif calcZeroTau && isempty(tau1S) 
       tauArrayUsecTotal=[0;tauArrayUsecTotal];
       tau1L=[0;tau1L];      
       tau3L=[0;tau3L];
   end 

    


NptsTCF = numel(tau1L);
avgFileRunTime=0;
 



fprintf('     tau1: %d points [%d - %d] usec\r',numel(tau1L),tau1L(1),tau1L(end));
fprintf('     tau3: %d points [%d - %d] usec\r',numel(tau3L),tau3L(1),tau3L(end));


%%
%--------------------------------------------------------------------------
% Extra information
%--------------------------------------------------------------------------

if testMode                                                        
    %     disp('Testmode: Will only iterate over 1 molecule')
else
    NfilesToIterate = Nfiles;
end


%---------------------
% -----------------------------------------------------
%*********** Choose Values for tau2 ******************
%--------------------------------------------------------------------------
% tau2ArrayUnits_us = [0 1 10 100 1e3 1e4 1e5 1e6];
% tau2ArrayUnits_us = [1e6 1e5 1e4 1e3 1e2 1 0];

tau2ArrayUnits = tau2ArrayUnits_us./usres;


fprintf(['The dimensions of the 4pt will be tau1 = %d:%d usec (%d pts)' ...
    'and tau3 = %d:%d usec (%d pts)\r\n'],tau1L(1),tau1L(end),length(tau1L),tau3L(1),tau3L(end),length(tau3L));


%--------------------------------------------------------------------------
% Vary tau2
%--------------------------------------------------------------------------


    tau2_units = tau2ArrayUnits;
    tau2ValUsec = tau2_units*usres;
    
    %--------------------------------------------------------------------------
    % Create a folder to store the output
    %----------------------------------------------------------------------
    tau2_subFolderName = ['tau2-' num2str(tau2ValUsec,'%06d') 'us'];
    
    if subsetOpt && ~calcZeroTau
    folderOutName = ['3ptTCFs_SL_' subsetName '_' constructName '_' quantityName  filesep() tau2_subFolderName];
    elseif subsetOpt && calcZeroTau
    folderOutName = ['3ptTCFs_SL_' subsetName '_' constructName '_' quantityName '_calcZeroTau' filesep() tau2_subFolderName];
    elseif ~subsetOpt && ~calcZeroTau
    folderOutName = ['3ptTCFs_SL' '_' constructName '_' quantityName filesep() tau2_subFolderName];
    elseif ~subsetOpt && calcZeroTau
    folderOutName = ['3ptTCFs_SL' '_' constructName '_' quantityName '_calcZeroTau' filesep() tau2_subFolderName];
    end
    
    if save_opt
        if exist(folderOutName,'dir') ~= 7
            mkdir(folderOutName);
            disp(['Making a folder "' folderOutName '"  to hold the four point TCFs']);
        end
    end
%%
if checkTaus
tauCheckArray=zeros(numel(tau1L),numel(tau3L));

%     loops to check the population fo the tau-arrays when calculating the
%     "short x short" sqaure and the "short x long_UB" rectangle
for i=1:numel(tau1S)
    for j=1:numel(tau3S)
        tauCheckArray(i,j)=tauCheckArray(i,j)+1;
    end 
end 

for i=1:numel(tau1S)
    for j=(numel(tau3S)+1):numel(tau3L)
      tauCheckArray(i,j)=tauCheckArray(i,j)+1;
        
    end    
end 

for i=(numel(tau1S)+1):numel(tau1L)
    for j=1:numel(tau3S)
      tauCheckArray(i,j)=tauCheckArray(i,j)+1;
        
    end    
end 

if usres==1000
    for i=1:numel(tau1L)
    for j=1:numel(tau3L)
        tauCheckArray(i,j)=tauCheckArray(i,j)+1;
    end 
    end 
end

figure(1)
surf(tau1L,tau3L,tauCheckArray);
logx;
logy;

end

   
%%
    parfor i = 1:NfilesToIterate
        tic
        
%         check to see if the current trace is a member of the assigned
%         filesList
        %         Nfiles_soFar = Nfiles_soFar + 1;
        fileName = fileNames(i).name;
%         scanIDFormat = 'YYYYMMDD-###';
      if ~strcmp(quantityName,'FRET')
        scanID = fileName([1:15]);
      elseif strcmp(quantityName,'FRET')
        scanID = fileName([1:12]);  
      end
        
%         if certain workers end up with many continues and others the
%         actual files, there can be worker 'idleness' which reduces the
%         overall efficency of the parfor loop. If this becomes a
%         signficant proiblem, a list of the intended filenames can be
%         generated prior to the parfor loop and the iterations done over
%         those elements, loading in the file by name for each iteration. 
        if subsetOpt
          if ~ismember(scanID, filesList)
            continue;
          end 
        end

%         define the unit tau variables within the parfor loop 0 since
%         definitiion outside the loop leaves them as a 'broadcast'
%         variable whic reuires constant communciation between workers and
%         client, possibly leading to signifcant overhead in parallel. 
        if usres==1000 || usres==10000
        tau1LUnits = tau1L/usres;
        tau3LUnits = tau3L/usres;
        elseif usres==10 || usres==100
        tau1LUnits = tau1L/usres;
        tau3LUnits = tau3L/usres;
        tau1SUnits = tau1S/usres;
        tau3SUnits = tau3S/usres;    
        elseif strcmp(quantityName,'FRET')
        tau1LUnits = tauArrayUsecTotal/usres;   
        tau3LUnits = tauArrayUsecTotal/usres;
        end
        
    newFileEndName = 'fourptTCF.mat';
%         steps = numel(tau1ArrayUnits)*numel(tau3ArrayUnits);
%         step = 0;
        
            curFile=load(fileName);
            
            if  strcmp(quantityName,'P1PhaseFact')            
            f = curFile.P1PhaseFact;
            fAvg = nanmean(f);
            
            elseif strcmp(quantityName,'P1Amp')
            f = curFile.P1PhaseFact;
            f = abs(f);
            fAvg = nanmean(f);
            
            elseif strcmp(quantityName,'P1PhaseFactN')
            f = curFile.P1PhaseFactN;          
            fAvg = nanmean(f);
            
            elseif strcmp(quantityName,'P1AmpAvg')
            f = curFile.P1PhaseFactN;
            f = abs(f);
            fAvg = nanmean(f);   
            
            elseif strcmp(quantityName,'FRET')
            traceData = load(fileName,'counts1','counts2','res','scanID');
            don = traceData.counts1;
            acc = traceData.counts2;
            
            FRET = acc./(acc + don);
            f = FRET;
            fAvg=nanmean(f); 
            
            elseif strcmp(quantityName, 'AvgAmpCrossRate')
            f = curFile.P1PhaseFactN;
            f = abs(f);
            fAvg = nanmean(f); 
            g = curFile.counts;
            g = g/(usres*1e-6); 
            gAvg = nanmean(g); 
                
            end
        
        foutName = [scanID '_' num2str(usres,'%06d') 'us_tau2eq' num2str(tau2ValUsec,'%06d') 'us_fourptTCF.mat'];
        %         foutName = [fileName([1:end-length(original_fileEndName)]) newFileEndName]
        file_path = [folderOutName  filesep() foutName];
        
% %         perform a check here to see if the file is already in the
%           output folder, if skipMode is on, then skip it
        if  exist(file_path,'file') == 2 && skipMode == 1
            fprintf('%d/%d: Already detected %s in %s. Skipping.\n',i,Nfiles,foutName,folderOutName);
%             continue; 
%             load(file_path,'fourptTCF','NpairsTau_fourptTCF','tau1arrayUsec','tau2ValUsec','tau3arrayUsec','scanID','res');
%             tau1arraySec = tau1arrayUsec*1e-6;
%             tau3arraySec = tau3arrayUsec*1e-6;
            
        else
            fprintf('%d/%d: Now computing the 4pt TCF of %s. Saving as %s in %s.\n',i,Nfiles,scanID,foutName,folderOutName);
            
         
            
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
            fourptTCF = zeros(numel(tau1LUnits),numel(tau3LUnits));
            NpairsTau_fourptTCF = zeros(size(fourptTCF));
            
      if ~strcmp(quantityName,'FRET')
        if usres==10 || usres==100 
            
%             FIRST CALCULATE THE INNER SQUARE OF THE 4PT TCF (TAU1S X
%             TAU3S) THEN CALCULATE THE SIDES
            %--------------------------------------------------------------------------
            % Vary Tau 1
            %--------------------------------------------------------------------------
            for tau1_idx=1:numel(tau1S)
               
                tau1_units = tau1SUnits(tau1_idx);
                tau1ValUsec = tau1SUnits(tau1_idx)*res*1e6;
                if verbose_mode
                    disp(['          Tau1 = ' num2str(tau1ValUsec*1e-6) ' seconds']);
                end
                
                %--------------------------------------------------------------------------
                % Vary Tau 3
                %--------------------------------------------------------------------------
              for tau3_idx=1:numel(tau3S)
                    tau3_units = tau3SUnits(tau3_idx);
                    tau3ValUsec = tau3SUnits(tau3_idx)*res*1e6;
                   
    %                    
     %--------------------------------------------------------------------------
     % Step through the trajectories one element at a time.
     %--------------------------------------------------------------------------
     %                     
      for t = 1:NptsRaw - tau1_units - tau3_units
          %--------------------------------------------------
          % Speed up the calculation by checking for NaN's
          % before making a product
          %--------------------------------------------------
          
        if ~crossCorr
         if isnan(f(t)) 
            continue;
          elseif isnan(f(t+tau1_units))
            continue
          elseif isnan(f(t+tau1_units + tau3_units))
            continue
         end
         
          %------------------------------------------------------------------
          % If fourptTCF_temp ~= NaN, add it to the fourptTC
          %------------------------------------------------------------------
          %Not sure about the conj here - for 3-pt 01/17/2023
         fourptTCF_temp = (f(t)-fAvg)*((f(t+tau1_units)-fAvg))*((f(t + tau1_units + tau3_units)-fAvg));
         fourptTCF(tau1_idx,tau3_idx) = fourptTCF(tau1_idx,tau3_idx) + fourptTCF_temp;
         NpairsTau_fourptTCF(tau1_idx,tau3_idx) =  NpairsTau_fourptTCF(tau1_idx,tau3_idx) + 1;  
         
        elseif crossCorr
            
        if isnan(f(t))  
            continue;
          elseif isnan(g(t+tau1_units))
            continue          
          elseif isnan(g(t+tau1_units + tau3_units))
            continue
         end
         
          %------------------------------------------------------------------
          % If fourptTCF_temp ~= NaN, add it to the fourptTC
          %------------------------------------------------------------------
         fourptTCF_temp = (f(t)-fAvg)*((g(t+tau1_units)-gAvg))*((g(t + tau1_units + tau3_units)-gAvg));
         fourptTCF(tau1_idx,tau3_idx) = fourptTCF(tau1_idx,tau3_idx) + fourptTCF_temp;
         NpairsTau_fourptTCF(tau1_idx,tau3_idx) =  NpairsTau_fourptTCF(tau1_idx,tau3_idx) + 1;     
        
        end
        
                  
       end
                    
                    
                end
            end
            
            
%             calulate the FIRST "side"
         for tau1_idx=1:numel(tau1S)
               
                tau1_units = tau1SUnits(tau1_idx);
                tau1ValUsec = tau1SUnits(tau1_idx)*res*1e6;
                if verbose_mode
                    disp(['          Tau1 = ' num2str(tau1ValUsec*1e-6) ' seconds']);
                end
                
                %--------------------------------------------------------------------------
                % Vary Tau 3
                %--------------------------------------------------------------------------
              for tau3_idx=(numel(tau3S)+1):numel(tau3L)
                    tau3_units = tau3LUnits(tau3_idx);
                    tau3ValUsec = tau3LUnits(tau3_idx)*res*1e6;
                    
                    
     %--------------------------------------------------------------------------
     % Step through the trajectories one element at a time.
     %--------------------------------------------------------------------------
     %                     
      for t = 1:NptsRaw - tau1_units - tau3_units
          %--------------------------------------------------
          % Speed up the calculation by checking for NaN's
          % before making a product
          %--------------------------------------------------
          
        if ~crossCorr
         if isnan(f(t)) 
            continue;
          elseif isnan(f(t+tau1_units))
            continue          
          elseif isnan(f(t+tau1_units + tau3_units))
            continue
         end
         
          %------------------------------------------------------------------
          % If fourptTCF_temp ~= NaN, add it to the fourptTC
          %------------------------------------------------------------------
         fourptTCF_temp = (f(t)-fAvg)*((f(t+tau1_units)-fAvg))*((f(t + tau1_units + tau3_units)-fAvg));
         fourptTCF(tau1_idx,tau3_idx) = fourptTCF(tau1_idx,tau3_idx) + fourptTCF_temp;
         NpairsTau_fourptTCF(tau1_idx,tau3_idx) =  NpairsTau_fourptTCF(tau1_idx,tau3_idx) + 1;  
         
        elseif crossCorr
            
            if isnan(f(t)) 
            continue;
          elseif isnan(g(t+tau1_units))
            continue         
          elseif isnan(g(t+tau1_units + tau3_units))
            continue
         end
         
          %------------------------------------------------------------------
          % If fourptTCF_temp ~= NaN, add it to the fourptTC
          %------------------------------------------------------------------
         fourptTCF_temp = (f(t)-fAvg)*((g(t+tau1_units)-gAvg))*((g(t + tau1_units + tau3_units)-gAvg));
         fourptTCF(tau1_idx,tau3_idx) = fourptTCF(tau1_idx,tau3_idx) + fourptTCF_temp;
         NpairsTau_fourptTCF(tau1_idx,tau3_idx) =  NpairsTau_fourptTCF(tau1_idx,tau3_idx) + 1;       
        
        end
                    
                    
       end
     end
          
    %             calulate the SECOND "side"
         for tau1_idx=(numel(tau1S)+1):numel(tau1L)
             
                tau1_units = tau1LUnits(tau1_idx);
                tau1ValUsec = tau1LUnits(tau1_idx)*res*1e6;
                if verbose_mode
                    disp(['          Tau1 = ' num2str(tau1ValUsec*1e-6) ' seconds']);
                end
                
                %--------------------------------------------------------------------------
                % Vary Tau 3
                %--------------------------------------------------------------------------
              for tau3_idx=1:numel(tau3S)
                    tau3_units = tau3SUnits(tau3_idx);
                    tau3ValUsec = tau3SUnits(tau3_idx)*res*1e6;
                    
                    
     %--------------------------------------------------------------------------
     % Step through the trajectories one element at a time.
     %--------------------------------------------------------------------------
     %                     
      for t = 1:NptsRaw - tau1_units  - tau3_units
          %--------------------------------------------------
          % Speed up the calculation by checking for NaN's
          % before making a product
          %--------------------------------------------------
        if ~crossCorr
         if isnan(f(t)) 
            continue;
          elseif isnan(f(t+tau1_units))
            continue
          elseif isnan(f(t+tau1_units + tau3_units))
            continue
         end
         
          %------------------------------------------------------------------
          % If fourptTCF_temp ~= NaN, add it to the fourptTC
          %------------------------------------------------------------------
         fourptTCF_temp = (f(t)-fAvg)*((f(t+tau1_units)-fAvg))*((f(t + tau1_units + tau3_units)-fAvg));
         fourptTCF(tau1_idx,tau3_idx) = fourptTCF(tau1_idx,tau3_idx) + fourptTCF_temp;
         NpairsTau_fourptTCF(tau1_idx,tau3_idx) =  NpairsTau_fourptTCF(tau1_idx,tau3_idx) + 1;   
         
        elseif crossCorr
            
          if isnan(f(t)) 
            continue;
          elseif isnan(g(t+tau1_units))
            continue
          elseif isnan(g(t+tau1_units + tau3_units))
            continue
         end
         
          %------------------------------------------------------------------
          % If fourptTCF_temp ~= NaN, add it to the fourptTC
          %------------------------------------------------------------------
         fourptTCF_temp = (f(t)-fAvg)*((g(t+tau1_units)-gAvg))*((g(t + tau1_units + tau3_units)-gAvg));
         fourptTCF(tau1_idx,tau3_idx) = fourptTCF(tau1_idx,tau3_idx) + fourptTCF_temp;
         NpairsTau_fourptTCF(tau1_idx,tau3_idx) =  NpairsTau_fourptTCF(tau1_idx,tau3_idx) + 1;
         
        end
                  
         end
                    
                    
       end
     end   
            
         
        end
        end
      
    if ~strcmp(quantityName,'FRET')       
     if usres==1000 || usres==10000
          %--------------------------------------------------------------------------
          % Vary Tau 1
          %--------------------------------------------------------------------------
        for tau1_idx=1:numel(tau1L)
               
          tau1_units = tau1LUnits(tau1_idx);
          tau1ValUsec = tau1LUnits(tau1_idx)*res*1e6;
               
          if verbose_mode
                    disp(['          Tau1 = ' num2str(tau1ValUsec*1e-6) ' seconds']);
                end
                
                %--------------------------------------------------------------------------
                % Vary Tau 3
                %--------------------------------------------------------------------------
          for tau3_idx=1:numel(tau3L)
                    tau3_units = tau3LUnits(tau3_idx);
                    tau3ValUsec = tau3LUnits(tau3_idx)*res*1e6;
                   
    %                    
     %--------------------------------------------------------------------------
     % Step through the trajectories one element at a time.
     %--------------------------------------------------------------------------
     %                     
      for t = 1:NptsRaw - tau1_units - tau3_units
          %--------------------------------------------------
          % Speed up the calculation by checking for NaN's
          % before making a product
          %--------------------------------------------------
          
        if ~crossCorr
            
         if isnan(f(t)) 
            continue;
          elseif isnan(f(t+tau1_units))
            continue
          elseif isnan(f(t+tau1_units + tau3_units))
            continue
         end
         
          %------------------------------------------------------------------
          % If fourptTCF_temp ~= NaN, add it to the fourptTC
          %------------------------------------------------------------------
         fourptTCF_temp = (f(t)-fAvg)*((f(t+tau1_units)-fAvg))*((f(t + tau1_units + tau3_units)-fAvg));
         fourptTCF(tau1_idx,tau3_idx) = fourptTCF(tau1_idx,tau3_idx) + fourptTCF_temp;
         NpairsTau_fourptTCF(tau1_idx,tau3_idx) =  NpairsTau_fourptTCF(tau1_idx,tau3_idx) + 1;    
         
        elseif crossCorr
            
            if isnan(f(t)) 
            continue;
          elseif isnan(g(t+tau1_units))
            continue
          elseif isnan(g(t+tau1_units + tau3_units))
            continue
         end
         
          %------------------------------------------------------------------
          % If fourptTCF_temp ~= NaN, add it to the fourptTC
          %------------------------------------------------------------------
         fourptTCF_temp = (f(t)-fAvg)*((g(t+tau1_units)-gAvg))*((g(t + tau1_units + tau3_units)-gAvg));
         fourptTCF(tau1_idx,tau3_idx) = fourptTCF(tau1_idx,tau3_idx) + fourptTCF_temp;
         NpairsTau_fourptTCF(tau1_idx,tau3_idx) =  NpairsTau_fourptTCF(tau1_idx,tau3_idx) + 1;    
        
        end
                      
            
       end
                    
                   
       end
      end
               
     end
    end
    
    if strcmp(quantityName,'FRET')
          %--------------------------------------------------------------------------
          % Vary Tau 1
          %--------------------------------------------------------------------------
        for tau1_idx=1:numel(tau1LUnits)
               
          tau1_units = tau1LUnits(tau1_idx);
          tau1ValUsec = tau1LUnits(tau1_idx)*res*1e6;
               
          if verbose_mode
                    disp(['          Tau1 = ' num2str(tau1ValUsec*1e-6) ' seconds']);
                end
                
                %--------------------------------------------------------------------------
                % Vary Tau 3
                %--------------------------------------------------------------------------
          for tau3_idx=1:numel(tau3LUnits)
                    tau3_units = tau3LUnits(tau3_idx);
                    tau3ValUsec = tau3LUnits(tau3_idx)*res*1e6;
                   
    %                    
     %--------------------------------------------------------------------------
     % Step through the trajectories one element at a time.
     %--------------------------------------------------------------------------
     %                     
      for t = 1:NptsRaw - tau1_units - tau3_units
          %--------------------------------------------------
          % Speed up the calculation by checking for NaN's
          % before making a product
          %--------------------------------------------------
          
       if ~crossCorr
         if isnan(f(t)) 
            continue;
          elseif isnan(f(t+tau1_units))
            continue
          elseif isnan(f(t+tau1_units + tau3_units))
            continue
         end
         
          %------------------------------------------------------------------
          % If fourptTCF_temp ~= NaN, add it to the fourptTC
          %------------------------------------------------------------------
         fourptTCF_temp = (f(t)-fAvg)*((f(t+tau1_units)-fAvg))*((f(t + tau1_units + tau3_units)-fAvg));
         fourptTCF(tau1_idx,tau3_idx) = fourptTCF(tau1_idx,tau3_idx) + fourptTCF_temp;
         NpairsTau_fourptTCF(tau1_idx,tau3_idx) =  NpairsTau_fourptTCF(tau1_idx,tau3_idx) + 1; 
         
       elseif crossCorr
           
         if isnan(f(t)) 
            continue;
          elseif isnan(g(t+tau1_units))
            continue
          elseif isnan(g(t+tau1_units + tau3_units))
            continue
         end
         
          %------------------------------------------------------------------
          % If fourptTCF_temp ~= NaN, add it to the fourptTC
          %------------------------------------------------------------------
         fourptTCF_temp = (f(t)-fAvg)*((g(t+tau1_units)-gAvg))*((g(t + tau1_units + tau3_units)-gAvg));
         fourptTCF(tau1_idx,tau3_idx) = fourptTCF(tau1_idx,tau3_idx) + fourptTCF_temp;
         NpairsTau_fourptTCF(tau1_idx,tau3_idx) =  NpairsTau_fourptTCF(tau1_idx,tau3_idx) + 1;   
       
       end
                      
            
       end
                    
                   
                end
            end
               
     end
   
           
     
            %--------------------------------------------------------------------------
            % Normalize the Arrays
            %---------------------------------------
            fourptTCF = fourptTCF./NpairsTau_fourptTCF;
            fourptTCF(isnan(fourptTCF))=0; 
            
            %/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\
            % bf4ptTCF.m (END)
            %/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\
            if save_opt
%                 foutName = [scanID '_' num2str(usres,'%06d') 'us_tau2eq' num2str(tau2ValUsec,'%06d') 'us_fourptTCF.mat'];
               fname=[folderOutName filesep() foutName]; 
                parsaveC4_v2(fname,fourptTCF,NpairsTau_fourptTCF,tauArrayUsecTotal,tau1S,tau2ValUsec,tau3S,tau3L,tau1L,scanID,res);               
                disp(['     Saving variables in ' folderOutName ' as ' foutName]);
            end
            
            disp('Done Calculating the 4 pt TCF');
            
           end 
            
            toc
%             fileTime=toc;
%             avgFileRunTime=avgFileRunTime + fileTime;
        end
            
           
    end 
            totalRunTime= toc
            avgFileRunTime=totalRunTime/NfilesToIterate