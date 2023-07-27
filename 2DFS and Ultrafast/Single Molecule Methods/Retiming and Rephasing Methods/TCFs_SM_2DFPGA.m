% GOAL: 
% 1. bring in a 2D SM data set. Retime the t21 and t43 axes, such that
% the corrected time axes refer to the true inter/intra pulse delays. 
% 2. With the corrected time axes, go back into the 2D SM raw data, form the
% phaseFactor trajectry for each of the delays with T21,T43<30fs at some resolution (due to
% the 14-15fs pulse duration, notign beyond 30fs should  not have significant
% overlap at all). 
% 3. Perform a TCF analysis of the raw data across all delays within the
% 30fs window (along both axes). 
% 4. save the processed phaseFactor (at the set resolution) to a file along
% with the relevant TCFs of that trajectory (phase, amp, phaseFactor, counts), 
% labeled with the interpulse delays. The delays will have to fall within an interval due to the
% individual timing corrections almost never being the same. 
% 5. let the interpulse delay set the output folder, within each folder
% there shoukd be sub directory for "Trajectories" and "TCFs" with each TCF
% saved to its own folder for the use of wAVG plottign schemes. 

%%
% load in the data, get the timing corrections, build a "trueTimes" matrix
% to match the size of the underlying raw data (photon streams) 
fileFolder = 'C:\Users\Jack\Documents\Marcus Lab\2D_SingleMolecule\Sm_avgFolders\TimeLapses\20191217_mol1_1hr\20191217-164312-2DFPGA';
outFolderGen= 'C:\Users\Jack\Documents\Marcus Lab\2D_SingleMolecule\Sm_avgFolders\TCFs';
constructName='E3cE4'; 
% scanID for the current 2d data set under operation. May need to change
% depending on the structure of the directory where data is stored 
scanID=fileFolder(96:end); 


Tdes=0.6;
intTime=42;
plotOpt=0;
printOpt=0;
scanNum=1;

% time in fs beyond which no TCFs will be calculated for the interpulse
% pair delay 
cutoffTime=30; 
% resolution in usec for the trajCalc to run
res=500;

% generate a TauArray for the TCF calculations
NptsTCFDesired = 450;
NptsToAdd = NptsTCFDesired - 10;
tauArrayUsec = [1:9, round(logspace(1, 7,NptsToAdd))]';
tauArrayUsec = tauArrayUsec(tauArrayUsec>=res);   
 
% This resulting Tau array will be equal at each index to the number of
% bins which seperate the two points in the TCF. 
% IN OTHER WORDS: THIS IS THE BINS TAU ARRAY, DENOTING THE SPACING OF BINS
% AT THE CURRENT RESOLTUTION TCF
tauArray = unique(floor(tauArrayUsec/(res)));
% The tauArray will be in sec, running from 0 to 3 seconds.
tauArray = [0; tauArray];

% need 11 spaces for 10 intervals
outDelays=linspace(0,30,11); 

% generate output folders for the TCFs to be saved to. Set ranges for the
% "trueTime" delay, within which the TCFs resulting will be saved and
% processed together. 
outFolder1= ([outFolderGen filesep() constructName '_' num2str(res) 'usec_2Pt_TCFs'  filesep() 'Linear TCFs' filesep() 'Tau=0-3fs']);
outFolder2= ([outFolderGen filesep() constructName '_' num2str(res) 'usec_2Pt_TCFs'  filesep() 'Linear TCFs' filesep() 'Tau=3-6fs']);
outFolder3= ([outFolderGen filesep() constructName '_' num2str(res) 'usec_2Pt_TCFs'  filesep() 'Linear TCFs' filesep() 'Tau=6-9fs']);
outFolder4= ([outFolderGen filesep() constructName '_' num2str(res) 'usec_2Pt_TCFs'  filesep() 'Linear TCFs' filesep() 'Tau=9-12fs']);
outFolder5= ([outFolderGen filesep() constructName '_' num2str(res) 'usec_2Pt_TCFs'  filesep() 'Linear TCFs' filesep() 'Tau=12-15fs']);
outFolder6= ([outFolderGen filesep() constructName '_' num2str(res) 'usec_2Pt_TCFs'  filesep() 'Linear TCFs' filesep() 'Tau=15-18fs']);
outFolder7= ([outFolderGen filesep() constructName '_' num2str(res) 'usec_2Pt_TCFs'  filesep() 'Linear TCFs' filesep() 'Tau=18-21fs']);
outFolder8= ([outFolderGen filesep() constructName '_' num2str(res) 'usec_2Pt_TCFs'  filesep() 'Linear TCFs' filesep() 'Tau=21-24fs']);
outFolder9= ([outFolderGen filesep() constructName '_' num2str(res) 'usec_2Pt_TCFs'  filesep() 'Linear TCFs' filesep() 'Tau=24-27fs']);
outFolder10= ([outFolderGen filesep() constructName '_' num2str(res) 'usec_2Pt_TCFs'  filesep() 'Linear TCFs' filesep() 'Tau=27-30fs']);

% create a structure to hold the output folder names for use within the
% loop
foldTable= struct('f1',outFolder1,'f2',outFolder2,'f3',outFolder3,'f4',outFolder4,'f5',outFolder5,'f6',outFolder6,'f7',outFolder7,'f8',outFolder8,'f9',outFolder9,'f10',outFolder10); 

% switch the t21 axis time to be "positive" (negative from stage moving
% opposit that of t43 time axis). put axes into fs from ps
[XmatTemp, YmatTemp, tb1, tb2,photGrid] = linear2DFPGA_RT_Lapse(fileFolder, intTime);
    
[~, ~, Xshifts, Yshifts, Xangles, Yangles, ~, ~, ~, ~, NIFFT,Xzero] = interpLinear_Lapse(plotOpt,XmatTemp, YmatTemp, tb1, tb2, Tdes,scanNum,photGrid,printOpt);

tb1=tb1*(-1)*1e3;
tb2=tb2*1e3;

trueTimesX= repmat(tb1,length(tb2),1) - Xshifts(:,1); 
trueTimesY= repmat(tb2,length(tb1),1) - Yshifts(:,1);
trueTimesY=trueTimesY.'; 
 
%%

% loop over the data folder, load in each integration period 1 by 1 and
% perform a TCF of the time resolved data. save to an output folder based
% on the trueTime (must be below 30fs) 

xl = length(tb1);
yl = length(tb2);
photGrid= zeros(xl, yl);        


        for yi = 1:yl
            for xi = 1:xl
%                 evaluate the times from the trueTimes, only proceed with
%                 the calculation if the times are within the allowed range
%                 (the truetimes matrices are both 9x10 so t43 x t21 in
%                 dimension) 
                disp(['Currently on step x=' num2str(xi) ' y=' num2str(yi) ' grid size=' num2str(xl) 'x' num2str(yl)]); 
                
%                 NOTE: 11/16/2022 - think about if this conditonal is true
%                 or not needed to maintain the grid
                if trueTimesX(yi,xi)>=0 && trueTimesX(yi,xi)<=cutoffTime 
%                     generate a label for the tru time zero
                  timeZero=trueTimesX(yi,xi);
                  if timeZero>=10
                  label=num2str(timeZero);
                  label=([label(1:2) '_' label(4:end)]);
                  else
                  label=num2str(timeZero);
                  label=([label(1) '_' label(3:end)]);
                  end
                    
                timeID = fullfile(fileFolder, [num2str(xi-1,'%02.0f') num2str(yi-1, '%02.0f') 'time.bin']);
                p1file = fullfile(fileFolder, [num2str(xi-1,'%02.0f') num2str(yi-1, '%02.0f') 'p1.bin']);
                timeFile = fopen(timeID);
                p1ID = fopen(p1file);
                time = fread(timeFile,Inf,'uint64=>uint64',0,'s');
%                 times are converted to seconds here 
                times = ((double (time-time(1))./(8e7))+1);
                times=times-times(1);
                timesUsec=times*1e6; 
                p1 = fread(p1ID,Inf,'float64=>double',0,'s');
                p1r = (2.*pi).*p1;
                [intWind]=find(times<=intTime);
%                At this point isolate the phase and time lists from zero
%                up to the integration time set
                p1=p1(1:intWind(end));
                photGrid(xi,yi)=numel(p1); 
                
%                 Here: insert a function call to both a) calculate the
%                 trace at a resolution, and then b) calculate the TCF of
%                 the proccesed trace (trajectory) 

%                 generate the proper outout folder from the general table
%                 above. Then save all 4 TCFs to their respective folders 
                curList=find(outDelays<=trueTimesX(yi,xi));
                curOutFold= getfield(foldTable,['f' num2str(curList(end))]);
                
                [trajAvgP1, trajRawP1, trajLength ] = trajCalcFPGA(p1r, timesUsec, res);
                
%                 phase TCF
                [ tcfArr, tcfNpairs ] = TCFcalc(angle(trajAvgP1), tauArray); 
                curOutPhase=[curOutFold  filesep() 'P1' filesep() 'Phase'];                
                if exist(curOutPhase, 'dir')~= 7
                    mkdir(curOutPhase);
                end    
                varOfArray=nanvar(angle(trajAvgP1));
                save([curOutPhase filesep() scanID '_' label 'fs'], 'tcfArr','tcfNpairs','trajLength','trajAvgP1','trajRawP1','tauArray','tauArrayUsec','varOfArray');

%                 amp TCF
                [ tcfArr, tcfNpairs ] = TCFcalc( abs(trajAvgP1), tauArray ); 
                 curOutAmp=[curOutFold  filesep() 'P1' filesep() 'Amp'];
                if exist(curOutAmp, 'dir')~= 7
                    mkdir(curOutAmp);
                end
                varOfArray=nanvar(abs(trajAvgP1));
                save([curOutAmp filesep() scanID '_' label 'fs'], 'tcfArr', 'tcfNpairs','trajLength','trajAvgP1','trajRawP1','tauArray','tauArrayUsec','varOfArray');  

                 
%                 phase Factor avg TCF
                [ tcfArr, tcfNpairs] = TCFcalc((trajAvgP1), tauArray );
                 curOutPhsFctAvg=[curOutFold  filesep() 'P1' filesep() 'PhaseFactAvg'];
                 if exist(curOutPhsFctAvg, 'dir')~= 7
                    mkdir(curOutPhsFctAvg);
                 end
                 varOfArray=nanvar((trajAvgP1));
                 save([curOutPhsFctAvg filesep() scanID '_' label 'fs'], 'tcfArr', 'tcfNpairs','trajLength','trajAvgP1','trajRawP1','tauArray','tauArrayUsec','varOfArray');  

%                  phase Factor raw TCF
                [ tcfArr, tcfNpairs ] = TCFcalc((trajRawP1), tauArray );                
                curOutPhsFctRaw=[curOutFold filesep() 'P1' filesep() 'PhaseFactRaw' ];  
                 if exist(curOutPhsFctRaw, 'dir')~= 7
                    mkdir(curOutPhsFctRaw);
                 end 
                varOfArray=nanvar((trajRawP1));
                save([curOutPhsFctRaw filesep() scanID '_' label 'fs'], 'tcfArr', 'tcfNpairs','trajLength','trajAvgP1','trajRawP1','tauArray','tauArrayUsec','varOfArray');  

        
                end 
                
                
                if trueTimesY(yi,xi)>=0 && trueTimesY(yi,xi)<=cutoffTime                    
%                   generate a label for the true time zero
                  timeZero=trueTimesY(yi,xi);
                  if timeZero>=10
                  label=num2str(timeZero);
                  label=([label(1:2) '_' label(4:end)]);
                  else
                  label=num2str(timeZero);
                  label=([label(1) '_' label(3:end)]);
                  end
                  
                  timeID = fullfile(fileFolder, [num2str(xi-1,'%02.0f') num2str(yi-1, '%02.0f') 'time.bin']);
                  p2file = fullfile(fileFolder, [num2str(xi-1,'%02.0f') num2str(yi-1, '%02.0f') 'p2.bin']);
                  timeFile = fopen(timeID);
                p2ID = fopen(p2file);
                time = fread(timeFile,Inf,'uint64=>uint64',0,'s');
                times = ((double (time-time(1))./(8e7))+1);
                times=times-times(1);
                timesUsec=times*1e6; 
                p2 = fread(p2ID,Inf,'float64=>double',0,'s');
                p2r = (2.*pi).*p2;
                [intWind]=find(times<=intTime);
%                At this point isolate the phase and time lists from zero
%                up to the integration time set
                p2=p2(1:intWind(end)); 
                photGrid(xi,yi)=numel(p2); 
                
%                  Here: insert a function call to both a) calculate the
%                 trace at a resolution, and then b) calculate the TCF of
%                 the proccesed trace (trajectory) 

%                 generate the proper outout folder from the general table
%                 above. Then save all 4 TCFs to their respective folders 
                curList=find(outDelays<=trueTimesY(yi,xi));
                curOutFold= getfield(foldTable,['f' num2str(curList(end))]);
                
                [trajAvgP2, trajRawP2, trajLength ] = trajCalcFPGA(p2r, timesUsec, res);
                
%                 phase TCF
                [ tcfArr, tcfNpairs ] = TCFcalc( angle(trajAvgP2), tauArray ); 
                curOutPhase=[curOutFold  filesep() 'P2' filesep() 'Phase'];                
                if exist(curOutPhase, 'dir')~= 7
                    mkdir(curOutPhase);
                end
                varOfArray=nanvar(angle(trajAvgP2));
                save([curOutPhase filesep() scanID '_' label 'fs'], 'tcfArr','tcfNpairs','trajLength','trajAvgP2','trajRawP2','tauArray','tauArrayUsec','varOfArray');

%                 amp TCF
                [ tcfArr, tcfNpairs ] = TCFcalc( abs(trajAvgP2), tauArray ); 
                 curOutAmp=[curOutFold  filesep() 'P2' filesep() 'Amp'];
                if exist(curOutAmp, 'dir')~= 7
                    mkdir(curOutAmp);
                end
                varOfArray=nanvar(abs(trajAvgP2));
                save([curOutAmp filesep() scanID '_' label 'fs'], 'tcfArr','tcfNpairs','trajLength','trajAvgP2','trajRawP2','tauArray','tauArrayUsec','varOfArray');  

                 
%                 phase Factor avg TCF
                [ tcfArr, tcfNpairs ] = TCFcalc((trajAvgP2), tauArray );
                 curOutPhsFctAvg=[curOutFold  filesep() 'P2' filesep() 'PhaseFactAvg'];
                 if exist(curOutPhsFctAvg, 'dir')~= 7
                    mkdir(curOutPhsFctAvg);
                 end
                 varOfArray=nanvar((trajAvgP2));
                 save([curOutPhsFctAvg filesep() scanID '_' label 'fs'], 'tcfArr', 'tcfNpairs','trajLength','trajAvgP2','trajRawP2','tauArray','tauArrayUsec','varOfArray');  

%                  phase Factor raw TCF
                [ tcfArr, tcfNpairs ] = TCFcalc((trajRawP2), tauArray );                
                curOutPhsFctRaw=[curOutFold filesep() 'P2' filesep() 'PhaseFactRaw' ];  
                 if exist(curOutPhsFctRaw, 'dir')~= 7
                    mkdir(curOutPhsFctRaw);
                 end
                 varOfArray=nanvar((trajRawP2));
                save([curOutPhsFctRaw filesep() scanID '_' label 'fs'], 'tcfArr', 'tcfNpairs','trajLength','trajAvgP2','trajRawP2','tauArray','tauArrayUsec','varOfArray');  

        
                    
                end 
                
                fclose('all');
                
            end
            
        end
%               
      
