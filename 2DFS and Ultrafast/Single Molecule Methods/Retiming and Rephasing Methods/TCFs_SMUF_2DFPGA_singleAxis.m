
%AUTHOR: Jack Maurer 
% Goal:  1). bring in 2dfpga data sets, of 1XN
%        2). retime the data, cut it, sort it, calc the trajs and TCFs
%        3). calc the wAVG

%%
fileFolder = 'C:\Users\Ryzen 5\Dropbox\SM_UF_Data\1DTimePhaseTags_2DFPGA\E3E4\03082023\';
outFolderGen= 'C:\Users\Ryzen 5\Dropbox\SM_UF_Data\1DTimePhaseTags_2DFPGA\Cy3Films\OutputData\20230322\ChosenMolecules';
ConstructName='E3E4'; 

cd(fileFolder); 
autoCutAll=1;
plotTimeDomain=0;
doubleWideFolders=1; 

if doubleWideFolders
outFolderGen= 'C:\Users\Ryzen 5\Dropbox\SM_UF_Data\1DTimePhaseTags_2DFPGA\E3E4\OutputData\doubleWide\03082023\ChosenMolecules';
end

filesList = dir('*2DFPGA'); 
% scanNum=11;
numFiles = numel(filesList); 
curPath = pwd;
zeroPad=10; 
splinePoints=200; 
usres=60000;

msres=usres/1000;

c0 = 0.000299792458; % speed of light mm/fs
lam_mono = 550; % monochromater wavelength in nanometers
Fmono = 1e7 / lam_mono;

%get the list of files NOT to cut
[elim_IDlist] = snapShot(fileFolder,usres);
%%
% for i=1:1
    for i=1:length(filesList)

    
 
        fileFolder= [curPath filesep() filesList(i).name]; 
        scanID = filesList(i).name(1:15);

        tb1 = dlmread(fullfile(fileFolder, 'timebase1.txt'));
        tb2 = dlmread(fullfile(fileFolder, 'timebase2.txt'));
        xl = length(tb1);
        yl = length(tb2);
        
        if i==1
            OGTimeInc = (tb2(2)-tb2(1))*1000;
            OGTimeInc=round(OGTimeInc);
            if doubleWideFolders
            halfBinWidth=OGTimeInc/2;
            else
             halfBinWidth=OGTimeInc/4;   
            end
            arbTimeAxis=linspace(1,length(tb2),length(tb2));
            arbTimeAxis=[0,arbTimeAxis*OGTimeInc]; 
            lowerBounds=arbTimeAxis-halfBinWidth; 
            upperBounds=arbTimeAxis+halfBinWidth;
            globBounds=[];
            foldTable=[];
            
            for j=1:length(lowerBounds)
                
                if doubleWideFolders
                globBounds=[globBounds;[lowerBounds(j), upperBounds(j)]];
                outFolder1 = ([outFolderGen filesep() ConstructName '_' num2str(usres) 'usec'  filesep() 'Tau=' num2str(lowerBounds(j)) '-' num2str(upperBounds(j)) 'fs']);
                foldTable=[foldTable, struct('f',outFolder1)];
                else 
                    
                if j==length(lowerBounds)
                globBounds=[globBounds;[lowerBounds(j), upperBounds(j)]]; 
                else
                globBounds=[globBounds;[lowerBounds(j), upperBounds(j)] ; [upperBounds(j), lowerBounds(j+1)]];
                end
                
                if j==length(lowerBounds)
                outFolder1 = ([outFolderGen filesep() ConstructName '_' num2str(usres) 'usec'  filesep() 'Tau=' num2str(lowerBounds(j)) '-' num2str(upperBounds(j)) 'fs']);
                foldTable=[foldTable, struct('f',outFolder1)];
                else
                outFolder1 = ([outFolderGen filesep() ConstructName '_' num2str(usres) 'usec'  filesep() 'Tau=' num2str(lowerBounds(j)) '-' num2str(upperBounds(j)) 'fs']);
                outFolder2 = ([outFolderGen filesep() ConstructName '_' num2str(usres) 'usec'  filesep() 'Tau=' num2str(upperBounds(j)) '-' num2str(lowerBounds(j+1)) 'fs']);
                foldTable=[foldTable, struct('f',outFolder1),struct('f',outFolder2)];

                end
                
                end
                
            end
            % create a structure to hold the output folder names for use within the
            % loop
%             foldTable= struct('f1',outFolder1,'f2',outFolder2,'f3',outFolder3,'f4',outFolder4,'f5',outFolder5,'f6',outFolder6,'f7',outFolder7,'f8',outFolder8,'f9',outFolder9,'f10',outFolder10);             
        end
        

%       zero the remaining steps to start from step 1=0fs
%       ( need the -1 since the X-values are in the
%       -x direction )
        tb1=-1*(tb1 - tb1(1));
        tb2=(tb2 - tb2(1));
        tb1FS=tb1*1000;
        tb2FS=tb2*1000;
        
        T = 1000*(mean(tb1(2:end) - tb1(1:end-1))); % Uses stage positions directly, could replace by 0.0004/c0
        NFFT = 2^(nextpow2(numel(tb1))+1);            
        Fs = 10/(c0*T);
        Faxis = Fs/2*linspace(-1,1,NFFT);
        
        
        XmatTemp = zeros(xl, yl);
        YmatTemp = zeros(xl, yl);
%         
%         if j==1
%            Xmat= zeros(xl, yl);
%            Ymat= zeros(xl, yl);
%         end 
        
%%


        for yi = 1:yl
            for xi = 1:xl
                p1file = fullfile(fileFolder, [num2str(xi-1,'%02.0f') num2str(yi-1, '%02.0f') 'p1.bin']);
                p2file = fullfile(fileFolder, [num2str(xi-1,'%02.0f') num2str(yi-1, '%02.0f') 'p2.bin']);
                timefile = fullfile(fileFolder, [num2str(xi-1,'%02.0f') num2str(yi-1, '%02.0f') 'time.bin']);
                p1ID = fopen(p1file);
                p2ID = fopen(p2file);
                timeID = fopen(timefile);
                p1 = fread(p1ID,Inf,'float64=>double',0,'s');
                p2 = fread(p2ID,Inf,'float64=>double',0,'s');
                time = fread(timeID,Inf,'uint64=>uint64',0,'s');
                fclose('all');
%                 This will form a matrix of the 5 (x) and 8 (y) Khz components at
%                 every time step. The linear Abs Spectra can be obtained
%                 by transformation of the rows/columns. The 5Khz will only
%                 be defined along the rows and the 8Khz only defined
%                 alogn the columns, as the fringe pattern of the 5Khz is
%                 static along columns as is the 8Khz fringe pattern along
%                 the rows. 
%         
%                 XmatTemp(xi, yi) = mean(exp(2i*pi*(p1))*exp(-1i*myPhase));
%                 YmatTemp(xi, yi) = mean(exp(2i*pi*(p2))*exp(-1i*myPhase));
                
                   XmatTemp(xi, yi) = mean(exp(2i*pi*(p1)));
                YmatTemp(xi, yi) = mean(exp(2i*pi*(p2)));
            end
        end
        
        if plotTimeDomain
        figure()
         plot((tb2FS),real(XmatTemp))
        hold on
        plot((tb2FS),imag(XmatTemp))
        plot((tb2FS),abs(XmatTemp))
        end
        %%
        
%         define a spline axis
        splineAxis=linspace(tb2FS(1),tb2FS(end),splinePoints);
        mySplineData=spline(tb2FS,abs(XmatTemp),splineAxis);
        [val,idx]= max(mySplineData);
        tZeroDelta=splineAxis(idx); 
        
        correctTimeAxis=tb2FS-tZeroDelta;
        
        %now that the time corrections are obtained, process each file if
        %its not in the ElimList
        for yi = 1:yl
            for xi = 1:xl
                if sum(contains(elim_IDlist,[scanID '_' num2str(xi*yi)]))>0
                    continue;
                end
                p1file = fullfile(fileFolder, [num2str(xi-1,'%02.0f') num2str(yi-1, '%02.0f') 'p1.bin']);
                p2file = fullfile(fileFolder, [num2str(xi-1,'%02.0f') num2str(yi-1, '%02.0f') 'p2.bin']);
                timefile = fullfile(fileFolder, [num2str(xi-1,'%02.0f') num2str(yi-1, '%02.0f') 'time.bin']);
                p1ID = fopen(p1file);
                p2ID = fopen(p2file);
                timeID = fopen(timefile);
                p1 = fread(p1ID,Inf,'float64=>double',0,'s');
                p2 = fread(p2ID,Inf,'float64=>double',0,'s');
                time = fread(timeID,Inf,'uint64=>uint64',0,'s');
                fclose('all');
                
                times = ((double (time-time(1))./(8e7))+1);
                % take out the first time from all later times, making event 1 occur at
                % t=0sec, then convert the times to usec in another array
                times=times-times(1);
                timesUsec=times*1e6; 
                % scan length is in seconds, scanLengths is just a holder for all scan
                % times at a glance
                scanLength= times(end);
                
                numWindows=floor((scanLength*1e6)/usres); 
                plotTime=linspace(0,scanLength,numWindows); 
                countsTraj = zeros(1,numWindows);
                PhaseFctTrajP1 = zeros(1,numWindows);
                PhaseFctTrajP2 = zeros(1,numWindows);
                testTraj = zeros(1,numWindows);
                
                
                      XquadP1 = mean(cos(p1.*2*pi));
                      YquadP1 = mean(sin(p1.*2*pi));                   
                      shiftAngleP1=atan2(YquadP1,XquadP1);
                      p1 = (p1.*2*pi)-(shiftAngleP1);
                      
                        for n=1:numWindows
                           Phases1= p1(find(timesUsec<=n*usres & timesUsec>=(n-1)*usres ));                    
                        % create the 30ms trajectories for plotting purposes. 
                           countsTraj(n)= numel(Phases1);
                           PhaseFctTrajP1(n)= mean(exp(1i*(Phases1)));
                           testTraj(n)= mean(exp(1i*(mod(Phases1,2*pi))));
                        end 
                        user_choice = '';
                        traceNumber = 1;
                        
                        %start user input block here
                        if autoCutAll
                            
                            traceFileName= [scanID '-' num2str(yi) '-cut.mat']; 
                            curT21= abs(correctTimeAxis(yi));
                            curSort=globBounds-curT21;
                            [idx,bool]=find(curSort(:,1)<0);
                            properIdx=idx(end);
                            outputFolderName=foldTable(properIdx).f;
                            
                            if exist(outputFolderName, 'dir')~= 7
                                 mkdir(outputFolderName);
                            end                          
                            fileOutName = [outputFolderName filesep() traceFileName];
                            
                            disp(['     Saving All of ' traceFileName]);
                            firstFrame=1;
                            lastFrame=numel(plotTime); 
                            save(fileOutName,'timesUsec','scanLength','shiftAngleP1','p1','firstFrame','lastFrame');
                        
                        else
                            
                           while strcmp(user_choice,'q') ~= 1
                               
                        %     if traceNumber > length(traceFiles)
                        %         disp('*****Reached the end of the line. Restarting at molNumber = 1;');
                        %         traceNumber = 1;
                        %     elseif traceNumber < 1
                        %         disp(['*****Reached the beggining of the line. Restarting at molNumber = ' num2str(length(traceFiles)) ]);
                        %         traceNumber = length(traceFiles);
                        %     end
                        
                            traceFileName= [scanID '-' num2str(yi) '-cut.mat']; 
                            curT21= abs(correctTimeAxis(yi));
                            curSort=globBounds-curT21;
                            [idx,bool]=find(curSort(:,1)<0);
                            properIdx=idx(end);
                            outputFolderName=foldTable(properIdx).f;
                            if exist(outputFolderName, 'dir')~= 7
                                 mkdir(outputFolderName);
                             end
                            fileOutName = [outputFolderName filesep() traceFileName];

                        %   
                            %----------| Fluorescence Intensity |-----------------------------
                            figure(1);
                            subplot(3,1,1);
                            plot(plotTime,countsTraj,'g');
                            title([ConstructName '  "' scanID '"  {\color{red} ' num2str(msres) ' ms}'],'fontsize',17);

                            ylabel(['{\color{green}Cy3}'],'FontSize',15);
                            xlabel('Time (sec) ','fontsize',15);
                            grid on;
                            ylim([0,inf]);
                            drawnow();

                            %----------| Phase Trajectory |----------------------------------
                            subplot(3,1,2);
                        %     plot(plotTime,angle(PhaseFctTrajP2),'r',plotTime,angle(PhaseFctTrajP1),'b');
                            plot(plotTime,angle(PhaseFctTrajP1),'b');
                            xlabel('Time (sec) ','fontsize',15);
                            ylabel('Phase: {\color{red}P2} {\color{blue}P1}','fontsize',15);
                            ylim([-pi,pi]);
                            grid on;
                            drawnow();

                            %----------| Amplitude Trajectory |----------------------------------
                            subplot(3,1,3);
                        %     plot(plotTime,abs(PhaseFctTrajP2),'r',plotTime,abs(PhaseFctTrajP1),'b');
                            plot(plotTime,abs(PhaseFctTrajP1),'b'); 
                            xlabel('Time (sec) ','fontsize',15);
                            ylabel('Amplitude: {\color{red}P2} {\color{blue}P1}','fontsize',15);
                            ylim([0,1]);
                            grid on;
                            drawnow();
                            set(findobj('type','axes'),'fontsize',15);


                            %------------------------------------------------------------------
                            % Program has detected the trace file in chosen molecules folder
                            %------------------------------------------------------------------
                            if exist(fileOutName,'file') == 2

                        %         if the file already has been cut and processed, overlya the
                        %         results of the old cut for comparison 
                                hold on;        
                                disp(['     Already Detected ' traceFileName ' in ' outputFolderName]);      
                                load(fileOutName, 'firstFrame','lastFrame'); 
                                %----------| Fluorescence Intensity |-----------------------------


                                subplot(3,1,1);
                                hold on;
                                plot(plotTime(firstFrame:lastFrame),countsTraj(firstFrame:lastFrame),'r');
                                hold off;

                                %----------| Phase Trajectory |----------------------------------
                                subplot(3,1,2);
                                hold on;
                        %         plot(plotTime(firstFrame:lastFrame),angle(PhaseFctTrajP1(firstFrame:lastFrame)),'k',plotTime(firstFrame:lastFrame),angle(PhaseFctTrajP2(firstFrame:lastFrame)),'y');
                                plot(plotTime(firstFrame:lastFrame),angle(PhaseFctTrajP1(firstFrame:lastFrame)),'k');       
                                hold off;

                                %----------| Amp Trajectory |----------------------------------
                                subplot(3,1,3);
                                hold on;
                        %         plot(plotTime(firstFrame:lastFrame),abs(PhaseFctTrajP1(firstFrame:lastFrame)),'k',plotTime(firstFrame:lastFrame),abs(PhaseFctTrajP2(firstFrame:lastFrame)),'y');
                                plot(plotTime(firstFrame:lastFrame),abs(PhaseFctTrajP1(firstFrame:lastFrame)),'k');        
                                hold off;
                            end

                            %-----------| USER INTERFACE |----------------------------------
                            prompt = ['     ' scanID ' : refresh/all/cut/next?/delete [enter/a/c/n/d]?: '];
                        %    num2str(traceNumber) '/' num2str(NtraceFiles)  **old input above**
                            user_choice = input(prompt,'s');
                            if isempty(user_choice)
                                user_choice = 's';
                            end
                            switch user_choice
                                case 'n'
                                    break;
                                case 'a'%save the entire trace
                                    disp(['     Saving All of ' traceFileName]);
                                    firstFrame=1;
                                    lastFrame=numel(plotTime); 
                                    save(fileOutName,'timesUsec','scanLength','shiftAngleP1','p1','firstFrame','lastFrame');
%                                     NcutFiles = NcutFiles + 1;
                                    traceNumber = traceNumber - 1;
                                case 'c' %Curtrail (cut) the trace at specified time
                                    disp('     Click on the photobleach time or press [enter] for all');
                                    [endtime,~] = ginput;
                                    if isempty(endtime)
                                        endtime = plotTime(end);
                                    end
                                    if numel(endtime) == 1
                        %                 first frame and lastFrame are the indices of plotTime
                        %                 (incremented at res) which are closest to the posistion
                        %                 of the click
                                        firstFrame = 1;
                                        lastFrame = find(plotTime <= endtime,1,'last');
                                    elseif numel(endtime) == 2
                                        firstFrame = find(plotTime >= endtime(1),1,'first');
                                        lastFrame = find(plotTime <= endtime(2),1,'last');
                                    elseif numel(endtime) > 2
                                        disp('Too many inputs. Click twice only.');
                                        traceNumber = traceNumber - 1;
                                        continue
                                    end
                        %             determine the nearest time to the "cut time" in the raw data
                        %             files, this will then be what is saved to the output folder.
                        %             plotTime is in seconds.
                                    cut_indices = find(timesUsec<=plotTime(lastFrame)*1e6 & timesUsec>=plotTime(firstFrame)*1e6 );
                                    timesUsec = timesUsec([cut_indices]);                                   
                                    p1 = p1([cut_indices]);


                                    disp(['     Saving frames ' num2str(firstFrame) ' - ' num2str(lastFrame) ' of ' traceFileName]);
                                    save(fileOutName,'timesUsec','scanLength','shiftAngleP1','p1','firstFrame','lastFrame');
%                                     NcutFiles = NcutFiles + 1;
                                    traceNumber = traceNumber - 1;
                                case 'd'
                                    delete(fileOutName);
%                                     NcutFiles = NcutFiles - 1;
                                    traceNumber = traceNumber - 1;
                                case 'b'
                                    traceNumber = traceNumber - 2;
                                case 'g'
                                    skipNum = input('Skip how many molecules?');
                                    traceNumber = traceNumber + skipNum -1;
                        %         case 'l'+
                        %             if showLabView_mode
                        %                 showLabView_mode = 0;
                        %             else
                        %                 showLabView_mode = 1;
                        %             end
                        %             traceNumber = traceNumber - 1;
                            end


                            traceNumber = traceNumber + 1;
                           end
                        end
                         %end of the user choice block
                    
                    
            end
        end
        
end 
%%

% % START: alternate retime based on FT interpolation - NOT GREAT for
% sparsely sampled time domain

% Tdes=0.25;
% Ts=abs(tb2(2)-tb2(1))*1000;
% NFFT = 2^(nextpow2(length(tb2))); 
% NIFFT = 2^(nextpow2((Ts/Tdes)*NFFT)); 
% Npad = abs(NIFFT-NFFT);
% f = ceil((-NFFT/2:NFFT/2-1))/(Ts*NFFT); % freq axis in PHz
% windowX = delayedGaussian(tb2*1000, 0, 15);
% 
% %     this is the case of rows for FFT
% XLin = XmatTemp(1,:).*windowX;
% 
% xPhases=zeros(1,size(XmatTemp,1));
% 
% for j=1:length(xPhases)
%    xPhases(j)=angle(XLin(j,1));
% end
% 
% % multiply the first entry in the matrix by 1/2 to make the transform
% % symmetric (see weird NMR reference...)
% XLin(:,1) = XLin(:,1)*0.5; 
% 
% XFLin = fftshift(fft(XLin, NFFT, 2)); % Signal in the frequency domain, FFT shift applied
% 
% 
% figure()
% plot(f, abs(XFLin.'));
% % plot(f, real(SFLin.'), f, imag(SFLin.'), f, abs(SFLin.'));
% title('X Linear Abs vectors - raw/NOT-retimed');
%  
% 
% % try dropping the zeros at the edge of the freq domain, OR the wrap around
% % point where freq goes + to - 
% XFLinzeros=[real(XFLin) , zeros(size(XFLin,1),Npad)];
% unshiftedXLinNew= circshift(XFLinzeros, -floor(NFFT/2),2);
% 
% RunshXLin = unshiftedXLinNew;  
% 
% XnLin =(ifft(RunshXLin, NIFFT, 2));
% 
% % block that may be neccesary to zero out the newly added space due to
% % interpolation in the time domain
% XnLin(:,((ceil(NIFFT/2)+1):end))=0;
% 
% Tn = (Ts*NFFT)/NIFFT;
% tn = ((0:NIFFT-1))*Tn;
% % tnX=(-(NIFFT/length(tx)):((NIFFT/length(tx))*(length(tx)-1))-1)*Tn; 
% tnX= tn-Ts; 
% 
% % rescaling and plotting matrices for X
% XnLin(:,1) = XnLin(:,1)*2; 
% XLin(:,1) = XLin(:,1)*2; 
% XnLinPlot= XnLin.';
% XLinPlot= XLin.';
% 
% XLinPlot=XLinPlot*NFFT/NIFFT;
%  
% figure()
% plot(tn, abs(XnLinPlot));
% title('X Linear Scans -Abs- Interpolated');
% xlabel('Time (fs)'); 
% ylabel('Magnitude (arb)'); 
% % plot(tn, real(SnLinPlot), tn, imag(SnLinPlot))
% hold on
% ax = gca;
% ax.ColorOrderIndex = 1;
% % plot(tb2*1000', abs(XmatTemp),'o');
% hold off
% 
% % create a vector for the maxima in the interpolated X and Y linear data,
% % make a vector of the indices as well as the value 
% [txShifts,IxShifts]=max(abs(XnLin.'));
% 
% % the maxtrix for X and Y shifts contaisn the raw time by which the scans
% % are off from time zero and the index of the maximum in the interpolated
% % time domain of order NIFFT. IN THE CASE OF THE Y SHIFTS, AS LONG AS THE
% % TY AXIS STARTS FROM ZERO IT SHOULD BE THE iTH+1 ELEMENT SINCE tnY(1)=0
% % Xshifts=[tnX([IxShifts])',IxShifts'];
% % Xzero=find(tnX==0);
% tn(IxShifts)
