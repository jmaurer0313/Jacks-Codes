% Goal: Create a trace cutting routine that maintains the "click for cut"
% functionality of the original, but goes straight from the raw data to a
% cut, processed and saved file for direct use in XYQuadCalc/TCF scripts.
% This will replace traceCutter->tracetimeGetter->chosenAPD2Mat work flow
% for FPGA data. Could be implemented for Cliff data too (in the Case of
% FRET measurements) 

% ConstructName='(+1)Dimer_t21=0'; 
ConstructName='(+15)Dimer';
AutoCutTraces=0;

rawDir=pwd; 
traceFiles = dir(['*Stream']);
phaseMap=linspace(0,1,64); 

startTrace=1;

if isempty(traceFiles) == 1
    error(['Cannot find any tracefiles with format "-Stream". Exiting Program.']);
    
end
NtraceFiles = numel(traceFiles);

% set a usres for the cutting, by default it will be 30000 *30ms) for
% traceCutting
usres = 30000;
% usres = 100000;

msres = usres*10^-3;%('i.e. 1000 usec = 1 msec')
res = usres*10^-6;%(i.e. 1 million useconds turn into 1 second)
disp(['Setting the resolution equal to ' num2str(usres) ' microseconds (' num2str(msres) ' msec) (' num2str(res) ' sec)']);

% get the upper level dir to create a chosen moilecules folder inm, that
% way the raw data and chosen molecules are at the same level in the dir
% structure
cd .. 
upperDir=pwd;
cd(rawDir); 

% create an output folder for the saved data traces after cutting 
outputFolderName = [upperDir filesep() 'ChosenMolecules' filesep() ConstructName];
if exist(outputFolderName,'dir') == 0
    disp(['Making a folder called ' outputFolderName ' to hold the data.']);
    mkdir(outputFolderName);
elseif exist(outputFolderName,'dir') == 7
    disp(['Detected ' outputFolderName ' folder Already.']);
end

% SO NOW, open the stream files one by one, read them in, plot the data
% using a 3 panel subplot (counts, Phase, Amplitude) at 30ms resolution.
% Make the "cuts" to the plotted data, store those times, save the
% processed data equal to the subset designated by the cut times for each
% file. 
NcutFiles=0;
for i=startTrace:numel(traceFiles)
%     
%     
%     numel(filesList)
%       for f=3:3

 disp(['Currently on File ' num2str(i) ' of ' num2str(numel(traceFiles))]);
% mainDir and Fileslisr will ne offset by 2 in the loop to account for the
% '.' and '..' directories in mainDir
targetFolder=[rawDir,'\',traceFiles(i).name,'\'];
scanID = traceFiles(i).name(1:15);

% open the files inside of the the target folder
timeID = fopen([targetFolder 'time.bin']);
p1ID = fopen([targetFolder 'p1.bin']);
p2ID = fopen([targetFolder 'p2.bin']);

% read in the files designated at the target locations
time = fread(timeID,Inf,'uint64=>uint64',0,'s');
p1 = fread(p1ID,Inf,'float64=>double',0,'s');
p2 = fread(p2ID,Inf,'float64=>double',0,'s');

fclose('all'); 
    
% intialize an array of the order p1/p2 to store the mapped phases in
p1Map=zeros(1,numel(p1));
p2Map=zeros(1,numel(p2));


% times here gets converted to seconds
times = ((double (time-time(1))./(8e7))+1);

% take out the first time from all later times, making event 1 occur at
% t=0sec, then convert the times to usec in another array
times=times-times(1);
timesUsec=times*1e6; 

% scan length is in seconds, scanLengths is just a holder for all scan
% times at a glance
scanLength= times(end);
% scanLengths=[scanLengths ; scanLength];

% map the p1 and p2 lists onto 64 possible values, these can be converted
% to radians by taking n/64*2pi. Im going to just do p1 for now, since we
% have only a single ref (so p1 and p2 are totally equal)

    for j=1:numel(p1)
    
%         this current style of mapping will place the phases in the
%         "nearest" bin. To perfectly replicate what Cliffs electronics do
%         we'll need to have it always sorted into the closest bin ROUNDING
%         DOWN. leaving it for now cause this is just a check. 
        curMap1=abs(phaseMap-p1(j));
        [~, Imin]=min(curMap1); 
        p1Map(j)=((Imin/64)*2*pi); 
        
        curMap2=abs(phaseMap-p2(j));
        [~, Imin]=min(curMap2); 
        p2Map(j)=((Imin/64)*2*pi); 
        
    end

%     create a display trace at 30ms resolution for cutting purposes.
     numWindows=floor((scanLength*1e6)/usres); 
     plotTime=linspace(0,scanLength,numWindows); 
    countsTraj = zeros(1,numWindows);
    PhaseFctTrajP1 = zeros(1,numWindows);
    PhaseFctTrajP2 = zeros(1,numWindows);
    
%          *****THIS IS THE NEW REPHASE METHOD*******
% This method takes the X and Y quadrature for the ENTIRE TRACE, which
% should reflect the overall average phase as a value away from zero phase.
% By subtracting off this angle from all phases within the list the
% distribution of phases will be centered at zero on the interval of -pi to
% pi whicch atan2 returns. This minimizes "toggling" erros across -Theta to
% Theta in neighboring bins. 
   XquadP1 = mean(cos(p1Map));
   YquadP1 = mean(sin(p1Map));
   
   XquadP2 = mean(cos(p2Map));
   YquadP2 = mean(sin(p2Map));
   
   shiftAngleP1=(pi/32)*round(atan2(YquadP1,XquadP1)*(32/(pi)));
   shiftAngleP2=(pi/32)*round(atan2(YquadP2,XquadP2)*(32/(pi)));
   
   p1Map = p1Map-(shiftAngleP1);
   p2Map = p2Map-(shiftAngleP2); 

   
   for n=1:numWindows
   Phases1= p1Map(find(timesUsec<=n*usres & timesUsec>=(n-1)*usres ));
   Phases2= p2Map(find(timesUsec<=n*usres & timesUsec>=(n-1)*usres ));
%    TotalPhases =[Phases1;Phases2];
    p1Map=p1Map.'; 
    p2Map=p2Map.'; 
    
% create the 30ms trajectories for plotting purposes. 
   countsTraj(n)= numel(Phases1);
   PhaseFctTrajP1(n)= mean(exp(1i*(Phases1)));
   PhaseFctTrajP2(n)= mean(exp(1i*(Phases2)));
   end 
   
   
%    open a figure to plot the trajectory in
figure(1);
%     set(gcf,'Position',[1 462 1500 900]);
    set(gcf,'Position',[50 110 1200 650]);
    set(gcf,'Color','w');

   
%    Begin plotting and cutting the traces based on user input

user_choice = '';
traceNumber = 1;

if AutoCutTraces
            traceFileName= [scanID '-cut.mat']; 
            fileOutName = [outputFolderName filesep() traceFileName];
            firstFrame=1;
            lastFrame=numel(plotTime); 
            save(fileOutName,'p1Map','p2Map','timesUsec','scanLength','shiftAngleP1','shiftAngleP2','p1','p2','firstFrame','lastFrame');
            NcutFiles = NcutFiles + 1;
            traceNumber = traceNumber - 1;
else

   while strcmp(user_choice,'q') ~= 1
%     if traceNumber > length(traceFiles)
%         disp('*****Reached the end of the line. Restarting at molNumber = 1;');
%         traceNumber = 1;
%     elseif traceNumber < 1
%         disp(['*****Reached the beggining of the line. Restarting at molNumber = ' num2str(length(traceFiles)) ]);
%         traceNumber = length(traceFiles);
%     end
    traceFileName= [scanID '-cut.mat']; 
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
            save(fileOutName,'p1Map','p2Map','timesUsec','scanLength','shiftAngleP1','shiftAngleP2','p1','p2','firstFrame','lastFrame');
            NcutFiles = NcutFiles + 1;
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
            p1Map = p1Map([cut_indices]);
            p2Map = p2Map([cut_indices]);
            p1 = p1([cut_indices]);
            p2 = p2([cut_indices]);
           
            
            disp(['     Saving frames ' num2str(firstFrame) ' - ' num2str(lastFrame) ' of ' traceFileName]);
            save(fileOutName,'p1Map','p2Map','timesUsec','scanLength','shiftAngleP1','shiftAngleP2','p1','p2','firstFrame','lastFrame');
            NcutFiles = NcutFiles + 1;
            traceNumber = traceNumber - 1;
        case 'd'
            delete(fileOutName);
            NcutFiles = NcutFiles - 1;
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
   
end
disp([' Cut a total of ' num2str(NcutFiles) ' files during this session.']);