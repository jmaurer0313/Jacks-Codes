% Goal: read in the stream files from the polarization experiment performed
% with the FPGA, convert the numbers 0->1 to 64 bins (on 0:2pi) and then
% histogram the reuslts.

% NOTE: must be in the directory containing the stream files to operate on
% (ideally they are moved outside of "streams" to prevent operation on the
% entire folder, which is mostly not polarization data).

traceFolder= pwd;
filesList= dir('*-Stream');

scanLengths=[];

phaseMap=linspace(0,1,64); 
histEdges=linspace(0,2*pi,64); 

outfolderName= ('FPGA_PolarizationScans');
if exist(outfolderName, 'dir')~=7
    mkdir(outfolderName)
end

% generate a general number of rows/columns for a subplot which handles an
% arbitrary number of scans
numPlots=numel(filesList);
dimRow= floor(sqrt(numPlots)); 
dimCol= floor(sqrt(numPlots));

if ((dimRow*dimCol)<numPlots)
    dimRow=dimRow+1;
    if (dimRow*dimCol)<numPlots
        dimCol=dimCol+1;
    end
end 

for i=1:numel(filesList)
%     
%     
%     numel(filesList)
%       for f=3:3

 disp(['Currently on File ' num2str(i) ' of ' num2str(numel(filesList))]);
% mainDir and Fileslisr will ne offset by 2 in the loop to account for the
% '.' and '..' directories in mainDir
targetFolder=[traceFolder,'\',filesList(i).name,'\'];
scanID = filesList(i).name(1:15);

% open the files inside of the the target folder
timeID = fopen([targetFolder 'time.bin']);
p1ID = fopen([targetFolder 'p1.bin']);
p2ID = fopen([targetFolder 'p2.bin']);

% read in the files designated at the target locations
time = fread(timeID,Inf,'uint64=>uint64',0,'s');
p1 = fread(p1ID,Inf,'float64=>double',0,'s');
p2 = fread(p2ID,Inf,'float64=>double',0,'s');
    
% intialize an array of the order p1/p2 to store the mapped phases in
p1Map=zeros(1,numel(p1));
% p2Map=zeros(1,numel(p2));


% times here gets converted to seconds
times = ((double (time-time(1))./(8e7))+1);

% take out the first time from all later times, making event 1 occur at
% t=0sec, then convert the times to usec in another array
times=times-times(1);
timesUsec=times*1e6; 

% scan length is in minutes, scanLengths is just a holder for all scan
% times at a glance
scanLength= times(end)/60;
scanLengths=[scanLengths ; scanLength];

% map the p1 and p2 lists onto 64 possible values, these can be converted
% to radians by taking n/64*2pi. Im going to just do p1 for now, since we
% have only a single ref (so p1 and p2 are totally equal)

    for j=1:numel(p1)
    
%         this current style of mapping will place the phases in the
%         "nearest" bin. To perfectly replicate what Cliffs electronics do
%         we'll need to have it always sorted into the closest bin ROUNDING
%         DOWN. leaving it for now cause this is just a check. 
        curMap=abs(phaseMap-p1(j));
        [minVal, Imin]=min(curMap); 
        p1Map(j)=((Imin/64)*2*pi); 
        
    end
% save the processed mapped files to an output folder, can access the raw
% data there for additonal plotting/export
     foutName=[outfolderName filesep() scanID '_rawData'];
     save(foutName, 'p1Map','p1','times','timesUsec','scanLength');    

%      subplot figure to watch the data as it comes off
     figure(1)
     suptitle('Raw phase histograms using FPGA'); 
     subplot(dimRow,dimCol,i);
     histogram(p1Map,histEdges);
     xlim([0,2*pi]);
     title(scanID, 'FontSize',7); 
     xlabel('Phase(rads)'); 
     ylabel('Occurences'); 
     drawnow(); 
     
     

end 
