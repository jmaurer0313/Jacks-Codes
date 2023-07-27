% check and histogram the time elapsed between photon arrivals.
disp(['Your current directory is' string(pwd)]);
traceFolder= pwd;
filesList= dir('*-Stream');
mainDir=dir;  

for f=1:1
%     numel(filesList)
%       for f=3:3
binNum1=1; 
 disp(['Currently on File ' num2str(f) ' of ' num2str(numel(filesList))]);
% mainDir and Fileslisr will ne offset by 2 in the loop to account for the
% '.' and '..' directories in mainDir
targetFolder=[traceFolder,'\',filesList(f).name,'\'];
scanID = filesList(f).name(1:15);

timeID = fopen([targetFolder 'time.bin']);
p1ID = fopen([targetFolder 'p1.bin']);
p2ID = fopen([targetFolder 'p2.bin']);

time = fread(timeID,Inf,'uint64=>uint64',0,'s');
p1 = fread(p1ID,Inf,'float64=>double',0,'s');
p2 = fread(p2ID,Inf,'float64=>double',0,'s');
    

% Since the sorting script starts at binNum=1, all times must be 1 greater
% to avoid indexing error
times = ((double (time-time(1))./(8e7))+1);
%  if f==2
%      times=times(1:round(0.8*length(time))
p1r = (2.*pi).*p1;
p2r = (2.*pi).*p2;

deltaList=zeros(1,numel(times)-1);

for n= 1:(numel(times)-1)
    
    deltaT= times(n+1)-times(n); 
    deltaList(n)=deltaT;  
    
end



end 

histogram(deltaList, 10000); 
    
    