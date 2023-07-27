disp(['Your current directory is' string(pwd)]);
traceFolder= pwd;
filesList= dir('*-Stream');
mainDir=dir;  

nH = 15;
harmsTotal = zeros(1, nH+1);
countsTotal=0;
stdMat=zeros(numel(filesList), nH+1); 
counts=[];
for f=1:numel(filesList)
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
% p2r = (2.*pi).*p2;


harms = zeros(1, nH);
counts = [counts;numel(p1r)];
for j=0:nH
    harms(:, j+1)= mean(exp(1i*j*p1r))*numel(p1r);
    stdMat(f,j+1)= abs(mean(exp(1i*j*p1r)));
    
end

harmsTotal = harmsTotal + harms;
countsTotal = countsTotal+numel(p1r); 




%     std of the indiv file means weighted by the number of photons from that file
    


% phi1 = mean(exp(1i*p1r));
% phi2 = mean(exp(1i*2*p1r));
% phi3 = mean(exp(1i*3*p1r));
% phi4 = mean(exp(1i*4*p1r));

% bar(0:nH, abs(harms))
% hold on




end 
stdDev = std(stdMat,counts); 
harmsTotal = harmsTotal./(countsTotal*numel(filesList));
harmRatio = abs(harmsTotal(3))/abs(harmsTotal(2))
bar(0:nH, abs(harmsTotal));
hold on
plot(0:nH, stdDev);    
hold off
disp(['Ratio of 2f/1f is ' num2str(harmRatio)]);
disp(['Std of 2f is ' num2str(stdDev(3))]);