disp(['Your current directory is' string(pwd)]);
traceFolder= pwd;
filesList= dir('*-Stream');
mainDir=dir;  

nH = 4;
harmsRPTotal = zeros(1, nH+1);
harmsNRPTotal = zeros(1, nH+1);
countsTotal=0;

harmsRHref = zeros(numel(filesList), nH);
harmsLHref = zeros(numel(filesList), nH);

harmsRP = zeros(numel(filesList), nH);
harmsNRP = zeros(numel(filesList), nH);
% This size is such that every single harmonic of every single file has a
% standard dev defined for it
stdMatRP=zeros(numel(filesList), nH+1); 
stdMatNRP=zeros(numel(filesList), nH+1);
counts=[];
scanLengths=[];
figOneRef=figure();
figNonLin=figure();
figTwoRef=figure();

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

fclose('all');

% Since the sorting script starts at binNum=1, all times must be 1 greater
% to avoid indexing error
times = ((double (time-time(1))./(8e7))+1);
scanLength= times(end)/60;
scanLengths=[scanLengths ; scanLength];
%  if f==2
%      times=times(1:round(0.8*length(time))
p1r = (2.*pi).*p1;
p2r = (2.*pi).*p2;

RP= p2r - p1r;
NRP= p2r + p1r;

if numel(RP)~=numel(NRP)
    disp('Mismatch in length and number of elements in phase lists'); 
    break
end 

% sum is nonrephasing, difference is rephasing

% generate a matrix with a riw for each file, and N columns for each row to
% tabulate the N harmonics of both the RP and NRP phase combinations.
% Average the resulting rows afterward to preserve the individual file
% content.


% tracks the number of counts for each file so they can be weighted
% properly
countsRP = [counts;numel(RP)];
countsNRP = [counts;numel(NRP)];


Bsample = zeros(length(RP), nH);
for j=1:nH
    Bsample(:,j)=exp(1i*j*RP);
end
bfunc = @(Bsample) abs(mean(Bsample));
[bootstat,bootsam] = bootstrp(400,bfunc,Bsample);
absies = abs(harmsRP(5,2:end));
erries = std(bootstat);
percenterr = erries./absies*100;

% When calcualting the uncertainty in a SINGLE FILE at a SINGLE HARMONIC,
% it should be the stdDev(phaseList)/sqrt(Npoints) where Npoints is the
% number of points which goes into determining the mean/average phasefactor
for j=0:nH
    
    harmsRHref(f, j+1) = mean(exp(1i*j*p1r));
    harmsLHref(f, j+1) = mean(exp(1i*j*p2r));
    
    harmsRP(f, j+1)= mean(exp(1i*j*RP));
    harmsNRP(f, j+1)= mean(exp(1i*j*NRP));
    
    stdMatRP(f,j+1)= std(exp(1i*j*RP))/sqrt(numel(RP));
    stdMatNRP(f,j+1)= std(exp(1i*j*NRP))/sqrt(numel(NRP));
    
    
end

harmsRPTotal = harmsRPTotal + harmsRP(f,:)* numel(RP);
harmsNRPTotal = harmsNRPTotal + harmsNRP(f,:)* numel(NRP);



end

countsRPTotal = sum(countsRP);
countsNRPTotal = sum(countsNRP);

% stdDev = std(stdMat,counts); 

harmsRPTotal = harmsRPTotal./(countsRPTotal*numel(filesList));
harmsRPratio = abs(harmsRPTotal(2))/abs(harmsRPTotal(1));

harmsNRPTotal = harmsNRPTotal./(countsNRPTotal*numel(filesList));
harmsNRPratio = abs(harmsNRPTotal(2))/abs(harmsNRPTotal(1));

% figRP=figure();
% bar(0:nH, abs(harmsRPTotal));
% 
% figNRP=figure();
% bar(0:nH, abs(harmsNRPTotal));
% 
% % plot(0:nH, stdDev);    
% hold off

disp(['Ratio of RP 1f/0f is ' num2str(harmsRPratio)]);
disp(['Ratio of NRP 1f/0f is ' num2str(harmsNRPratio)]);

% generate a series of subplots. For each molecule display 
numPlots=numel(filesList);
    DimRow= floor(sqrt(numPlots)); 
    DimCol= floor(sqrt(numPlots)); 
    
    if ((DimRow*DimCol)< numPlots)
        DimRow=DimRow+1;
        if(DimRow*DimCol)< numPlots
            DimCol=DimCol+1; 
        end
    end
    
    
    twoRefLabel={'5kHz 1f' '8kHz 1f' 'RP 1f' 'NRP 1f'};
    oneRefLabel={'5kHz 1f' '8kHz 1f' '5kHz 2f' '8kHz 2f'}; 
    nonLinLabel={'RP 1f' 'NRP 1f' 'RP 2f' 'NRP 2f'};
    
    figure(figTwoRef);
    for k=1:numel(filesList)  
% subplot routine for the two ref combinations (one ref combos to follow) 

subplot(DimRow,DimCol,k);
    
%       plot(ticksFRET,FRETHist);
    data = [abs(harmsRHref(k,2)) abs(harmsLHref(k,2)) abs(harmsRP(k,2)) abs(harmsNRP(k,2))];
    bar([1 2 3 4],[abs(harmsRHref(k,2)) abs(harmsLHref(k,2)) abs(harmsRP(k,2)) abs(harmsNRP(k,2))]) ;
    hold on;
    set(gca, 'XTickLabel',twoRefLabel);
    ylabel('Absolute Magnitude','FontSize',8);
    ylim([0 1.1*max(data)]) 
    title(['Scan=' scanID '  Duration=' num2str(scanLengths(k)) 'min'],'FontSize',8); 
%     axis tight;    
    drawnow();
    sgtitle('Linear and Nonlinear First Harmonic Magnitude Comparisons');
    
    end 
    
    figure(figOneRef);
    for l=1:numel(filesList)  
% subplot routine for the two ref combinations (one ref combos to follow) 

subplot(DimRow,DimCol,l);
    
%       plot(ticksFRET,FRETHist);
    data = [abs(harmsRHref(l,2)) abs(harmsLHref(l,2)) abs(harmsRHref(l,3)) abs(harmsLHref(l,3))];
    bar([1 2 3 4],[abs(harmsRHref(l,2)) abs(harmsLHref(l,2)) abs(harmsRHref(l,3)) abs(harmsLHref(l,3))]) ;
    hold on;
    set(gca, 'XTickLabel',oneRefLabel);
    ylabel('Absolute Magnitude','FontSize',8);
    ylim([0 1.1*max(data)]) 
    title(['Scan=' scanID '  Duration=' num2str(scanLengths(l)) 'min'],'FontSize',8); 
%     axis tight;    
    drawnow();
    sgtitle('Linear Harmonic Magnitude Comparisons');
    
    end 
    
    figure(figNonLin);
    for z=1:numel(filesList)  
% subplot routine for the two ref combinations (one ref combos to follow) 

subplot(DimRow,DimCol,z);
    
%       plot(ticksFRET,FRETHist);
    data = [abs(harmsRP(z,2)) abs(harmsNRP(z,2)) abs(harmsRP(z,3)) abs(harmsNRP(z,3))];
    bar([1 2 3 4],[abs(harmsRP(z,2)) abs(harmsNRP(z,2)) abs(harmsRP(z,3)) abs(harmsNRP(z,3))]) ;
    hold on;
    set(gca, 'XTickLabel',nonLinLabel);
    ylabel('Absolute Magnitude','FontSize',8);
    ylim([0 1.1*max(data)]) 
    title(['Scan=' scanID '  Duration=' num2str(scanLengths(z)) 'min'],'FontSize',8); 
%     axis tight;    
    drawnow();
    sgtitle('Nonlinear Harmonic Magnitude Comparisons');
    
    end 
