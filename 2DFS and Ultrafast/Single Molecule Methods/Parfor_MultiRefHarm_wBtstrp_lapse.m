% GOAL: load in single integration periods from 2D data sets to perform
% bootstrapping error analysis on. 
 

% *********SET PARAMETERS*******

% provide directory containing data files and a integration time in seconds
% (will loop over all individual scans in the containing folder)
fileFolder= 'C:\Users\Jack\Documents\Marcus Lab\2D_SingleMolecule\Sm_avgFolders\20200303_mol1';
startDir=pwd;
% provide the coordinates for the delay space integration period of interest
t21step=2;
t43step=1;
% specify the number of harmonics and number of bootstraps to perform. 
nH = 3;
numBoots=20;
% choose either a number of photons OR integration time (set the booleans accordingly)
int=1;
nPhot=0;

nPhotMax=5e6; 
intTime=42; 
% ******************************

cd(fileFolder); 
filesList= dir('*2DFPGA'); 
numFiles=numel(filesList); 
curPath=pwd;


harmsXref = zeros(numel(filesList), nH+1);
harmsYref = zeros(numel(filesList), nH+1);

harmsRP = zeros(numel(filesList), nH+1);
harmsNRP = zeros(numel(filesList), nH+1);
% This size is such that every single harmonic of every single file has a
% standard dev defined for it (including DC)
errMatRPx=zeros(numel(filesList), nH+1); 
errMatNRPx=zeros(numel(filesList), nH+1);
errMatRPy=zeros(numel(filesList), nH+1); 
errMatNRPy=zeros(numel(filesList), nH+1);
errMatRPabs=zeros(numel(filesList), nH+1); 
errMatNRPabs=zeros(numel(filesList), nH+1);

errMatX=zeros(numel(filesList), nH+1); 
errMatY=zeros(numel(filesList), nH+1); 


scanIDList=[];
scanLengths=[];
t21ax=[];
t43ax=[];

figOneRef=figure();
figNonLin=figure();
figNonLinQuad=figure();
figTwoRef=figure();

% Generate a ScanID List prior to the parrallel loop to reference later on
for w=1:numel(filesList)
   
scanID = filesList(w).name(1:15);
scanIDList = [scanIDList; scanID];

end 


parfor f=1:numFiles
    disp(['Currently on File ' num2str(f) ' of ' num2str(numel(filesList))]);
    fileFolder= [curPath filesep() filesList(f).name]; 


        tb1 = dlmread(fullfile(fileFolder, 'timebase1.txt'));
        tb2 = dlmread(fullfile(fileFolder, 'timebase2.txt'));
        xl = length(tb1);
        yl = length(tb2);
        
        xi=t21step;
        yi=t43step;
        
        numPhots=0;
              
                timeID = fullfile(fileFolder, [num2str(xi-1,'%02.0f') num2str(yi-1, '%02.0f') 'time.bin']);
                p1file = fullfile(fileFolder, [num2str(xi-1,'%02.0f') num2str(yi-1, '%02.0f') 'p1.bin']);
                p2file = fullfile(fileFolder, [num2str(xi-1,'%02.0f') num2str(yi-1, '%02.0f') 'p2.bin']);
                timeFile = fopen(timeID);
                p1ID = fopen(p1file);
                p2ID = fopen(p2file);
                
                if int
                time = fread(timeFile,Inf,'uint64=>uint64',0,'s');
                times = ((double (time-time(1))./(8e7))+1);
                times=times-times(1);
                p1 = fread(p1ID,Inf,'float64=>double',0,'s');
                p2 = fread(p2ID,Inf,'float64=>double',0,'s');
                [intWind]=find(times<=intTime);
%                At this point isolate the phase and time lists from zero
%                up to the integration time set
                scanLengths=[scanLengths;(times(intWind(end)))]; 
                p1=p1(1:intWind(end));
                p2=p2(1:intWind(end)); 
                numPhots = numPhots + numel(p1);
                p1r = (2.*pi).*p1;
                p2r = (2.*pi).*p2;
                RP= p2r - p1r;
                NRP= p2r + p1r;                
                fclose('all');
                
                end
                
                if nPhot && int
                    disp('Select either a number of photons OR an integration time (not both) - Cancel execution now'); 
                end 
                
                if nPhot
                    time = fread(timeID,nPhotMax,'uint64=>uint64',0,'s');
                    p1 = fread(p1ID,nPhotMax,'float64=>double',0,'s');
                    p2 = fread(p2ID,nPhotMax,'float64=>double',0,'s');
                    times = ((double (time-time(1))./(8e7))+1);
                    scanLength= times(end);
                    scanLengths=[scanLengths ; scanLength];
                    numPhots = numPhots + numel(p1);
                    p1r = (2.*pi).*p1;
                    p2r = (2.*pi).*p2;
                    RP= p2r - p1r;
                    NRP= p2r + p1r;                
                    fclose('all');
                end
                
                if numel(RP)~=numel(NRP)
                disp('Mismatch in length and number of elements in phase lists');
                end
                
                % Return time as a femtosecond column vector
                tb1 = tb1*1e3;
                tb2 = tb2*1e3;   
                
                if f==1
                t21ax=[t21ax;tb1];
                t43ax=[t43ax;tb2]; 
                end 
                
                % When calcualting the uncertainty in a SINGLE FILE at a SINGLE HARMONIC,
% it should be the stdDev(phaseList)/sqrt(Npoints) where Npoints is the
% number of points which goes into determining the mean/average phasefactor

BsampleRPx = zeros(length(RP), nH+1);
BsampleNRPx = zeros(length(NRP), nH+1);
BsampleRPabs = zeros(length(RP), nH+1);
BsampleNRPabs = zeros(length(NRP), nH+1);
BsampleRPy = zeros(length(RP), nH+1);
BsampleNRPy = zeros(length(NRP), nH+1);
BsampleX = zeros(length(p1r), nH+1);
BsampleY = zeros(length(p2r), nH+1);

  harmsXtemp= zeros(1, nH+1);
  harmsYtemp= zeros(1, nH+1);
  harmsRPtemp= zeros(1, nH+1);
  harmsNRPtemp= zeros(1, nH+1);


for j=0:nH
   
    
    BsampleRPabs(:,j+1)=abs(exp(1i*(j)*RP));
    BsampleNRPabs(:,j+1)=abs(exp(1i*(j)*NRP));
    BsampleRPx(:,j+1)=real(exp(1i*(j)*RP));
    BsampleNRPx(:,j+1)=real(exp(1i*(j)*NRP));
    BsampleRPy(:,j+1)=imag(exp(1i*(j)*RP));
    BsampleNRPy(:,j+1)=imag(exp(1i*(j)*NRP));
    BsampleX(:,j+1) = exp(1i*j*p1r);
    BsampleY(:,j+1) = exp(1i*j*p2r); 
    
    harmsXtemp(1, j+1) = mean(exp(1i*j*p1r));
    harmsYtemp(1, j+1) = mean(exp(1i*j*p2r));
    harmsRPtemp(1, j+1) = mean(exp(1i*j*RP));
    harmsNRPtemp(1, j+1) = mean(exp(1i*j*NRP));
       
    
end

    harmsXref(f,:) = harmsXtemp;
    harmsYref(f,:) = harmsYtemp;
    
    harmsRP(f,:) = harmsRPtemp;
    harmsNRP(f,:) = harmsNRPtemp;

% Bootstrap for RPx
bfuncRPx = @(BsampleRPx) (mean(BsampleRPx));
bootstatRPx = bootstrp(numBoots,bfuncRPx,BsampleRPx);
realRP = real(harmsRP(f,:));
errRPx = std(bootstatRPx);
perErrRPx = errRPx(2:end)./realRP(2:end)*100;

% Bootstrap for NRPx
bfuncNRPx = @(BsampleNRPx) (mean(BsampleNRPx));
bootstatNRPx = bootstrp(numBoots,bfuncNRPx,BsampleNRPx);
realNRP = real(harmsNRP(f,:));
errNRPx = std(bootstatNRPx);
perErrNRPx = errNRPx(2:end)./realNRP(2:end)*100;

% Bootstrap for RPy
bfuncRPy = @(BsampleRPy) (mean(BsampleRPy));
bootstatRPy = bootstrp(numBoots,bfuncRPy,BsampleRPy);
imagRP = imag(harmsRP(f,:));
errRPy = std(bootstatRPy);
perErrRPy = errRPy(2:end)./imagRP(2:end)*100;

% Bootstrap for NRPy
bfuncNRPy = @(BsampleNRPy) (mean(BsampleNRPy));
bootstatNRPy = bootstrp(numBoots,bfuncNRPy,BsampleNRPy);
imagNRP = imag(harmsNRP(f,:));
errNRPy = std(bootstatNRPy);
perErrNRPy = errNRPy(2:end)./imagNRP(2:end)*100;

% Bootstrap for RP abs
bfuncRPabs = @(BsampleRPabs) (mean(BsampleRPabs));
bootstatRPy = bootstrp(numBoots,bfuncRPabs,BsampleRPabs);
absRP = abs(harmsRP(f,:));
errRPabs = std(bootstatRPy);
perErrRPabs = errRPabs(2:end)./absRP(2:end)*100;

% Bootstrap for NRP abs
bfuncNRPabs = @(BsampleNRPabs) (mean(BsampleNRPabs));
bootstatNRPabs = bootstrp(numBoots,bfuncNRPabs,BsampleNRPabs);
absNRP = abs(harmsNRP(f,:));
errNRPabs = std(bootstatNRPabs);
perErrNRPabs = errNRPabs(2:end)./absNRP(2:end)*100;

% Bootstrap for X
bfuncX = @(BsampleX) abs(mean(BsampleX));
bootstatX = bootstrp(numBoots,bfuncX,BsampleX);
absX = abs(harmsXref(f,:));
errX = std(bootstatX);
perErrX = errX(2:end)./absX(2:end)*100;

% Bootstrap for Y
bfuncY = @(BsampleY) abs(mean(BsampleY));
bootstatY = bootstrp(numBoots,bfuncY,BsampleY);
absY = abs(harmsYref(f,:));
errY = std(bootstatY);
perErrY = errY(2:end)./absY(2:end)*100;


% harmsRPTotal = harmsRPTotal + harmsRP(f,:)* numel(RP);
% harmsNRPTotal = harmsNRPTotal + harmsNRP(f,:)* numel(NRP);

    errMatRPx(f,:)= errRPx;
    errMatNRPx(f,:)= errNRPx;
    errMatRPy(f,:)= errRPy;
    errMatNRPy(f,:)= errNRPy;
    errMatRPabs(f,:)= errRPabs;
    errMatNRPabs(f,:)= errNRPabs;
    
    errMatX(f,:)= errX;
    errMatY(f,:)= errY;

end

% generate a series of subplots. For each molecule scan display 
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
    nonLinQuadLabel={'Re RP 1f' 'Re NRP 1f' 'Im RP 1f' 'Im NRP 1f'};
    
    figure(figTwoRef);
    for k=1:numel(filesList)  
% subplot routine for the two ref combinations (one ref combos to follow) 

subplot(DimRow,DimCol,k);
    
%       plot(ticksFRET,FRETHist);
    data = [abs(harmsXref(k,2)) abs(harmsYref(k,2)) abs(harmsRP(k,2)) abs(harmsNRP(k,2))];
    errhigh= [errMatX(k,2) errMatY(k,2) errMatRPabs(k,2) errMatNRPabs(k,2)]; 
    errlow= [errMatX(k,2) errMatY(k,2) errMatRPabs(k,2) errMatNRPabs(k,2)];
    bar([1 2 3 4],[abs(harmsXref(k,2)) abs(harmsYref(k,2)) abs(harmsRP(k,2)) abs(harmsNRP(k,2))])
  
    hold on;
    
    er = errorbar([1 2 3 4],data,errlow,errhigh);    
    er.Color = [0 0 0]; 
    er.LineWidth= 2; 
    er.LineStyle = 'none';

    set(gca, 'XTickLabel',twoRefLabel);
    ylabel('Absolute Magnitude','FontSize',8);
    ylim([0 1.1*max(data)]) 
    title(['Scan=' scanIDList(k,:) '  Duration=' num2str(scanLengths(k)) 'sec'],'FontSize',8); 
%     axis tight;    
    drawnow();
    sgtitle(['Linear and Nonlinear First Harmonic Magnitude Comparisons - t21=' num2str(t21ax(t21step)) ' t43=' num2str(t43ax(t43step))]);
    
    end 
    
    figure(figOneRef);
    for l=1:numel(filesList)  
% subplot routine for the two ref combinations (one ref combos to follow) 

subplot(DimRow,DimCol,l);
    
%       plot(ticksFRET,FRETHist);
    data = [abs(harmsXref(l,2)) abs(harmsYref(l,2)) abs(harmsXref(l,3)) abs(harmsYref(l,3))];
    errhigh= [errMatX(l,2) errMatY(l,2) errMatX(l,3) errMatY(l,3)]; 
    errlow= [errMatX(l,2) errMatY(l,2) errMatX(l,3) errMatY(l,3)];
    bar([1 2 3 4],[abs(harmsXref(l,2)) abs(harmsYref(l,2)) abs(harmsXref(l,3)) abs(harmsYref(l,3))]) ;
    hold on;
    
    er = errorbar([1 2 3 4],data,errlow,errhigh);    
    er.Color = [0 0 0]; 
     er.LineWidth= 2;
    er.LineStyle = 'none';
    
    set(gca, 'XTickLabel',oneRefLabel);
    ylabel('Absolute Magnitude','FontSize',8);
    ylim([0 1.1*max(data)]) 
    title(['Scan=' scanIDList(l) '  Duration=' num2str(scanLengths(l)) 'sec'],'FontSize',8); 
%     axis tight;    
    drawnow();
    sgtitle(['Linear Harmonic Magnitude Comparisons - t21=' num2str(t21ax(t21step)) ' t43=' num2str(t43ax(t43step))]);
    
    end 
    
    figure(figNonLin);
    for z=1:numel(filesList)  
% subplot routine for the two ref combinations (one ref combos to follow) 

subplot(DimRow,DimCol,z);
    
%       plot(ticksFRET,FRETHist);
    data = [abs(harmsRP(z,2)) abs(harmsNRP(z,2)) abs(harmsRP(z,3)) abs(harmsNRP(z,3))];
    errhigh= [errMatRPabs(z,2) errMatNRPabs(z,2) errMatRPabs(z,3) errMatNRPabs(z,3)]; 
    errlow= [errMatRPabs(z,2) errMatNRPabs(z,2) errMatRPabs(z,3) errMatNRPabs(z,3)];
    bar([1 2 3 4],[abs(harmsRP(z,2)) abs(harmsNRP(z,2)) abs(harmsRP(z,3)) abs(harmsNRP(z,3))]) ;
    hold on;
    
    er = errorbar([1 2 3 4],data,errlow,errhigh);    
    er.Color = [0 0 0]; 
    er.LineWidth= 2;
    er.LineStyle = 'none';
    
    upLim=1.1*max(data); 
    
    set(gca, 'XTickLabel',nonLinLabel);
    ylabel('Absolute Magnitude','FontSize',8);
    ylim([0 upLim]) 
    title(['Scan=' scanIDList(z) '  Duration=' num2str(scanLengths(z)) 'sec'],'FontSize',8); 
%     axis tight;    
    drawnow();
    sgtitle(['Nonlinear Harmonic Magnitude Comparisons - t21=' num2str(t21ax(t21step)) ' t43=' num2str(t43ax(t43step))]);
    
    end 
    
    figure(figNonLinQuad);
    for m=1:numel(filesList)  
% subplot routine for the two ref combinations (one ref combos to follow) 

subplot(DimRow,DimCol,m);
    
%       plot(ticksFRET,FRETHist);
    data = [real(harmsRP(m,2)) real(harmsNRP(m,2)) imag(harmsRP(m,2)) imag(harmsNRP(m,2))];
    errhigh= [errMatRPx(m,2) errMatNRPx(m,2) errMatRPy(m,2) errMatNRPy(m,2)]; 
    errlow= [errMatRPx(m,2) errMatNRPx(m,2) errMatRPy(m,2) errMatNRPy(m,2)];
    bar([1 2 3 4],[real(harmsRP(m,2)) real(harmsNRP(m,2)) imag(harmsRP(m,2)) imag(harmsNRP(m,2))]) ;
    hold on;
    
    er = errorbar([1 2 3 4],data,errlow,errhigh);    
    er.Color = [0 0 0]; 
    er.LineWidth= 2;
    er.LineStyle = 'none';
    
    upLim=max(errhigh)+0.1*max(errhigh);
    lowLim=min(errlow)+0.1*min(errlow); 
    
    set(gca, 'XTickLabel',nonLinQuadLabel);
    ylabel('Absolute Magnitude','FontSize',8);
    ylim([-inf inf]) 
    title(['Scan=' scanIDList(m) '  Duration=' num2str(scanLengths(m)) 'sec'],'FontSize',8); 
%     axis tight;    
    drawnow();
    sgtitle(['Nonlinear Quadrature Comparisons - t21=' num2str(t21ax(t21step)) ' t43=' num2str(t43ax(t43step))]);
    
    end                     
                    
              
                        
           cd(startDir); 