% GOAL: load in single integration periods from 2D data sets to perform
% bootstrapping error analysis on. 

% UPDATE: 12/7/2020: this version of the program is NOT in parrallel. IT
% instead can handle a list of N photons and will calcualte all relevant
% quantities over that list of N photons with associated errors 
 

% *********SET PARAMETERS*******

% provide directory containing data files and a integration time in seconds
% (will loop over all individual scans in the containing folder)
fileFolder= 'C:\Users\Jack\Documents\Marcus Lab\2D_SingleMolecule\Sm_avgFolders\20200303_mol1';
startDir=pwd;
% provide the coordinates for the delay space integration period of interest
t21step=2;
t43step=1;
% specify the number of harmonics and number of bootstraps to perform. 
nH = 2;
numBoots=200;
allFigs=0;
% choose either a number of photons OR integration time (set the booleans accordingly)
int=0;
nPhot=1;

% nPhotMax=5e6; 
% create a vector of number of photons to be analyzed and specify the
% harmonic to be examined
curHarm=1; 
nPhotsMax=[1e3 5e3 1e4 3e4 5e4 7.5e4 1e5];
intTime=42; 
% ******************************

cd(fileFolder); 
filesList= dir('*2DFPGA'); 
numFiles=numel(filesList); 
curPath=pwd;

% for the photon list, a magntidue of the RP/NRP signal and its abs error
% is needed (all for the 1f right now)
RPlistX=zeros(numFiles,numel(nPhotsMax)); 
NRPlistX=zeros(numFiles,numel(nPhotsMax)); 
erRPlistX=zeros(numFiles,numel(nPhotsMax)); 
erNRPlistX=zeros(numFiles,numel(nPhotsMax)); 

RPlistY=zeros(numFiles,numel(nPhotsMax)); 
NRPlistY=zeros(numFiles,numel(nPhotsMax)); 
erRPlistY=zeros(numFiles,numel(nPhotsMax)); 
erNRPlistY=zeros(numFiles,numel(nPhotsMax)); 

RPlistAbs=zeros(numFiles,numel(nPhotsMax)); 
NRPlistAbs=zeros(numFiles,numel(nPhotsMax)); 
erRPlistAbs=zeros(numFiles,numel(nPhotsMax)); 
erNRPlistAbs=zeros(numFiles,numel(nPhotsMax));

xSiglistX=zeros(numFiles,numel(nPhotsMax)); 
erxSigListX=zeros(numFiles,numel(nPhotsMax)); 
xSiglistY=zeros(numFiles,numel(nPhotsMax)); 
erxSigListy=zeros(numFiles,numel(nPhotsMax)); 
 

if allFigs
figOneRef=figure();
figNonLin=figure();
figNonLinQuad=figure();
figTwoRef=figure();
end
figNPhotReal=figure(); 
figNPhotImag=figure();
figNPhotAbs=figure();


for k=1:numel(nPhotsMax)
curNPhot=nPhotsMax(k); 
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



% Generate a ScanID List prior to the parrallel loop to reference later on
for w=1:numel(filesList)
   
scanID = filesList(w).name(1:15);
scanIDList = [scanIDList; scanID];

end 


for f=1:1
%     numFiles
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
                
%                 if int
%                 time = fread(timeFile,Inf,'uint64=>uint64',0,'s');
%                 times = ((double (time-time(1))./(8e7))+1);
%                 times=times-times(1);
%                 p1 = fread(p1ID,Inf,'float64=>double',0,'s');
%                 p2 = fread(p2ID,Inf,'float64=>double',0,'s');
%                 [intWind]=find(times<=intTime);
% %                At this point isolate the phase and time lists from zero
% %                up to the integration time set
%                 scanLengths=[scanLengths;(times(intWind(end)))]; 
%                 p1=p1(1:intWind(end));
%                 p2=p2(1:intWind(end)); 
%                 numPhots = numPhots + numel(p1);
%                 p1r = (2.*pi).*p1;
%                 p2r = (2.*pi).*p2;
%                 RP= p2r - p1r;
%                 NRP= p2r + p1r;                
%                 fclose('all');
%                 
%                 end
                
%                 if nPhot && int
%                     disp('Select either a number of photons OR an integration time (not both) - Cancel execution now'); 
%                 end 
                
%                 if nPhot
%                     time = fread(timeID,curNPhot,'uint64=>uint64',0,'s');
%                     p1 = fread(p1ID,curNPhot,'float64=>double',0,'s');
%                     p2 = fread(p2ID,curNPhot,'float64=>double',0,'s');
%                     times = ((double (time-time(1))./(8e7))+1);
                time = fread(timeFile,Inf,'uint64=>uint64',0,'s');
                times = ((double (time-time(1))./(8e7))+1);
                times=times-times(1);
                p1 = fread(p1ID,Inf,'float64=>double',0,'s');
                p2 = fread(p2ID,Inf,'float64=>double',0,'s');
                
%                 store the total mean/phase of the signal
                    p1rTot = (2.*pi).*p1;
                    p2rTot = (2.*pi).*p2;
                    RPTot= p2rTot - p1rTot;
                    NRPTot= p2rTot + p1rTot; 
                    
                    RPphase=angle(mean(exp(-1i*RPTot)));
                    NRPphase=angle(mean(exp(-1i*NRPTot)));
                    
                times=times(1:curNPhot);
                p1=p1(1:curNPhot);
                p2=p2(1:curNPhot);
                    scanLength= times(end);
                    scanLengths=[scanLengths ; scanLength];
                    numPhots = numPhots + numel(p1);
                    p1r = (2.*pi).*p1;
                    p2r = (2.*pi).*p2;
                    RP= p2r - p1r;
                    NRP= p2r + p1r;                
                    fclose('all');
%                 end
                
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
BsampleRP = zeros(length(RP), nH+1);
BsampleNRP = zeros(length(NRP), nH+1);
BsampleRPy = zeros(length(RP), nH+1);
BsampleNRPy = zeros(length(NRP), nH+1);
BsampleX = zeros(length(p1r), nH+1);
BsampleY = zeros(length(p2r), nH+1);

  harmsXtemp= zeros(1, nH+1);
  harmsYtemp= zeros(1, nH+1);
  harmsRPtemp= zeros(1, nH+1);
  harmsNRPtemp= zeros(1, nH+1);


for j=0:nH
   
    
    BsampleRP(:,j+1)=(exp(-1i*(j)*RP));
    BsampleNRP(:,j+1)=(exp(-1i*(j)*NRP));
    BsampleRPx(:,j+1)=real(exp(-1i*(j)*RP));
    BsampleNRPx(:,j+1)=real(exp(-1i*(j)*NRP));
    BsampleRPy(:,j+1)=imag(exp(-1i*(j)*RP));
    BsampleNRPy(:,j+1)=imag(exp(-1i*(j)*NRP));
    BsampleX(:,j+1) = exp(-1i*j*p1r);
    BsampleY(:,j+1) = exp(-1i*j*p2r); 
    
    harmsXtemp(1, j+1) = mean(exp(-1i*j*p1r));
    harmsYtemp(1, j+1) = mean(exp(-1i*j*p2r));
    harmsRPtemp(1, j+1) = mean(exp(-1i*j*RP));
    harmsNRPtemp(1, j+1) = mean(exp(-1i*j*NRP));
       
    
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
bfuncRPabs = @(BsampleRPabs) (mean(BsampleRP));
bootstatRPabs = bootstrp(numBoots,bfuncRPabs,BsampleRP);
absRP = abs(harmsRP(f,:));
errRPabs = std(bootstatRPabs);
perErrRPabs = errRPabs(2:end)./absRP(2:end)*100;

% Bootstrap for NRP abs
bfuncNRPabs = @(BsampleNRPabs) (mean(BsampleNRP));
bootstatNRPabs = bootstrp(numBoots,bfuncNRPabs,BsampleNRP);
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

% bootstrap using the Tiemo method to be able to examine the Abs (or real
% phased to have no imag component)
[xerrRP, yerrRP, abserrRP, phierrRP, rawbootstrapRP, fftbootstrapRP] = bootstrapFourierStatistics_phased(numBoots, nH, RP, RPphase);
[xerrNRP, yerrNRP, abserrNRP, phierrNRP, rawbootstrapNRP, fftbootstrapNRP] = bootstrapFourierStatistics_phased(numBoots, nH, NRP, NRPphase);



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

    
    RPlistX(f,k)= real(harmsRP(f,curHarm+1));
    NRPlistX(f,k)= real(harmsNRP(f,curHarm+1));
    erRPlistX(f,k)=xerrRP(curHarm+1); 
    erNRPlistX(f,k)= xerrNRP(curHarm+1); 

    RPlistY(f,k)= imag(harmsRP(f,curHarm+1));
    NRPlistY(f,k)= imag(harmsNRP(f,curHarm+1));
    erRPlistY(f,k)=yerrRP(curHarm+1); 
    erNRPlistY(f,k)= yerrNRP(curHarm+1);
    
    RPlistAbs(f,k)= abs(harmsRP(f,curHarm+1));
    NRPlistAbs(f,k)= abs(harmsNRP(f,curHarm+1));
    erRPlistAbs(f,k)=abserrRP(curHarm+1); 
    erNRPlistAbs(f,k)= abserrNRP(curHarm+1);
    
%    xSiglistY(f,k)= imag(harmsXref(f,curHarm+1));
%    xSiglistX(f,k)= real(harmsXref(f,curHarm+1));
%    erxSiglistY(f,k)=errRPy(curHarm+1); 
%    erxSiglistY(f,k)= errNRPy(curHarm+1); 
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
    
    if allFigs
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
    title(['Scan=' scanIDList(m) '  Duration=' num2str(scanLengths(m)) 'min'],'FontSize',8); 
%     axis tight;    
    drawnow();
    sgtitle(['Nonlinear Quadrature Comparisons - t21=' num2str(t21ax(t21step)) ' t43=' num2str(t43ax(t43step))]);
    
    end
    
    end
    
end 

figure(figNPhotReal)
% fig1 = figure(1);
% clf(1)
fig1pos = get(figNPhotReal, 'Position');
set(figNPhotReal, 'PaperUnits','centimeters')
set(figNPhotReal, 'PaperSize', [18, 18/fig1pos(3)*fig1pos(4)])
% num of files
for z=1:1
%     size(RPlistX,1)
%     num of elements in photon list
    for n=1:size(RPlistX,2)
hold on
scatter(sqrt(nPhotsMax(n)),(abs(real(mean(exp(-1i*RPTot))*exp(1i*RPphase)))/abs(erRPlistX(z,n)))/abs(mean(exp(-1i*RPTot))),'Filled','b');
scatter(sqrt(nPhotsMax(n)),(abs(real(mean(exp(-1i*NRPTot))*exp(1i*NRPphase)))/abs(erNRPlistX(z,n)))/abs(mean(exp(-1i*NRPTot))),'Filled','r');

% erRP = errorbar(nPhotsMax(n),RPlistX(z,n),erRPlistX(z,n),erRPlistX(z,n)); 
% erNRP = errorbar(nPhotsMax(n),RPlistY(z,n),erRPlistY(z,n),erRPlistY(z,n)); 
% 
% erRP.Color = [0 0 0]; 
% erRP.LineWidth= 2;
% erRP.LineStyle = 'none';
% 
% erNRP.Color = [0 0 0]; 
% erNRP.LineWidth= 2;
% erNRP.LineStyle = 'none';



    end 
    line(sqrt(nPhotsMax),sqrt(nPhotsMax/2)); 
end 

title('Sideband Real Part Scaling with N Photons');
xlabel('Number of Photons');
ylabel('Magnitude (arb)');
% xlim([0.0 1.1]); 
% ylim([0 0.3]); 
grid on
legend('RP Re', 'NRP Re');
logx;
logy;
% if printOpt
% print -dpdf -r300 -bestfit LinPreAtt.pdf
% end 


figure(figNPhotImag)
% fig1 = figure(1);
% clf(1)
fig2pos = get(figNPhotImag, 'Position');
set(figNPhotImag, 'PaperUnits','centimeters')
set(figNPhotImag, 'PaperSize', [18, 18/fig2pos(3)*fig2pos(4)])
% num of files
for z=1:1
%     size(NRPlistX,1)
%     num of elements in photon list
    for n=1:size(NRPlistX,2)
hold on
% if you phase these to be purely imag, you dont need to take the imag part
% any longer
scatter(sqrt(nPhotsMax(n)),(abs((mean(exp(-1i*RPTot))*exp(1i*(pi-RPphase))))/abs(erRPlistY(z,n)))/abs(mean(exp(-1i*RPTot))),'Filled','b');
scatter(sqrt(nPhotsMax(n)),(abs((mean(exp(-1i*NRPTot))*exp(1i*(pi-NRPphase))))/abs(erNRPlistY(z,n)))/abs(mean(exp(-1i*NRPTot))),'Filled','r');


% erRP = errorbar(sqrt(nPhotsMax(n)),abs(NRPlistX(z,n)),abs(erNRPlistX(z,n)),abs(erNRPlistX(z,n))); 
% erNRP = errorbar(sqrt(nPhotsMax(n)),abs(NRPlistY(z,n)),abs(erNRPlistY(z,n)),abs(erNRPlistY(z,n))); 

% erRP.Color = [0 0 0]; 
% erRP.LineWidth= 2;
% erRP.LineStyle = 'none';
% 
% erNRP.Color = [0 0 0]; 
% erNRP.LineWidth= 2;
% erNRP.LineStyle = 'none';

    end 
    line(sqrt(nPhotsMax),sqrt(nPhotsMax/2)); 
    ylim([1 1000]);
end 

title('Sideband Imag Part Scaling with N Photons');
xlabel('Number of Photons');
ylabel('Magnitude (arb)');
% xlim([0.0 1.1]); 
% ylim([0 0.3]); 
grid on
legend('RP Im', 'NRP Im'); 
logx;
logy;                   


figure(figNPhotAbs)
% fig1 = figure(1);
% clf(1)
fig3pos = get(figNPhotAbs, 'Position');
set(figNPhotAbs, 'PaperUnits','centimeters')
set(figNPhotAbs, 'PaperSize', [18, 18/fig3pos(3)*fig3pos(4)])
% num of files
for z=1:1
%     size(NRPlistX,1)
%     num of elements in photon list
    for n=1:size(RPlistAbs,2)
hold on
scatter(sqrt(nPhotsMax(n)),(abs(mean(exp(-1i*RPTot)))/abs(erRPlistAbs(z,n)))/abs(mean(exp(-1i*RPTot))),'Filled','b');
scatter(sqrt(nPhotsMax(n)),(abs(mean(exp(-1i*NRPTot)))/abs(erNRPlistAbs(z,n)))/abs(mean(exp(-1i*NRPTot))),'Filled','r');


% erRP = errorbar(sqrt(nPhotsMax(n)),abs(NRPlistX(z,n)),abs(erNRPlistX(z,n)),abs(erNRPlistX(z,n))); 
% erNRP = errorbar(sqrt(nPhotsMax(n)),abs(NRPlistY(z,n)),abs(erNRPlistY(z,n)),abs(erNRPlistY(z,n))); 

% erRP.Color = [0 0 0]; 
% erRP.LineWidth= 2;
% erRP.LineStyle = 'none';
% 
% erNRP.Color = [0 0 0]; 
% erNRP.LineWidth= 2;
% erNRP.LineStyle = 'none';

    end 
    line(sqrt(nPhotsMax),sqrt(nPhotsMax/2)); 
    ylim([1 1000]);
    
end 

title('Sideband Abs Scaling with N Photons');
xlabel('Number of Photons');
ylabel('Magnitude (arb)');
% xlim([0.0 1.1]); 
% ylim([0 0.3]); 
grid on
legend('RP Abs', 'NRP Abs'); 
logx;
logy;                   
              
                        
           cd(startDir); 