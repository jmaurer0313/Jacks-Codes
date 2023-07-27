% === file structure of Ek.dat and Speck.dat ===
% x=1\t Intensity\t Phase\t Real()\t Imaginary()\t\n
% :
% :
% :
% x=N\t Intensity\t Phase\t Real()\t Imaginary()\t\n
% === EOF ===
% where x:= time in Ek.dat and x:= lambda in Speck.dat
clear all 
close all
% ****************PROPER FILE NAME REQUIRED FOR ANALYSIS ******************
% mySpeckMat = dlmread('FROG_TRACE_05052022_trace3_realigned_wGlass.bin.Speck.dat');
% myEkMat = dlmread('FROG_TRACE_05052022_trace3_realigned_wGlass.bin.Ek.dat'); 

% mySpeckMat = dlmread('FROGtrace_04112022_NOPA_w1inchGlass_PrismComp_trace1.bin.Speck.dat');
% myEkMat = dlmread('FROGtrace_04112022_NOPA_w1inchGlass_PrismComp_trace1.bin.Ek.dat');

%% Claire Data
%--------------------------------------------------------------------------
% without glass
%--------------------------------------------------------------------------
% mySpeckMat = dlmread('20220517_FROGtrace_Pharos_noGlass_185370_v1_4p5e4iter.bin.Speck.dat');
% myEkMat = dlmread('20220517_FROGtrace_Pharos_noGlass_185370_v1_4p5e4iter.bin.Ek.dat');

% mySpeckMat = dlmread('20220517_FROGtrace_Pharos_noGlass_185370_v2_3p6e4iter.bin.Speck.dat');
% myEkMat = dlmread('20220517_FROGtrace_Pharos_noGlass_185370_v2_3p6e4iter.bin.Ek.dat');

% need updated v3!! ... don't use this one
% mySpeckMat = dlmread('20220517_FROGtrace_Pharos_noGlass_185370_v3_2e4iter.bin.Speck.dat');
% myEkMat = dlmread('20220517_FROGtrace_Pharos_noGlass_185370_v3_2e4iter.bin.Ek.dat');

% mySpeckMat = dlmread('20220517_FROGtrace_Pharos_noGlass_185370_v4_6e4iter.bin.Speck.dat');
% myEkMat = dlmread('20220517_FROGtrace_Pharos_noGlass_185370_v4_6e4iter.bin.Ek.dat');


%--------------------------------------------------------------------------
% with glass
%--------------------------------------------------------------------------
% mySpeckMat = dlmread('20220517_FROGtrace_Pharos_1inGlass_185370_v1_3p5e4iter.bin.Speck.dat');
% myEkMat = dlmread('20220517_FROGtrace_Pharos_1inGlass_185370_v1_3p5e4iter.bin.Ek.dat');

% mySpeckMat = dlmread('20220517_FROGtrace_Pharos_1inGlass_185370_v2_3p5e4iter.bin.Speck.dat');
% myEkMat = dlmread('20220517_FROGtrace_Pharos_1inGlass_185370_v2_3p5e4iter.bin.Ek.dat');

% mySpeckMat = dlmread('20220517_FROGtrace_Pharos_1inGlass_185370_v3-4e4iter.bin.Speck.dat');
% myEkMat = dlmread('20220517_FROGtrace_Pharos_1inGlass_185370_v3-4e4iter.bin.Ek.dat');
% 
% mySpeckMat = dlmread('20220517_FROGtrace_Pharos_1inGlass_185370_v4_3e4iter.bin.Speck.dat');
% myEkMat = dlmread('20220517_FROGtrace_Pharos_1inGlass_185370_v4_3e4iter.bin.Ek.dat');

% mySpeckMat = dlmread('20220517_FROGtrace_Pharos_1inGlass_185370_v5_2e5iter.bin.Speck.dat');
% myEkMat = dlmread('20220517_FROGtrace_Pharos_1inGlass_185370_v5_2e5iter.bin.Ek.dat');

%% Claire Data compiled
%--------------------------------------------------------------------------
% without glass
%--------------------------------------------------------------------------
% cd('/Users/clairealbrecht/Dropbox/Claire_Dropbox/FROG/FROGtraces_Claire/20220517/Analyzed_noGlass')
% fileNamesPrefix = {'20220517_FROGtrace_Pharos_noGlass_185370_v1_4p5e4iter';...
%     '20220517_FROGtrace_Pharos_noGlass_185370_v2_3p6e4iter';...
%     '20220517_FROGtrace_Pharos_noGlass_185370_v4_6e4iter'};

%--------------------------------------------------------------------------
% with glass
%--------------------------------------------------------------------------    
cd('C:\Users\Jack\Documents\Marcus Lab\2DFS\FROG\FROGtraces\NOPA_wGlassSet')
fileNamesPrefix = {'FROG_TRACE_05052022_trace3_realigned_wGlass';...
    'FROG_TRACE_05052022_trace4_realigned_wGlass';...
    'FROGtrace_NOPA_w1inchGlass_noCompress_04182022_trace1_80pPow';...
    'FROGtrace_NOPA_w1inchGlass_noCompress_04182022_trace2_80pPow'};

%%
% ************************************************************************

% flipTime=1;
% if flipTime
%  myEkMat(:,2)= fliplr(myEkMat(:,2));   
%  myEkMat(:,3)= flipud(myEkMat(:,3));   
%  
%  mySpeckMat(:,3) = flipud(mySpeckMat(:,3));
% end 

% *********SET CUTOFFS AND CENTER WAVELENGTH ********8
% cutoff wavelengths in nm
lambdaLB = 530;
lambdaUB = 580; 
centerLambda = 555;

% lambdaLB = 1000;
% lambdaUB = 1035;
% centerLambda = 1020;
% lambdaLB = 1004;
% lambdaUB = 1035;
% centerLambda = 1018;

applyHigherTerms=1;
splineLength=1000;
% **************************************************

freqUB= 3e8/(lambdaLB*1e-9);
freqLB= 3e8/(lambdaUB*1e-9);
centerFreq = 3e8/(centerLambda*1e-9); 


%% Loop over the files with the same condition

mag2ndOrder = zeros(length(fileNamesPrefix),4);
mag3rdOrder = zeros(length(fileNamesPrefix),4);

Avg2ndOrderCurve_3rdDeg = zeros(1,splineLength-2);
Avg3rdOrderCurve_3rdDeg = zeros(1,splineLength-3);

Avg2ndOrderCurve_4thDeg = zeros(1,splineLength-2);
Avg3rdOrderCurve_4thDeg = zeros(1,splineLength-3);

Avg2ndOrderCurve_5thDeg = zeros(1,splineLength-2);
Avg3rdOrderCurve_5thDeg = zeros(1,splineLength-3);

Avg2ndOrderCurve_6thDeg = zeros(1,splineLength-2);
Avg3rdOrderCurve_6thDeg = zeros(1,splineLength-3);

for i = 1:length(fileNamesPrefix)
    mySpeckMat = dlmread([fileNamesPrefix{i} '.bin.Speck.dat']);
    myEkMat = dlmread([fileNamesPrefix{i} '.bin.Ek.dat']);
    
    
% first reconstruct the plots from Frogger() as a soft check on errors made
% in the data files
% time domain plot
figure(1)
plot(myEkMat(:,1),(myEkMat(:,2)*(max(abs(myEkMat(:,3)))/max(myEkMat(:,2)))) )
hold on
plot(myEkMat(:,1),(myEkMat(:,3)+abs(min(myEkMat(:,3)))) )
legend('Intensity', 'Phase'); 
title('Time domain representation of the Field'); 

% Freq domain plot
figure(2)
plot(mySpeckMat(:,1),mySpeckMat(:,2)*(max(abs(mySpeckMat(:,3)))/max(mySpeckMat(:,2))))
hold on
plot(mySpeckMat(:,1),(mySpeckMat(:,3)+abs(min(mySpeckMat(:,3)))) )
legend('Intensity', 'Phase'); 
title('Frequency domain representation of the Field'); 

%%
% next, the relevant region of the spectral phase should be isolated, so
% that only the region within the bandwidth of the laser is considered -
% since the spectral phase outside that region is ill conditioned and
% unphysical - there seems to be a constant offset in the reported
% frequency compared with the known bandwidth center, could consider
% adjusting but will keep for now since only the magnitude and slope of the
% spectral phase matters 

% convert spectral axis to frequency 
specAxisTotal = 3e8./(mySpeckMat(:,1)*1e-9);
specPhaseTotal = mySpeckMat(:,3);

% fit a gaussian to the bandwidth to establish the center frequency
[fitobj, goodness, output, convmsg]= fit(specAxisTotal,mySpeckMat(:,2),'gauss1');
centerLambda=fitobj.b1;

% indices are reversed due to the wavelength to freq conversion
firstIndList = find(specAxisTotal>=freqLB);
secondIndList = find(specAxisTotal<=freqUB);

% specCutAxis = specAxisTotal(firstIndList(1):secondIndList(end)); 
% specPhaseCut= specPhaseTotal(firstIndList(1):secondIndList(end));
specCutAxis = specAxisTotal(secondIndList(1):firstIndList(end)); 
specPhaseCut= specPhaseTotal(secondIndList(1):firstIndList(end));

specCutAxisNoCarrier = (2*pi*specCutAxis) -(2*pi*centerFreq); 

splineAxis=linspace(specCutAxis(1),specCutAxis(end),splineLength); 
specSpline= spline(specCutAxis, specPhaseCut,splineAxis);
lambdaSplineAxis= (3e8./(splineAxis*1e-9));


splineAxisNoCarrier=linspace(specCutAxisNoCarrier(1),specCutAxisNoCarrier(end),splineLength); 
specSplineNoCarrier= spline(specCutAxisNoCarrier, specPhaseCut,splineAxisNoCarrier);  


%% CA: adding the derivatives of raw data plot back in from v1

% now with the relevant region isolated from the spectral phase, take the
% 1st, 2nd and 3rd derivatives to isolate the relevant dispersion terms
diffArr= splineAxis-circshift(splineAxis,1,2);
h = mode(diffArr(2:end));       % step size
X = splineAxis;                 % domain
f = specSpline;                 % range
firstOrd = diff(f)/h;           % first derivative
secondOrd = diff(firstOrd)/h;   % second derivative
thirdOrd = diff(secondOrd)/h;   % third derivative

% % now make a plot of the successive derivatives over the domain
% figure(33)
% clf
% % figure(5)
% plot(splineAxisNoCarrier(1:length(firstOrd))+(2*pi*centerFreq), firstOrd*(1e15), ...
%     splineAxisNoCarrier(1:length(secondOrd))+(2*pi*centerFreq), (1/2)*secondOrd*(1e15^2), ...
%     splineAxisNoCarrier(1:length(thirdOrd))+(2*pi*centerFreq), (1/6)*thirdOrd*(1e15^3))
% hold on
% plot(splineAxisNoCarrier+(2*pi*centerFreq),zeros(1,length(splineAxisNoCarrier)),'--')
% legend('First Order' ,'Second Order' ,'Third Order');
% 
% 
% % plot(splineAxis(1:length(firstOrd)), firstOrd, ...
% %     splineAxis(1:length(secondOrd)), secondOrd, ...
% %     splineAxis(1:length(thirdOrd)), thirdOrd)
% % hold on
% % plot(splineAxis,zeros(1,length(splineAxis)),'--')
% % legend('First Order' ,'Second Order' ,'Third Order');
% title('Spectral Phase Derivatives - Disperion Orders'); 




%%

% apply a 6th-3rd degree polynomial fit to the no carrier axis data
    [poly6, poly6Struc] = polyfit(splineAxisNoCarrier,specSplineNoCarrier, 6);
    coeff6 = polyfit(splineAxisNoCarrier,specSplineNoCarrier, 6);
    y6 = polyval(poly6, splineAxisNoCarrier);

    [poly5, poly5Struc] = polyfit(splineAxisNoCarrier,specSplineNoCarrier, 5);
    coeff5 = polyfit(splineAxisNoCarrier,specSplineNoCarrier, 5);
    y5 = polyval(poly5, splineAxisNoCarrier);
    
    [poly4, poly4Struc] = polyfit(splineAxisNoCarrier,specSplineNoCarrier, 4);
    coeff4 = polyfit(splineAxisNoCarrier,specSplineNoCarrier, 4);
    y4 = polyval(poly4, splineAxisNoCarrier);

    [poly3, poly3Struc] = polyfit(splineAxisNoCarrier,specSplineNoCarrier, 3);
    coeff3 = polyfit(splineAxisNoCarrier,specSplineNoCarrier, 3);
    y3 = polyval(poly3, splineAxisNoCarrier);
    
%     plot an overlay of the raw data with te various degrees of polynomial
%     fit
figure(3)
plot(specCutAxisNoCarrier,specPhaseCut,'LineWidth',3)
hold on
plot(splineAxisNoCarrier,polyval(poly3, splineAxisNoCarrier),'--','LineWidth',1.5)
plot(splineAxisNoCarrier,polyval(poly4, splineAxisNoCarrier),'--','LineWidth',1.5)
plot(splineAxisNoCarrier,polyval(poly5, splineAxisNoCarrier),'--','LineWidth',1.5)
plot(splineAxisNoCarrier,polyval(poly6, splineAxisNoCarrier),'--','LineWidth',1.5)
legend('Raw Data', 'Poly3 Fit','Poly4 Fit','Poly5 Fit','Poly6 Fit');
title('Overlay of all Polynomial Fits'); 
xlabel('Frequency (w-w0)');
ylabel('Magnitude (arb)'); 

%     collect the no carrier (NC) terms based on degree of fit and order of
%     dispersion 
    poly6_2ndOrdNC = poly6(5)*(1e15^2);
    poly6_3rdOrdNC = poly6(4)*(1e15^3);
    poly6_4thOrdNC = poly6(3)*(1e15^4);
    poly6_5thOrdNC = poly6(2)*(1e15^5);
    poly6_6thOrdNC = poly6(1)*(1e15^6);
    
    poly5_2ndOrdNC = poly5(4)*(1e15^2);
    poly5_3rdOrdNC = poly5(3)*(1e15^3);
    poly5_4thOrdNC = poly5(2)*(1e15^4);
    poly5_5thOrdNC = poly5(1)*(1e15^5);
    
    poly4_2ndOrdNC = poly4(3)*(1e15^2);
    poly4_3rdOrdNC = poly4(2)*(1e15^3);
    poly4_4thOrdNC = poly4(1)*(1e15^4);
    
    poly3_2ndOrdNC = poly3(2)*(1e15^2);
    poly3_3rdOrdNC = poly3(1)*(1e15^3);
    
%     units = [1e15^2, 1e15^3, 1e15^4, 1e15^5, 1e15^6];
%     poly6_wUnits = flip(poly6(1:5)) .* units(1:5);
%     poly5_wUnits = flip(poly5(1:4)) .* units(1:4);
%     poly4_wUnits = flip(poly4(1:3)) .* units(1:3);
%     poly3_wUnits = flip(poly3(1:2)) .* units(1:2);
    
    mag2ndOrder(i,:) = [poly3_2ndOrdNC, poly4_2ndOrdNC,poly5_2ndOrdNC,poly6_2ndOrdNC];%,poly6_2ndOrdNC];
    mag3rdOrder(i,:) = [poly3_3rdOrdNC, poly4_3rdOrdNC,poly5_3rdOrdNC,poly6_3rdOrdNC];%,poly6_3rdOrdNC];

    figure(4)
    plot(poly6_2ndOrdNC,poly6_3rdOrdNC,'o','LineWidth',2);
    hold on
    plot(poly5_2ndOrdNC,poly5_3rdOrdNC,'o','LineWidth',2);
    plot(poly4_2ndOrdNC,poly4_3rdOrdNC,'o','LineWidth',2);
    plot(poly3_2ndOrdNC,poly3_3rdOrdNC,'o','LineWidth',2);
    legend('6th deg fit', '5th deg fit','4th deg fit','3rd deg fit'); 
    xlabel('2nd Order Magnitude (fs^2)');
    ylabel('3rd Order Magnitude (fs^3)');
    title('GDD and TOD Coefficents Obtained from Nth deg polyFit');

diffArr= splineAxisNoCarrier-circshift(splineAxisNoCarrier,1,2);
h = mode(diffArr(2:end));       % step size
X = splineAxisNoCarrier;    % domain
f = y3;      % range
firstOrd = diff(f)/h;   % first derivative
secondOrd = diff(firstOrd)/h;   % second derivative
thirdOrd = diff(secondOrd)/h;   % second derivative

Avg2ndOrderCurve_3rdDeg = Avg2ndOrderCurve_3rdDeg + abs(secondOrd);
Avg3rdOrderCurve_3rdDeg = Avg3rdOrderCurve_3rdDeg + abs(thirdOrd);

% now make a plot of the successive derivatives over the domain using the
% 3rd order polynomial fit
figure(5)
plot(splineAxisNoCarrier(1:length(firstOrd))+(2*pi*centerFreq), firstOrd*(1e15), ...
    splineAxisNoCarrier(1:length(secondOrd))+(2*pi*centerFreq), (1/2)*secondOrd*(1e15^2), ...
    splineAxisNoCarrier(1:length(thirdOrd))+(2*pi*centerFreq), (1/6)*thirdOrd*(1e15^3))
hold on
plot(splineAxisNoCarrier+(2*pi*centerFreq),zeros(1,length(splineAxisNoCarrier)),'--')
legend('First Order' ,'Second Order' ,'Third Order');
title('Spectral Phase Derivatives - 3rd deg PolyFit - Disperion Orders'); 
ylabel('Magnitude (fs, fs^2, fs^3)');
xlabel('Frequency (Hz)') ;

if applyHigherTerms
%%
% take the derivatives of the 4th order poly fit with proper magnitudes 
diffArr= splineAxisNoCarrier-circshift(splineAxisNoCarrier,1,2);
h = mode(diffArr(2:end));       % step size
X = splineAxisNoCarrier;    % domain
f = polyval(poly4(2:end),splineAxisNoCarrier);      % range
% for higher order polynomials, isolation of the proper GD, GDD, TOD terms
% requires more derivatives to eliminate the higher ordered terms brought in
% by poly
firstOrd = diff(f)/h;   % first derivative
secondOrd = diff(firstOrd)/h;   % second derivative
thirdOrd = diff(secondOrd)/h;   % second derivative

Avg2ndOrderCurve_4thDeg = Avg2ndOrderCurve_4thDeg + abs(secondOrd);
Avg3rdOrderCurve_4thDeg = Avg3rdOrderCurve_4thDeg + abs(thirdOrd);

% now make a plot of the successive derivatives over the domain using the
% 3rd order polynomial fit
figure(6)
plot(splineAxisNoCarrier(1:length(firstOrd))+(2*pi*centerFreq), firstOrd*(1e15), ...
    splineAxisNoCarrier(1:length(secondOrd))+(2*pi*centerFreq), (1/2)*secondOrd*(1e15^2), ...
    splineAxisNoCarrier(1:length(thirdOrd))+(2*pi*centerFreq), (1/6)*thirdOrd*(1e15^3))
hold on
plot(splineAxisNoCarrier+(2*pi*centerFreq),zeros(1,length(splineAxisNoCarrier)),'--')
legend('First Order' ,'Second Order' ,'Third Order');
title('Spectral Phase Derivatives -4th deg PolyFit - Disperion Orders');
ylabel('Magnitude (fs, fs^2, fs^3)');
xlabel('Frequency (Hz)') ;

% take the derivatives of the 5th order poly fit with proper magnitudes 
diffArr= splineAxisNoCarrier-circshift(splineAxisNoCarrier,1,2);
h = mode(diffArr(2:end));       % step size
X = splineAxisNoCarrier;    % domain
f = polyval(poly5(3:end),splineAxisNoCarrier);      % range
% for higher order polynomials, isolation of the proper GD, GDD, TOD terms
% requires more derivatives to eliminate the higher ordered terms brought in
% by poly
firstOrd = diff(f)/h;   % first derivative
secondOrd = diff(firstOrd)/h;   % second derivative
thirdOrd = diff(secondOrd)/h;   % second derivative

Avg2ndOrderCurve_5thDeg = Avg2ndOrderCurve_5thDeg + abs(secondOrd);
Avg3rdOrderCurve_5thDeg = Avg3rdOrderCurve_5thDeg + abs(thirdOrd);

% now make a plot of the successive derivatives over the domain using the
% 3rd order polynomial fit
figure(7)
plot(splineAxisNoCarrier(1:length(firstOrd))+(2*pi*centerFreq), firstOrd*(1e15), ...
    splineAxisNoCarrier(1:length(secondOrd))+(2*pi*centerFreq), (1/2)*secondOrd*(1e15^2), ...
    splineAxisNoCarrier(1:length(thirdOrd))+(2*pi*centerFreq), (1/6)*thirdOrd*(1e15^3))
hold on
plot(splineAxisNoCarrier+(2*pi*centerFreq),zeros(1,length(splineAxisNoCarrier)),'--')
legend('First Order' ,'Second Order' ,'Third Order');
title('Spectral Phase Derivatives -5th deg PolyFit - Disperion Orders');
ylabel('Magnitude (fs, fs^2, fs^3)');
xlabel('Frequency (Hz)') ;


% take the derivatives of the 6th order poly fit with proper magnitudes 
diffArr= splineAxisNoCarrier-circshift(splineAxisNoCarrier,1,2);
h = mode(diffArr(2:end));       % step size
X = splineAxisNoCarrier;    % domain
f = polyval(poly6(4:end),splineAxisNoCarrier);     % range
% for higher order polynomials, isolation of the proper GD, GDD, TOD terms
% requires more derivatives to eliminate the higher ordered terms brought in
% by poly
firstOrd = diff(f)/h;   % first derivative
secondOrd = diff(firstOrd)/h;   % second derivative
thirdOrd = diff(secondOrd)/h;   % second derivative

Avg2ndOrderCurve_6thDeg = Avg2ndOrderCurve_6thDeg + abs(secondOrd);
Avg3rdOrderCurve_6thDeg = Avg3rdOrderCurve_6thDeg + abs(thirdOrd);

% now make a plot of the successive derivatives over the domain using the
% 3rd order polynomial fit
figure(8)
plot(splineAxisNoCarrier(1:length(firstOrd))+(2*pi*centerFreq), firstOrd*(1e15), ...
    splineAxisNoCarrier(1:length(secondOrd))+(2*pi*centerFreq), (1/2)*secondOrd*(1e15^2), ...
    splineAxisNoCarrier(1:length(thirdOrd))+(2*pi*centerFreq), (1/6)*thirdOrd*(1e15^3))
hold on
plot(splineAxisNoCarrier+(2*pi*centerFreq),zeros(1,length(splineAxisNoCarrier)),'--')
legend('First Order' ,'Second Order' ,'Third Order');
title('Spectral Phase Derivatives -6th deg PolyFit - Disperion Orders');
ylabel('Magnitude (fs, fs^2, fs^3)');
xlabel('Frequency (Hz)') ;
% 
end



end

Avg2ndOrderCurve_3rdDeg = Avg2ndOrderCurve_3rdDeg./length(fileNamesPrefix);
Avg3rdOrderCurve_3rdDeg = Avg3rdOrderCurve_3rdDeg./length(fileNamesPrefix);

Avg2ndOrderCurve_4thDeg = Avg2ndOrderCurve_4thDeg./length(fileNamesPrefix);
Avg3rdOrderCurve_4thDeg = Avg3rdOrderCurve_4thDeg./length(fileNamesPrefix);

Avg2ndOrderCurve_5thDeg = Avg2ndOrderCurve_5thDeg./length(fileNamesPrefix);
Avg3rdOrderCurve_5thDeg = Avg3rdOrderCurve_5thDeg./length(fileNamesPrefix);

Avg2ndOrderCurve_6thDeg = Avg2ndOrderCurve_6thDeg./length(fileNamesPrefix);
Avg3rdOrderCurve_6thDeg = Avg3rdOrderCurve_6thDeg./length(fileNamesPrefix);


%%

polyfitorders = [3,4,5,6]; %[2,3,4,5,6];
figure(9)
clf
lgd = [];
for i = 1:length(fileNamesPrefix)
scatter(polyfitorders, abs(mag2ndOrder(i,:)))
lgd_temp = ['scan: ',num2str(i)];
lgd = [lgd; lgd_temp];
hold on
end
scatter(polyfitorders, mean(abs(mag2ndOrder)),'k','marker','x','linewidth',2)
text(polyfitorders(1)*(1+0.01), mean(abs(mag2ndOrder(:,1))),['{\mu} = ' num2str(round(mean(abs(mag2ndOrder(:,1))),2))])
text(polyfitorders(2)*(1+0.01), mean(abs(mag2ndOrder(:,2))),['{\mu} = ' num2str(round(mean(abs(mag2ndOrder(:,2))),2))])
text(polyfitorders(3)*(1+0.01), mean(abs(mag2ndOrder(:,3))),['{\mu} = ' num2str(round(mean(abs(mag2ndOrder(:,3))),2))])
text(polyfitorders(4)*(1+0.01), mean(abs(mag2ndOrder(:,4))),['{\mu} = ' num2str(round(mean(abs(mag2ndOrder(:,4))),2))])
xticks(polyfitorders)
legend([lgd; 'mean   '])
set(gcf,'color','w')
xlabel('Degree of polynomial fit')
ylabel('Magnitude of second order dispersion')
title('Comparing GDD across multiple scans')

figure(10)
clf
lgd = [];
for i = 1:length(fileNamesPrefix)
scatter(polyfitorders, abs(mag3rdOrder(i,:)))
lgd_temp = ['scan: ',num2str(i)];
lgd = [lgd; lgd_temp];
hold on
end
scatter(polyfitorders, mean(abs(mag3rdOrder)),'k','marker','x','linewidth',2)
text(polyfitorders(1)*(1+0.01), mean(abs(mag3rdOrder(:,1))),['{\mu} = ' num2str(round(mean(abs(mag3rdOrder(:,1))),2))])
text(polyfitorders(2)*(1+0.01), mean(abs(mag3rdOrder(:,2))),['{\mu} = ' num2str(round(mean(abs(mag3rdOrder(:,2))),2))])
text(polyfitorders(3)*(1+0.01), mean(abs(mag3rdOrder(:,3))),['{\mu} = ' num2str(round(mean(abs(mag3rdOrder(:,3))),2))])
text(polyfitorders(4)*(1+0.01), mean(abs(mag3rdOrder(:,4))),['{\mu} = ' num2str(round(mean(abs(mag3rdOrder(:,4))),2))])
xticks(polyfitorders)
legend([lgd; 'mean   '])
set(gcf,'color','w')
xlabel('Degree of polynomial fit')
ylabel('Magnitude of third order dispersion')
title('Comparing TOD across multiple scans')
 
figure(11)
plot(splineAxisNoCarrier(1:length(Avg2ndOrderCurve_3rdDeg))+(2*pi*centerFreq), (1/2)*Avg2ndOrderCurve_3rdDeg*(1e15^2), ...
    splineAxisNoCarrier(1:length(Avg2ndOrderCurve_4thDeg))+(2*pi*centerFreq), (1/2)*Avg2ndOrderCurve_4thDeg*(1e15^2), ...
    splineAxisNoCarrier(1:length(Avg2ndOrderCurve_5thDeg))+(2*pi*centerFreq), (1/2)*Avg2ndOrderCurve_5thDeg*(1e15^2), ...
    splineAxisNoCarrier(1:length(Avg2ndOrderCurve_6thDeg))+(2*pi*centerFreq), (1/2)*Avg2ndOrderCurve_6thDeg*(1e15^2));
hold on
plot(splineAxisNoCarrier+(2*pi*centerFreq),zeros(1,length(splineAxisNoCarrier)),'--')
legend('3rd deg' ,'4th deg', '5th deg','6th deg');
title('Average GDD Curves - various polynomial degrees');
ylabel('Magnitude (fs^2)');
xlabel('Frequency (Hz)') ;

figure(12)
plot(splineAxisNoCarrier(1:length(Avg3rdOrderCurve_3rdDeg))+(2*pi*centerFreq), (1/6)*Avg3rdOrderCurve_3rdDeg*(1e15^3), ...
    splineAxisNoCarrier(1:length(Avg3rdOrderCurve_4thDeg))+(2*pi*centerFreq), (1/6)*Avg3rdOrderCurve_4thDeg*(1e15^3), ...
    splineAxisNoCarrier(1:length(Avg3rdOrderCurve_5thDeg))+(2*pi*centerFreq), (1/6)*Avg3rdOrderCurve_5thDeg*(1e15^3), ...
    splineAxisNoCarrier(1:length(Avg3rdOrderCurve_6thDeg))+(2*pi*centerFreq), (1/6)*Avg3rdOrderCurve_6thDeg*(1e15^3));
hold on
plot(splineAxisNoCarrier+(2*pi*centerFreq),zeros(1,length(splineAxisNoCarrier)),'--')
legend('3rd deg' ,'4th deg', '5th deg','6th deg');
title('Average TOD Curves - various polynomial degrees');
ylabel('Magnitude (fs^3)');
xlabel('Frequency (Hz)') ;

figure(13)
plot(lambdaSplineAxis(1:length(Avg2ndOrderCurve_3rdDeg)), (1/2)*Avg2ndOrderCurve_3rdDeg*(1e15^2), ...
    lambdaSplineAxis(1:length(Avg2ndOrderCurve_4thDeg)), (1/2)*Avg2ndOrderCurve_4thDeg*(1e15^2), ...
    lambdaSplineAxis(1:length(Avg2ndOrderCurve_5thDeg)), (1/2)*Avg2ndOrderCurve_5thDeg*(1e15^2), ...
    lambdaSplineAxis(1:length(Avg2ndOrderCurve_6thDeg)), (1/2)*Avg2ndOrderCurve_6thDeg*(1e15^2));
hold on
plot(lambdaSplineAxis,zeros(1,length(splineAxisNoCarrier)),'--')
legend('3rd deg' ,'4th deg', '5th deg','6th deg');
title('Average GDD Curves - various polynomial degrees');
ylabel('Magnitude (fs^2)');
xlabel('Frequency (Hz)') ;
