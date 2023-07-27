% === file structure of Ek.dat and Speck.dat ===
% x=1\t Intensity\t Phase\t Real()\t Imaginary()\t\n
% :
% :
% :
% x=N\t Intensity\t Phase\t Real()\t Imaginary()\t\n
% === EOF ===
% where x:= time in Ek.dat and x:= lambda in Speck.dat

% ****************PROPER FILE NAME REQUIRED FOR ANALYSIS ******************
% mySpeckMat = dlmread('FROGtrace_NOPA_w1inchGlass_noCompress_04182022_trace2_80pPow.bin.Speck.dat');
% myEkMat = dlmread('FROGtrace_NOPA_w1inchGlass_noCompress_04182022_trace2_80pPow.bin.Ek.dat'); 

mySpeckMat = dlmread('FROGtrace_NOPA_w1inchGlass_noCompress_04182022_trace1_80pPow.bin.Speck.dat');
myEkMat = dlmread('FROGtrace_NOPA_w1inchGlass_noCompress_04182022_trace1_80pPow.bin.Ek.dat');

% mySpeckMat = dlmread('FROG_TRACE_05052022_trace3_realigned_wGlass.bin.Speck.dat');
% myEkMat = dlmread('FROG_TRACE_05052022_trace3_realigned_wGlass.bin.Ek.dat');
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

applyHigherTerms=1;
% **************************************************

freqUB= 3e8/(lambdaLB*1e-9);
freqLB= 3e8/(lambdaUB*1e-9);
centerFreq = 3e8/(centerLambda*1e-9); 


% first reconstruct the plots from Frogger() as a soft check on errors made
% in the data files
%%
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

splineAxis=linspace(specCutAxis(1),specCutAxis(end),1000); 
specSpline= spline(specCutAxis, specPhaseCut,splineAxis);

splineAxisNoCarrier=linspace(specCutAxisNoCarrier(1),specCutAxisNoCarrier(end),1000); 
specSplineNoCarrier= spline(specCutAxisNoCarrier, specPhaseCut,splineAxisNoCarrier);  

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
    
    poly5_2ndOrdNC = poly5(4)*(1e15^2);
    poly5_3rdOrdNC = poly5(3)*(1e15^3);
    
    poly4_2ndOrdNC = poly4(3)*(1e15^2);
    poly4_3rdOrdNC = poly4(2)*(1e15^3);
    
    poly3_2ndOrdNC = poly3(2)*(1e15^2);
    poly3_3rdOrdNC = poly3(1)*(1e15^3);
    
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