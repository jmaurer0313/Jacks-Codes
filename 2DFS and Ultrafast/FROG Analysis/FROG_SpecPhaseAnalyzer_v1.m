% === file structure of Ek.dat and Speck.dat ===
% x=1\t Intensity\t Phase\t Real()\t Imaginary()\t\n
% :
% :
% :
% x=N\t Intensity\t Phase\t Real()\t Imaginary()\t\n
% === EOF ===
% where x:= time in Ek.dat and x:= lambda in Speck.dat
% mySpeckMat = dlmread('FROGtrace_NOPA_w1inchGlass_noCompress_04182022_trace2_80pPow.bin.Speck.dat');
% myEkMat = dlmread('FROGtrace_NOPA_w1inchGlass_noCompress_04182022_trace2_80pPow.bin.Ek.dat'); 

% mySpeckMat = dlmread('FROGtrace_NOPA_w1inchGlass_noCompress_04052022_trace1_85pPow_FINAL_v2.bin.Speck.dat');
% myEkMat = dlmread('FROGtrace_NOPA_w1inchGlass_noCompress_04052022_trace1_85pPow_FINAL_v2.bin.Ek.dat');

mySpeckMat = dlmread('FROG_pharos_noGlass_04202022_trace2.bin.Speck.dat');
myEkMat = dlmread('FROG_pharos_noGlass_04202022_trace2.bin.Ek.dat');

% flipTime=1;
% if flipTime
%  myEkMat(:,2)= fliplr(myEkMat(:,2));   
%  myEkMat(:,3)= flipud(myEkMat(:,3));   
%  
%  mySpeckMat(:,3) = flipud(mySpeckMat(:,3));
% end 

% lambdaLB = 530;
% lambdaUB = 580; 
lambdaLB = 1000;
lambdaUB = 1035;

applyHigherTerms=0;
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


specAxisTotal = mySpeckMat(:,1);
specPhaseTotal = mySpeckMat(:,3);

% fit a gaussian to the bandwidth to establish the center frequency
[fitobj, goodness, output, convmsg]= fit(specAxisTotal,mySpeckMat(:,2),'gauss1');
centerLambda=fitobj.b1;

firstIndList = find(specAxisTotal>=lambdaLB);
secondIndList = find(specAxisTotal<=lambdaUB);

specCutAxis = specAxisTotal(firstIndList(1):secondIndList(end)); 
specPhaseCut= specPhaseTotal(firstIndList(1):secondIndList(end));

splineAxis=linspace(specCutAxis(1),specCutAxis(end),1000); 
specSpline= spline(specCutAxis, specPhaseCut,splineAxis);


% now with the relevant region isolated from the spectral phase, take the
% 1st, 2nd and 3rd derivatives to isolate the relevant dispersion terms
diffArr= splineAxis-circshift(splineAxis,1,2);
h = mode(diffArr(2:end));       % step size
X = splineAxis;    % domain
f = specSpline;      % range
firstOrd = diff(f)/h;   % first derivative
secondOrd = diff(firstOrd)/h;   % second derivative
thirdOrd = diff(secondOrd)/h;   % second derivative

% now make a plot of the successive derivatives over the domain
figure(3)
plot(splineAxis(1:length(firstOrd)), firstOrd, ...
    splineAxis(1:length(secondOrd)), secondOrd, ...
    splineAxis(1:length(thirdOrd)), thirdOrd)
hold on
plot(splineAxis,zeros(1,length(splineAxis)),'--')
legend('First Order' ,'Second Order' ,'Third Order');
title('Spectral Phase Derivatives - Disperion Orders'); 

%%

% now apply the deviations from the center wavelength of the bandwidth to
% 1st, 2nd, 3rd order (to try and fall in line with the dispersion expansion
% equation)
figure(4)
plot(splineAxis(1:length(firstOrd)), firstOrd.*(splineAxis(1:length(firstOrd))-centerLambda), ...
    splineAxis(1:length(secondOrd)), secondOrd.*((splineAxis(1:length(secondOrd))-centerLambda).^2), ...
    splineAxis(1:length(thirdOrd)), thirdOrd.*((splineAxis(1:length(thirdOrd))-centerLambda).^3))
hold on
plot(splineAxis,zeros(1,length(splineAxis)),'--')
legend('First Order' ,'Second Order' ,'Third Order');
title('Spectral Phase Derivatives - Disperion Orders with Magnitudes'); 


%%
% now use a polynomial fit to the spectral phase curve and take its
% derivatives rather than going from raw data - take the fit up to poly3 to
% capture whatever contribution from third order might be present
    [poly3, poly3Struc] = polyfit(splineAxis,specSpline, 3);
    coeff = polyfit(splineAxis,specSpline, 3);
    y3 = polyval(poly3, splineAxis);
    
 % now with the relevant region isolated from the spectral phase, take the
% 1st, 2nd and 3rd derivatives to isolate the relevant dispersion terms
diffArr= splineAxis-circshift(splineAxis,1,2);
h = mode(diffArr(2:end));       % step size
X = splineAxis;    % domain
f = y3;      % range
firstOrd = diff(f)/h;   % first derivative
secondOrd = diff(firstOrd)/h;   % second derivative
thirdOrd = diff(secondOrd)/h;   % second derivative

% now make a plot of the successive derivatives over the domain
figure(5)
plot(splineAxis(1:length(firstOrd)), firstOrd, ...
    splineAxis(1:length(secondOrd)), secondOrd, ...
    splineAxis(1:length(thirdOrd)), thirdOrd)
hold on
plot(splineAxis,zeros(1,length(splineAxis)),'--')
legend('First Order' ,'Second Order' ,'Third Order');
title('Spectral Phase Derivatives -PolyFit - Disperion Orders');    

if applyHigherTerms
%%
% now apply the deviations from the center wavelength of the bandwidth to
% 1st, 2nd, 3rd order (to try and fall in line with the dispersion expansion
% equation)
figure(6)
plot(splineAxis(1:length(firstOrd)), firstOrd.*(splineAxis(1:length(firstOrd))-centerLambda), ...
    splineAxis(1:length(secondOrd)), secondOrd.*((splineAxis(1:length(secondOrd))-centerLambda).^2), ...
    splineAxis(1:length(thirdOrd)), thirdOrd.*((splineAxis(1:length(thirdOrd))-centerLambda).^3))
hold on
plot(splineAxis,zeros(1,length(splineAxis)),'--')
legend('First Order' ,'Second Order' ,'Third Order');
title('Spectral Phase Derivatives - Disperion Orders with Magnitudes'); 

%%
%%
% now use a polynomial fit to the spectral phase curve and take its
% derivatives rather than going from raw data - take the fit up to poly4 to
% capture whatever contribution from third order and beyond might be present
    [poly4, poly4Struc] = polyfit(splineAxis,specSpline, 4);
    coeff = polyfit(splineAxis,specSpline, 4);
    y4 = polyval(poly4, splineAxis);
    
 % now with the relevant region isolated from the spectral phase, take the
% 1st, 2nd and 3rd derivatives to isolate the relevant dispersion terms
diffArr= splineAxis-circshift(splineAxis,1,2);
h = mode(diffArr(2:end));       % step size
X = splineAxis;    % domain
f = y4;      % range
firstOrd = diff(f)/h;   % first derivative
secondOrd = diff(firstOrd)/h;   % second derivative
thirdOrd = diff(secondOrd)/h;   % second derivative

% now make a plot of the successive derivatives over the domain
figure(7)
plot(splineAxis(1:length(firstOrd)), firstOrd, ...
    splineAxis(1:length(secondOrd)), secondOrd, ...
    splineAxis(1:length(thirdOrd)), thirdOrd)
hold on
plot(splineAxis,zeros(1,length(splineAxis)),'--')
legend('First Order' ,'Second Order' ,'Third Order');
title('Spectral Phase Derivatives -PolyFit - Disperion Orders');  

%%
%%
% % now use a polynomial fit to the spectral phase curve and take its
% % derivatives rather than going from raw data - take the fit up to poly2 to
% % capture contributions up to 2nd order
%     [poly2, poly2Struc] = polyfit(splineAxis,specSpline, 2);
%     coeff = polyfit(splineAxis,specSpline, 2);
%     y2 = polyval(poly2, splineAxis);
%     
%  % now with the relevant region isolated from the spectral phase, take the
% % 1st, 2nd and 3rd derivatives to isolate the relevant dispersion terms
% diffArr= splineAxis-circshift(splineAxis,1,2);
% h = mode(diffArr(2:end));       % step size
% X = splineAxis;    % domain
% f = y2;      % range
% firstOrd = diff(f)/h;   % first derivative
% secondOrd = diff(firstOrd)/h;   % second derivative
% thirdOrd = diff(secondOrd)/h;   % second derivative
% 
% % now make a plot of the successive derivatives over the domain
% figure(8)
% plot(splineAxis(1:length(firstOrd)), firstOrd, ...
%     splineAxis(1:length(secondOrd)), secondOrd, ...
%     splineAxis(1:length(thirdOrd)), thirdOrd)
% hold on
% plot(splineAxis,zeros(1,length(splineAxis)),'--')
% legend('First Order' ,'Second Order' ,'Third Order');
% title('Spectral Phase Derivatives -PolyFit - Disperion Orders');  
% 
end 