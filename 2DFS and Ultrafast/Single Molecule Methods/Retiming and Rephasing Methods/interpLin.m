% GOAL: take in the X and Y linear data, as well as the tb1 nd tb2 time
% bases and a Tdes (desired time resolution post interpolation). Window, transform, interpllate, back transform and then find
% maxima for all the runs of X and Y. Store these in a vector of
% "timeshifts" for each dimension. In the case of Y, where we may want to
% perform the time shifts as one average, implement an outliers function to
% toss any values which are well outside the others.

function[Xshifts,Yshifts, tn, NIFFT] = interpLinear(XLin, YLin, tX, ty, Tdes) 

if Tdes<1e-3
    disp('Desired time resolution is too small');
    return; 
end 

if abs(Tdes)>abs(tx(1)-tx(2))
    disp('Desired time resolution must be less than input time resolution');
    return; 
end 

% Ts must be in fs for the rest of the script to port over easily. The
% incoming time base is in ps (from the labview) so it must be scaled to fs
if abs(tx(1)-tx(2))==abs(ty(1)-ty(2))
Ts=abs(tx(1)-tx(2))*1e3; 
else
    disp('Step sizes do not match between X and Y');
    return;
end 
% the X time base comes in from the matlab with opposite sign than the Y
% due to the motion of the stages in the +/- diretion, correct it with a
% sign change here
tx=(-1*tx);
c0 = 299.792458; % speed of light in nm/fs
NFFT = 2^(nextpow2(length(tx))); 
Npad = 2^nextpow2((Ts/Tdes)*NFFT); 

% Introduce a windowing function to be applied to the linear 
% data in order to eliminate any ringing in the interploated time domain as
% a result of nonzero signal at the edge of the original raw data in time.
% IMPORTANTLY: the window must be applied to the data prior to the first
% FFT, otherwise the result of padding by a small number of zeros
% NFFT-length(S) will give a sharp edge in the intial time domain, which
% rings in the frequency domain. THE VALUES CHOSEN HERE ARE MEANT TO
% INTRODUCE A GENTLE WINDOW WHICH STARTS VERY EARLY BUT DECAYS QUIT SLOWLY.
% GIVEN OUR COARSELY SAMPLE DATA THIS APPEARS TO THE ONLY SCHEME WHICH
% WORKS EFFECTIVELY
windowX = delayedGaussian(tx, 0, 50);
windowY = delayedGaussian(ty, 0, 50);

%     this is the case of rows for FFT
XLin = XLin.*windowX; 
XLin(:,1) = XLin(:,1)*0.5; 
XLinPlot= XLin.';

%     this is the case of columns for FFT
YLin = YLin.*windowY.';
YLin(1,:) = YLin(1,:)*0.5;  
YLinPlot= YLin;

f = ceil((-NFFT/2:NFFT/2-1))/(Ts*NFFT); % freq axis in PHz

figure(1)
plot(ty, abs(YLinPlot));
title('Y Time Domain Data (Abs)'); 

figure(2)
plot(tx, abs(XLinPlot));    
title('X Time Domain Data (Abs)'); 


XFLin = fftshift(fft(XLin, NFFT, 2)); % Signal in the frequency domain, FFT shift applied
YFLin = fftshift(fft(YLin, NFFT, 1)); % Signal in the frequency domain, FFT shift applied

figure(3)
plot(f, abs(YFLin));
% plot(f, real(SFLin), f, imag(SFLin), f, abs(SFLin));
title('Y Linear Abs vectors');


figure(4)
plot(f, abs(XFLin.'));
% plot(f, real(SFLin.'), f, imag(SFLin.'), f, abs(SFLin.'));
title('X Linear Abs vectors');

% try dropping the zeros at the edge of the freq domain, OR the wrap around
% point where freq goes + to - 
YFLinzeros=[real(YFLin) ; zeros(Npad,size(YFLin,2))];
XFLinzeros=[real(XFLin) , zeros(size(XFLin,1),Npad)];

% Shift the newly zero-padded array by the amount which it was off from
% zero to start, then calculate the proper complex vector for each row/col
% (depending on Dim) from the real valued vectors in the freq domain.
% hilbert operates on columns when a matirx is passed.

unshiftedYLinNew= circshift(YFLinzeros, -floor(NFFT/2),1);
% real and complex unshifted matrices
RunshYLin = unshiftedYLinNew;
CunshYLin = complex(RunshYLin, imag(hilbert(-RunshYLin)));


unshiftedXLinNew= circshift(XFLinzeros, -floor(NFFT/2),2);
% real and complex unshifted matrices, transposed here to handle the hlbert
% operation properly
RunshXLin = unshiftedXLinNew.';
CunshXLin = complex(RunshXLin, imag(hilbert(-RunshXLin)));
CunshXLin=CunshXLin.'; 


% determine the new numbe of points and their values in the time domain for
% the interpolated data
if length(CunshXLin)==length(CunshYLin)
NIFFT = length(CunshXLin);
else
    disp('NIFFT does not match between the X and Y dimension');
    return;
end 
% new interpolated time axes with NIFFT points
Tn = (Ts*NFFT)/NIFFT;
tn = ((0:NIFFT-1))*Tn;

XnLin =(ifft(CunshXLin, NIFFT, 2));
YnLin =(ifft(CunshYLin, NIFFT, 1));


% rescaling and plotting matrices for X
XnLin(:,1) = XnLin(:,1)*2; 
XLin(:,1) = XLin(:,1)*2; 
XnLinPlot= XnLin.';
XLinPlot= XLin.';

% rescaling and plotting matrices for Y
YnLin(1,:) = YnLin(1,:)*2;
YLin(1,:) = YLin(1,:)*2;
YnLinPlot= YnLin;
YLinPlot= YLin;


XLinPlot=XLinPlot*NFFT/NIFFT;
YLinPlot=YLinPlot*NFFT/NIFFT;
 

figure(5)
plot(tn, abs(XnLinPlot));
% plot(tn, real(SnLinPlot), tn, imag(SnLinPlot))
hold on
ax = gca;
ax.ColorOrderIndex = 1;
plot(tx, abs(XLinPlot),'o');
hold off

figure(6)
plot(tn, abs(YnLinPlot));
% plot(tn, real(SnLinPlot), tn, imag(SnLinPlot))
hold on
ax = gca;
ax.ColorOrderIndex = 1;
plot(ty, abs(YLinPlot),'o');
hold off


end