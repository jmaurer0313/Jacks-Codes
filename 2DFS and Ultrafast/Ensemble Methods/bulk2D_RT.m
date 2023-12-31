%AUTHOR: Jack Maurer 
% code to retime and then rephase the bulk 2D data using the nonlinear
% signals post interpolation 

% bring in the data and time axes 
c0 = 0.000299792458; % speed of light mm/fs
base_folder = 'C:\Users\Jack\Documents\Marcus Lab\Data\2DFS\PopStudy\20201020-191604';
[dmat, smat, t21ax, t43ax] = bulk2Dread(base_folder);

% this is the desired time resoltuion after interpolation. this will set
% the number of zeros used in the frequency domain to acheive the proper
% spacing in time upon back transformation (typically Ts is ~2.66fs for a
% bulk 2D scan to start, so this is roughly 10 fold interpolation)
Tdes=0.6;
plotOpt=1; 


if length(t21ax)==length(t43ax) && (t43ax(2)==t21ax(2))
    disp('Times equal in length and step size');
    Ts=t43ax(2)-t43ax(1); 
else
    disp('Time axes unbalanced and possibly different in resolution, must append to generate Ts');  
end 
% time axes come in terms of fs, so they do NOT need ot be rescaled
% (unlike FPGA data, which gets reported in picoseconds
NFFT = 2^(nextpow2(length(t21ax))+2); 
NIFFT = 2^(nextpow2((Ts/Tdes)*NFFT)); 
Npad = abs(NIFFT-NFFT);

% introduce windows for each dimension to ensure the data goes to zero
% prior to transformation (two windows will be totally equal if data is
% symmetric in resolution and length). the set (taxis^2, 55, 10) has a 65fs
% time domain going to 1e-7~1e-8 by the end -- should be enough
[wX, wY] = meshgrid(t21ax,t43ax);
window = delayedGaussian(sqrt(wX.^2 + wY.^2), 50, 10);

smatW=smat.*window;
dmatW=dmat.*window;

lam_mono = 553; % monochromater wavelength in nanometers
Fmono = 1e7 / lam_mono;
% T = 7/6/1e3*2/c0;
Fs = 10/(c0*Ts);
Faxis = Fs/2*linspace(-1,1,NFFT);

f = ceil((-NFFT/2:NFFT/2-1))/(Ts*NFFT); % freq axis in PHz
fWN=f*0.000029979; 
% plot the raw data with the window applied to ensure it looks right prior
% to transform.
if plotOpt
 figure(1)
 subplot(2,2,1)
contourf(t21ax, t43ax, real(smatW),10)
title('NRP real - windowed')

axis equal
subplot(2,2,2)
contourf(t21ax, t43ax, imag(smatW),10)
title('NRP imag - windowed')
axis equal
subplot(2,2,3)
contourf(t21ax, t43ax, real(dmatW),10)
title('RP real - windowed')
axis equal
subplot(2,2,4)
contourf(t21ax, t43ax, imag(dmatW),10)
title('RP imag - windowed')
axis equal
end 

if plotOpt
 figure(2)
 hold on
 for i=1:8
 plot(t43ax, abs(smatW(i,:)));
 end 
end 

% first task is to get the timing corrections along each dimension, store
% them, generate A NEW window which is corrected by the shifts, window the
% raw data again. Transform, interpolate THEN zero out up to the true time
% zero, circshift the extra zeros out to the edge of the time domain. 
%%
% **********************************************************************
% TIMING CORRECTIONS ARE GOING TO BE BASED ON THE DMAT (RP) SIGNAL SINCE
% IT'S ABSOLUTE VALUE IS TYPICALLY QUITE GAUSSIAN AROUND TIME ZERO.
% *********************************************************************

%  try phasing the individual rows here and storing their values in a vector
% for reverse correction after interpolation
xRPphases=zeros(1,size(dmatW,1));
xNRPphases=zeros(1,size(smatW,1)); 

yRPphases=zeros(1,size(dmatW,2));
yNRPphases=zeros(1,size(smatW,2)); 

for m=1:length(xRPphases)
    xRPphases(m)=angle(dmatW(m,1)); 
    xNRPphases(m)=angle(smatW(m,1)); 
end 

for n=1:length(xRPphases)
    yRPphases(n)=angle(dmatW(1,n)); 
    yNRPphases(n)=angle(smatW(1,n)); 
end 

% the "x" phases are going to the rows (so should be a column vector to
% multiply them) "y" phases go to the columns, so a row vector is needed
xRPphases=xRPphases.';
xNRPphases=xNRPphases.';

%"x" dimension (should be rows) 
smatWPx=smatW.*(exp(1i*(-1*xNRPphases))); 
dmatWPx=dmatW.*(exp(1i*(-1*xRPphases)));

%"y" dimension (should be columns) 
smatWPy=smatW.*(exp(1i*(-1*yNRPphases))); 
dmatWPy=dmatW.*(exp(1i*(-1*yRPphases)));

%%
smatWPx(:,1) = smatWPx(:,1)*0.5;
dmatWPx(:,1) = dmatWPx(:,1)*0.5;

smatWPy(1,:) = smatWPy(1,:)*0.5;
dmatWPy(1,:) = dmatWPy(1,:)*0.5;

% FIRST TRANSFORM FROM TIME to OMEGA
smatWPxFT = fftshift(fft(smatWPx, NFFT, 2),2); % Signal in the frequency domain, FFT shift applied 
dmatWPxFT = fftshift(fft(dmatWPx, NFFT, 2),2); % Signal in the frequency domain, FFT shift applied 

% FIRST TRANSFORM FROM TIME to OMEGA
smatWPyFT = fftshift(fft(smatWPy, NFFT, 1),1); % Signal in the frequency domain, FFT shift applied 
dmatWPyFT = fftshift(fft(dmatWPy, NFFT, 1),1); % Signal in the frequency domain, FFT shift applied 


% try using a gaussian window to help cut down the noise ih the FFT, such
% that the zero padding produces a very clean result
arbAxis=linspace(0,NFFT,NFFT); 
gauss = @(x, a) sqrt(a/pi)*exp(-a*x.^2); % general gaussian of std = 1/sqrt(2a);
width=0.002;
centerX=56;
centerY=73;
centerYs=56; 

FFTwindowx=(gauss(arbAxis-centerX, width)); 
FFTwindowy=(gauss(arbAxis-centerY, width)); 
FFTwindowYs=(gauss(arbAxis-centerYs, width)); 

FFTwindowx=FFTwindowx*(max(max(abs(dmatWPxFT)))/max(FFTwindowx));
FFTwindowy=FFTwindowy*(max(max(abs(dmatWPyFT)))/max(FFTwindowy));
FFTwindowYs=FFTwindowYs*(max(max(abs(dmatWPyFT)))/max(FFTwindowYs));

% THIS WINDOW FUNCTION MUST BE CHANGED FOR: A)ALTERNATE LASER BANDWITH
% B)DIFFERING VALUES OF NFFT, SUCH THAT THE NUMBER OF POINTS IN THE FFT
% DIFFERS 

% possible that the scaling which the gaussian window performs miht distort
% one dimension versus the other. May need to resort to a equal amplitude
% variable center gaussian to achieve best outcome


if plotOpt
 figure(3)
 hold on
 for i=1:3
 plot(real(dmatWPxFT(i,:)));
 plot(real(smatWPxFT(i,:)));
 plot(FFTwindowx); 
 title('X dim');
 
 end 
 
 figure(4)
 hold on
 for i=1:3
 plot(real(dmatWPyFT(:,i)));
 plot(real(smatWPyFT(:,i)));
 plot(FFTwindowy.'); 
 plot(FFTwindowYs.'); 
 title('Y dim');
 
 end 
 
 figure(5)
 hold on
 for i=1:3
 plot(real(smatWPxFT(i,:)).*FFTwindowx);
 plot(real(dmatWPxFT(i,:)).*FFTwindowx);
%  plot(FFTwindow); 
 end 
 
 figure(6)
 hold on
 for i=1:3
 plot(real(smatWPyFT(:,i)).*(FFTwindowYs.'));
 plot(real(dmatWPyFT(:,i)).*(FFTwindowy.'));
%  plot(FFTwindow); 
 end 
 
end 
% 
for i=1:size(dmatWPxFT,1)
 smatWPxFT(i,:)=smatWPxFT(i,:).*FFTwindowx;
 dmatWPxFT(i,:)=dmatWPxFT(i,:).*FFTwindowx;
 
 smatWPyFT(:,i)=smatWPyFT(:,i).*(FFTwindowYs.');
 dmatWPyFT(:,i)=dmatWPyFT(:,i).*(FFTwindowy.');
 
 end 
%%

% ZERO PADDING
% try dropping the zeros at the edge of the freq domain after having windowed, OR the wrap around
% point where freq goes + to - 
smatWPxFTz=[real(smatWPxFT) , zeros(size(smatWPxFT,1),Npad)];
dmatWPxFTz=[real(dmatWPxFT) , zeros(size(dmatWPxFT,1),Npad)];

smatWPyFTz=[real(smatWPyFT) ; zeros(Npad, size(smatWPyFT,2))];
dmatWPyFTz=[real(dmatWPyFT) ; zeros(Npad, size(dmatWPyFT,2))];

% Shift the newly zero-padded array by the amount which it was off from
% zero to start, then calculate the proper complex vector for each row/col
% (depending on Dim) from the real valued vectors in the freq domain (since phase was zero going in).
smatWPxFTz= circshift(smatWPxFTz, -floor(NFFT/2),2);
dmatWPxFTz= circshift(dmatWPxFTz, -floor(NFFT/2),2);

smatWPyFTz= circshift(smatWPyFTz, -floor(NFFT/2),1);
dmatWPyFTz= circshift(dmatWPyFTz, -floor(NFFT/2),1);

dmatIntPx = (ifft(dmatWPxFTz, NIFFT, 2));
smatIntPx = (ifft(smatWPxFTz, NIFFT, 2));

dmatIntPy = (ifft(dmatWPyFTz, NIFFT, 1));
smatIntPy = (ifft(smatWPyFTz, NIFFT, 1));

dmatIntPx(:,((ceil(NIFFT/2)+1):end))=0;
smatIntPx(:,((ceil(NIFFT/2)+1):end))=0;

dmatIntPy(((ceil(NIFFT/2)+1):end),:)=0;
smatIntPy(((ceil(NIFFT/2)+1):end),:)=0;

% MIGHT NEED A FACTOR OF 2 TWICE HERE, ONCE FOR THE 0.5 ENTERING THE FFT
% AND ONCE FOR CHOPPING THE NEGATIVE FREQS OFF THE IFFT RESULT
dmatIntPx=dmatIntPx.*2;
smatIntPx=smatIntPx.*2;

dmatIntPy=dmatIntPy.*2;
smatIntPy=smatIntPy.*2;

% reapply the phases post-interpolation 
smatIntPx=smatIntPx.*(exp(1i*(1*xNRPphases))); 
dmatIntPx=dmatIntPx.*(exp(1i*(1*xRPphases)));

smatIntPy=smatIntPy.*(exp(1i*(1*yNRPphases))); 
dmatIntPy=dmatIntPy.*(exp(1i*(1*yRPphases)));

% these are the new time axes for the interpolated 2D data BUT TAKEN FROM
% ZERO. Until the X-axis is time corrected, plotting against these values
% will not be indivative of the true value of the signal,until the timing
% is corrected (they are reflective of the original times written to file)
Tsn = (Ts*NFFT)/NIFFT;
tn = ((0:NIFFT-1))*Tsn;
tnX=tn;
tnY=tn;

% the scaling for this plot is off due to the window applied to the FFT. to
% check the interpolation is done correctly, you cna simpyl turn off the
% application of the window (with the NIFFT/NFFT factor included below) and
% the coarse/interpolated time domain will overlay
if plotOpt
 figure(7)
 hold on
 for i=1:3
%  plot(t43ax, abs(dmatW(i,:)),'o');
 plot(tn, abs(dmatIntPx(i,:)*(NIFFT/NFFT)));
 plot(tn, abs(smatIntPx(i,:)*(NIFFT/NFFT)));
 title('Dmat and Smat X time Dom interp');
 end 
 
 figure(8)
 hold on
 for i=1:3
%  plot(t43ax, abs(dmatW(i,:)),'o');
 plot(tn, abs(dmatIntPy(:,i)*(NIFFT/NFFT)));
 plot(tn, abs(smatIntPy(:,i)*(NIFFT/NFFT)));
 title('Dmat and Smat Y time Dom interp');
 end 
 
end 


% now have the time domain interpolated in the x-dimension (rows), need to
% store a vector of the shifts to be applied
[dmatxShifts,dmatxISs]=max(abs(dmatIntPx.'));
[dmatyShifts,dmatyISs]=max(abs(dmatIntPy));
% [smatxShifts,smatxISs]=max(abs(dmatIntPx.'));

% these vectors now stores the times and number of indices by which each
% row/column needs to be shifted to acheive the theoretical time zero along
% each dimension
t43shifts=[tn([dmatxISs])',dmatxISs'];
t21shifts=[tn([dmatyISs])',dmatyISs'];

% Here the "xshifts" are the t43 stage, which increments only once per
% every run of t21, it's mis-timing

% apply the time shifts based on the RP matrix (dmat) which is more
% gaussian in the absolute than the corrsponding NRP (smat)

% ****NOTE: THE X DIMENSION MUST BE HANDLED IN TOTALITY PRIOR TO FFT AND
% INTERPOLATION OF THE Y DIMENSION (OR VICE VERSA) FOR THE SIZE AND
% ORDERING OF THE MATRICES TO MATCH THE NUMBER OF CORRECTIONS IN THE TIME
% LIST --- BEST PATH FORWARD IS TO GET THESE TIMIGN COREECTIONS AND THEN
% RE-DO THE FFT WITH THESE NEW TIMING CORRECTIONS IN MIND. THE PREVIOUS
% CODE RELIED ON KNOWLEDGE OF THE TIMING CORRECTIONS PRIOR TO FFT

% CREATE A NEW WINDOW THAT ACCOUNTS FOR THE TIMING CORRECTIONS TO BE
% APPLIED

% these are rows of the proper time axes (fot the t21 stage, which resets
% at each the end up each N step run, more prone to error, should contain a
% variety of timing corrections)
wt21= repmat(t21ax',size(smat,2),1) - repmat(t21shifts(:,1),1,size(smat,1));

% ONLY THE AVERAGE SHIFT WILL BE APPLIED IN T43 SINCE THERE IS TRULY ONLY A
% SINGLE TIMIGN ERROR GIVEN THAT THE STAGE MAKES A TOTAL OF n STEPS, RATHER
% THAN NXN STEPS AS IN THE CASE OF T21 STAGE
avgShift= mean(t43shifts(:,1));
avgShiftI= round(mean(t43shifts(:,2)));
avgShiftVec=ones(1,size(t43shifts,1))*avgShift; 
% these are columns of the proper time axes (corresponding to the t43 stage
% which steps just once for each run of the t21 stage)
% wt43= repmat(t43ax,1,size(smat,1)) - repmat(t43shifts(:,1)',size(smat,2),1);
wt43= repmat(t43ax,1,size(smat,1)) - repmat(avgShiftVec,size(smat,2),1);

window2 = delayedGaussian(sqrt(wt21.^2 + wt43.^2), 50, 10);


% START A NEW ROUND OF WINDOWING AND TRANSFORMS TO ACHEIVE THE RETIMED
% SPECTRA

smatW2=smat.*window2;
dmatW2=dmat.*window2;

% ALL OPERATIONS FOR THE X DIMENSION WILL BE PERFORMED FIRST, THEN Y, IN
% ORDER TO RETAIN A SINGLE MATRIX A THE END WHICH IS TIME CORERECTED IN
% BOTH AXES.

% the phases recorded for the first point in the time series are the same
% as the intial assesment (have not been altered by ore than 1e-15rad after all
% opertations)

%"x" dimension (should be rows) 
smatWPx2=smatW2.*(exp(1i*(-1*xNRPphases))); 
dmatWPx2=dmatW2.*(exp(1i*(-1*xRPphases)));

smatWPx2(:,1) = smatWPx2(:,1)*0.5;
dmatWPx2(:,1) = dmatWPx2(:,1)*0.5;

% FIRST TRANSFORM FROM TIME to OMEGA (x dim)
smatWPxFT2 = fftshift(fft(smatWPx2, NFFT, 2),2); % Signal in the frequency domain, FFT shift applied 
dmatWPxFT2 = fftshift(fft(dmatWPx2, NFFT, 2),2); % Signal in the frequency domain, FFT shift applied 

if plotOpt
 figure(9)
 hold on
 for i=1:3
 plot(real(dmatWPxFT2(i,:)));
 plot(real(smatWPxFT2(i,:)));
 plot(FFTwindowx); 
 title('X dim');
 
 end 
end 

 
if plotOpt
 figure(10)
 hold on
 for i=1:3
 plot(real(dmatWPxFT2(i,:)).*FFTwindowx);
 plot(real(smatWPxFT2(i,:)).*FFTwindowx);
 
 title('X dim');
 
 end 
end

for i=1:size(dmatWPxFT,1)
 smatWPxFT2(i,:)=smatWPxFT2(i,:).*FFTwindowx;
 dmatWPxFT2(i,:)=dmatWPxFT2(i,:).*FFTwindowx; 
end 

% ZERO PADDING
% try dropping the zeros at the edge of the freq domain after having windowed, OR the wrap around
% point where freq goes + to - 
smatWPxFTz2=[real(smatWPxFT2) , zeros(size(smatWPxFT2,1),Npad)];
dmatWPxFTz2=[real(dmatWPxFT2) , zeros(size(dmatWPxFT2,1),Npad)];

smatWPxFTz2= circshift(smatWPxFTz2, -floor(NFFT/2),2);
dmatWPxFTz2= circshift(dmatWPxFTz2, -floor(NFFT/2),2);

dmatIntPx2 = (ifft(dmatWPxFTz2, NIFFT, 2));
smatIntPx2 = (ifft(smatWPxFTz2, NIFFT, 2));

dmatIntPx2(:,((ceil(NIFFT/2)+1):end))=0;
smatIntPx2(:,((ceil(NIFFT/2)+1):end))=0;

dmatIntPx2=dmatIntPx2.*2;
smatIntPx2=smatIntPx2.*2;

% reapply the phases post-interpolation 
smatIntPx2=smatIntPx2.*(exp(1i*(1*xNRPphases))); 
dmatIntPx2=dmatIntPx2.*(exp(1i*(1*xRPphases)));

% the x dimension is now ready to be retimed post inerpolation based on the
% coreections obtained from the first round of analysis (not yet clear
% which set of shifts is the 'correct' one. will be obvious after applying
% correction.
for i=1:size(smatIntPx2,1)
    if (t21shifts(i,2))~=0
    smatIntPx2(i,(1:t21shifts(i,2)-1))=0; 
    smatIntPx2(i,:)=circshift((smatIntPx2(i,:)),-(t21shifts(i,2)-1),2);
%     smatIntPx2(i,:)=circshift((smatIntPx2(i,:)),-(t21shifts(i,2)),2);  

    
    dmatIntPx2(i,(1:t21shifts(i,2)-1))=0; 
    dmatIntPx2(i,:)=circshift((dmatIntPx2(i,:)),-(t21shifts(i,2)-1),2);  
    end
end

if plotOpt
 figure(11)
 hold on
 for i=1:3
 plot(abs(dmatIntPx2(i,:)));
 plot(abs(smatIntPx2(i,:)));
 
 title('X dim retimed');
 
 end 
end

% assign the newly timed x dim matirx as the starting point for the y dim
smatWPy2=smatIntPx2;
dmatWPy2=dmatIntPx2;

% need to get new phases for the now longer dimension along Y (due to x
% interpolation)
yRPphases2=zeros(1,size(dmatWPy2,2));
yNRPphases2=zeros(1,size(smatWPy2,2)); 

for n=1:size(smatWPy2,2)
    yRPphases2(n)=angle(dmatWPy2(1,n)); 
    yNRPphases2(n)=angle(smatWPy2(1,n)); 
end 

% pre phase the y dimension
%"y" dimension (should be columns) 
smatWPy2=smatWPy2.*(exp(1i*(-1*yNRPphases2))); 
dmatWPy2=dmatWPy2.*(exp(1i*(-1*yRPphases2)));

if plotOpt
 figure(12)
 hold on
 for i=1:3
 plot(t43ax,real(smatWPy2(:,i)));
 plot(t43ax,real(dmatWPy2(:,i)));
 title('Y dim real time domain post X interp');
 end
end 

smatWPy2(1,:) = smatWPy2(1,:)*0.5;
dmatWPy2(1,:) = dmatWPy2(1,:)*0.5;

% FIRST TRANSFORM FROM TIME to OMEGA (y dim)
smatWPyFT2 = fftshift(fft(smatWPy2, NFFT, 1),1); % Signal in the frequency domain, FFT shift applied 
dmatWPyFT2 = fftshift(fft(dmatWPy2, NFFT, 1),1); % Signal in the frequency domain, FFT shift applied 

if plotOpt
 figure(13)
 hold on
 for i=1:3
 plot(real(dmatWPyFT2(:,i)));
 plot(real(smatWPyFT2(:,i)));
 plot(FFTwindowy.'); 
 plot(FFTwindowYs.'); 
 title('Y dim FFT post interp');
 end

 figure(14)
 hold on
 for i=1:3
 plot(real(dmatWPyFT2(:,i)).*FFTwindowy.');
 plot(real(smatWPyFT2(:,i)).*FFTwindowYs.');
%  plot(FFTwindowy.'); 
%  plot(FFTwindowYs.'); 
 title('Y dim FFT post interp - windowed');
 end
end

smatWPyFTz2=[real(smatWPyFT2) ; zeros(Npad, size(smatWPyFT2,2))];
dmatWPyFTz2=[real(dmatWPyFT2) ; zeros(Npad, size(dmatWPyFT2,2))];


smatWPyFTz2= circshift(smatWPyFTz2, -floor(NFFT/2),1);
dmatWPyFTz2= circshift(dmatWPyFTz2, -floor(NFFT/2),1);

dmatIntPy2 = (ifft(dmatWPyFTz2, NIFFT, 1));
smatIntPy2 = (ifft(smatWPyFTz2, NIFFT, 1));

dmatIntPy2(((ceil(NIFFT/2)+1):end),:)=0;
smatIntPy2(((ceil(NIFFT/2)+1):end),:)=0;

dmatIntPy2=dmatIntPy2.*2;
smatIntPy2=smatIntPy2.*2;

% reapply the phases post-interpolation 
smatIntPy2=smatIntPy2.*(exp(1i*(1*yNRPphases2))); 
dmatIntPy2=dmatIntPy2.*(exp(1i*(1*yRPphases2)));

% retime the y dimension
for i=1:size(smatIntPy2,2)
    if (avgShiftI)~=0
    smatIntPy2((1:(avgShiftI-1)),i)=0; 
    smatIntPy2(:,i)=circshift((smatIntPy2(:,i)),-(avgShiftI-1),1);
%     smatIntPx2(i,:)=circshift((smatIntPx2(i,:)),-(t21shifts(i,2)),2);  

    
    dmatIntPy2((1:(avgShiftI-1)),i)=0; 
    dmatIntPy2(:,i)=circshift((dmatIntPy2(:,i)),-(avgShiftI-1),1); 
    end
end

% now phase by the new time zero
smatIntPy2 = exp(1i*(-angle(smatIntPy2(1,1))))*smatIntPy2;
dmatIntPy2 = exp(1i*(-angle(dmatIntPy2(1,1))))*dmatIntPy2;

smatIntPy2= smatIntPy2.*(NIFFT/NFFT);
dmatIntPy2= dmatIntPy2.*(NIFFT/NFFT);

if plotOpt
    figure(15)
    subplot(2,2,1)
    contourf(tn, tn, real(smatIntPy2),10)
    title('NRP real - retimed')
    axis equal
    
    subplot(2,2,2)
    contourf(tn, tn, imag(smatIntPy2),10)
    title('NRP imag - retimed')
    axis equal
    
    subplot(2,2,3)
    contourf(tn, tn, real(dmatIntPy2),10)
    title('RP real - retimed')
    axis equal
    
    subplot(2,2,4)
    contourf(tn, tn, imag(dmatIntPy2),10)
    title('RP imag - retimed')
    axis equal
    
end 



