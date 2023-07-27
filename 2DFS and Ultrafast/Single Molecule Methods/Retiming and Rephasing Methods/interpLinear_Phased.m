% GOAL: take in the X and Y linear data, as well as the tb1 nd tb2 time
% bases and a Tdes (desired time resolution post interpolation). Window, transform, interpllate, back transform and then find
% maxima for all the runs of X and Y. Store these in a vector of
% "timeshifts" for each dimension. In the case of Y, where we may want to
% perform the time shifts as one average, implement an outliers function to
% toss any values which are well outside the others.

function[XnLin, YnLin, Xshifts, Yshifts, Xangles, Yangles, tnX, tnY, tx, ty, NIFFT,Xzero] = interpLinear_Phased(plotOpt,XLin, YLin, tx, ty, Tdes) 

if Tdes<1e-3
    disp('Desired time resolution is too small');
    return; 
end 

if abs(Tdes)>abs(tx(1)-tx(2))*1e3
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
tx=tx*1e3;
ty=ty*1e3;
c0 = 299.792458; % speed of light in nm/fs
NFFT = 2^(nextpow2(length(tx))); 
NIFFT = 2^(nextpow2((Ts/Tdes)*NFFT)); 
Npad = abs(NIFFT-NFFT);

% Introduce a windowing function to be applied to the linear 
% data in order to eliminate any ringing in the interploated time domain as
% a result of nonzero signal at the edge of the original raw data in time.
% IMPORTANTLY: the window must be applied to the data prior to the first
% FFT, otherwise the result of padding by a small number of zeros
% NFFT-length(S) will give a sharp edge in the intial time domain, which
% rings in the frequency domain. THE VALUES CHOSEN HERE ARE MEANT TO
% INTRODUCE A GENTLE WINDOW WHICH STARTS VERY EARLY BUT DECAYS QUITE SLOWLY.
% GIVEN OUR COARSELY SAMPLE DATA THIS APPEARS TO THE ONLY SCHEME WHICH
% WORKS EFFECTIVELY
windowX = delayedGaussian(tx, 0, 45);
windowY = delayedGaussian(ty, 0, 45);

%     this is the case of rows for FFT
XLin = XLin.*windowX;
XLinPlot= XLin.';



%     this is the case of columns for FFT
YLin = YLin.*windowY.';
YLinPlot= YLin;



f = ceil((-NFFT/2:NFFT/2-1))/(Ts*NFFT); % freq axis in PHz

if plotOpt
figure(1)
plot(ty, abs(YLinPlot));
title('Y Time Domain Data (Abs)- Raw'); 
xlabel('Time (fs)');
ylabel('Magnitude (arb)'); 

figure(2)
plot(tx, abs(XLinPlot));    
title('X Time Domain Data (Abs) -Raw');
xlabel('Time (fs)');
ylabel('Magnitude (arb)'); 

end

% phase the X and Y linear signals prior to transformation. Storing them in
% a vector
xPhases=zeros(1,size(XmatTemp,1));
yPhases=zeros(1,size(YmatTemp,2));

for j=1:length(xPhases)
   xPhases(j)=angle(XLin(j,1));
end

for n=1:length(yPhases)
   yPhases(n)=angle(YLin(1,n));
end

% multiply the first entry in the matrix by 1/2 to make the transform
% symmetric (see weird NMR reference...)
XLin(:,1) = XLin(:,1)*0.5; 
YLin(1,:) = YLin(1,:)*0.5;  

XFLin = fftshift(fft(XLin, NFFT, 2)); % Signal in the frequency domain, FFT shift applied
YFLin = fftshift(fft(YLin, NFFT, 1)); % Signal in the frequency domain, FFT shift applied

if plotOpt
figure(3)
plot(f, abs(YFLin));
% plot(f, real(SFLin), f, imag(SFLin), f, abs(SFLin));
title('Y Linear Abs vectors - raw/NOT-retimed');


figure(4)
plot(f, abs(XFLin.'));
% plot(f, real(SFLin.'), f, imag(SFLin.'), f, abs(SFLin.'));
title('X Linear Abs vectors - raw/NOT-retimed');
end 

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
% NIFFT = size(CunshXLin,2);
else
    disp('NIFFT does not match between the X and Y dimension');
    return;
end 
% new interpolated time axes with NIFFT points
Tn = (Ts*NFFT)/NIFFT;
tn = ((0:NIFFT-1))*Tn;
% tnX=(-(NIFFT/length(tx)):((NIFFT/length(tx))*(length(tx)-1))-1)*Tn; 
tnX= tn-Ts; 
tnY=tn; 

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
 
if plotOpt
figure(5)
plot(tnX, abs(XnLinPlot));
title('X Linear Scans -Abs- Interpolated');
xlabel('Time (fs)'); 
ylabel('Magnitude (arb)'); 
% plot(tn, real(SnLinPlot), tn, imag(SnLinPlot))
hold on
ax = gca;
ax.ColorOrderIndex = 1;
plot(tx', abs(XLinPlot),'o');
hold off

figure(6)
plot(tn, abs(YnLinPlot));
title('Y Linear Scans -Abs- Interpolated');
xlabel('Time (fs)'); 
ylabel('Magnitude (arb)'); 
% plot(tn, real(SnLinPlot), tn, imag(SnLinPlot))
hold on
ax = gca;
ax.ColorOrderIndex = 1;
plot(ty, abs(YLinPlot),'o');
hold off
end

% create a vector for the maxima in the interpolated X and Y linear data,
% make a vector of the indices as well as the value 
[txShifts,IxShifts]=max(abs(XnLin.'));
[tyShifts,IyShifts]=max(abs(YnLin));

% the maxtrix for X and Y shifts contaisn the raw time by which the scans
% are off from time zero and the index of the maximum in the interpolated
% time domain of order NIFFT. IN THE CASE OF THE Y SHIFTS, AS LONG AS THE
% TY AXIS STARTS FROM ZERO IT SHOULD BE THE iTH+1 ELEMENT SINCE tnY(1)=0
Xshifts=[tnX([IxShifts])',IxShifts'];
Yshifts=[tnY([IyShifts]+1)',IyShifts'];
Xzero=find(tnX==0);

% create a vector of phase angles for the Xmaxs and Ymaxs seperately
Xangles=zeros(1,size(Xshifts,1));
Yangles=zeros(1,size(Yshifts,1));

for v=1:size(Xshifts,1)
    Xangles(v)=angle(XnLin(v,Xshifts(v,2))); 
end 

for t=1:size(Yshifts,1)
    Yangles(t)=angle(YnLin(Yshifts(t,2),t)); 
end 

% in the case of the Yshifts, throw out any outliers based on the 3 scaled
% mediam absolute deviations which matlab uses within the function call.
% replace the elements which are outliers in the time and index list by the
% average of the remaining elements
tavg=0;
tcount=0;

IAvg=0;
Icount=0;

TF= isoutlier(Yshifts(:,1));
% define the average of the remaining elements for the times
for i=1:numel(Yshifts(:,1))
if ~TF(i)
tavg=tavg+Yshifts(i,1);
tcount=tcount+1;
end
end
tavg=tavg/tcount;

% now go back through the list and replace the outliers with the average
for j=1:numel(Yshifts(:,1))
if TF(j)
Yshifts(j,1)=tavg;
end
end

% define the average of the remaining elements for the indices
for m=1:numel(Yshifts(:,2))
if ~TF(m)
IAvg=IAvg+Yshifts(m,2);
Icount=Icount+1;
end
end
IAvg=ceil(IAvg/Icount);

for z=1:numel(Yshifts(:,2))
if TF(z)
Yshifts(z,2)=IAvg;
end
end

% Now zero out and circshift the linear data according to the shift indices
% FOR THE X DATA: The time axis reported will no longer be valid for this
% shifted linear data, but is the proper time axis for the 2D data to be
% windowed over. For the instance of plotting the Xlin data post timing
% correction, the domain should either not be speicifed or a new time axis
% should be created to match that of tnY
for i=1:size(XnLin,1)
%     if the difference btw the zero time index on the x axis and the shift is zero, 
%  then the x row is already centered at time zero.
    if Xshifts(i,2)~=0
    XnLin(i,(1:Xshifts(i,2)-1))=0; 
    XnLin(i,:)=circshift((XnLin(i,:)),-(Xshifts(i,2)-1),2);
%     else
%     XnLin(i,(1:Xzero-1))=0; 
%     XnLin(i,:)=circshift((XnLin(i,:)),-(Xzero-1),2);   
    end
end 

for j=1:size(YnLin,2)
    if Yshifts(j,2)~=0
    YnLin((1:Yshifts(j,2)-1),j)=0; 
    YnLin(:,j)=circshift((YnLin(:,j)),-(Yshifts(j,2)-1),1);
    end
end 

% Xphases=exp(-1i*angle(XnLin(:,1)));
% XnLin=XnLin.*Xphases; 
if plotOpt
figure(8)
plot(real(fftshift(fft(XnLin,2),2)).');
title('X freq domain post interp');
end

end