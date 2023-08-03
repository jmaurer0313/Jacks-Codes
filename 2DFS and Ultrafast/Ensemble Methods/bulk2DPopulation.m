%AUTHOR: Jack Maurer 

c0 = 0.000299792458; % speed of light mm/fs

base_folder = 'C:\Users\Jack\Documents\Marcus Lab\Data\2DFS\PopStudy\20201027\20201027-181815';
[dmat, smat, t21ax, t43ax] = bulk2Dread(base_folder);
retimed=1;

windowt23 = delayedGaussian(t43ax, 100, 80);
windowt23_2D= repmat(windowt23,1,20);

% window based on the post-freq filter size and length of the time domain
% *the time domain is lost in the filter, pure indices are returned. The
% data is still the length of the original time axis, but with the wrong
% units, therefore a window will need the same mumber of indices as the
% original data but with an onset and slope that matches the pure value of
% the indexes
windowt23_Postfilt = delayedGaussian(linspace(1,length(t43ax),length(t43ax)), round(length(t43ax)/2), round(0.5*(length(t43ax)/2)));

% % window along the z direction
% windowt23_2D=windowt23_2D.';

% Rough phasing (must be manually adjusted here) 
% dmat = exp(1i*(-angle(dmat(1,1))))*dmat;
% smat = exp(1i*(-angle(smat(1,1))))*smat;

% % window along the z direction
% smat=smat.*windowt23_2D; 
% dmat=dmat.*windowt23_2D; 

% block of code to interpolate the very sparsely sampled time domain using
% cubic spline interpolation (rather than Fourier) to set the time zero
% along the t21 dimension more sensitively. Once timing corrections are
% applied, reexamine the t23 dimension for evidence of beating. Corrections
% will be based on the dmat (RP) signal as in the past/previous retimign
% efforts.
%%
% first interpolate the abs(dmat) to get the timing corrections
spIn=0:1:t21ax(end);
dtOut=spline(t21ax,abs(dmat.'),spIn);

dmatOut=spline(t21ax,(dmat.'),spIn);
smatOut=spline(t21ax,(smat.'),spIn);

% creat a vector of timing corrections
[mVal, Im]=max(dtOut.'); 

Iadj=Im-min(Im);


% in order to retain as much of the time domainas possible (minimal
% 'zeroing' out of the pre-T0 portion) we will take the column with the
% minimum index corresponding to its maxima and shift all other columns
% back relative to that in principle this is no different than the standard
% approach, except that the individual shifts are offset by a constant
% value, which is the smallest shift in the overall set. 

% start with a simple circshift without zeroing out the edges
for j=1:numel(Iadj)
    if Iadj(j)>0
    dmatOut(j,:)=circshift(dmatOut(j,:),Iadj(j),2); 
    smatOut(j,:)=circshift(smatOut(j,:),Iadj(j),2); 
    end
end 

dmatOut=dmatOut.';
smatOut=smatOut.'; 

%%
 if ~retimed
figure(1)
hold on
title('Raw T23 signals -Bkgd Subtracted- Real parts');
xlabel('T23(fs)');
ylabel('Demod. Mag.(arb)'); 

for i=1:round(size(smat,1)/2)
    
plot(t43ax, real(smat(i,:))-mean(real(smat(i,:))));
plot(t43ax, real(dmat(i,:))-mean(real(dmat(i,:))));

end 

figure(2)
hold on
title('Raw T23 signal -Bkgd Subtracted- Imag Parts');
xlabel('T23(fs)');
ylabel('Demod. Mag.(arb)'); 

for i=1:round(size(smat,1)/2)
    
plot(t43ax, imag(dmat(i,:))-mean(imag(dmat(i,:))));
plot(t43ax, imag(smat(i,:))-mean(imag(smat(i,:))));

end 

figure(3)
subplot(2,2,1)
contourf(t43ax, t21ax, real(smat),10)
title('NRP Real')
xlabel('T23 (fs)')
ylabel('T21 (fs)') 
% axis equal
subplot(2,2,2)
contourf(t43ax, t21ax, imag(smat),10)
title('NRP Imag')
xlabel('T23 (fs)')
ylabel('T21 (fs)')
% axis equal
subplot(2,2,3)
contourf(t43ax, t21ax, real(dmat),10)
title('RP Real')
xlabel('T23 (fs)')
ylabel('T21 (fs)')
% axis equal
subplot(2,2,4)
contourf(t43ax, t21ax, imag(dmat),10)
title('RP Imag')
xlabel('T23 (fs)')
ylabel('T21 (fs)')
% axis equal



% impose a filter onto the signals with variable range and then examine the
% resulting filtered signal 

% a period of 100fs corresponds to a frequnecy of 1.0e13

% this is the filter pass interval in Hz
% filterInt= [1e6 1e19];
% filterInt= [0.3e12 60e12];
filterInt= [1.5e12 20e12];
% filterInt= [25e12 60e12];
% convert t43ax to seconds
timeAxis=t43ax.*1e-15; 

for j=1:size(smat,1)
    smat(j,:)=(smat(j,:)*angle(smat(j,1)));
end 

for j=1:size(dmat,1)
    dmat(j,:)=(dmat(j,:)*angle(dmat(j,1)));
end 

for j=1:size(smat,1)
%     size(smat,1)
% smatRIn= timeseries((smat(j,:)-mean(smat(j,:))),timeAxis);
smatRIn= timeseries((real(smat(j,:))),timeAxis);
smatROut= idealfilter(smatRIn,filterInt,'pass');

dmatRIn= timeseries(real(dmat(j,:)),timeAxis);
dmatROut= idealfilter(dmatRIn,filterInt,'pass');

smatImIn= timeseries(imag(smat(j,:)),timeAxis);
smatImOut= idealfilter(smatImIn,filterInt,'pass');

dmatImIn= timeseries(imag(dmat(j,:)),timeAxis);
dmatImOut= idealfilter(dmatImIn,filterInt,'pass');

figure(5)
hold on
plot((smatROut));
title('NRP Filtered Time domain Signal-Bandpass: 10-2000 cm^{-1} - Real Part'); 
xlabel('T23 (sec)'); 

figure(6)
hold on
plot((dmatROut));
title('RP Filtered Time domain Signal-Bandpass: 10-2000 cm^{-1} - Real Part'); 
xlabel('T23 (sec)'); 

figure(9)
hold on
plot((smatImOut));
title('NRP Filtered Time domain Signal-Bandpass: 10-2000 cm^{-1} - Imag Part'); 
xlabel('T23 (sec)'); 

figure(10)
hold on
plot((dmatImOut));
title('RP Filtered Time domain Signal-Bandpass: 10-2000 cm^{-1} - Imag Part'); 
xlabel('T23 (sec)'); 


T = mean(t43ax(2:end) - t43ax(1:end-1)); % Uses stage positions directly, could replace by 0.0004/c0
NFFT = 1024;
Fs = 10/(c0*T);
% this frequency axis should be in wavenumbers (assuming the time steps are
% in fs
Faxis = Fs/2*linspace(-1,1,NFFT);

% get the data from time series back out of the object
smatReFilt= getdatasamples(smatROut,[1:length(smat)]);
smatReFilt=squeeze(smatReFilt(1,1,:));
sFToutRe= fftshift(fft(smatReFilt,NFFT,1),1);

% get the data from time series back out of the object
dmatReFilt= getdatasamples(dmatROut,[1:length(dmat)]);
dmatReFilt=squeeze(dmatReFilt(1,1,:));
dFToutRe= fftshift(fft(dmatReFilt,NFFT,1),1);




figure(7)
hold on
title('NRP FT of Filtered Signal - Bandpass: 10-2000 cm^{-1}');
plot(Faxis,abs(sFToutRe))
xlim([0 inf]);
xlabel('Wavenumbers');
ylabel('Mag (arb)'); 

figure(8)
hold on
title('RP FT of Filtered Signal - Bandpass: 10-2000 cm^{-1}');
plot(Faxis,abs(dFToutRe))
xlim([0 inf]);
xlabel('Wavenumbers');
ylabel('Mag (arb)'); 

end 

%*************Begin the retimed section of the code - same but retimed************* 
 else
   
     
     figure(1)
hold on
title('Raw T23 signals -Bkgd Subtracted- Real parts');
xlabel('T23(fs)');
ylabel('Demod. Mag.(arb)'); 

for i=1:round(size(smat,1)/2)
    
plot(t43ax, real(smatOut(i,:))-mean(real(smatOut(i,:))));
plot(t43ax, real(dmatOut(i,:))-mean(real(dmatOut(i,:))));

end 

figure(2)
hold on
title('Raw T23 signal -Bkgd Subtracted- Imag Parts');
xlabel('T23(fs)');
ylabel('Demod. Mag.(arb)'); 

for i=1:round(size(smat,1)/2)
    
plot(t43ax, imag(dmatOut(i,:))-mean(imag(dmatOut(i,:))));
plot(t43ax, imag(smatOut(i,:))-mean(imag(smatOut(i,:))));

end 

figure(3)
subplot(2,2,1)
contourf(t43ax, spIn, real(smatOut),10)
title('NRP Real')
xlabel('T23 (fs)')
ylabel('T21 (fs)') 
% axis equal
subplot(2,2,2)
contourf(t43ax, spIn, imag(smatOut),10)
title('NRP Imag')
xlabel('T23 (fs)')
ylabel('T21 (fs)')
% axis equal
subplot(2,2,3)
contourf(t43ax, spIn, real(dmatOut),10)
title('RP Real')
xlabel('T23 (fs)')
ylabel('T21 (fs)')
% axis equal
subplot(2,2,4)
contourf(t43ax, spIn, imag(dmatOut),10)
title('RP Imag')
xlabel('T23 (fs)')
ylabel('T21 (fs)')
% axis equal



% impose a filter onto the signals with variable range and then examine the
% resulting filtered signal 

% a period of 100fs corresponds to a frequnecy of 1.0e13

% this is the filter pass interval in Hz
% filterInt= [1e6 1e19];
% filterInt= [0.3e12 60e12];

%Current Filter (appears to be 50cm-1 to ~700cm-1
filterInt= [1.5e12 20e12];

% wider test filter from 10-2000cm-1
filterInt= [0.3e12 60e12];

% filterInt= [25e12 60e12];
% convert t43ax to seconds
timeAxis=t43ax.*1e-15; 

for j=1:size(smatOut,1)
    smatOut(j,:)=(smatOut(j,:)*angle(smatOut(j,1)));
end 

for j=1:size(dmatOut,1)
    dmatOut(j,:)=(dmatOut(j,:)*angle(dmatOut(j,1)));
end 

for j=1:5
%     size(smatOut,1)
%     size(smat,1)
% smatRIn= timeseries((smat(j,:)-mean(smat(j,:))),timeAxis);
smatRIn= timeseries((smatOut(j,:)),timeAxis);
smatROut= idealfilter(smatRIn,filterInt,'pass');

dmatRIn= timeseries((dmatOut(j,:)),timeAxis);
dmatROut= idealfilter(dmatRIn,filterInt,'pass');

smatImIn= timeseries(imag(smatOut(j,:)),timeAxis);
smatImOut= idealfilter(smatImIn,filterInt,'pass');

dmatImIn= timeseries(imag(dmatOut(j,:)),timeAxis);
dmatImOut= idealfilter(dmatImIn,filterInt,'pass');

figure(5)
hold on
plot((smatROut));
title('NRP Filtered Time domain Signal-Bandpass: 10-2000 cm^{-1} - Real Part'); 
xlabel('T23 (sec)'); 

figure(6)
hold on
plot((dmatROut));
title('RP Filtered Time domain Signal-Bandpass: 10-2000 cm^{-1} - Real Part'); 
xlabel('T23 (sec)'); 

figure(7)
hold on
plot((smatImOut));
title('NRP Filtered Time domain Signal-Bandpass: 10-2000 cm^{-1} - Imag Part'); 
xlabel('T23 (sec)'); 

figure(8)
hold on
plot((dmatImOut));
title('RP Filtered Time domain Signal-Bandpass: 10-2000 cm^{-1} - Imag Part'); 
xlabel('T23 (sec)'); 


T = mean(t43ax(2:end) - t43ax(1:end-1)); % Uses stage positions directly, could replace by 0.0004/c0
NFFT = 1024;
Fs = 10/(c0*T);
% this frequency axis should be in wavenumbers (assuming the time steps are
% in fs
Faxis = Fs/2*linspace(-1,1,NFFT);

%******REAL PARTS FILTERED******
% get the data from time series back out of the object
smatReFilt= getdatasamples(smatROut,[1:length(smat)]);
smatReFilt=squeeze(smatReFilt(1,1,:));
sFToutRe= fftshift(fft(smatReFilt,NFFT,1),1);

% get the data from time series back out of the object
dmatReFilt= getdatasamples(dmatROut,[1:length(dmat)]);
dmatReFilt=squeeze(dmatReFilt(1,1,:));
dFToutRe= fftshift(fft(dmatReFilt,NFFT,1),1);

% *****IMAG PARTS FILTERED ******
% get the data from time series back out of the object
smatImFilt= getdatasamples(smatImOut,[1:length(smat)]);
smatImFilt=squeeze(smatImFilt(1,1,:));
sFToutIm= fftshift(fft(smatImFilt,NFFT,1),1);

% get the data from time series back out of the object
dmatImFilt= getdatasamples(dmatImOut,[1:length(dmat)]);
dmatImFilt=squeeze(dmatImFilt(1,1,:));
dFToutIm= fftshift(fft(dmatImFilt,NFFT,1),1);

figure(9)
hold on
title('NRP Real FT of Filtered Signal - Bandpass: 10-2000 cm^{-1}');
plot(Faxis,abs(sFToutRe))
xlim([0 inf]);
xlabel('Wavenumbers');
ylabel('Mag (arb)'); 

figure(10)
hold on
title('RP Real FT of Filtered Signal - Bandpass: 10-2000 cm^{-1}');
plot(Faxis,abs(dFToutRe))
xlim([0 inf]);
xlabel('Wavenumbers');
ylabel('Mag (arb)'); 

figure(11)
hold on
title('NRP Imag FT of Filtered Signal - Bandpass: 10-2000 cm^{-1}');
plot(Faxis,abs(sFToutIm))
xlim([0 inf]);
xlabel('Wavenumbers');
ylabel('Mag (arb)'); 

figure(12)
hold on
title('RP Imag FT of Filtered Signal - Bandpass: 10-2000 cm^{-1}');
plot(Faxis,abs(dFToutIm))
xlim([0 inf]);
xlabel('Wavenumbers');
ylabel('Mag (arb)'); 

%Now perform the same transform of the filtered signal but window it
%properly - first by typical 1 sided window, but then also by using a
%double sided window (flat topped two sided gaussian) 
windowt23_Postfilt=reshape(windowt23_Postfilt,[size(smat,2),1]);
%******REAL PARTS FILTERED******
% get the data from time series back out of the object
smatReFilt= getdatasamples(smatROut,[1:length(smat)]);
smatReFilt=squeeze(smatReFilt(1,1,:));
sFToutReWind= fftshift(fft(smatReFilt.*windowt23_Postfilt,NFFT,1),1);

% get the data from time series back out of the object
dmatReFilt= getdatasamples(dmatROut,[1:length(dmat)]);
dmatReFilt=squeeze(dmatReFilt(1,1,:));
dFToutReWind= fftshift(fft(dmatReFilt.*windowt23_Postfilt,NFFT,1),1);

% *****IMAG PARTS FILTERED ******
% get the data from time series back out of the object
smatImFilt= getdatasamples(smatImOut,[1:length(smat)]);
smatImFilt=squeeze(smatImFilt(1,1,:));
sFToutImWind= fftshift(fft(smatImFilt.*windowt23_Postfilt,NFFT,1),1);

% get the data from time series back out of the object
dmatImFilt= getdatasamples(dmatImOut,[1:length(dmat)]);
dmatImFilt=squeeze(dmatImFilt(1,1,:));
dFToutImWind= fftshift(fft(dmatImFilt.*windowt23_Postfilt,NFFT,1),1);

figure(13)
hold on
title('NRP Real FT of Filtered Signal + window - Bandpass: 10-2000 cm^{-1}');
plot(Faxis,abs(sFToutReWind))
xlim([0 inf]);
xlabel('Wavenumbers');
ylabel('Mag (arb)'); 

figure(14)
hold on
title('RP Real FT of Filtered Signal + window - Bandpass: 10-2000 cm^{-1}');
plot(Faxis,abs(dFToutReWind))
xlim([0 inf]);
xlabel('Wavenumbers');
ylabel('Mag (arb)'); 

figure(15)
hold on
title('NRP Imag FT of Filtered Signal + window - Bandpass: 10-2000 cm^{-1}');
plot(Faxis,abs(sFToutImWind))
xlim([0 inf]);
xlabel('Wavenumbers');
ylabel('Mag (arb)'); 

figure(16)
hold on
title('RP Imag FT of Filtered Signal + window - Bandpass: 10-2000 cm^{-1}');
plot(Faxis,abs(dFToutImWind))
xlim([0 inf]);
xlabel('Wavenumbers');
ylabel('Mag (arb)'); 



end 
     
 end 
