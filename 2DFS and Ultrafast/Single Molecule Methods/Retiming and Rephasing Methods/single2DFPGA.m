%constants
c0 = 0.000299792458; % speed of light mm/fs
% set(groot, 'defaultFigureRenderer', 'Painters')


% base_folder = 'D:\Data\2DFPGA\Sm_avgFolders\20200212_mol2'; % good mono
base_folder = 'C:\Users\Jack\Documents\Marcus Lab\2D_SingleMolecule\Sm_avgFolders\TimeLapses\20200311_mol3_2hr';
name =[ base_folder((length(base_folder)-12):length(base_folder)-5),' ', base_folder((length(base_folder)-3):end)];

% Set and integration time and decide whether or not to adjust the
% integration time during averaging. 
Lapse=0;
intTime=30; 

pythonmap = load('C:\Users\Jack\Documents\Marcus Lab\2D_SingleMolecule\Sm_avgFolders\test_data\pythonmap.mat'); 
pythMap= pythonmap.newmap; 

if Lapse
[dmat, smat, t21ax, t43ax] = avg2DFPGA_Lapse(intTime, base_folder);
else
[dmat, smat, t21ax, t43ax] = avg2DFPGA(base_folder);
end 

retime_from_file = 0;
no_retime = 0; %always off for now

prefix = 'a';
x_pref = 'RHS3';
y_pref = 'LHS1';


dropmat = [1,0]; % functionality placeholder, any value other than [1,0] will break code.

% t21mat = dlmread(fullfile(base_folder, ['t21' prefix])); %figure out how to handle case where nothing need be dropped.
t21ax(dropmat(1))=[]; % Does anything 

% t21ax = t21mat(:,1);
% t43ax = dlmread(fullfile(base_foldr, ['t43' prefix]));

% Retiming logic
if retime_from_file&&~no_retime
    t21_0 = fitAC(fullfile(base_folder, x_pref), 0.5, 0.5); 
    t43_0 = fitAC(fullfile(base_folder, y_pref), 0.5, 0.5);
elseif no_retime
    t21_0=t21ax(1); % assume correct time zero
    t43_0=t43ax(1); % assume correct time zero
else
    t21_0=0; % put in correct time here otherwise make sure you avoid this case
    t43_0=0;
end

% Calculate Absolute Time axes
t21ax = -1*(t21ax - t21_0);
t43ax = (t43ax - t43_0);

% % Load Data Matrices
% dreal = dlmread(fullfile(base_folder, ['x3' prefix]));
% dimag = dlmread(fullfile(base_folder, ['y3' prefix]));
% sreal = dlmread(fullfile(base_folder, ['x13' prefix]));
% simag = dlmread(fullfile(base_folder, ['y13' prefix]));

dmat(dropmat(1),:)=[];
smat(dropmat(1),:)=[];


% dmat = complex(dreal,dimag)';
% smat = complex(sreal,simag)';
dAngle=angle(dmat(1,1));
sAngle=angle(smat(1,1));

% Rough phasing
dmat = exp(1i*(-angle(dmat(1,1))))*dmat;
smat = exp(1i*(-angle(smat(1,1))))*smat;

figure(1)
sgtitle(name);
subplot(2,2,1)
contourf(t21ax, t43ax, real(smat),10)
title('NRPreal')
axis equal
subplot(2,2,2)
contourf(t21ax, t43ax, imag(smat),10)
title('NRPimag')
axis equal
subplot(2,2,3)
contourf(t21ax, t43ax, real(dmat),10)
title('RPreal')
axis equal
subplot(2,2,4)
contourf(t21ax, t43ax, imag(dmat),10)
title('RPimag')
axis equal


% Window the thing
[wX, wY] = meshgrid(t21ax,t43ax);
window = delayedGaussian(sqrt(wX.^2 + wY.^2), 40, 32);
figure(2)
contourf(t21ax, t43ax, window)
axis equal

% Plot windowed spectra
dmatW = dmat.*window;
smatW = smat.*window;

figure(3)
sgtitle(name)
subplot(2,2,1)
contourf(t21ax, t43ax, real(smatW),10)
title('NRPreal Windowed')
axis equal
subplot(2,2,2)
contourf(t21ax, t43ax, imag(smatW),10)
title('NRPimag Windowed')
axis equal
subplot(2,2,3)
contourf(t21ax, t43ax, real(dmatW),10)
title('RPreal Windowed')
axis equal
subplot(2,2,4)
contourf(t21ax, t43ax, imag(dmatW),10)
title('RPimag Windowed')
axis equal

T = mean(t21ax(2:end) - t21ax(1:end-1)); % Uses stage positions directly, could replace by 0.0004/c0
NFFT = 16;
lam_mono = 550; % monochromater wavelength in nanometers
Fmono = 1e7 / lam_mono;
% T = 7/6/1e3*2/c0;
Fs = 10/(c0*T);
Faxis = Fs/2*linspace(-1,1,NFFT);

% % Retime by zero padding in frequency first
% nret = 2^nextpow2(length(dmat));
% totalT = T*nret;
% required_resolution = 0.05; % in fs for retiming
% ninv = 2^nextpow2(nret*T/required_resolution); % zero_padding used for inverse FFT back to time.
% Tret = totalT/ninv; %actual retiming increment
% dmatTemp = ifftshift(fliplr(fftshift(fft2(dmat, nret, nret)))); % preparing for zero padding
% [~ , I] = max(dmatTemp(:));  % Find spectrum peak
% [I_row, I_col] = ind2sub(size(dmatTemp),I);  % Indeces of spectrum peak
% row_insert = mod(round(0.5*(abs(nret-2*I_row)))+I_row, nret);
% col_insert = mod(round(0.5*(abs(nret-2*I_col)))+I_col, nret);
% dmatTempP = zeros(ninv);
% dmatTempP(1:row_insert,1:col_insert) = dmatTemp(1:row_insert,1:col_insert);
% dmatTempP(end - row_insert + 1:end,1:col_insert) = dmatTemp(end-row_insert + 1:end,1:col_insert);
% dmatTempP(1:row_insert,end - col_insert + 1: end) = dmatTemp(1:row_insert,end - col_insert + 1: end);
% dmatTempP(end - row_insert + 1:end,end - col_insert + 1: end) = dmatTemp(end - row_insert + 1:end,end - col_insert + 1: end);
% dmatret = ifft(dmatTempP);
% 
% 
% Fsource = (10/(c0*T))*linspace(-1,1,nret);
% Fsourceshift = fftshift(Fsource);
% peakF = 5700; % spectrum peak minus monoch)
up = ones(size(dmatW));


freqDomRP=fft2(dmatW, NFFT,NFFT);
maxScalRP= max(max(abs(real(freqDomRP))));
normFreqRP= freqDomRP/maxScalRP; 
figure(4)
contourf(1e-3*(Faxis+Fmono), 1e-3*(Faxis+Fmono), fliplr(real(fftshift(normFreqRP))),15)
axis equal
colormap(pythMap);
colorbar;
caxis([-1 1])
title(['RP(\omega) real ' name])
xlabel('Wavenumbers (x10^3 cm^{-1})');
ylabel('Wavenumbers (x10^3 cm^{-1})');
xlim([17 20.5]); 
ylim([17 20.5]);
% alpha color
% alpha scaled


freqDomNRP=fft2(smatW, NFFT,NFFT);
maxScalNRP= max(max(abs(real(freqDomNRP))));
normFreqNRP= freqDomNRP/maxScalNRP; 
figure(5)
contourf(1e-3*(Faxis+Fmono), 1e-3*(Faxis+Fmono), (real(fftshift(normFreqNRP))),15)
axis equal
colormap(pythMap);
colorbar;
caxis([-1 1]);
title(['NRP(\omega) real ' name])
xlabel('Wavenumbers (x10^3 cm^{-1})');
ylabel('Wavenumbers (x10^3 cm^{-1})');
xlim([17 20.5]); 
ylim([17 20.5]);

figure(6)
contourf(1e-3*(Faxis+Fmono), 1e-3*(Faxis+Fmono), (imag(fftshift(normFreqNRP))),15)
axis equal
colormap(pythMap);
colorbar;
caxis([-1 1]);
title(['NRP(\omega) Imag ' name]);
xlabel('Wavenumbers (x10^3 cm^{-1})');
ylabel('Wavenumbers (x10^3 cm^{-1})');
xlim([17 20.5]); 
ylim([17 20.5]);

figure(7)
contourf(1e-3*(Faxis+Fmono), 1e-3*(Faxis+Fmono), fliplr((imag(fftshift(normFreqRP)))),15)
axis equal
colormap(pythMap);
colorbar;
caxis([-1 1]);
title(['RP(\omega) Imag ' name])
xlabel('Wavenumbers (x10^3 cm^{-1})');
ylabel('Wavenumbers (x10^3 cm^{-1})');
xlim([17 20.5]); 
ylim([17 20.5]);


% % Retime by zero padding in frequency first
% nret = (length(dmatW));
% totalT = T*(nret-1);
% required_resolution = 0.05; % in fs for retiming
% ninv = 2^nextpow2(totalT/required_resolution); % zero_padding used for inverse FFT back to time.
% Tret = totalT/(ninv-1); %actual retiming increment
% dmatinterp = interpft(interpft(dmatW, ninv, 1),ninv,2);
% smatinterp = interpft(interpft(smatW, ninv, 1),ninv,2);
% xinterp = linspace(t21ax(1),t21ax(end),ninv);
% yinterp = linspace(t43ax(1),t43ax(end),ninv);
% % Rewindow
% [wwX, wwY] = meshgrid(xinterp,yinterp);
% windoww = delayedGaussian(sqrt(wwX.^2 + wwY.^2), 40, 20);
% down1 = downsample(dmatinterp.*windoww, 64,0);
% down2 = downsample(down1.', 64, 0).';
% sdown1 = downsample(smatinterp.*windoww, 64,0);
% sdown2 = downsample(sdown1.', 64, 0).';
% 
% down2 = exp(1i*(-angle(down2(1,1))))*down2;
% sdown2 = exp(1i*(-angle(sdown2(1,1))))*sdown2;
% 
% figure(8)
% contourf(real(down2),15)
% axis equal
% figure(10)
% contourf(real(sdown2),15)
% axis equal
% 
% 
% figure(9)
% contourf(1e-3*(Faxis+Fmono), 1e-3*(Faxis+Fmono), fliplr(real(fftshift(fft2(down2, NFFT,NFFT)))),15)
% axis equal
% title('real FFT dmat retimed')
% 
% figure(10)
% contourf(1e-3*(Faxis+Fmono), 1e-3*(Faxis+Fmono), real(fftshift(fft2(smat, NFFT,NFFT))+ fliplr(fftshift(fft2(dmat, NFFT,NFFT)))),10)
% axis equal
% title(['Total real ' name])
% xlabel('Wavenumbers (x10^3 cm^{-1})');
% ylabel('Wavenumbers (x10^3 cm^{-1})');
% xlim([17 20.5]); 
% ylim([17 20.5]); 
% % plot(0:31, trapz(abs(down2)))
% 
% figure(11)
% contourf(1e-3*(Faxis+Fmono), 1e-3*(Faxis+Fmono), imag(fftshift(fft2(smat, NFFT,NFFT))+ fliplr(fftshift(fft2(dmat, NFFT,NFFT)))),10)
% axis equal
% title(['Total imag ' name])
% xlabel('Wavenumbers (x10^3 cm^{-1})');
% ylabel('Wavenumbers (x10^3 cm^{-1})');
% xlim([17 20.5]); 
% ylim([17 20.5]); 
% % plot(0:31, trapz(abs(down2)))
% 
% figure(12)
% contourf(1e-3*(Faxis+Fmono), 1e-3*(Faxis+Fmono), fliplr(real(fftshift(fft2(dmatW, NFFT,NFFT)))),15)
% axis equal
% title(['RP(\omega) real ' name])
% xlabel('Wavenumbers (x10^3 cm^{-1})');
% ylabel('Wavenumbers (x10^3 cm^{-1})');
% xlim([17 20.5]); 
% ylim([17 20.5]); 
% 
% figure(13)
% contourf(1e-3*(Faxis+Fmono), 1e-3*(Faxis+Fmono), fliplr(imag(fftshift(fft2(dmatW, NFFT,NFFT)))),15)
% axis equal
% title(['RP(\omega) imag ' name])
% xlabel('Wavenumbers (x10^3 cm^{-1})');
% ylabel('Wavenumbers (x10^3 cm^{-1})');
% xlim([17 20.5]); 
% ylim([17 20.5]); 
% 
% figure(14)
% contourf(1e-3*(Faxis+Fmono), 1e-3*(Faxis+Fmono), (real(fftshift(fft2(smatW, NFFT,NFFT)))),15)
% axis equal
% title(['NRP(\omega) real ' name])
% xlabel('Wavenumbers (x10^3 cm^{-1})');
% ylabel('Wavenumbers (x10^3 cm^{-1})');
% xlim([17 20.5]); 
% ylim([17 20.5]); 
% 
% figure(15)
% contourf(1e-3*(Faxis+Fmono), 1e-3*(Faxis+Fmono), (imag(fftshift(fft2(smatW, NFFT,NFFT)))),15)
% axis equal
% title(['NRP(\omega) imag ' name])
% xlabel('Wavenumbers (x10^3 cm^{-1})');
% ylabel('Wavenumbers (x10^3 cm^{-1})');
% xlim([17 20.5]); 
% ylim([17 20.5]); 


