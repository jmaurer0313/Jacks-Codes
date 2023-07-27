% single2DFPGA_RT (same as original single2DFPGA but for the new retiming
% script. 

%constants
c0 = 0.000299792458; % speed of light mm/fs
% set(groot, 'defaultFigureRenderer', 'Painters')
%*****OPTIONS FOR HIGH QUALITY FIGURES*****
% fontname = 'Times New Roman'; 
% set(groot,'DefaultAxesFontName',fontname,'DefaultTextFontName',fontname);
% set(groot,'defaultAxesFontSize',12);
% set(groot,'defaultFigureRenderer', 'Painters');

% base_folder = 'D:\Data\2DFPGA\Sm_avgFolders\20200212_mol2'; % good mono
% good Dimer Path
% C:\Users\Jack\Documents\Marcus Lab\2D_SingleMolecule\Sm_avgFolders\TimeLapses\20191217_mol1_3hr
% good mono path
% 'C:\Users\Jack\Documents\Marcus Lab\2D_SingleMolecule\Sm_avgFolders\TimeLapses\20200311_mol3_2hr';


base_folder = 'C:\Users\Jack\Documents\Marcus Lab\2D_SingleMolecule\Sm_avgFolders\TimeLapses\20200311_mol3_2hr';

% normal SM avg folder name
% name =[ base_folder((length(base_folder)-12):length(base_folder)-5),' ', base_folder((length(base_folder)-3):end)];

% time lapse enabled name
name =[ base_folder((length(base_folder)-12):length(base_folder)-9),' ',base_folder((length(base_folder)-7):length(base_folder)-4), ' ', base_folder((length(base_folder)-2):end)]

% Set and integration time and decide whether or not to adjust the
% integration time during averaging. Set whether you want to INDIVIDUALLY
% retime and rephase the data 
Lapse=0;
intTime=42; 
retiming=1;
rephase=1;
printOpt=0;



% Set the plotOpt to toggle whether or not to see the underlying plots in
% all nested function calls
plotOpt=0;
axLim=60;

% Load the python color map
pythonmap = load('C:\Users\Jack\Documents\Marcus Lab\2D_SingleMolecule\Sm_avgFolders\test_data\pythonmap.mat'); 
pythMap= pythonmap.newmap; 

% generate the proper colormap for matching purposes to the bulk data (this
% is the code to attempt a close match. Above code utilzies the actual RGB
% values from Dylan
% mycolors= jet;
% mymap= mycolors(8:56,:); 

if Lapse
 
[dmat, smat, t21ax, t43ax, window2Dfinal,numPhots] = avg2DFPGA_RT_Lapse(intTime, plotOpt, base_folder,retiming,rephase);

else
[dmat, smat, t21ax, t43ax, window2Dfinal,numPhots] = avg2DFPGA_RT(plotOpt, base_folder,retiming,rephase);

end 
% [dmat, smat, t21ax, t43ax, window2Dfinal] = avg2DFPGA_RT(plotOpt, base_folder);

if ~rephase
% Rough phasing
dmat = exp(1i*(-angle(dmat(1,1))))*dmat;
smat = exp(1i*(-angle(smat(1,1))))*smat;
end

% convert the number of photons to a string which can be used as scientific
% notation
powTen= numel(num2str(numPhots))-1;
redPhot= numPhots/(10^powTen);
photStr= num2str(redPhot);
photLabel = [' N=' photStr(1:5) ,'E',num2str(powTen)]; 

name=[name, photLabel]; 

prefix = 'a';
x_pref = 'RHS3';
y_pref = 'LHS1';



figure(1)
sgtitle(name);
subplot(2,2,1)
contourf(t21ax, t43ax, real(smat),10)
title('NRPreal')
% xlim([0 50]);
% ylim([0 50]);
axis equal

subplot(2,2,2)
contourf(t21ax, t43ax, imag(smat),10)
title('NRPimag')
% xlim([0 50]);
% ylim([0 50]);
axis equal

subplot(2,2,3)
contourf(t21ax, t43ax, real(dmat),10)
title('RPreal')
% xlim([0 50]);
% ylim([0 50]);
axis equal

subplot(2,2,4)
contourf(t21ax, t43ax, imag(dmat),10)
title('RPimag')
% xlim([0 50]);
% ylim([0 50]);
axis equal


% Window the thing
% [wX, wY] = meshgrid(t21ax,t43ax);
% window = delayedGaussian(sqrt(wX.^2 + wY.^2), 50, 20);
% figure(2)
% contourf(t21ax, t43ax, window)
% axis equal

% Plot windowed spectra
dmatW = dmat.*window2Dfinal;
smatW = smat.*window2Dfinal;

figure(3)
sgtitle(name);
subplot(2,2,1)
contourf(t21ax(1:axLim), t43ax(1:axLim), real(smatW((1:axLim),(1:axLim))))
title('Retimed NRPreal Windowed')
axis equal

subplot(2,2,2)
contourf(t21ax(1:axLim), t43ax(1:axLim), imag(smatW((1:axLim),(1:axLim))))
title('Retimed NRPimag Windowed')
% axis([0 50 0 50]);
% ylim([0 50]);
axis equal

subplot(2,2,3)
contourf(t21ax(1:axLim), t43ax(1:axLim), real(dmatW((1:axLim),(1:axLim))))
% xlim([0 50]);
% ylim([0 50]);
title('Retimed RPreal Windowed')
axis equal

subplot(2,2,4)
contourf(t21ax(1:axLim), t43ax(1:axLim), imag(dmatW((1:axLim),(1:axLim))))
% xlim([0 50]);
% ylim([0 50]);
title('Retimed RPimag Windowed')
axis equal

T = t21ax(2)-t21ax(1); % Uses interpolated time domain axis to calculate
NFFT = length(t21ax);
lam_mono = 550; % monochromater wavelength in nanometers
Fmono = 1e7 / lam_mono; % convert to wavenumbers
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
% up = ones(size(dmatW));



freqDom=fft2(dmatW, NFFT,NFFT);
maxScal= max(max(abs(real(freqDom))));
normFreq= freqDom/maxScal; 

fig4 = figure(4);
clf(4)
fig4pos = get(fig4, 'Position');
set(fig4, 'PaperUnits','centimeters')
set(fig4, 'PaperSize', [18, 18/fig4pos(3)*fig4pos(4)])
% figure(4)
contourf(1e-3*(Faxis+Fmono), 1e-3*(Faxis+Fmono), fliplr(real(fftshift(normFreq))),15)
axis equal
colormap(pythMap);
colorbar;
caxis([-1 1]);
if ~retiming
    title(['NOT-Retimed RP(\omega) real ' name]) 
else
    title(['Retimed RP(\omega) real ' name]) 
end 
title(['Retimed RP(\omega) real ' name])
xlabel('Wavenumbers (x10^3 cm^{-1})');
ylabel('Wavenumbers (x10^3 cm^{-1})');
xlim([17 20.5]); 
ylim([17 20.5]); 
if printOpt
print -dpdf -r300 -bestfit RPreal.pdf
end


freqDom=fft2(smatW, NFFT,NFFT);
maxScal= max(max(abs(real(freqDom))));
normFreq= freqDom/maxScal; 
% figure(5)
fig5 = figure(5);
clf(5)
fig5pos = get(fig5, 'Position');
set(fig5, 'PaperUnits','centimeters')
set(fig5, 'PaperSize', [18, 18/fig5pos(3)*fig5pos(4)])
contourf(1e-3*(Faxis+Fmono), 1e-3*(Faxis+Fmono), (real(fftshift(normFreq))),15)
axis equal
colormap(pythMap);
colorbar;
caxis([-1 1]);
if ~retiming
    title(['NOT-Retimed NRP(\omega) real ' name]) 
else
    title(['Retimed NRP(\omega) real ' name]) 
end 
% title(['Retimed NRP(\omega) real ' name])
xlabel('Wavenumbers (x10^3 cm^{-1})');
ylabel('Wavenumbers (x10^3 cm^{-1})');
xlim([17 20.5]); 
ylim([17 20.5]); 
if printOpt
print -dpdf -r300 -bestfit NRPreal.pdf
end



freqDom=fft2(dmatW, NFFT,NFFT);
maxScal= max(max(abs(imag(freqDom))));
normFreq= freqDom/maxScal; 
% figure(6)
fig6 = figure(6);
clf(6)
fig6pos = get(fig6, 'Position');
set(fig6, 'PaperUnits','centimeters')
set(fig6, 'PaperSize', [18, 18/fig6pos(3)*fig6pos(4)])
contourf(1e-3*(Faxis+Fmono), 1e-3*(Faxis+Fmono),(fliplr(imag(fftshift(normFreq)))),15)
axis equal
colormap(pythMap);
colorbar;
caxis([-1 1]);
if ~retiming
    title(['NOT-Retimed RP(\omega) Imag ' name]) 
else
    title(['Retimed RP(\omega) Imag ' name]) 
end 
% title(['Retimed RP(\omega) Imag ' name])
xlabel('Wavenumbers (x10^3 cm^{-1})');
ylabel('Wavenumbers (x10^3 cm^{-1})');
xlim([17 20.5]); 
ylim([17 20.5]); 
if printOpt
print -dpdf -r300 -bestfit RPimag.pdf
end


freqDom=fft2(smatW, NFFT,NFFT);
maxScal= max(max(abs(imag(freqDom))));
normFreq= freqDom/maxScal;
% figure(7)
fig7 = figure(7);
clf(7)
fig7pos = get(fig7, 'Position');
set(fig7, 'PaperUnits','centimeters')
set(fig7, 'PaperSize', [18, 18/fig7pos(3)*fig7pos(4)])
contourf(1e-3*(Faxis+Fmono), 1e-3*(Faxis+Fmono), (imag(fftshift(normFreq))),15)
axis equal
colormap(pythMap);
colorbar;
caxis([-1 1]);
if ~retiming
    title(['NOT-Retimed NRP(\omega) Imag ' name]) 
else
    title(['Retimed NRP(\omega) Imag ' name]) 
end 
% title(['Retimed NRP(\omega) Imag '  name])
xlabel('Wavenumbers (x10^3 cm^{-1})');
ylabel('Wavenumbers (x10^3 cm^{-1})');
xlim([17 20.5]); 
ylim([17 20.5]); 
if printOpt
print -dpdf -r300 -bestfit NRPimag.pdf
end


% figure(9)
fig9 = figure(9);
clf(9)
fig9pos = get(fig9, 'Position');
set(fig9, 'PaperUnits','centimeters')
set(fig9, 'PaperSize', [18, 18/fig9pos(3)*fig9pos(4)])
contourf(1e-3*(Faxis+Fmono), 1e-3*(Faxis+Fmono), real(fftshift(fft2(smatW, NFFT,NFFT))+ fliplr(fftshift(fft2(dmatW, NFFT,NFFT))))/max(max(real(fftshift(fft2(smatW, NFFT,NFFT))+ fliplr(fftshift(fft2(dmatW, NFFT,NFFT)))))),15)
axis equal
colormap(pythMap);
colorbar;
caxis([-1 1]);
if ~retiming
    title(['NOT-Retimed Total real ' name]) 
else
    title(['Retimed Total real ' name]) 
end 
% title(['Retimed Total real ' name])
xlabel('Wavenumbers (x10^3 cm^{-1})');
ylabel('Wavenumbers (x10^3 cm^{-1})');
xlim([17 20.5]); 
ylim([17 20.5]);
if printOpt
print -dpdf -r300 -bestfit TotReal.pdf
end


% figure(10)
fig10 = figure(10);
clf(10)
fig10pos = get(fig10, 'Position');
set(fig10, 'PaperUnits','centimeters')
set(fig10, 'PaperSize', [18, 18/fig10pos(3)*fig10pos(4)])
contourf(1e-3*(Faxis+Fmono), 1e-3*(Faxis+Fmono), imag(fftshift(fft2(smatW, NFFT,NFFT))+ fliplr(fftshift(fft2(dmatW, NFFT,NFFT))))/max(max(imag(fftshift(fft2(smatW, NFFT,NFFT))+ fliplr(fftshift(fft2(dmatW, NFFT,NFFT)))))),15)
axis equal
colormap(pythMap);
colorbar;
caxis([-1 1]);
if ~retiming
    title(['NOT-Retimed Total Imag ' name]) 
else
    title(['Retimed Total Imag ' name]) 
end 
% title(['Retimed Total Imag ' name])
xlabel('Wavenumbers (x10^3 cm^{-1})');
ylabel('Wavenumbers (x10^3 cm^{-1})');
xlim([17 20.5]); 
ylim([17 20.5]);
if printOpt
print -dpdf -r300 -bestfit TotImag.pdf
end

% plot(0:31, trapz(abs(down2)))
% % 
% % % Retime by zero padding in frequency first
% % nret = (length(dmatW));
% % totalT = T*(nret-1);
% % required_resolution = 0.05; % in fs for retiming
% % ninv = 2^nextpow2(totalT/required_resolution); % zero_padding used for inverse FFT back to time.
% % Tret = totalT/(ninv-1); %actual retiming increment
% % dmatinterp = interpft(interpft(dmatW, ninv, 1),ninv,2);
% % smatinterp = interpft(interpft(smatW, ninv, 1),ninv,2);
% % xinterp = linspace(t21ax(1),t21ax(end),ninv);
% % yinterp = linspace(t43ax(1),t43ax(end),ninv);
% % % Rewindow
% % [wwX, wwY] = meshgrid(xinterp,yinterp);
% % windoww = delayedGaussian(sqrt(wwX.^2 + wwY.^2), 40, 20);
% % down1 = downsample(dmatinterp.*windoww, 64,0);
% % down2 = downsample(down1.', 64, 0).';
% % sdown1 = downsample(smatinterp.*windoww, 64,0);
% % sdown2 = downsample(sdown1.', 64, 0).';
% % 
% % down2 = exp(1i*(-angle(down2(1,1))))*down2;
% % sdown2 = exp(1i*(-angle(sdown2(1,1))))*sdown2;

% figure(7)
% contourf(real(down2),15)
% axis equal
% figure(10)
% contourf(real(sdown2),15)
% axis equal
% 
% 
% figure(8)
% contourf(1e-3*(Faxis+Fmono), 1e-3*(Faxis+Fmono), fliplr(real(fftshift(fft2(down2, NFFT,NFFT)))),15)
% axis equal
% title('real FFT dmat retimed')
% 
% figure(9)
% contourf(1e-3*(Faxis+Fmono), 1e-3*(Faxis+Fmono), real(fftshift(fft2(smat, NFFT,NFFT))+ fliplr(fftshift(fft2(dmat, NFFT,NFFT)))),10)
% axis equal
% title('Total real')
% xlabel('Wavenumbers (x10^3 cm^{-1})');
% ylabel('Wavenumbers (x10^3 cm^{-1})');
% xlim([17 20.5]); 
% ylim([17 20.5]); 
% % plot(0:31, trapz(abs(down2)))
% 
% figure(6)
% contourf(1e-3*(Faxis+Fmono), 1e-3*(Faxis+Fmono), imag(fftshift(fft2(smat, NFFT,NFFT))+ fliplr(fftshift(fft2(dmat, NFFT,NFFT)))),10)
% axis equal
% title('Total imag')
% xlabel('Wavenumbers (x10^3 cm^{-1})');
% ylabel('Wavenumbers (x10^3 cm^{-1})');
% xlim([17 20.5]); 
% ylim([17 20.5]); 
% % plot(0:31, trapz(abs(down2)))
% 
% figure(11)
% contourf(1e-3*(Faxis+Fmono), 1e-3*(Faxis+Fmono), fliplr(real(fftshift(fft2(dmatW, NFFT,NFFT)))),15)
% axis equal
% title('RP(\omega) real')
% xlabel('Wavenumbers (x10^3 cm^{-1})');
% ylabel('Wavenumbers (x10^3 cm^{-1})');
% xlim([17 20.5]); 
% ylim([17 20.5]); 
% 
% figure(12)
% contourf(1e-3*(Faxis+Fmono), 1e-3*(Faxis+Fmono), fliplr(imag(fftshift(fft2(dmatW, NFFT,NFFT)))),15)
% axis equal
% title('RP(\omega) imag')
% xlabel('Wavenumbers (x10^3 cm^{-1})');
% ylabel('Wavenumbers (x10^3 cm^{-1})');
% xlim([17 20.5]); 
% ylim([17 20.5]); 
% 
% figure(13)
% contourf(1e-3*(Faxis+Fmono), 1e-3*(Faxis+Fmono), (real(fftshift(fft2(smatW, NFFT,NFFT)))),15)
% axis equal
% title('NRP(\omega) real')
% xlabel('Wavenumbers (x10^3 cm^{-1})');
% ylabel('Wavenumbers (x10^3 cm^{-1})');
% xlim([17 20.5]); 
% ylim([17 20.5]); 
% 
% figure(14)
% contourf(1e-3*(Faxis+Fmono), 1e-3*(Faxis+Fmono), (imag(fftshift(fft2(smatW, NFFT,NFFT)))),15)
% axis equal
% title('NRP(\omega) imag')
% xlabel('Wavenumbers (x10^3 cm^{-1})');
% ylabel('Wavenumbers (x10^3 cm^{-1})');
% xlim([17 20.5]); 
% ylim([17 20.5]); 
% 
% 
