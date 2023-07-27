%set up edges, "centers" is basically the phase axis of 2n+1 odd numbers.
% close all
% JCP size: one column = 8.5cm, two column = 17cm width
fontname = 'Times New Roman'; 
set(groot,'DefaultAxesFontName',fontname,'DefaultTextFontName',fontname);
set(groot,'defaultAxesFontSize',12)
nboot = 2000;
timerange = [-100 100];
% sntimerange = [-55,55];
% snindex = time > sntimerange(1) & time < sntimerange(2);
waverange = [1.5 2.18]*1e4;
% if exist('fitobject', 'var')
%     normfact = fitobject.a1;
% else
%     normfact = 0.21668;
% end
normfact = 1;

n = 14; % 2n+1 total histogram bins, produces one bin centered at zero, necessary.
N = 2*n +1;
edges = linspace(-1, 1, 2*n +2);
edgesl = edges(1:end-1);
edgesr = edges(2:end);
centers = 0.5*(edgesl + edgesr);
% load actual histogram data
load('C:\Users\Jack\Documents\Marcus Lab\2D_SingleMolecule\Auto\20200310-192133-auto\auto.mat');
% outpath = '\ProcessedData\20190207-FPGA-20200306\17k-16.5\';
% outprefix = '205002-17k-16.5';

% Changed Amrs save path so it doesnt look for the folder above and just
% write to the individual files
outpath = '\ProcessedData';
outprefix = '';
% figout = ['D:\Users\amroa\Dropbox\FPGA SNR analysis\Matlab\1DFPGAscan\ProcessedData\Molecule\20191007-201138\' outprefix '\Fig'];

time = time*1000 ; %ps to fs
% - 3.1200
% time = time*1000 - 5.62675126822481; %ps to fs
pL = length(photons);
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

%% Histogram plotting section
[tmin, Imin] = min(abs(time));
% Imin = floor(length(time)*.5); % Which stage position to plot?
h = photons{Imin};
h = (h - 180)/180;
[y, ~] = histcounts(h, edges);
% Calculate zero phase
f10 = sum(1/length(photons{Imin})*exp(-1i*pi*(photons{Imin}-180)/180));
f20 = sum(1/length(photons{Imin})*exp(-2i*pi*(photons{Imin}-180)/180));

% histogfig = figure;
fig1 = figure(1);
set(fig1,'renderer','Painters', 'Units', 'centimeters')
bar(centers, y, 'FaceColor',[0.8500 0.3250 0.0980],'FaceAlpha', 0.8, 'EdgeColor',[0 0 0],'LineWidth',1.25)
ax1 = gca;
title('Photon Detection-Phase Histogram','FontSize', 14);
xlabel('detection phase (\pi rad)', 'FontSize', 14)
ylabel('counts', 'FontSize', 14)
ax1.XLimSpec = 'Tight';
ax1.LineWidth = 1.5;
ax1.TickLength = [.02 .02];
fig1pos = get(fig1, 'Position');
set(fig1, 'PaperUnits','centimeters')
set(fig1, 'PaperSize', [8.5, 8.5/fig1pos(3)*fig1pos(4)])

print -dpdf -r300 -bestfit FPGA-hist.pdf
% export_fig test.pdf -q101
% print(fig1, '-dpdf', 'D:\Users\amroa\Desktop\testp.pdf'); 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%% FFT time zero to show harmonics

ys = ifftshift(y);

%pad and transform
padding = 4;
M = N*padding;
beta = (-(M/2):(M/2-1))/padding;
transformed = fft(repmat(ys, 1, padding))/M;
transformed(1) = 0.5*transformed(1);
first = transformed(padding+1);
Y = fftshift(transformed);
reconstructed = 2*(abs(first)*cos(centers*pi + angle(first))+abs(transformed(1)));
% plot(centers, reconstructed, centers, y);

% fftfig = figure;
fig2 = figure(2);
ax2 = gca;
plot(beta, real(Y), beta, imag(Y), beta, abs(Y), '--k', 'LineWidth',1.25)
title 'Fourier transform of tiled histogram';
xlabel('harmonic', 'FontSize', 14)
ylabel('counts', 'FontSize', 14)
legend('Re','Im','Abs', 'FontSize', 12)
ax2.XLimSpec = 'Tight';
ax2.LineWidth = 1.5;
ax2.TickLength = [.02 .02];
fig2pos = get(fig2, 'Position');
set(fig2, 'PaperUnits','centimeters')
set(fig2, 'PaperSize', [8.5, 8.5/fig2pos(3)*fig2pos(4)])

print -dpdf -r300 -bestfit FPGA-FFT.pdf


%% analyze many
yz = zeros(pL, length(centers));
for i = 1:pL
%     [yi, ~] = histcounts(photons{i}, edges);
    [yi, ~] = histcounts((photons{i}-180)/180, edges);
    yz(i, :) = yi;
end

% bootstrapped Tiemo™ replacement
Error = zeros(pL,4);
f1tiemo = zeros(pL,1);
f2tiemo = zeros(pL,1);
fboots = zeros(pL, nboot); % for frequency domain error bars

parfor i = 1:pL  % here's where you can turn off parallel
    [xerr, yerr, abserr, phierr, rawbootstrap, fftbootstrap] = bootstrapFourierStatistics(nboot, 1, pi*(photons{i}-180 + angle(f10)*180/pi)/180);
    f1tiemo(i) = sum(1/length(photons{i})*exp(-1i*pi*(photons{i}-180 + angle(f10)*180/pi)/180));
    f2tiemo(i) = sum(1/length(photons{i})*exp(-2i*pi*(photons{i}-180 + angle(f20)*180/pi)/180));
    fboots(i, :) = fftbootstrap(:,2).';
    Error(i,:) = [xerr(2) yerr(2) abserr(2) phierr(2)];
end

transformedz = fft(repmat(ifftshift(yz, 2), 1, padding),[],2)/M;
% step axis
steps = (0:499)*50;
c0 = 299.792458; %nm/fs
% time = steps*2/c0;
steps = steps/1e6;
dc = transformedz(:, 0*padding +1);
% f1 = transformedz(:, 1*padding +1);
f1 = f1tiemo;
f2 = f2tiemo;
% f2 = transformedz(:, 2*padding +1);
f3 = transformedz(:, 3*padding +1);

f1normphased = f1/normfact;

ErrorMeans = NaN(2,4);
Errorsmoothed = movmean(Error, 3, 1);
ErrorMeans(1,:) = Errorsmoothed(Imin, :);
ErrorMeans(2,:) = ErrorMeans(1,:)./normfact;
errbars = NaN(size(f1));

% time = (steps)*2000000/c0-177/2;
% time = linspace(-1,1,199)*99;

fig3 = figure(3);
clf(3)
fig3pos = get(fig3, 'Position');
set(fig3, 'PaperUnits','centimeters')
set(fig3, 'PaperSize', [18, 18/fig3pos(3)*fig3pos(4)])
% plot(time, (real(f1)), time, imag(f1), time, abs(f1), '--')
% errorbar(time, real(f1), Error(:,1), 'LineWidth', 1.25, 'color', [0.8500 0.3250 0.0980])
% hold on
% errorbar(time, imag(f1), Error(:,2), 'LineWidth', 1.25, 'color', [0 0.4470 0.7410])
% errorbar(time, abs(f1), Error(:,3), '--', 'LineWidth', 1.25, 'color', [0.4660 0.6740 0.1880])
% hold off
errbars(Imin) = ErrorMeans(1,1);
errorbar(time, (real(f1)), errbars, '-o', 'LineWidth', 1.25, 'color', [0.8500 0.3250 0.0980])
hold on
errbars(Imin) = ErrorMeans(1,2);
errorbar(time, -1.*imag(f1), errbars, '-o', 'LineWidth', 1.25, 'color', [0 0.4470 0.7410])
errbars(Imin) = NaN;
errbars(Imin+2) = ErrorMeans(1,3);
errorbar(time, (abs(f1)), errbars, '--', 'LineWidth', 1.25, 'color', [0.4660 0.6740 0.1880])
hold off
ax3 = gca;
title 'First Harmonic';
xlabel('delay (fs)')
ylabel('f1 (au)')
lgd3 = legend('Re','Im','Abs');
lgd3.Units = 'centimeters';
lgd3.FontSize = 10;
ax3.XLimSpec = 'Tight';
ax3.XLim = timerange;
ax3.LineWidth = 1.5;
ax3.TickLength = [.02 .02];

dim = [.2 .56 .3 .3];
str = ['S/N X = ' num2str(max(real(f1))/ErrorMeans(1,1)) newline 'S/N Y = ' num2str(max(real(f1))/ErrorMeans(1,2))];
annotation('textbox',dim,'String',str,'FitBoxToText','on', 'Units', 'normalized');

print -dpdf -r300 -bestfit FPGA-f1-err.pdf


%% All error time plot
fig4 = figure(4);
clf(4)
fig4pos = get(fig4, 'Position');
set(fig4, 'PaperUnits','centimeters')
set(fig4, 'PaperSize', [20, 20/fig4pos(3)*fig4pos(4)])
plot(time, (real(f1)), time, imag(f1), time, abs(f1), '--')
errorbar(time, real(f1), Error(:,1), '-o', 'LineWidth', 1.25, 'color', [0.8500 0.3250 0.0980])
hold on
errorbar(time, imag(f1), Error(:,2), '-o', 'LineWidth', 1.25, 'color', [0 0.4470 0.7410])
errorbar(time, abs(f1), Error(:,3), 'o--', 'LineWidth', 1.25, 'color', [0.4660 0.6740 0.1880])
hold off
% errbars(Imin) = ErrorMeans(1,1);
% errorbar(time, (real(f1)), errbars, 'LineWidth', 1.25, 'color', [0.8500 0.3250 0.0980])
% hold on
% errbars(Imin) = ErrorMeans(1,2);
% errorbar(time, imag(f1), errbars, 'LineWidth', 1.25, 'color', [0 0.4470 0.7410])
% errbars(Imin) = NaN;
% errbars(Imin+2) = ErrorMeans(1,3);
% errorbar(time, (abs(f1)), errbars, '--', 'LineWidth', 1.25, 'color', [0.4660 0.6740 0.1880])
hold off
ax4 = gca;
title 'First Harmonic';
xlabel('delay (fs)')
ylabel('f1 (au)')
lgd4 = legend('Re','Im','Abs');
lgd4.Units = 'centimeters';
lgd4.FontSize = 10;
ax4.XLimSpec = 'Tight';
ax4.XLim = timerange;
ax4.LineWidth = 1.5;
ax4.TickLength = [.02 .02];

print -dpdf -r300 -bestfit FPGA-f1-fullerr.pdf


%% Plot F2
fig5 = figure(5);
fig5pos = get(fig5, 'Position');
set(fig5, 'PaperUnits','centimeters')
set(fig5, 'PaperSize', [12, 12/fig5pos(3)*fig5pos(4)])
% plot(time, (real(f2)), time, imag(f2), time, abs(f2), '--')
plot(time, real(f2), 'LineWidth', 1.25, 'color', [0.8500 0.3250 0.0980])
hold on
plot(time, imag(f2), 'LineWidth', 1.25, 'color', [0 0.4470 0.7410])
plot(time, abs(f2), '--', 'LineWidth', 1.25, 'color', [0.4660 0.6740 0.1880])
hold off
% errbars(Imin) = ErrorMeans(1,1);
% errorbar(time, (real(f1)), errbars, 'LineWidth', 1.25, 'color', [0.8500 0.3250 0.0980])
% hold on
% errbars(Imin) = ErrorMeans(1,2);
% errorbar(time, imag(f1), errbars, 'LineWidth', 1.25, 'color', [0 0.4470 0.7410])
% errbars(Imin) = NaN;
% errbars(Imin+2) = ErrorMeans(1,3);
% errorbar(time, (abs(f1)), errbars, '--', 'LineWidth', 1.25, 'color', [0.4660 0.6740 0.1880])
% hold off
ax5 = gca;
title 'Second Harmonic';
xlabel('delay (fs)')
ylabel('f2 (au)')
lgd5 = legend('Re','Im','Abs');
lgd5.Units = 'centimeters';
lgd5.FontSize = 10;
ax5.XLimSpec = 'Tight';
ax5.XLim = timerange;
ax5.LineWidth = 1.5;
ax5.TickLength = [.02 .02];

print -dpdf -r300 -bestfit FPGA-f2-noerr.pdf



% figure(4)
% plot(time, (real(f2)), time, imag(f2), time, real(f3))
% title 'Second Harmonic';
% xlabel('delay (fs)')
% ylabel('counts')
% legend('real','imaginary','real f3')
% 
%
% figure(5)
% plot(time, (abs(dc-147)), time, abs(f1))
% title 'DC Component';
% xlabel('delay (fs)')
% ylabel('counts')
% legend('dc-147','abs(f1)')
% 
% CHANGE
% figure(6)
% Ts = time(end) - time(end-1); %when time available.
% freq = ((-length(time)/2+1):(length(time)/2))/(Ts*length(time));
% waven = freq./(c0*1e-7) + 1e7/560;
% F2 = fftshift(abs(fft(f2)))/length(time);
% expected = 1e7/266 - 2e7/560;
% plot(waven, F2,[expected expected], [0 1.1*max(F2)], '--')
% 
% % hold on
% % line()
% % hold off
% title 'Fourier transform of Second Harmonic';
% xlabel('1/\lambda (cm^{-1})')
% ylabel('counts')
% legend('fft(f2)', 'expected')

%% Fourier Transform

fig7 = figure(7);
clf(7)

L = length(time);
% NFFT=2^(nextpow2(L)+1);
NFFT=L;

% Ts = 50*2/c0; %for mms
Ts = mean(time(2:end) - time(1:end-1)); %when time available.
freq = ceil((-NFFT/2:NFFT/2-1))/(Ts*NFFT);
% depending on whether an X or Y interferometer is beign run, the sign of
% Freq must change... might be due to times coming in +/- between each case
waven = freq./(c0*1e-7) + 1e7/550;
if exist('waven', 'file')==0
    save('waven.mat','waven')
end
padded = zeros(NFFT,1);
phased = f1.*exp(-1i*angle(f1(Imin)));
padded(1:(L-Imin+1))=phased(Imin:L);
padded(1)=0.5*padded(1);
% padded1=padded;
% padded(1:Imin)=0.5*(padded(1:Imin)+conj(flipud(phased(1:Imin))));
F1 = fftshift(fft(padded, NFFT))/L;
% F11 = fftshift(fft(padded1, NFFT))/L;
% expected = 1e7/532;
% plot(waven, abs(F1),[expected expected], [0 1.1*max(abs(F1))], '--')

% Start FFT Error
paddedn = zeros(NFFT,nboot);
phasedn = fboots;
paddedn(1:(L-Imin+1),:)=phasedn(Imin:L,:);
paddedn(1,:)=0.5*paddedn(1,:);
F1n = fftshift(fft(paddedn, NFFT,1),1)/L; % remember to normalize by max(real(F1))
FError =  std(real(F1n),0,2) + 1i.*std(imag(F1n),0,2);
FErrorab = mean(std(abs(F1n),0, 2));
% End FFT error

[wmin, IFmax] = max(real(F1));
% IFmax = 141;

errbarsF = NaN(size(waven));

errbarsF(IFmax) = mean(real(FError))/max(real(F1));
errorbar(waven, (real(F1))/max(real(F1)), errbarsF, '-o','LineWidth', 1.25, 'color', [0.8500 0.3250 0.0980])
% newly added abs of the plot
hold on
errbarsF(IFmax) = mean(abs(FError))/max(abs(F1));
errorbar(waven, (abs(F1))/max(abs(F1)), errbarsF, '-o','LineWidth', 1.25, 'color', [1 1 0])

% depending on whether an X or Y interferometer is beign run, the sign of
% imag component must change... might be due to times coming in +/- between each case

errbarsF(IFmax) = mean(imag(FError))/max(real(F1));
errorbar(waven, imag(F1)/max(real(F1)), errbarsF, '-o', 'LineWidth', 1.25, 'color', [0 0.4470 0.7410])
errbarsF(IFmax) = NaN;
% plot(waven, imag(F1)/max(real(F1)), 'LineWidth', 1.25, 'color', [0 0.4470 0.7410])
hold off
ax7 = gca;
title 'Fourier transform of First Harmonic';
xlabel('1/\lambda (cm^{-1})')
ylabel('normalized counts')
% xlim([17000 20500]);
if exist('product', 'var')
    hilbilly = hilbert(product/max(product));
    hold on
    plot(waven, real(hilbilly), '--', 'LineWidth', 1.25, 'color', [100 61 97]/256)
    plot(waven, -1*imag(hilbilly), '--', 'LineWidth', 1.25, 'color', [0.0660 0.6740 0.1880])
    chH = get(gca,'Children');
%     set(gca,'Children',flipud(chH))
    hold off
    lgd7 = legend('Re', 'Im', 'Re(exp.)');
else
    lgd7 = legend('Re','Abs', 'Im');
end
lgd7.Units = 'centimeters';
lgd7.FontSize = 10;
ax7.XLimSpec = 'Tight';
ax7.LineWidth = 1.5;
ax7.TickLength = [.02 .02];
% ax7.XLim = waverange;
ax7.XLim = [17000 20500];
fig7pos = get(fig7, 'Position');
set(fig7, 'PaperUnits','centimeters')
set(fig7, 'PaperSize', [18, 18/fig7pos(3)*fig7pos(4)])

dim = [.18 .55 .3 .3];
str = ['S/N = ' num2str(max(real(F1))/mean([mean(real(FError)) mean(imag(FError))]))];
annotation('textbox',dim,'String',str,'FitBoxToText','on', 'Units', 'normalized');

print -dpdf -r300 -bestfit FPGA-F1FT.pdf

% This added for F2, 20191021
fig8 = figure(8);

padded2 = zeros(NFFT,1);
phased2 = f2.*exp(-1i*angle(f2(Imin)));
padded2(1:(L-Imin+1))=phased2(Imin:L);
% padded2(1:Imin)=0.5*(padded2(1:Imin)+conj(flipud(phased2(1:Imin))));
% padded3 = conj(flipud(phased2(1:Imin)));
F2 = fftshift(fft(padded2, NFFT))/L;
% F11 = fftshift(fft(padded1, NFFT))/L;
expected = 2e7/526-1e7/550;
% plot(waven, abs(F1),[expected expected], [0 1.1*max(abs(F1))], '--')

plot(waven, abs(F2)/max(abs(F2)),'LineWidth', 1.25, 'color', [0.8500 0.3250 0.0980])
hold on
% plot(waven, imag(F2)/max(abs(F2)), 'LineWidth', 1.25, 'color', [0 0.4470 0.7410])
plot([expected expected], [0 1.1], '--', 'LineWidth', 1.25, 'color', [0 0.4470 0.7410])
hold off
ax8 = gca;
title 'Fourier transform of Second Harmonic';
xlabel('1/\lambda (cm^{-1})')
ylabel('normalized counts')

lgd8 = legend('Abs(F2)', '526-2F');
lgd8.Units = 'centimeters';
lgd8.FontSize = 10;

ax8.XLimSpec = 'Tight';
ax8.LineWidth = 1.5;
ax8.TickLength = [.02 .02];
% ax8.XLim = waverange;
fig8pos = get(fig8, 'Position');
set(fig8, 'PaperUnits','centimeters')
set(fig8, 'PaperSize', [12, 12/fig7pos(3)*fig7pos(4)])

print -dpdf -r300 -bestfit FPGA-F2FT.pdf

% Angle error away from Imin
Errorforangle = Errorsmoothed(:,4); % make sure to plot the right hand side and check that it is well behaved.
Errorforangle(Imin-16:Imin+16) = [];
angerrbad = mean(Errorforangle);

% writemat = [time.' real(f1tiemo) Error(:,1) imag(f1tiemo) Error(:,2) abs(f1tiemo) Error(:,3) angle(f1tiemo) Error(:,4)];
writemat = [time.' real(f1) imag(f1) abs(f1) angle(f1)];
writematn = [time.' real(f1normphased) imag(f1normphased) abs(f1normphased) angle(f1normphased)];
writematw = [waven.' real(F1/max(abs(F1))) imag(F1/max(abs(F1))) abs(F1/max(abs(F1))) angle(F1/max(abs(F1)))];
writematw2 = [waven.' real(F2/max(abs(F2))) imag(F2/max(abs(F2))) abs(F2/max(abs(F2))) angle(F2/max(abs(F2)))];
dlmwrite([pwd outpath outprefix '-f1t.txt'],writemat, 'newline', 'pc', 'precision', 10);
dlmwrite([pwd outpath outprefix '-f1tnormphased.txt'],writematn, 'newline', 'pc', 'precision', 10);
dlmwrite([pwd outpath outprefix '-F1w.txt'],writematw,'newline', 'pc', 'precision', 10);
% dlmwrite([pwd outpath outprefix '-F2w.txt'],writematw2,'newline', 'pc', 'precision', 10);
dlmwrite([pwd outpath outprefix '-err.txt'],Error,'newline', 'pc', 'precision', 10);


%% time correction
figcorr = figure(9);
options = fitoptions('gauss2');
options.lower = [0,-Inf,10,-Inf,-50,100];
options.upper = [Inf,Inf,Inf,Inf,50,Inf];
fitobject = fit(time.', abs(f1), 'gauss2', options);
plot(fitobject, time, abs(f1))
timecorr = fitobject.b1;
display(time(Imin))

%% amp correction
figampcorr = figure(10);
options2 = fitoptions('gauss1');
fitobject2 = fit(time.', real(f1), 'gauss1', options2);
plot(fitobject2, time, real(f1))