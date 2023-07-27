%set up edges, "centers" is basically the phase axis of 2n+1 odd numbers.
% close all
% JCP size: one column = 8.5cm, two column = 17cm width
fontname = 'Times New Roman'; 
set(groot,'DefaultAxesFontName',fontname,'DefaultTextFontName',fontname);
set(groot,'defaultAxesFontSize',12)
nboot = 1000;
nbootsc = 1000;
timerange = [-100 100];
% sntimerange = [-55,55];
% snindex = time > sntimerange(1) & time < sntimerange(2);
waverange = [1.5 2.18]*1e4;

% if exist('fitobject', 'var')
%     normfact = fitobject.a1;
% else
%     normfact = 1;
% end

normfact = 1;

n = 14; % 2n+1 total histogram bins, produces one bin centered at zero, necessary.
N = 2*n +1;
edges = linspace(-1, 1, 2*n +2);
edgesl = edges(1:end-1);
edgesr = edges(2:end);
centers = 0.5*(edgesl + edgesr);
% load actual histogram data
load('D:\Users\amroa\Dropbox\FPGA SNR analysis\Data\Auto\20190211-205227-auto\auto.mat');
outpath = '\ProcessedData\20190207-FPGA-20200306\17k-16.5\';
outprefix = '205227-17k-16.5';
% figout = 'D:\Users\amroa\Dropbox\FPGA SNR analysis\Matlab\1DFPGAscan\ProcessedData\20190207-FPGA\17k-16.5';

time = time*1000; %ps to fs
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
stdX = std(real(exp(-1i*pi*(photons{Imin}-180 + angle(f10)*180/pi)/180)));
stdY = std(imag(exp(-1i*pi*(photons{Imin}-180 + angle(f10)*180/pi)/180)));
stdout = [stdX, stdY];
meanX = mean(real(exp(-1i*pi*(photons{Imin}-180 + angle(f10)*180/pi)/180)));
meanY = mean(imag(exp(-1i*pi*(photons{Imin}-180 + angle(f10)*180/pi)/180)));
meanout = [meanX, meanY];


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

for i = 1:pL
    [xerr, yerr, abserr, phierr, rawbootstrap, fftbootstrap] = bootstrapFourierStatistics(nboot, 1, pi*(photons{i}-180 + angle(f10)*180/pi)/180);
    f1tiemo(i) = sum(1/length(photons{i})*exp(-1i*pi*(photons{i}-180 + angle(f10)*180/pi)/180));
    Error(i,:) = [xerr(2) yerr(2) abserr(2) phierr(2)];
end


%scaling plot from long run
nscalingpts = 100;
minphot = 10;
maxphot = length(photons{Imin});
sqmaxphot = floor(sqrt(maxphot));
sqminphot = floor(sqrt(minphot));
scaling = @(x) x.^2;
linear = linspace(sqminphot, sqmaxphot, nscalingpts);
scalingpts = floor(scaling(linear));
extras = ceil([1943,804,315,133,36.8888888900000,17490]);
scalingpts = [scalingpts extras];
scalingpts = sort(scalingpts);
Nsc = length(scalingpts);

Errorsc = zeros(Nsc,4);
f1sc = zeros(Nsc,1);
f1bt = zeros(nbootsc,1);
photonsmin = photons{Imin};

for j = 1:Nsc
    scalingpt = scalingpts(j);    
    for i=1:nbootsc
        currentphot = randsample(photonsmin, scalingpt, true);       
        f1bt(i) = sum(1/length(currentphot)*exp(-1i*pi*(currentphot-180 + angle(f10)*180/pi)/180));
    end
    Errorsc(j,:) = [std(real(f1bt)) std(imag(f1bt)) std(abs(f1bt)) std(angle(f1bt))];
end


transformedz = fft(repmat(ifftshift(yz, 2), 1, padding),[],2)/M;
% step axis
steps = (0:499)*50;
c0 = 299.792458; %nm/fs
% time = steps*2/c0;
steps = steps/1e6;
dc = transformedz(:, 0*padding +1);
f1 = transformedz(:, 1*padding +1);
rate = round(length(photons{Imin})/540);
vis = 2*abs(f1(Imin))/dc(Imin);
f1 = f1tiemo;
f2 = transformedz(:, 2*padding +1);
f3 = transformedz(:, 3*padding +1);

f1normphased = f1/normfact;

ErrorMeans = NaN(2,4);
ErrorMeans(1,:) = Error(Imin, :);
ErrorMeans(2,:) = ErrorMeans(1,:)./normfact;
errbars = NaN(size(f1));

% time = (steps)*2000000/c0-177/2;
% time = linspace(-1,1,199)*99;

fig3 = figure(3);
fig3pos = get(fig3, 'Position');
set(fig3, 'PaperUnits','centimeters')
set(fig3, 'PaperSize', [12, 12/fig3pos(3)*fig3pos(4)])
% plot(time, (real(f1)), time, imag(f1), time, abs(f1), '--')
% errorbar(time, real(f1), Error(:,1), 'LineWidth', 1.25, 'color', [0.8500 0.3250 0.0980])
% hold on
% errorbar(time, imag(f1), Error(:,2), 'LineWidth', 1.25, 'color', [0 0.4470 0.7410])
% errorbar(time, abs(f1), Error(:,3), '--', 'LineWidth', 1.25, 'color', [0.4660 0.6740 0.1880])
% hold off
errbars(Imin) = ErrorMeans(1,1);
errorbar(time, (real(f1)), errbars, 'LineWidth', 1.25, 'color', [0.8500 0.3250 0.0980])
hold on
errbars(Imin) = ErrorMeans(1,2);
errorbar(time, imag(f1), errbars, 'LineWidth', 1.25, 'color', [0 0.4470 0.7410])
errbars(Imin) = NaN;
errbars(Imin) = ErrorMeans(1,3);
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
% ax3.XLim = timerange;
ax3.LineWidth = 1.5;
ax3.TickLength = [.02 .02];

print -dpdf -r300 -bestfit FPGA-f1-err.pdf





fig4 = figure(4);
fig4pos = get(fig4, 'Position');
set(fig4, 'PaperUnits','centimeters')
set(fig4, 'PaperSize', [12, 12/fig4pos(3)*fig4pos(4)])
plot(time, (real(f1)), time, imag(f1), time, abs(f1), '--')
errorbar(time, real(f1), Error(:,1), 'LineWidth', 1.25, 'color', [0.8500 0.3250 0.0980])
hold on
errorbar(time, imag(f1), Error(:,2), 'LineWidth', 1.25, 'color', [0 0.4470 0.7410])
errorbar(time, abs(f1), Error(:,3), '--', 'LineWidth', 1.25, 'color', [0.4660 0.6740 0.1880])
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
% ax4.XLim = timerange;
ax4.LineWidth = 1.5;
ax4.TickLength = [.02 .02];

print -dpdf -r300 -bestfit FPGA-f1-fullerr.pdf



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
% % CHANGE
% figure(6)
% Ts = 50*2/c0;
% freq = ((-length(time)/2):(length(time)/2)-1)/(Ts*length(time));
% waven = freq./(c0*1e-7);
% F2 = fftshift(abs(fft(f2)))/length(time);
% expected = -1e7/266 + 2e7/560;
% plot(waven, F2,[expected expected], [0 1.1*max(F2)], '--')
% 
% % hold on
% % line()
% % hold off
% title 'Fourier transform of Second Harmonic';
% xlabel('1/\lambda (cm^{-1})')
% ylabel('counts')
% legend('fft(f2)', 'expected')

fig5 = figure(5);
fig5pos = get(fig5, 'Position');
set(fig5, 'PaperUnits','centimeters')
set(fig5, 'PaperSize', [18, 18/fig5pos(3)*fig5pos(4)])
loglog(scalingpts, Errorsc(:, 1), scalingpts, Errorsc(:,2), scalingpts, Errorsc(:,3), scalingpts, Errorsc(:,4))
ax5 = gca;
title 'Bootstrapping error as a function of sample size';
xlabel('sample size (photons)')
ylabel('error (au)')
lgd5 = legend('Re','Im','Abs', 'Angle');
lgd5.Units = 'centimeters';
lgd5.FontSize = 10;
ax5.XLimSpec = 'Tight';
% ax5.XLim = timerange;
ax5.LineWidth = 1.5;
ax5.TickLength = [.02 .02];

print -dpdf -r300 -bestfit FPGA-f1-errsc.pdf

% fig7 = figure(7);
% L = length(time);
% NFFT=2^(nextpow2(L)+1);
% % Ts = 50*2/c0; %for mms
% Ts = time(end) - time(end-1); %when time available.
% freq = ((-NFFT/2+1):(NFFT/2))/(Ts*NFFT);
% waven = freq./(c0*1e-7) + 1e7/560;
% if exist('waven', 'file')==0
%     save('waven.mat','waven')
% end
% padded = zeros(NFFT,1);
% phased = f1.*exp(-1i*angle(f1(Imin)));
% padded(1:(L-Imin+1))=phased(Imin:L);
% % padded1=padded;
% % padded1(1:Imin)=0.5*(padded1(1:Imin)+conj(flipud(phased(1:Imin))));
% F1 = fftshift(fft(padded, NFFT))/L;
% % F11 = fftshift(fft(padded1, NFFT))/L;
% expected = 1e7/532;
% % plot(waven, abs(F1),[expected expected], [0 1.1*max(abs(F1))], '--')
% 
% plot(waven, real(F1)/max(real(F1)),'LineWidth', 1.25, 'color', [0.8500 0.3250 0.0980])
% hold on
% plot(waven, imag(F1)/max(real(F1)), 'LineWidth', 1.25, 'color', [0 0.4470 0.7410])
% hold off
% ax7 = gca;
% title 'Fourier transform of First Harmonic';
% xlabel('1/\lambda (cm^{-1})')
% ylabel('normalized counts')
% if exist('product', 'var')
%     hilbilly = hilbert(product/max(product));
%     hold on
%     plot(waven, real(hilbilly), '--', 'LineWidth', 1.25, 'color', [100 61 97]/256)
% %     plot(waven, -1*imag(hilbilly), '--', 'LineWidth', 1.25, 'color', [0.0660 0.6740 0.1880])
%     chH = get(gca,'Children');
% %     set(gca,'Children',flipud(chH))
%     hold off
%     lgd7 = legend('Re', 'Im', 'Re(exp.)');
% else
%     lgd7 = legend('Re', 'Im');
% end
% lgd7.Units = 'centimeters';
% lgd7.FontSize = 10;
% ax7.XLimSpec = 'Tight';
% ax7.LineWidth = 1.5;
% ax7.TickLength = [.02 .02];
% ax7.XLim = waverange;
% fig1pos = get(fig1, 'Position');
% set(fig7, 'PaperUnits','centimeters')
% set(fig7, 'PaperSize', [12, 12/fig1pos(3)*fig1pos(4)])
% 
% print -dpdf -r300 -bestfit FPGA-F1FT.pdf


writemat = [scalingpts.' Errorsc];
dlmwrite([pwd outpath outprefix '-errsc.txt'],writemat,'newline', 'pc', 'precision', 10);
save([pwd outpath outprefix '-errsc.mat'], 'scalingpts','Errorsc');


% % time correction
% figcorr = figure(8);
% options = fitoptions('gauss2');
% fitobject = fit(time.', abs(f1), 'gauss2', fitoptions);
% plot(fitobject, time, abs(f1))
% timecorr = fitobject.b1;
% display(time(Imin))