%% Constants
c0 = 299.792458; % speed of light in nm/fs
mono = 560; % monochromater wavelength in nm
f0 = c0/532-c0/mono; % frequency in PHz, green.
% Ts = 1/(2.5*f0); % experimental sampling increment in fs
Ts = 7; % experimental sampling (step size) increment in fs
t = (0:15)*Ts; % time axis
sigma = 50; % target is 30 fs decay
a = 1/(2*(sigma.^2)); % input width to gaussian function
deltat= 2; % mistiming in fs, this will be the parameter of 
% interest for resetting T0 for all underlying scans
NFFT = 2^(nextpow2(length(t))+2);

%% Reference function
gauss = @(x, a) sqrt(a/pi)*exp(-a*x.^2); % general gaussian of std = 1/sqrt(2a);
window = delayedGaussian(t, 0, 150);

S = (gauss(t, a).*exp(2i*pi*f0*t)).*window;% signal in time, a gaussian with an oscillating complex exponential
Sp = (gauss(t-65, a).*exp(2i*pi*f0*t)).*window; % Optional phase
% shift applied

% S= XmatTemp(1,:); 
% S = ifftshift(S);
figure(1)
plot(t, real(S), t, imag(S),t,abs(S))
title('Signal with no timing error - Real/Imag');
xlabel('Time (fs)');
ylabel('Magnitude (arb)');

figure(6)
plot(t, real(S), t, real(Sp))
hold on
plot(t, abs(S),'--', t, abs(Sp),'--')
title('Real parts of Timed and Mistimed signals'); 
xlabel('Time (fs)');
ylabel('Magnitude (arb)');
% plot(t, real(S), t, imag(S),t,abs(S), t, abs(S.*window))
S(1) = 0.5*S(1); % weird but necessary, some NMR reference somewhere.
Sp(1)=0.5*Sp(1); 



figure(2)
f = ceil((-NFFT/2:NFFT/2-1))/(Ts*NFFT); % freq axis in PHz
f = f/c0*1e7; % optional conversion into wavenumbers.
SFuns = fft(S, NFFT);
SF = fftshift(fft(S, NFFT)); % Signal in the frequency domain, FFT shift applied 
SpF = fftshift(fft(Sp, NFFT)); % Signal in the frequency domain, FFT shift applied 
plot(f, real(SF), f, imag(SF), f, abs(SF), f, abs(SFuns))

figure(3)
plot(f+mono, real(SF), f+mono, real(SpF))
title('Real part of Fourier Transform');
xlabel('Frequency (1/cm)');
ylabel('Magnitude (arb)'); 
xlim([-500,2500]);
% hold on
% plot(f, imag(hilbert(-real(fft(S, NFFT)))))
% hold off

% The circshift bit (where to drop the zeros)
forshift = real(ifftshift(SF));
% try dropping the zeros at the edge of the freq domain, OR the wrap around
% point where freq goes + to -
Npad = 1024;
% zeroes must be 
SFzeros=[real(SF) zeros(1,Npad)];
SpFzeros=[real(SpF) zeros(1,Npad)];

unshiftedNew= circshift(SFzeros, -floor(NFFT/2));
unshiftedNewP= circshift(SpFzeros, -floor(NFFT/2));


Runsh = unshiftedNew;
RunshP = unshiftedNewP;
NIFFT = length(Runsh);

Sn=ifft(Runsh, NIFFT); 
SnP=ifft(RunshP, NIFFT); 

Sn=Sn*2;
SnP=SnP*2;

Sn((ceil(NIFFT/2)+1):end)=0; 
SnP((ceil(NIFFT/2)+1):end)=0; 
% Cunsh = complex(Runsh, imag(hilbert(-Runsh)));

Tsn = (Ts*NFFT)/NIFFT;
tn = ((0:NIFFT-1))*Tsn;


% Currently the insertion point for zeroes is set by an inspection of the
% maximum distance from the peak of SF. This can be done in a equivalent
% way using the minimum of SF rather than needing to choose a numerical
% constant for each scan.

% [M, Ipad] = min(abs(ifftshift(f)-insertpt));
% [M, Ipad] = min(abs(ifftshift(SF)));
% shifted = circshift(forshift, length(f)-Ipad);
% padded = [shifted zeros(1,Npad)];
% unshifted = circshift(padded, -length(f)+Ipad);

% Runsh = unshifted;
% Cunsh = complex(Runsh, imag(hilbert(-Runsh)));
% NIFFT = length(Cunsh);
% Tsn = (Ts*NFFT)/NIFFT;
% tn = ((0:NIFFT-1))*Tsn;
% tn = ceil((-NIFFT/2:NIFFT/2-1))*Tsn;

% Sn = (ifft(Cunsh, NIFFT));
% Sn(1) = 2*Sn(1);

S(1) = 2*S(1); % undo.
S = S*NFFT/NIFFT;

Sp(1) = 2*Sp(1); % undo.
Sp = Sp*NFFT/NIFFT;


figure(5)
% plot(1:NIFFT, real(Cunsh), 1:NIFFT, imag(Cunsh), 1:NIFFT, abs(Cunsh))
plot(tn, real(Sn), tn, imag(Sn), tn, abs(Sn))
hold on
ax = gca;
ax.ColorOrderIndex = 1;
plot(t, real(S),'o', t, imag(S),'o');
title('Interpolated Retimed Data'); 
xlabel('Time (fs)');
ylabel('Magnitude (arb)');
xlim([0,250]); 
hold off

figure(4)
plot(tn, real(SnP), tn, imag(SnP), tn, abs(SnP))
hold on
ax = gca;
ax.ColorOrderIndex = 1;
plot(t, real(Sp),'o', t, imag(Sp),'o');
title('Interpolated Mistimed Data'); 
xlabel('Time (fs)');
ylabel('Magnitude (arb)');
xlim([0,250]); 

hold off



%% Function which needs to be retimed

% Sdelta = gauss(t-deltat, a).*exp(2i*pi*f0*(t-deltat)); % signal to retime
% Sdelta(1) = 0.5*Sdelta(1); % weird but necessary, some NMR reference somewhere.
% % Sdelta = ifftshift(Sdelta); % in case of symmetric
% 
% figure(3)
% plot(t, real(Sdelta), t, imag(Sdelta), t, abs(Sdelta))
% SdF = fftshift(fft(Sdelta, NFFT));
% SdF = SdF.*exp(2i*pi*deltat*(f0));
% figure(4)
% % hold on
% plot(f, real(SdF), f, imag(SdF), f, abs(SdF), f, real(SF))
% % hold off
% figure(5)
% Sr = (ifft(ifftshift(real(SdF))));
% % Sr= Sr(1:length(t));
% plot(abs(fftshift(Sr)))