%AUTHOR: Jack Maurer 
% bulk read of the linear data with plot

base_folder='C:\Users\PMECS\Documents\AutoCorrData\Cy3DNA\20201020-161919';
time = dlmread(fullfile(base_folder, 'timeDomain'));
xQuad= dlmread(fullfile(base_folder, 'XQuad'));
yQuad= dlmread(fullfile(base_folder, 'YQuad'));

linSig= xQuad + 1i*yQuad; 

figure(1)
hold on
plot(time, abs(linSig)); 
plot(time, real(linSig)); 
plot(time, imag(linSig)); 
title('Linear Scan - time domain'); 
xlabel('Time (fs)'); 
ylabel('Magnittude (arb)'); 
legend('Abs','Real','Imag'); 
hold off

% small block to correct timing (in fs)

% deltaT is the distance (in time) off from zero as determined from the
% gauss curve fit to the linear scan/autocorr data.

deltaT= 10; 
(deltaT*1e-15)*(2.98e8)*(1e3)/2



