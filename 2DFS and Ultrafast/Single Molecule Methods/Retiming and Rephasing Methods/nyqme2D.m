function [nyquistnm, nyquistfs,steps] = nyqme2D(maxwave,maxtime, mononm)
%NYQME [nyquistnm,steps] = nyqme(maxwave(cm^-1),maxtime(fs))
%   Returns [nyquistnm, steps] where the first is the Nyquist limiting
%   step size in nm and the latter is the total number of such steps necessary
%   to sample maxwave, the frequency in wavenumbers, and maxtime, the total
%   length of the scan in time. 

c0 = 299.792458; %nm/fs
stepnm = 10; % minimum step size in nm
% maxwave = 2000; % 1/lambda in wavenumbers
% maxtime = 500; % in fs
lam = 1e7/abs(maxwave-1e7/mononm);

nyquistnm = 0.25*lam; % minimum step: 1/2 for Nyq, 1/2 for double pass, nm.
nyquistnm = stepnm*floor(nyquistnm./stepnm);
% nyquistmm = nyquistnm*1e-6;

nyquistfs = 2*nyquistnm/c0;
steps = ceil(maxtime/nyquistfs)+1;

end

