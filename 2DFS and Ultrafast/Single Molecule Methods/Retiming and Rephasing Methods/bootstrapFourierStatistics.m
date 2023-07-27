%% bootstrapFourierStatistics (bsFS) bootstraps the input data, and 
% calculates the positive Fourier Transform of the data. This assumes input
% data is real. If it is complex, then negative frequencies must also be
% calculated, and this will take some reworking.
%
% Inputs:  
%    data is a vector of values
%    nboot is the number of bootstrap samples to run. 1000 is a good default.
%    nfreq is the number of frequencies to calculate.
% Outputs:
%    xerr, yerr, abserr, phierr: these each return vectors of length 
%    nfreq+1. The nth entry corresponds to the error of the n-1st fourier
%    component. The 1st entry should always be 0, since the size of your
%    data doesn't change.
%    rawbootstrap: This is the raw output from bootstrapping. It is a 
%    nboot*length(data) size matrix.
%    fftbootstrap: This is a nboot*(nfreq+1) sized matric containg the 
%    positive fourier transform of each bootstrapping iteration
%
%function [xerr, yerr, abserr, phierr, rawbootstrap, fftbootstrap]
%   = bootstrapFourierStatistics(nboot, nfreq, data)

%%%
%Last edited March 5th 2019  Tiemo Landes
%%%

%% In place Fourier Transform Method. For a discrete run of data summing
%% exp(i*data) is the same as taking the Fourier Transform of the
%% histogrammed data modulo binning.

function [xerr, yerr, abserr, phierr, rawbootstrap, fftbootstrap] ...
    = bootstrapFourierStatistics(nboot, nfreq, data)

    %rawbootstrap is a nboot*length(data) matrix
    N = length(data);
    rawbootstrap = bootstrp(nboot, @(x) x, data);
    fftbootstrap = zeros(nboot, nfreq+1);
    for k = 0:(nfreq)
        %Calculate the kth fourier component of the data
        fftbootstrap(:,k+1) = 1/N * sum(exp(-1i*k*rawbootstrap),2);
    end
    xerr = std(real(fftbootstrap),0,1);
    yerr = std(imag(fftbootstrap),0,1);
    abserr = std(abs(fftbootstrap),0,1);
    phierr = std(angle(fftbootstrap),0,1);

end
% Histogram Method
%
% function [xerr, yerr, abserr, phierr, rawbootstrap, histbootstrap, ...
%     ffthist] = bootstrapFourierStatistics(nboot, nhist, data)
% 
%     %rawbootstrap is a nboot*length(data) matrix
%     edges = linspace(0, 2*pi, nhist+1);
%     rawbootstrap = bootstrp(nboot, @(x) x, data);
%     histbootstrap = histc(rawbootstrap, edges, 1);
%     histbootstrap = histbootstrap(1:(end-1),:);
%     ffthist = fft(histbootstrap,nhist,1);
%     xerr = std(real(ffthist),0,2);
%     yerr = std(imag(ffthist),0,2);
%     abserr = std(abs(ffthist),0,2);
%     phierr = std(angle(ffthist),0,2);
% 
% end