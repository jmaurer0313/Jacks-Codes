function w = delayedGaussian(x,c,s)
%DELAYEDGAUSSIAN Delayed gaussian window.
%   DELAYEDGAUSSIAN(x, c, s) returns an gaussian with transition width s=
%   3*sigma starting with cutoff c. For x<=c, w is one.
% 
%   Author(s): A. Tamimi

w = ones(size(x));
shifted = x-c;
index = shifted > 0;
w(index) = exp(-4.5*(shifted(index).^2)/(s^2));

end