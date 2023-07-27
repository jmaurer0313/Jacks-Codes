function [] = txt2frg(filePath,outputName)
%TXT2FRG Summary of this function goes here
% Load the txt file with path filePath and save it as the .frg file for binner
% === file structure of .txt ===
% width(temporal domain) of FROG trace in pixel
% height(spectral domain) of FROG trace in pixel
% temporal_calibration (fs/px)
% spectral_calibration (nm/px)
% central_wavelength (nm)
% FROG_TRACE_RAW_DATA

if nargin==1
    [~, outputName, ~] = fileparts(filePath);
    outputName = convertStringsToChars(outputName);
end
fileID = fopen(filePath,'r');
traceData = fscanf(fileID,'%f');
numDelayPoints = traceData(1);
numSpectralPoints = traceData(2);
delayStep = traceData(3);
specResolution = traceData(4);
centerWavelength = traceData(5);
trace = traceData(6:end);
trace = reshape(trace,numSpectralPoints,numDelayPoints); 

% now cut the trace down to be just the relevant region of the FROG
% spectrum - to make the matrix symmetric for binner 
% lam = (-numSpectralPoints/2:numSpectralPoints/2-1) * specResolution + centerWavelength;
% trace = trace(49:348,:);

header = [numDelayPoints;numSpectralPoints;delayStep;specResolution;centerWavelength];
save([outputName '.frg'],'header','-ASCII')
save([outputName '.frg'],'trace','-APPEND','-ASCII')

disp(['saving the frog trace as ' [outputName '.frg']])