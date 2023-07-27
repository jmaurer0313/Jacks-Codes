% divde up an interval with multplie decades in time such that there is an
% equal number of bins in each decade, and then fill those bins with
% occurneces according to a random array that is indexed mapped onto the
% properly spaced-evenly binned array.

function [valuesArray,finalLogArray] = logRandUniform(lowerBound, upperBound,pointsPerDecade,totalRands)
% first construct an array containing an equal number of bins for each decade in time, 
% based on on end point values passed in 
% lowerBound=6e-6;
% upperBound=3e-1;
% pointsPerDecade=100;
% totalRands=100000; 
logArray=[];
lowExp=floor(log10(lowerBound));
highExp=ceil(log10(upperBound)); 
lowLim= 1*10^(lowExp); 
highLim= 1*10^(highExp); 
numDecades=log10(highLim/lowLim); 

if numDecades==1
  pointsPerDecade=pointsPerDecade*10;  
end

for i=1:numDecades
   logArray=[logArray linspace(lowLim*(10^(i-1)),0.99*lowLim*(10^i),pointsPerDecade)];    
end

% now with the logArray in hand - which is an equal number of bins between
% each decade in time - generate random numbers on the order of 0->length(logArray) and 
% then map the random numbers onto the index of logArray and select numbers
% in that way (after cutting the edges of logArray to match the bounds
% passed in)
finalLogArray=logArray(logArray>=lowerBound);
finalLogArray=finalLogArray(finalLogArray<=upperBound); 
randIndexArray=round(rand(1,totalRands)*(length(finalLogArray)-1))+1; 

valuesArray=finalLogArray(randIndexArray);

% LA = log(A); LB = log(B);
% myArray= exp(LA + (LB-LA) * rand(1,n));

end