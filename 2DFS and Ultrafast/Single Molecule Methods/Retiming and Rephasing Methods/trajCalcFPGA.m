% GOAL: function to take in a photon stream with a time and phase list from
% the FPGA. Calculate a phase factor trajectory (averaged and unaveraged)
% at the assigned resolution (in usec). 

function [trajAvg, trajRaw, trajLength ] = trajCalcFPGA(phases, timesUsec, res)

% it will be assumed that the phase list coming in has been scaled to a
% range of 0 to 2Pi
  XquadP1 = mean(cos(phases));
  YquadP1 = mean(sin(phases));
   
%    could be something about going from 0->2pi to -angle -> angle, if
%    issues arise could be from this subtraction 
   shiftAngle1=atan2(YquadP1,XquadP1);
   
   phases = (phases)-(shiftAngle1);
   
   curTrajP1= [timesUsec, phases];
   binNum1=1; 
%    binSize is assigned in usec
    binSize=res;
   trajLength=(max(timesUsec)-min(timesUsec));
   Bins=floor(trajLength/binSize); 
%     trajAvg= zeros(1,Bins);
    trajRaw= zeros(1,Bins);
    counts=zeros(1,Bins);
  
% set the temp holders to zero for the start of a new file/trajectory
    tempCounts=0;
    temptrajRaw=0;
%     temptrajAvg=0;

    for j=1:(numel(timesUsec))
%         disp(['Current Ch1 iteration is ' int2str(j)]);
        if(curTrajP1(j,1)<(binNum1*binSize))
        tempCounts=tempCounts+1;
        temptrajRaw=temptrajRaw + exp(1i*(curTrajP1(j,2)));
%          temptrajAvg=temptrajAvg + exp(1i*(curTrajP1Map(j,2)));
           
        else          
                      
           counts(binNum1) = tempCounts;
           trajRaw(binNum1)= temptrajRaw;        

%            take the bin num to the proper bin which the current time value falls 
%            in rather than update the counter by 1
            binNum1= floor(curTrajP1(j,1)/binSize)+1;
        
        tempCounts=1;
        temptrajRaw= exp(1i*(curTrajP1(j,2)));        
      
                                
        end
        
    end
    
    trajAvg=trajRaw./counts; 



end 