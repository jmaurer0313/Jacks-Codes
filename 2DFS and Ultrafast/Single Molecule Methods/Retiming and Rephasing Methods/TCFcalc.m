function [ avgTCF, numPairs ] = TCFcalc( F, tauArray )

avg= nanmean(F);
TCFarr= zeros(1,numel(tauArray));
numPairs= zeros(1,numel(tauArray));
     for j = 1:numel(tauArray)
         curTau = tauArray(j);
         for t=1:numel(F)-curTau
           if isnan(F(t))||isnan(F(t+curTau))
           continue
           else
             TCFtemp = (conj(F(t))-conj(avg))*(F(t+curTau)-avg);
             TCFarr(j) = TCFarr(j) + TCFtemp; 
             numPairs(j) = numPairs(j)+1; 
           end
         end
     end
avgTCF = TCFarr./numPairs;

end

