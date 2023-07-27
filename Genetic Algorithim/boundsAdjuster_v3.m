function [boundsArrayAdj,atEdges]= boundsAdjuster_v3(boundsArray,windowMag,paramAdjs,absLow,absHigh,Nparam,Nstates,fixedDecades)

atEdges=zeros(1,(Nparam-2*Nstates));
% UPDATE: 03-16-2022 - adding sigma of the gaussians to the fitting
% routine, will require a rework of Nparam vs Nstates and so on. mainly
% that the Ntimes=Nparam - 2*Nstates, since there is now an average and a
% sigma for every state
if sum(paramAdjs)~=0
for m=1:(Nparam-2*Nstates)
    
    if ~fixedDecades
    intInterval=(boundsArray(m,2)-boundsArray(m,1));
    
    boundsArray(m,1) = boundsArray(m,1)+(windowMag*(intInterval)*paramAdjs(m)); 
    boundsArray(m,2) = boundsArray(m,2)+(windowMag*(intInterval)*paramAdjs(m)); 
    
    if boundsArray(m,1)<absLow
        boundsArray(m,1)=absLow;
        atEdges(m)=1; 
    end
        if abs(boundsArray(m,2)-boundsArray(m,1))<intInterval || boundsArray(m,2)<boundsArray(m,1)
          boundsArray(m,2)=boundsArray(m,1)+intInterval;  
         
        end
    
    
    if boundsArray(m,2)>absHigh
       boundsArray(m,2)=absHigh;
       atEdges(m)=1; 
    end
       if abs(boundsArray(m,2)-boundsArray(m,1))<intInterval || boundsArray(m,2)<boundsArray(m,1)
          boundsArray(m,1)=boundsArray(m,2)-intInterval;   
          
        end
    
    elseif fixedDecades
        
        intInterval=(boundsArray(m,2)-boundsArray(m,1));
        lowExp=floor(log10(boundsArray(m,1)));
        highExp=ceil(log10(boundsArray(m,2))); 
        lowLim= 1*10^(lowExp); 
        highLim= 1*10^(highExp); 
        numDecades=log10(highLim/lowLim);    
%         lowLimMult=boundsArray(m,1)/lowLim;
%         highLimMult=boundsArray(m,2)/highLim; 
        
        boundsArray(m,1) = boundsArray(m,1)+(windowMag*(intInterval)*paramAdjs(m)); 
        

        if boundsArray(m,1)<absLow
            boundsArray(m,1)=absLow;
            atEdges(m)=1; 
        end
        
        boundsArray(m,2) = boundsArray(m,1)*10^(numDecades); 
        
        if boundsArray(m,2)>absHigh
            boundsArray(m,2)=absHigh;
            boundsArray(m,1)=absHigh*10^(-numDecades);
            atEdges(m)=1; 
        end
         
%         if boundsArray(m,2)<boundsArray(m,1)
%           boundsArray(m,1)=absHigh*10^(-numDecades);          
%         end    
      
        
    end
end
end
boundsArrayAdj=boundsArray; 

end 