% avg2DFPGA_RT

function [RPmatOut, NRPmatOut, tbx, tby, window2Dfinal,numPhots] = avg2DFPGA_RT(plotOpt, homeFolder,retiming,rephase)
startFolder= pwd;
cd(homeFolder); 
filesList= dir('*2DFPGA'); 
numFiles=numel(filesList); 
curPath=pwd;

Tdes=0.6; 

for s=1:numFiles
    
    
    fileFolder= [curPath filesep() filesList(s).name]; 
    
    [XmatTemp, YmatTemp, tb1, tb2] = linear2DFPGA_RT(fileFolder);
    
    [~, ~, Xshifts, Yshifts, Xangles, Yangles, ~, ~, ~, ~, NIFFT,Xzero] = interpLinear(plotOpt,XmatTemp, YmatTemp, tb1, tb2, Tdes);

    [RPmat, NRPmat, tb1, tb2,~,~,numPhots] = ind2DFPGA_RT(fileFolder);
    
    [ NRPtimed , RPtimed , tnX2 , tnY2 , window2Dfinal ] = interp2D_Phased( plotOpt, NRPmat, RPmat, tb1, tb2, NIFFT , Xshifts, Yshifts, Xzero, retiming, rephase);
  
    if s==1
            RPmatOut= RPtimed;
            NRPmatOut= NRPtimed;
    else
        RPmatOut=RPmatOut+RPtimed;
        NRPmatOut=NRPmatOut+NRPtimed; 
    end 
     

     
end 
 RPmatOut= RPmatOut./numFiles;
 NRPmatOut= NRPmatOut./numFiles;

tbx=tnX2;
tby=tnY2;
 


cd(startFolder); 
end 