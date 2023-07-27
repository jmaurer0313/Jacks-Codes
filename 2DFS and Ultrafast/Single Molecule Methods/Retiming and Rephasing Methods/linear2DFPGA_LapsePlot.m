% GOAL: create a variable int time linear plotting code which doesn't both
% to calculate the nonlinear contributions

fileFolder = 'C:\Users\Jack\Documents\Marcus Lab\2D_SingleMolecule\Sm_avgFolders\TimeLapses\20191217_mol1_1hr\20191217-164312-2DFPGA';
Tdes=0.6;
intTime= 42;
plotOpt=1;
printOpt=0;
scanNum=1;


[XmatTemp, YmatTemp, tb1, tb2,photGrid] = linear2DFPGA_RT_Lapse(fileFolder, intTime);
    
[~, ~, Xshifts, Yshifts, Xangles, Yangles, ~, ~, ~, ~, NIFFT,Xzero] = interpLinear_Lapse(plotOpt,XmatTemp, YmatTemp, tb1, tb2, Tdes,scanNum,photGrid,printOpt);