% GOAL: loop over the output folders from the GenAlg (later on
% incorporating patSea results) to make a table of the names, top 3 chi2,
% Nparams, etc. For a single model at a time have a variable option to view
% the histogram of the parameters in the best fits for the top 3 + higher
% chi2s to get a sense of the distribution of guesses within the allowed
% intervals of model_builder() 
%%
% paste a directory here to be looped over (must be an "outputs" folder
% from the GenAlg space)
% startFold= 'C:\Users\Ryzen 5\Dropbox\TempHolder_C3_250usecC2\GenAlg_Data\E3E4\100mM_NaCl\genAlgData_UNorm_wCtrlwWeights_10ms_c2_250usec_v2\Outputs';
startFold= 'C:\Users\Ryzen 5\Dropbox\TempHolder_C3_250usecC2\GenAlg_Data\E5E6\20mM_NaCl_0mM_MgCl2\genAlgData_UNorm_wCtrlwWeights_10ms_c2_250usec_v2\Outputs';
% startFold= 'C:\Users\Ryzen 5\Dropbox\TempHolder_C3_250usecC2\GenAlg_Data\E3E4\300mM_NaCl\genAlgData_UNorm_wCtrlwWeights_10ms_c2_250usec_v2\Outputs';
% startFold= 'C:\Users\Ryzen 5\Dropbox\TempHolder_C3_250usecC2\GenAlg_Data\E5E6\20mM_NaCl\genAlgData_UNorm_wCtrlwWeights_10ms_c2_250usec_v2\Outputs';
% startFold= 'C:\Users\Ryzen 5\Dropbox\TempHolder_C3_modGenPS\TalapasDownloads\GenAlg_Data\E3E4\100mM_NaCl\genAlgData_UNorm_wCtrlwWeights_10ms_c2_10000usec_v2\Outputs';
% startFold= 'C:\Users\Ryzen 5\Dropbox\TempHolder_C3_250usecC2\GenAlg_Data\E5E6\20mM_NaCl_0mM_MgCl2\genAlgData_UNorm_wCtrlwWeights_10ms_c2_250usec_v4\Outputs';

% startFold= 'C:\Users\Ryzen 5\Dropbox\TempHolder_C3_250usecC2\GenAlg_Data\(+15)Dimer\100mM_NaCl\genAlgData_UNorm_wCtrlwWeights_10ms_c2_250usec_v2\Outputs'
% startFold= 'C:\Users\Ryzen 5\Dropbox\TempHolder_C3_250usecC2\GenAlg_Data\(+15)Dimer\100mM_NaCl\genAlgData_UNorm_wCtrlwWeights_10ms_c2_250usec_v2\Outputs';
% startFold= 'C:\Users\Ryzen 5\Dropbox\TempHolder_C3_250usecC2\GenAlg_Data\E5E6\20mM_NaCl_0mM_MgCl2\genAlgData_UNorm_wCtrlwWeights_10ms_c2_250usec_v2\Outputs';
% startFold= 'C:\Users\Ryzen 5\Dropbox\TempHolder_C3_modGenPS\E5E6\20mM_NaCl_0mM_MgCl2\Outputs';
% startFold= 'C:\Users\Jack\Dropbox\TempHolder_modGenPS\E5E6\100mM_NaCl\Outputs';
% startFold='C:\Users\Jack\Dropbox\TempHolder_modGenPS\TestFolder\Outputs';
% startFold= 'C:\Users\Jack\Dropbox\TempHolder_modGenPS_100msHist\weightScheme_1\E5E6\Outputs';
programName= 'GenAlg_Nstate_updated2';
addControlMode=1;
tau2=0; 
% options for post talapas plotter
numModels2Plot=8;
% Patsea models to include in weight calc
modelsToKeep=8;

genAlgMode=0;
printModels=1;

%toggles for hist stats on speoific models versus all the best PS outcomes
histThisModel=0; 
PSstatsHist=0;

% For the non Global output view, specify the suffix number of the output folder,
% within the PatSea version (weight scheme) to be examined
outputNumber=17;
% Set the number of states and the model sets you wish to view, if model
% set name is contained then all models within the set will be added to the
% sorted table
Nstates=4;

if Nstates==4
modelSetsToView= {'model_lin', 'model_loop1' 'model_loopN'}; 
else
modelSetsToView= {'model_lin', 'model_loop1'}; 
end


% specify the version number (v2 ,v3 , v4 etc) for PatSea - should reflect
% alternate weighting schemes used on Talapas when in altWeightScheme mode
PatSea_version_Num=6; 

% outputNumber=20;

% list of models for which a histogtam subplot will be produced to examine
% distribution of parameters 
% model_list={'linear1324' 'pyramid13' 'fully_connected_minus31' 'linear3124'};
% model_list={'pyramid42' 'pyramid41' 'linear3124' 'pyramid32'};
% model_list={'linear1243' 'pyramid43' 'fully_connected_minus43' 'linear2134' 'pyramid14'};
% model_list = {'chain2341' 'chain4312' 'chain4132_cxn43' 'chain4321' 'chain4231'};
% model_list={'chain3142' 'chain3142_cxn34' 'chain4312_cxn41' 'chain4321_cxn42'}; 
model_list = {};


% parameters specific to PatSea - Now that the outputs are numbered in the
% v2 folder
PatSea_v2_mode=1;

% Generate the model sets for N=X in order to view only certain sets
[model_lin,model_loop1,model_loopN,model_brch] = model_generator_v5(Nstates,'model_loopN');
linearNames=model_lin(1:size(model_lin,1),1);
loop1Names=model_loop1(1:size(model_loop1,1),1);
if Nstates>3
loopNNames=model_loopN(1:size(model_loopN,1),1);
end
namesToKeep=[];

if sum(contains(modelSetsToView, 'model_lin'))>0
    namesToKeep=[namesToKeep ; linearNames];
end 
if sum(contains(modelSetsToView, 'model_loop1'))>0
    namesToKeep=[namesToKeep ; loop1Names];
end
if sum(contains(modelSetsToView, 'model_loopN'))>0
    namesToKeep=[namesToKeep ; loopNNames];
end
    


c2Exp = @ (x, A, tau) A*exp(-x/tau);
c4Exp = @ (x, A, tau) 2*((A*exp(-x/tau)).*(A*exp(-x/tau)).');

% now cd to the start folder and loop over the contained directories to
% assembly the lists for table generation 
cd(startFold);

%%
outputList=dir('GenAlg*'); 
modelNamesGA={};
modelNamesPS={};
modelNamesPSGlob={};
modelNameOutTracker={}; 
modelNamesPSmissing={};
modelNamesPSmissingGlob={};
chi2PS=[];
chi2PSGlob=[];
outNumGlob=[]; 
outNumTotal=[];
lowChi2=[];
secLowChi2=[];
thirdLowChi2=[];
nParams=[]; 
totalIters=[];
iterOfLowChi2=[];
figIdx=[];

% holders for the weighted chi2 across al patSea runs
wChiC2=[];
wChiC4=[];
wChiHist=[]; 

% holders for PS model set sorted lists
chi2PS_ML=[];
modelNamesPS_ML={};
PSparamsTable=[];

wChiC2_ML=[];
wChiC4_ML=[];
wChiHist_ML=[]; 

%same for GA results
wChiC2_GA=[];
wChiC4_GA=[];
wChiHist_GA=[];

% holders for the individual surface chi2s
histChi2s=[];
c2Chi2s=[];
c4Chi2s=[];
varNames= {'Model Name' , 'Lowest Chi2' , '2nd Lowest Chi2' , '3rd Lowest Chi2' , ...
            'Num Params' , 'Num Iterations' , 'Iter of Low Chi2' } ; 
modelListCounter=0;

for i=1:length(outputList)
    for j=1:length(model_list)
        if strcmp(convertCharsToStrings(outputList(i).name), convertCharsToStrings([programName '_' model_list{j}]))
        params_cur =[];
        figIdx=j; 
        histThisModel=1;
        modelListCounter=modelListCounter+1;
        break; 
        else
        histThisModel=0; 
        end
    end
   
   cd([startFold filesep() outputList(i).name]); 
   if exist([startFold filesep() outputList(i).name filesep() 'lowestChiSquare'])
    modelNamesGA=[modelNamesGA ; outputList(i).name(length(programName)+1:end)];
   cd([startFold filesep() outputList(i).name filesep() 'lowestChiSquare']); 
   load([startFold filesep() outputList(i).name filesep() 'lowestChiSquare' filesep() 'fitInputData.mat']);
   load([startFold filesep() outputList(i).name filesep() 'lowestChiSquare' filesep() 'BestFitResults.mat']);
   %    calculate weight functions for the 3 surfaces to sort the models based
%    on best outcomes in the 3 surfaces individually 
   weightC2func = 1./(sqrt(C2_exp_x));% 2-pt. TCF weighting function
   weightC2func=weightC2func./(sum(weightC2func));
   weightHistFunc=ones(1,length(targetHistogram));
   wC4func = 1./(sqrt(C4_tau1range)).*(1./(sqrt(C4_tau1range)));
   wC4func=wC4func./(sum(sum(wC4func)));
   
%    calculate histogram chi2 for best genAlg outcome
    [~, ~, hist_sim, ~] = histMaker_Nstate_v3normCorrect(P, A, sigmas, Amp_bins);
%     histChi2s = [histChi2s ; mean(((hist_sim-reshape(targetHistogram,size(hist_sim))).^2).*weightHistFunc)];
    if size(hist_sim) == size(targetHistogram)
        histChi2s = [histChi2s ; mean(((hist_sim-targetHistogram).^2).*weightHistFunc)];
    elseif size(hist_sim') == size(targetHistogram)
        histChi2s = [histChi2s ; mean(((hist_sim'-targetHistogram).^2).*weightHistFunc)];
    end
    
%     calculate the c2 chi2 for best genAlg outcome
    c2Control = c2Exp(C2_exp_x, c2ctrlAmp, c2ctrlTau)*C2_exp_y(1);
    [C2_sim,~] = PA2C2(P,A,C2_exp_x,yoff,addControlMode,c2Control);
    c2Chi2s =[c2Chi2s ; mean(((C2_sim - reshape(C2_exp_y,1,length(C2_exp_y))).^2).*reshape(weightC2func,1,length(weightC2func)))];
   
%     lastly for the C4
    P_C4 = K2P(K, C4_tau1range);
    c4Control = c4Exp(C4_tau1range,c4ctrlAmp,c2ctrlTau)*C4_tau2eq0_exp(1,1);
    [C4_tau2eq0_sim,~] = PAK2C4(P_C4,A,K,C4_tau1range,tau2,zoff,addControlMode,c4Control);
    c4Chi2s = [c4Chi2s ; mean(mean(((C4_tau2eq0_sim - C4_tau2eq0_exp).^2).*wC4func))];
    
   numPlots=Nparam; 
   DimRow= floor(sqrt(numPlots)); 
   DimCol= floor(sqrt(numPlots)); 
   
   if ((DimRow*DimCol)< numPlots)
       DimRow=DimRow+1;
       if(DimRow*DimCol)< numPlots
           DimCol=DimCol+1;
       end
   end
%    params for the current model in the model select list
if histThisModel
   params_cur=[tijs , A , sigmas , c2ctrlAmp , c2ctrlTau, c4ctrlAmp]; 
   lowchi2Params=[tijs , A , sigmas , c2ctrlAmp , c2ctrlTau, c4ctrlAmp]; 
end
   lowChi2=[lowChi2 ; chisquared]; 
   nParams=[nParams ; Nparam]; 
   iterOfLowChi2=[iterOfLowChi2 ; iter]; 
   if exist('chisquared_Weighted_array','var') == 1
      wChiC2_GA=[wChiC2_GA,chisquared_Weighted_array(2)];
      wChiC4_GA=[wChiC4_GA,chisquared_Weighted_array(3)];
      wChiHist_GA=[wChiHist_GA,chisquared_Weighted_array(1)]; 
   end
   cd ..
    
   end
   
   if exist([startFold filesep() outputList(i).name filesep() 'secondLowestChiSquare'])
   cd([startFold filesep() outputList(i).name filesep() 'secondLowestChiSquare']); 
   load([startFold filesep() outputList(i).name filesep() 'secondLowestChiSquare' filesep() 'BestFitResults.mat']);
   secLowChi2=[secLowChi2; chisquared]; 
   if histThisModel
   secLowChi2Params=[tijs , A , sigmas , c2ctrlAmp , c2ctrlTau, c4ctrlAmp];
   params_cur=[params_cur ; tijs , A , sigmas , c2ctrlAmp , c2ctrlTau, c4ctrlAmp];
   end
   cd ..
   end
   
  if exist([startFold filesep() outputList(i).name filesep() 'thirdLowestChiSquare'])
   cd([startFold filesep() outputList(i).name filesep() 'thirdLowestChiSquare']); 
   load([startFold filesep() outputList(i).name filesep() 'thirdLowestChiSquare' filesep() 'BestFitResults.mat']);
   thirdLowChi2=[thirdLowChi2; chisquared]; 
   if histThisModel
   thirdLowChi2Params=[tijs , A , sigmas , c2ctrlAmp , c2ctrlTau, c4ctrlAmp];
   params_cur=[params_cur ; tijs , A , sigmas , c2ctrlAmp , c2ctrlTau, c4ctrlAmp]; 
   end  
   cd ..
  end
  
   if exist([startFold filesep() outputList(i).name filesep() 'higherChiSquares'])
   cd([startFold filesep() outputList(i).name filesep() 'higherChiSquares']);  
   higherChi2s=dir('fitResults_iter*'); 
   totalIters=[totalIters; length(higherChi2s)]; 
   if histThisModel
   for n=1:length(higherChi2s)
       load([startFold filesep() outputList(i).name filesep() 'higherChiSquares'  filesep() higherChi2s(n).name]);
       params_cur=[params_cur ; tijs , A , sigmas , c2ctrlAmp , c2ctrlTau, c4ctrlAmp ]; 
   end
   end
   
   cd ..
   end
   
   if PatSea_v2_mode==1
       
       if exist([startFold filesep() outputList(i).name filesep() 'PatSea_Nstate_polz_outputs_v' num2str(PatSea_version_Num) filesep() 'output' num2str(outputNumber) filesep() 'BestFitResults_PatSea_Nstate_polz.mat'])
        cd([startFold filesep() outputList(i).name filesep() 'PatSea_Nstate_polz_outputs_v' num2str(PatSea_version_Num) filesep() 'output' num2str(outputNumber)]); 
        load([startFold filesep() outputList(i).name filesep() 'PatSea_Nstate_polz_outputs_v' num2str(PatSea_version_Num) filesep() 'output' num2str(outputNumber) filesep() 'BestFitResults_PatSea_Nstate_polz.mat']);
        modelNamesPS=[modelNamesPS ; outputList(i).name]; 
        chi2PS= [chi2PS ; chisquared];  
        patSeaParams=[tijs , A , sigma_A , c2CtrlAmp , C2C4_ctrlTime, c4CtrlAmp ];
        PSparamsTable=[PSparamsTable ; patSeaParams(length(tijs)+1:end)];
        if exist('chisquared_Weighted_array','var') == 1
        wChiC2=[wChiC2,chisquared_Weighted_array(2)];
        wChiC4=[wChiC4,chisquared_Weighted_array(3)];
        wChiHist=[wChiHist,chisquared_Weighted_array(1)]; 
        end
        
        if sum(contains(namesToKeep,outputList(i).name(length(programName)+2:end)))>0
        chi2PS_ML=[chi2PS_ML ; chisquared];
        modelNamesPS_ML=[modelNamesPS_ML ; outputList(i).name];
            if exist('chisquared_Weighted_array','var') == 1
            wChiC2_ML=[wChiC2_ML,chisquared_Weighted_array(2)];
            wChiC4_ML=[wChiC4_ML,chisquared_Weighted_array(3)];
            wChiHist_ML=[wChiHist_ML,chisquared_Weighted_array(1)];
            end
        end
       else
        modelNamesPSmissing=[modelNamesPSmissing ; outputList(i).name]; 
        patSeaParams=[]; 
       end
       
   else
   
       if exist([startFold filesep() outputList(i).name filesep() 'PatSea_Nstate_polz_output'])
        cd([startFold filesep() outputList(i).name filesep() 'PatSea_Nstate_polz_output']); 
        load([startFold filesep() outputList(i).name filesep() 'PatSea_Nstate_polz_output' filesep() 'BestFitResults_PatSea_Nstate_polz.mat']);
        modelNamesPS=[modelNamesPS ; outputList(i).name]; 
        chi2PS= [chi2PS ; chisquared];  
        patSeaParams=[tijs , A , sigma_A , c2CtrlAmp , C2C4_ctrlTime, c4CtrlAmp ];
       else
        patSeaParams=[]; 
       end
   
   end
   
%    now add a block to tabulate the "global" lowest chi2 for a particualt
%    weighitng scheme after mulitple trials (possible bounds adjustment) -
%    this will be dictated by PatSea.... v2/v3/....vN
   
if PatSea_v2_mode==1
       
       if exist([startFold filesep() outputList(i).name filesep() 'PatSea_Nstate_polz_outputs_v' num2str(PatSea_version_Num)])
        cd([startFold filesep() outputList(i).name filesep() 'PatSea_Nstate_polz_outputs_v' num2str(PatSea_version_Num)]); 
        totalOutputs = length(dir('output*'));
        modelNameOutTracker=[modelNameOutTracker ; outputList(i).name];
        outNumTotal=[outNumTotal ; totalOutputs]; 
        for j=1:totalOutputs
         if exist([startFold filesep() outputList(i).name filesep() 'PatSea_Nstate_polz_outputs_v' num2str(PatSea_version_Num) filesep() 'output' num2str(j) filesep() 'BestFitResults_PatSea_Nstate_polz.mat']);
        load([startFold filesep() outputList(i).name filesep() 'PatSea_Nstate_polz_outputs_v' num2str(PatSea_version_Num) filesep() 'output' num2str(j) filesep() 'BestFitResults_PatSea_Nstate_polz.mat']);
        modelNamesPSGlob = [modelNamesPSGlob ; outputList(i).name]; 
        chi2PSGlob = [chi2PSGlob ; chisquared]; 
        outNumGlob = [outNumGlob ; j]; 
         end
        end
%         patSeaParams=[tijs , A , sigma_A , c2CtrlAmp , C2C4_ctrlTime, c4CtrlAmp ];
       
       end
end


   if histThisModel
   figure(figIdx);
   sgtitle(['Model Name = ' model_list{figIdx}]);
   for param_idx=1:Nparam
       subplot(DimRow,DimCol,param_idx)  
       if param_idx<= length(tijs)
       [~, hist_bins] = logRandUniform(boundsArray(param_idx,1),boundsArray(param_idx,2),50,1);
       else
       hist_bins=linspace(boundsArray(param_idx,1),boundsArray(param_idx,2),30);
       end
       histogram(params_cur(:,param_idx),'BinEdges',hist_bins);
       [N,~] = histcounts(params_cur(:,param_idx),hist_bins);
       if param_idx<= length(tijs)
       set(gca, 'xscale','log')
       end
       maxCountCur = max(N);
       hold on
       progressLine = line([lowchi2Params(param_idx) lowchi2Params(param_idx) ],[0  maxCountCur],...
           'Color','blue','LineStyle','--','LineWidth',1); 
       progressLine = line([secLowChi2Params(param_idx) secLowChi2Params(param_idx) ],[0  maxCountCur],...
           'Color','red','LineStyle','--','LineWidth',1); 
       progressLine = line([thirdLowChi2Params(param_idx) thirdLowChi2Params(param_idx) ],[0  maxCountCur],...
           'Color','green','LineStyle','--','LineWidth',1); 
       if ~isempty(patSeaParams)
       progressLine = line([patSeaParams(param_idx) patSeaParams(param_idx) ],[0  maxCountCur],...
           'Color','cyan','LineStyle','--','LineWidth',1); 
       end
       ylabel('Frequency');
       title(cellstr(param_strings(param_idx)));
   end
   end
   
   cd ..
   
   
   
   
end

% create tables for the GenAlg and PatSea outcomes
myTable = table(modelNamesGA,lowChi2,secLowChi2, thirdLowChi2,nParams,...
    iterOfLowChi2);
myTableSorted = sortrows(myTable,{'lowChi2','secLowChi2'});
myTablePS = table(modelNamesPS,chi2PS);
myTableSortedPS = sortrows(myTablePS,{'chi2PS'});

if PSstatsHist
%sorted list by chi2 of the PS params for subplot histograms
if Nstates==4
myTablePSparams=table(modelNamesPS,chi2PS,PSparamsTable(:,1),PSparamsTable(:,2), ...
     PSparamsTable(:,3),PSparamsTable(:,4),PSparamsTable(:,5),PSparamsTable(:,6), ...
     PSparamsTable(:,7), PSparamsTable(:,8), PSparamsTable(:,9), PSparamsTable(:,10), ...
     PSparamsTable(:,11)); 

elseif Nstates==5
 myTablePSparams=table(modelNamesPS,chi2PS,PSparamsTable(:,1),PSparamsTable(:,2), ...
     PSparamsTable(:,3),PSparamsTable(:,4),PSparamsTable(:,5),PSparamsTable(:,6), ...
     PSparamsTable(:,7), PSparamsTable(:,8), PSparamsTable(:,9), PSparamsTable(:,10), ...
     PSparamsTable(:,11), PSparamsTable(:,12), PSparamsTable(:,13));   
end 

myTablePSparamsSorted= sortrows(myTablePSparams,{'chi2PS'});


   NparamNoTijs= size(PSparamsTable,2);
   adjBoundsArr=boundsArray(length(boundsArray)-NparamNoTijs+1:end,:);
   adjParamStrings=param_strings(length(boundsArray)-NparamNoTijs+1:end);
   figure();
   sgtitle(['Top ' num2str(modelsToKeep) ' models Optimized PS param distributions']);
   for param_idx=1:NparamNoTijs
       subplot(DimRow,DimCol,param_idx)  
%        hist_bins=linspace(adjBoundsArr(param_idx,1),adjBoundsArr(param_idx,2),30);
       histogram(myTablePSparamsSorted{2:modelsToKeep+1,param_idx+2},30);
%        histogram(myTablePSparamsSorted{2:modelsToKeep+1,param_idx+2},'BinEdges',hist_bins);
       ylabel('Frequency');
       title(cellstr(adjParamStrings(param_idx)));
   end
   end




%sorted table for the model set list variant
myTablePS_ML = table(modelNamesPS_ML,chi2PS_ML);
myTableSortedPS_ML = sortrows(myTablePS_ML,{'chi2PS_ML'});

% create a table for the "global" outputs list
myTablePSGlobal=table(modelNamesPSGlob,chi2PSGlob,outNumGlob);
myTablePSGlobalSorted=sortrows(myTablePSGlobal,{'chi2PSGlob'}); 

% tables for the hist, c2 , c4 sorted GenAlg
myTableHist = table(modelNamesGA,histChi2s);
myTableC2 = table(modelNamesGA,c2Chi2s);
myTableC4 = table(modelNamesGA,c4Chi2s);

HistTableSorted=sortrows(myTableHist,{'histChi2s'});
C2TableSorted=sortrows(myTableC2,{'c2Chi2s'});
C4TableSorted=sortrows(myTableC4,{'c4Chi2s'});

myHistNames = HistTableSorted.modelNamesGA(:);
myC2Names = C2TableSorted.modelNamesGA(:);
myC4Names = C4TableSorted.modelNamesGA(:);
myHistChi2s = HistTableSorted.histChi2s(:);
myC2Chi2s = C2TableSorted.c2Chi2s(:);
myC4Chi2s = C4TableSorted.c4Chi2s(:);

%for PAtSea
avgWeightChiC2 = mean(wChiC2);
avgWeightChiC4 = mean(wChiC4);
avgWeightChiHist = mean(wChiHist);
PatSea_AvgWeights=[avgWeightChiHist; avgWeightChiC2 ; avgWeightChiC4];

if length(wChiHist)==length(modelNamesGA)
% tables for the hist, c2 , c4 sorted PatSea
wChiHist=reshape(wChiHist,length(wChiHist),1);
wChiC2=reshape(wChiC2,length(wChiC2),1);
wChiC4=reshape(wChiC4,length(wChiC4),1);
myTableHist_PS = table(modelNamesGA,wChiHist);
myTableC2_PS = table(modelNamesGA,wChiC2);
myTableC4_PS = table(modelNamesGA,wChiC4);

HistTableSorted_PS=sortrows(myTableHist_PS,{'wChiHist'});
C2TableSorted_PS=sortrows(myTableC2_PS,{'wChiC2'});
C4TableSorted_PS=sortrows(myTableC4_PS,{'wChiC4'});
end 

%For GA
avgWeightChiC2_GA = mean(wChiC2_GA);
avgWeightChiC4_GA = mean(wChiC4_GA);
avgWeightChiHist_GA = mean(wChiHist_GA);
GenAlg_AvgWeights=[avgWeightChiHist_GA; avgWeightChiC2_GA ; avgWeightChiC4_GA];


finalSurfaceTable= table(myHistNames,myHistChi2s, myC2Names,myC2Chi2s, myC4Names,myC4Chi2s);

% write those tables to excel files for convienent viewing 
cd(startFold);
writetable(myTableSorted,'myTableSorted.xls');
writetable(myTableSortedPS,'myTableSortedPS.xls');
writetable(finalSurfaceTable,'surfaceSortedTable.xls');

% section to pass back over the top X fits from PatSea, and record their
% chi2weighted arrays, as a better measure of how much scalrs need to be
% adjusted when the data is actually close to the model. 
if size(myTableSortedPS,1)>0
myListToKeep=myTableSortedPS((1:modelsToKeep),1); 
wChiC2_PSkeep=[];
wChiC4_PSkeep=[];
wChiHist_PSkeep=[]; 
topFitsKMask=zeros(Nstates,Nstates); 

for i=1:max(size(myListToKeep))
    curName=cell2str((myListToKeep.modelNamesPS(i,1)));
    curName=curName(3:end-3); 
    if exist([startFold filesep() curName filesep() 'PatSea_Nstate_polz_outputs_v' num2str(PatSea_version_Num) filesep() 'output' num2str(outputNumber)])
        cd([startFold filesep() curName filesep() 'PatSea_Nstate_polz_outputs_v' num2str(PatSea_version_Num) filesep() 'output' num2str(outputNumber)]); 
        load([startFold filesep() curName filesep() 'PatSea_Nstate_polz_outputs_v' num2str(PatSea_version_Num) filesep() 'output' num2str(outputNumber) filesep() 'BestFitResults_PatSea_Nstate_polz.mat']);
%         modelNamesPS=[modelNamesPS ; curName]; 
        topFitsKMask=topFitsKMask+(K>0); 
        if exist('chisquared_Weighted_array','var') == 1
        wChiC2_PSkeep=[wChiC2_PSkeep,chisquared_Weighted_array(2)];
        wChiC4_PSkeep=[wChiC4_PSkeep,chisquared_Weighted_array(3)];
        wChiHist_PSkeep=[wChiHist_PSkeep,chisquared_Weighted_array(1)]; 
        end
        
    end
end

%for PAtSea kept models
avgWeightChiC2_PSkeep = mean(wChiC2_PSkeep);
avgWeightChiC4_PSkeep = mean(wChiC4_PSkeep);
avgWeightChiHist_PSkeep = mean(wChiHist_PSkeep);
PatSea_AvgWeights_PSkeep=[avgWeightChiHist_PSkeep; avgWeightChiC2_PSkeep ; avgWeightChiC4_PSkeep];
end

if printModels
if genAlgMode
   tableInput = myTableSorted;
else
   tableInput = myTableSortedPS;
end
  [newChi2List] = postTalapasPlotter_v2(numModels2Plot, tableInput, startFold, programName,genAlgMode,PatSea_version_Num,outputNumber);

% postTalapasPlotter(numModels2Plot, tableInput, startFold, programName,genAlgMode,PatSea_version_Num,outputNumber)
end
% 
% figure()
% surf(topFitsKMask./modelsToKeep)
% caxis([0 1]);
% colorbar;
% colormap('jet');
% view(0,90);

figure()
myBar=bar3(topFitsKMask./modelsToKeep);
for k = 1:length(myBar)
zdata = myBar(k).ZData;
myBar(k).CData = zdata;
myBar(k).FaceColor = 'interp';
end 
% set(gca,'XTickLabel',linspace(min(intensity),max(intensity),10))
% set(gca,'XTickLabel',edges_counts(1:7))
% set(gca,'YTickLabel',edges_amp*10)
view([0 90])
% ylim([min(intensity),max(intensity)]);
colormap('jet');
caxis([0 1]);
colorbar;
ylabel(['State 1']);
xlabel(['State 2']);