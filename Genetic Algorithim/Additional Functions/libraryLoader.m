% AUTHORS: Claire & Jack
% CREATED: Summer 2022
% PURPOSE: interface model_generator_v5 with modelbuilder or GenAlg&PatSea


function [K, rates, tijs, Nstates, Nparam, param_strings, boundsArray, maxMutationArray,maxMutationCountsArray, minMutationCountsArray, atEdges, boundsArray_noAdj,model_name,dbCell,rates_str_woDB,rates_idx_woDB,rates_idx_DB] = libraryLoader(loadMode, model_set, modelNum , Nstates, tijs, c2ctrlAmp, c2ctrlTau, c4ctrlAmp, K, model_name,dbCell,rates_str_woDB,rates_idx_woDB,rates_idx_DB,Nparam,param_strings, boundsArray, maxMutationArray,maxMutationCountsArray, minMutationCountsArray, atEdges, boundsArray_noAdj,model_lin,model_loop1,model_loopN,model_brch)
% global Nstates Nparam 
% cc
% input parameters
% loadMode = 1;
% N = 4;
% tijs = rand(3*N,1); %ones(3*N,1);
% % model_lin, model_loop1, model_loopN, model_brch
% model_set = ['model_loop1'];
% modelNum = 1;



N = Nstates;

% if loadMode == 1
% libraryPath = '/Users/calbrecht/Dropbox/MATLAB_programs/claire_programs/from_Jack/updatedgenalgcodes';
% fName = 'modelLibraries';
% load([libraryPath filesep() fName filesep() 'lib_N' num2str(Nstates) '.mat'])
% end

idxs = reshape(linspace(1,N*N,N*N),[N,N]);
% idxs = reshape(linspace(1,Nstates*Nstates,Nstates*Nstates),[Nstates,Nstates]);

% input (OLD): (tijs, A, sigma_A, modelName,SelectedTimeRange, paramAdjs, c2ctrlAmp, c2ctrlTau, c4ctrlAmp)
% output (OLD): [K, A, rates, tijs, Nstates, Nparam, sigma_A, param_strings, boundsArray, maxMutationArray,maxMutationCountsArray, minMutationCountsArray, atEdges, boundsArray_noAdj]

% NEW:
% input: (tijs, c2ctrlAmp, c2ctrlTau, c4ctrlAmp)
% output: [K, rates, tijs, Nstates, Nparam, param_strings, boundsArray, maxMutationArray,maxMutationCountsArray, minMutationCountsArray, atEdges, boundsArray_noAdj]



if loadMode == 1
    if strcmp(model_set, 'model_lin')
        [model_name, K, rates_idx_woDB, rates_str_woDB, chain, bool_mat] = model_lin{modelNum,:};
    elseif strcmp(model_set, 'model_loop1')
        [model_name, K, rates_idx_woDB, rates_str_woDB, chain, bool_mat_loop, rates_idx_DB, rates_str_DB, dbCell, cxns_all, cxns_cur, numLoops] = model_loop1{modelNum,:};
    elseif strcmp(model_set, 'model_loopN')
        [model_name, K, rates_idx_woDB, rates_str_woDB, chain, bool_mat_loopN, rates_idx_DB, rates_str_DB, dbCell, cxns_all, cxns_cur, numLoops] = model_loopN{modelNum,:};
    elseif strcmp(model_set, 'model_brch')
        [model_name, K, rates_idx_woDB, rates_str_woDB, chain, bool_mat_brch, poss_brchs, brch_cxns] = model_brch{modelNum,:};
    end
    
    %%%% APPLY FILTERS TO LOADED MODELS? %%%%


    lowerBound_tijs=1e-6;
    upperBound_tijs=80;
    tijs_bound = [lowerBound_tijs, upperBound_tijs];
%     pointsPerDecade=25;
%     totalRands=length(rates_str_woDB); 
%     [valuesArray,finalLogArray] = logRandUniform(lowerBound_tijs, upperBound_tijs,pointsPerDecade,totalRands); 
%     tijs = valuesArray;
    if N == 3
        A1_bounds = [0.001,0.08];%low amp
        A2_bounds = [0.08,0.25];%2nd lowest
        A3_bounds = [0.25,0.4];%2nd high amp
        A_bounds = [A1_bounds; A2_bounds; A3_bounds];
        str_As = {'A1','A2','A3'};
        str_sigs = {'sigma_A1','sigma_A2','sigma_A3'};

    elseif N == 4
        A1_bounds = [0.001,0.055];%low amp
        A2_bounds = [0.055,0.12];%2nd lowest
        A3_bounds = [0.12,0.25];%2nd high amp
        A4_bounds = [0.25,0.4]; %1st high amp
        A_bounds = [A1_bounds; A2_bounds; A3_bounds; A4_bounds];
        str_As = {'A1','A2','A3','A4'};
        str_sigs = {'sigma_A1','sigma_A2','sigma_A3','sigma_A4'};
    %     str_sigs = ['sig1','sig2','sig3','sig4'];
    elseif N == 5
        A1_bounds = [0.001,0.055];%low amp
        A2_bounds = [0.06,0.14];%2nd lowest
        A3_bounds = [0.13,0.2];%2nd high amp
        A4_bounds = [0.25,0.4]; %1st high amp
        A5_bounds = [0.25,0.4]; 
        A_bounds = [A1_bounds; A2_bounds; A3_bounds; A4_bounds; A5_bounds];
        str_As = {'A1','A2','A3','A4', 'A5'};
        str_sigs = {'sigma_A1','sigma_A2','sigma_A3','sigma_A4', 'sigma_A5'};
    %     str_sigs = ['sig1','sig2','sig3','sig4'];
    elseif N == 6
        A1_bounds = [0.5,0.7];            % HIGH fret State
        A2_bounds = [0.25,0.65];        % Med FRET state
        A3_bounds = [0.25,0.65];        % Med FRET state
        A4_bounds = [0.2,0.5];       % Med FRET state
        A5_bounds = [0.2,0.5]; %[0.01,0.4]; %[0.2,0.4];          % Low FRET state
        A6_bounds = [0.01,0.4]; %[0.2,0.4];          % Low FRET state
        A_bounds = [A1_bounds; A2_bounds; A3_bounds; A4_bounds; A5_bounds; A6_bounds];
        str_As = {'A1','A2','A3','A4', 'A5','A6'};
        str_sigs = {'sigma_A1','sigma_A2','sigma_A3','sigma_A4', 'sigma_A5','sigma_A6'};
    end
    lowerBound_sigs = ones(N,1)*0.025;
    upperBound_sigs = ones(N,1)*0.12;
    sig_bounds = [lowerBound_sigs, upperBound_sigs];

    boundsArray = [repmat(tijs_bound,length(rates_idx_woDB),1); A_bounds; sig_bounds];
    boundsArray_noAdj = boundsArray;
    atEdges = zeros(size(tijs));

    Nparam = length(boundsArray);
    Nstates = N;
    maxMutationCountsArray = zeros(1,Nparam);
    minMutationCountsArray = zeros(1,Nparam);
    if length(tijs) > (Nparam - (2 * Nstates))
        tijs = tijs(1:(Nparam - (2 * Nstates)));
    end

    Max_mut_factor = .2;
    sigma_A_mutate = 0.05;
    maxMutationArray = [Max_mut_factor * boundsArray(1:Nparam-(2*Nstates),2) / 2; ...
                    0.1 * ( boundsArray(Nparam-(2*Nstates)+1:Nparam-(Nstates),2) - boundsArray(Nparam-(2*Nstates)+1:Nparam-(Nstates),1));...
                    sigma_A_mutate * (boundsArray(Nparam-(Nstates)+1:end,2) - boundsArray(Nparam-(Nstates)+1:end,1))];%...

    param_strings = [rates_str_woDB(:)', str_As{:}, str_sigs{:}];           
end

if strcmp(model_set, 'model_lin') || strcmp(model_set, 'model_brch')
     % do this everytime w or wo loadmode
    rates = 1./tijs;

    %     K(diag(idxs)) = 0;
    K = zeros(Nstates);
    K(rates_idx_woDB) = rates;
    diagK = -sum(K,1);
    K(diag(idxs)) = diagK;
else

    rates = 1./tijs;

%     K(diag(idxs)) = 0;
    K = zeros(Nstates);
    K(rates_idx_woDB) = rates;

    for q = 1:size(dbCell,1)
        for_prod_idx_arr = dbCell{q,1};
        rev_prod_idx_arr = dbCell{q,2};
        for i = 1:length(rates_idx_DB)
            % Need to update prod arrays in modelGen for multiple loops!
            if length(for_prod_idx_arr) > length(rev_prod_idx_arr)
                rate_DB = prod(K(for_prod_idx_arr)) / prod(K(rev_prod_idx_arr));
            else
                 rate_DB = prod(K(rev_prod_idx_arr)) / prod(K(for_prod_idx_arr));
            end
            K(rates_idx_DB) = rate_DB;
        end
    end
    
    diagK = -sum(K,1);
    K(diag(idxs)) = diagK;
    if sum(K,1) ~= zeros(1,N)
        error(['ERROR: columns of K do not sum to zero!'])
    end
    
end



% global varyCtrlMode NctrlParams ctrlMutateAmpc2 ctrlMutateTau ctrlMutateAmpc4% added to vary the control timescale

% if varyCtrlMode == 1
%     NctrlParams = 3;
%     ctrlMutate = 0.1; 
%     c2ctrlAmp = 9.25e-3;
%     c2ctrlTau = 2.41;
if loadMode == 1
    c2ctrlAmp_bounds = c2ctrlAmp * [0.25, 1.75]; %[1-0.1, 1+0.1];
    c2ctrlTau_bounds = c2ctrlTau * [0.5, 1.5];
    c4ctrlAmp_bounds = c4ctrlAmp * [0.25, 1.75]; %[1-0.5, 2];
    
  
     boundsArray = [boundsArray; c2ctrlAmp_bounds; c2ctrlTau_bounds; c4ctrlAmp_bounds];
     
     Nparam = length(boundsArray);
     
    ctrlMutateAmpc2 = 0.1;
    ctrlMutateTau = 0.1;
    ctrlMutateAmpc4 = 0.1;
    ctrlMMA = [ctrlMutateAmpc2 * (boundsArray(Nparam-2, 2) - boundsArray(Nparam-2, 1)); ctrlMutateTau * (boundsArray(Nparam-1, 2) - boundsArray(Nparam-1, 1)); ctrlMutateAmpc4 * (boundsArray(Nparam, 2) - boundsArray(Nparam, 1))];
     
     maxMutationArray = [maxMutationArray; ctrlMMA];
     
     param_strings = {param_strings{:}, 'c2ctrlAmp','c2ctrlTau', 'c4ctrlAmp'};
     
%      all_params = [all_params];% c2ctrlAmp; c2ctrlTau; c4ctrlAmp];
%      all_params_str = {all_params_str{:},'c2ctrlAmp','c2ctrlTau', 'c4ctrlAmp'};

end
% end
end

%      maxMutationArray = [maxMutationArray; ...
%          ctrlMutateAmpc2 * (boundsArray(Nparam-2, 2) - boundsArray(Nparam-2, 1)); ...
%          ctrlMutateTau * (boundsArray(Nparam-1, 2) - boundsArray(Nparam-1, 1)); ...
%          ctrlMutateAmpc4 * (boundsArray(Nparam, 2) - boundsArray(Nparam, 1))];
     
     
