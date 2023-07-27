%--------------------------------------------------------------------------
% AUTHOR: Claire Albrecht & Brett Israels
%
% NAME: histMaker_Nstate.m
%
% CREATED:  November 2019
%
% PURPOSE:  Create histograms for any model
%
% INPUT:    (1) Conditional Prpbability matrix (P, from k2P)
%           (2) The FRET values for the model  (A, a vector)
%
% OUTPUT:   (1) FRET Histogram
%
% MODIFICATIONs:
% BI 20200116 (1) Changed plot mode to plotInternalMode
% BI 20200116 (2) Added FRET_bins as an input
%--------------------------------------------------------------------------
function [Peq, denom_hist_sim, hist_sim, Peq_names] = histMaker_Nstate_v3normCorrect(P, A, sigma_A, FRET_bins)

verboseMode = 0;
clockMode = 0;

% disp('**** Now in histMaker_Nstate');
% Calculate equilibrium probabilities
Peq = diag(P(:,:,end));

if verboseMode == 1
    if (sum(Peq) - 1)  < 1e-9
        disp('Equilibrium probabilities sum to 1! (within 1e-9 error)')
    else
        disp('Equilibrium probabilities do NOT sum to 1! (within 1e-9 error)')
        disp('The difference from one is:')
        disp(sum(Peq)  - 1)
        
    end
end

NumStates = length(Peq);  % Number of states

%--------------------------------------------------------------------------
% Need a list of Peq names for the legend of these histograms!
%--------------------------------------------------------------------------
Peq_sym = sym('Peq_',[NumStates,1]);            % Define list of strings for the legend of the Peq fits
Peq_char = [];
Peq_name = [];
Peq_names = [];
for stateIDX = 1:NumStates
    Peq_char = char(Peq_sym(stateIDX));        % Turn syms to char in a list
    Peq_name = [Peq_name; Peq_char];    % concatenate the list of chars
    Peq_names = [Peq_names; strrep(Peq_name(stateIDX,:),'eq','^{eq}')];  % Make the eq a super script
end

%% Calculate total histogram distribution
if clockMode == 1
    tic
end
% clear hist_sim denom_hist_sim hist_sim_mat
% hist_sim_mat = [];
hist_sim = 0;
% hist_sim_mat = 0;

for stateIDX = 1:length(Peq)
    %     hist_sim_temp = Peq(stateIDX) * exp(-(((FRET_bins - A(stateIDX))/sigma_A(stateIDX)).^2));
    
%     disp(['Peq shape: ',num2str(size(Peq))])
%     disp(['A shape: ',num2str(size(A))])
%     disp(['sigma_A shape: ',num2str(size(sigma_A))])
    
    %Nov 14, 2020 Use the new definition of  Histogram which has correct
    %normalization
    hist_sim_temp = Peq(stateIDX)/(sqrt(2*pi)*sigma_A(stateIDX))*exp(-((FRET_bins-A(stateIDX))/(sqrt(2)*sigma_A(stateIDX))).^2);
    
    %     hist_sim_mat = [hist_sim_mat; hist_sim_temp];
    hist_sim = hist_sim + hist_sim_temp;
end

denom_hist_sim = sum(hist_sim(:));

% hist_sim = sum(hist_sim_mat)./denom_hist_sim;
hist_sim = hist_sim./denom_hist_sim;
if clockMode == 1
    elapsedTime = toc;
    disp(['Time to calculate calculate Histogram for an ' num2str(NumStates) ' state model: ' num2str(elapsedTime) ]);
end

hist_sim = hist_sim';


end
