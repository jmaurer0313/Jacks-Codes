%__________________________________________________________________________
% AUTHOR:   Brett Israels and Claire Albrecht
%
% CREATED:  November 2019
%
% PURPOSE:  Solve the master equation for a rate matrix K numerically
%
% INPUT     K (rate matrix)
%
% OUTPUT:   P (matrix of conditional probabilites)
%           P(j,i,t) is the probability of going from i -> j in time t
%
% MODLOG:  BI 20191112 Inverse done with '\' not inv() (faster)
%          BI Replaces k2P to better match notation of paper
%          CA created numerical option
%
%__________________________________________________________________________

%P is the numerical NxNxt matrix of conditional probabilities
%V is the mode matrix
%p i the NxN matrix of conditional probabilities in terms of 't's

% function [P, V, p] = K2P_num(K,time)   %%% For the symbolic mode %%%
function [P, V, K, time] = K2P(K,time)
clockMode = 0;
verboseMode = 0;
switch nargin
    case 0
        NumStates = 8;
        clockMode = 0;
        %                 Simulate a K matrix
        switch NumStates
            case 2
                k12 = 12; k21 = 21;
                K = [ -k12 , k21;...
                    k12, -k21];
                
            case 3
                
                %Cyclical mode
                k12 = 12; k13 = 13; k21 = 21; k31 = 31; k23 = 23;
                k32 = k12*k23*k31/(k13*k21);
                
                K = [(-k12 - k13), k21, k31;...
                    k12, (-k21 - k23 ), k32;...
                    k13, k23, (-k31-k32);];
            case 8
                
                %                 K = [...
                %                     -(12 + 13), 21, 31, 0, 0, 0, 0, 0;...
                %                     12, -(21 + 23 + 24 + 25), (12*23*31/(13*21)), 42, 52, 0, 0, 0;...
                %                     13, 23, -(31 + ((12*23*31/(13*21)))), 0, 0, 0, 0, 0;...
                %                     0, 24, 0, -42, 0, 0, 0, 0;...
                %                     0, 25, 0, 0, -(52 + 56), 65, 0, 0;...
                %                     0, 0, 0, 0, 56, -(65 + 67 + 68), 76, 86;...
                %                     0, 0, 0, 0, 0, 67, -(76+78), 87;...
                %                     0, 0, 0, 0, 0, 68, 78, -(86 + 87);...
                %                     ];
                [K,~,time,rates] = paramSim_8state();
        end
        
        Npts = 243;
        time = [0:9,logspace(1,log10(3e6),Npts)]/1e6;
        
    case 1 %If you give it 1 arguement , assume the gp32 = 0.
        time = time_sim;
end
%%
% Determine the number of states in the system
NumStates = length(K);

%--------------------------------------------------------------------------
% Calculate the Eigenvalues and Eigenvectors of K matrix
%--------------------------------------------------------------------------
if clockMode == 1
    tic
    %     disp(['N = ' num2str(N)]);
end
[Evec, Lam_unsorted] = eig(K);
[Lam,ind] = sort(diag(Lam_unsorted),'descend'); %Lam is numerical eigenvalue
Lam(1)=0;
V = Evec(:,ind);
% V is the mode matrix (Each column is an eigenvector)
if clockMode == 1
    elapsedTime = toc;
    disp(['     Time to Calculate the eigenvalues/eigenvectors of K = ' num2str(elapsedTime) ' seconds']);
    % % (8state: < 1 ms)
end

%--------------------------------------------------------------------------
% Find the inverse of the matrix of eigenvectors
%--------------------------------------------------------------------------
if clockMode == 1
    tic
end
V_inv = inv(V);
if clockMode == 1
    elapsedTime = toc;
    disp(['     Time to find inverse of V = ' num2str(elapsedTime) ' seconds']);
    % % (8state: < 0 µs)
end

%--------------------------------------------------------------------------
% Use identity matrix to define the Initial conditions (Boundary conditions)
%--------------------------------------------------------------------------
InitCond = eye(NumStates);%Make an  identity matrix

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% For numerical P(t) calculation:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%----------------------------------------------------------------------
% Calculate the exponential time dependence (using time vector)
%----------------------------------------------------------------------
if clockMode == 1
    tic
end


Lam = diag(Lam);
sizeLam = size(Lam);

if sizeLam(2) == 1
    exp_LamTime = zeros([NumStates NumStates length(time)]);
    for i = 1:length(time)
        LamTime = Lam .* time(i);
        exp_LamTime_temp = diag(exp(Lam .* time(i)));
        exp_LamTime(:,:,i) = exp_LamTime_temp;
        %                 exp_LamTime = [exp_LamTime, exp_LamTime_temp];
        %             exp_LamTime = reshape(exp_LamTime, [N, N, length(time)]);
    end
    
elseif sizeLam(2) == NumStates
    Lam = diag(Lam); % Make the diagonal into a vector (Nx1)
    exp_LamTime = zeros([NumStates NumStates length(time)]);
    for i = 1:length(time)
        LamTime = Lam .* time(i);
        exp_LamTime_temp = diag(exp(Lam .* time(i)));
        exp_LamTime(:,:,i) = exp_LamTime_temp;
    end
    
    if clockMode == 1
        elapsedTime = toc;
        disp(['Time to evaluate exp_LamTime matrix numerically = ' num2str(elapsedTime)]);
    end
end

%%
%--------------------------------------------------------------------------
% Now, numerically Calculate P(t):
%--------------------------------------------------------------------------
% Calculate P(t) using the NxNxlength(time) exp_LamTime matrix
P = zeros([NumStates NumStates length(time)]);
if clockMode == 1
    tic
end
for i = 1:length(time)
    P_temp = V * exp_LamTime(:,:,i) * V_inv * InitCond;
    P(:, :, i) = P_temp;
end
%     P = reshape(P, [N N length(time)]);
if clockMode == 1
    elapsedTime = toc;
    disp(['Time to evaluate P matrix numerically as function of t = ' num2str(elapsedTime)]);
end

% Peq = diag(P(:,:,end));
