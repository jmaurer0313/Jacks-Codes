% AUTHOR:   Claire Albrecht
% CREATED:  September 2021
%
% PURPOSE:  Define a function to build the model and re-define the model at every
%           necessary step in the code.
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% list of model names:
%--------------------------------------------------------------------------
% 1. ssDNA_3state_cyclical (for FRET data with 3 states... no bound adjuster array included)
% 2. fully_connected (for polz data)
% 3. linear1234 (for other permutations, change ordering of numbers after
% 4. linear3124
% 5. center_1
% 6. pyrimid12
% 7. outer_loop
% 8. fully_connected_minus14 (fully connected except for the 14 connection
% 9. outer_loop_plus13 (outer loop plus one cross connection) 
% Note: all are for 4 states unless states otherwise
%--------------------------------------------------------------------------
function [K, A, rates, tijs, Nstates, Nparam, param_strings, boundsArray, maxMutationArray,maxMutationCountsArray, minMutationCountsArray, atEdges, boundsArray_noAdj] = model_builder_Polz_v4_10ms(tijs, A, model_name,paramAdjs,c2ctrlAmp,c2ctrlTau,c4ctrlAmp)
%--------------------------------------------------------------------------
if isempty(paramAdjs)
    paramAdjs=zeros(1,length(tijs)); 
end
windowMag=0.2;
absLow=0.1e-6;
absHigh=50;
fixedDecades=1;



if strcmp(model_name, 'ssDNA_3state_cyclical') == 1
    %----------------------------------------------------------------------
    % Initialzing parameters
    %----------------------------------------------------------------------
    Nstates = 3;
    t12_bounds = [20e-6,1000e-3];  %Paramater #1 is high--> med
    t13_bounds = [1e-3,500e-3];    %Paramater #2 is high --> low
    t21_bounds = [20e-6,1000e-3];%Paramater #3 is med --> high
    t23_bounds = [1e-6,1000e-3];%Paramater #4 is med --> low
    t31_bounds = [20e-6,2000e-3];  %Paramater #5 is low --> Medium = [1e-3,10e-3];%Paramater #5 %t32 is low --> Medium
    % *t32 wll be determined by the other rates
    
    % A1_bounds = [0.65,0.85];%Paramater #6 % HIGH fret State
    A1_bounds = [0.04,0.08];%Paramater #6 % HIGH fret State
    A2_bounds = [0.08,0.25];%Paramater #7 % Med FRET state
    A3_bounds = [0.25,0.4];%Paramater #8 %Low FRET state
    
    boundsArray = [t12_bounds; t13_bounds; t21_bounds; t23_bounds; t31_bounds; A1_bounds; A2_bounds; A3_bounds];
    Nparam = length(boundsArray);
    
    maxMutationCountsArray = zeros(1,Nparam);
    minMutationCountsArray = zeros(1,Nparam);
    
    % Maximum mutation value
    Max_mut_factor = .2;
    t12_mutate = Max_mut_factor*t12_bounds(2)/2;%100e-6;
    t13_mutate = Max_mut_factor*t13_bounds(2)/2;%100e-6;
    t21_mutate = Max_mut_factor*t21_bounds(2)/2;%500e-6;
    t23_mutate = Max_mut_factor*t23_bounds(2)/2;%10e-6;
    t31_mutate = Max_mut_factor*t31_bounds(2)/2;%100e-6;
    A1_mutate = (A1_bounds(2) - A1_bounds(1))*0.1;%0.05;
    A2_mutate = (A2_bounds(2) - A2_bounds(1))*0.1;%0.05;0.05;
    A3_mutate = (A3_bounds(2) - A3_bounds(1))*0.1;%0.05;0.05;
    maxMutationArray = [t12_mutate;t13_mutate;t21_mutate;t23_mutate;t31_mutate;A1_mutate;A2_mutate;A3_mutate];
    
    sigma_A1 = 0.044;
    sigma_A2 = 0.077;
    sigma_A3 = 0.104;
    sigma_A = [sigma_A1; sigma_A2; sigma_A3];
    
    param_strings = {'t12','t13','t21','t23','t31','A1','A2','A3'};
    % Note: the order of the strings needs to be the same as the variable
    % list so we can extract the parameters properly.
    
    % tij_strings = param_strings(1:(Nparam - Nstates));
    % A_strings = param_strings((Nparam - Nstates + 1):Nparam);
    
    %----------------------------------------------------------------------
    % Update model
    %----------------------------------------------------------------------
    % [A1, A2, A3] = feval(@(x) x{:}, num2cell(As));
    if length(tijs) > (Nparam - 2*Nstates)
        tijs = tijs(1:(Nparam - 2*Nstates));
    end
    %     disp(['model_builder shape tijs: ',num2str(size(tijs))])
    rates = 1 ./ tijs;
    %     disp(['model_builder shape rates: ',num2str(size(rates))])
    [k12, k13, k21, k23, k31] = feval(@(x) x{:}, num2cell(reshape(rates,[1, length(rates)])));
    k32 = (k12 * k23 * k31)/(k13 * k21);
    
    K = [-(k12 + k13), k21,           k31;...
        k12,          -(k21 + k23),   k32;...
        k13,           k23,          -(k31 + k32)];

    %--------------------------------------------------------------------------
    
% elseif strcmp(model_name, 'test_4state') == 1 % this is the fully connected model
elseif strcmp(model_name, 'fully_connected') == 1 % this is the fully connected model
    Nstates = 4;
    t12_bounds = [2e-6,1000e-3];  %Paramater #1 is high--> med
    t13_bounds = [2e-3,1000e-3];    %Paramater #2 is high --> low
    t14_bounds = [2e-6,1000e-3];
    t21_bounds = [2e-6,1000e-3];%Paramater #3 is med --> high
    t23_bounds = [2e-6,1000e-3];%Paramater #4 is med --> low
    t24_bounds = [2e-6,1000e-3];
    t31_bounds = [2e-6,1000e-3];  %Paramater #5 is low --> Medium = [1e-3,10e-3];%Paramater #5 %t32 is low --> Medium
    t34_bounds = [2e-6,1000e-3];
    t43_bounds = [2e-6,1000e-3];
    % *t32 wll be determined by the other rates
    
    % A1_bounds = [0.65,0.85];%Paramater #6 % HIGH fret State
    A1_bounds = [0.04,0.055];%low amp
    A2_bounds = [0.08,0.13];%2nd lowest
    A3_bounds = [0.15,0.2];%2nd high amp
    A4_bounds = [0.25,0.32]; %1st high amp
    
   %Bounds for standard deviation
     sigma_A1_bounds = [0.028,0.042]; %Paramater #8 %Low Amp state standard deviation
     sigma_A2_bounds = [0.05,0.07];%Paramater #9 %Med Amp state standard deviation
     sigma_A3_bounds = [0.065,0.09];%Paramater #10 %Med2 Amp state standard deviation
     sigma_A4_bounds = [0.09,0.14];%Paramater #10 %High Amp state standard deviation
     
    boundsArray = [t12_bounds; t13_bounds; t14_bounds; t21_bounds; t23_bounds; t24_bounds; t31_bounds; t34_bounds; t43_bounds; ...
                 A1_bounds; A2_bounds; A3_bounds; A4_bounds; ...
                 sigma_A1_bounds; sigma_A2_bounds; sigma_A3_bounds ; sigma_A4_bounds];
             
    Nparam = length(boundsArray);
    boundsArray_noAdj = boundsArray;
    
    [boundsArrayAdj,atEdges]= boundsAdjuster_v3(boundsArray,windowMag,paramAdjs,absLow,absHigh,Nparam,Nstates,fixedDecades);
    boundsArray=boundsArrayAdj;
    
    maxMutationCountsArray = zeros(1,Nparam);
    minMutationCountsArray = zeros(1,Nparam);
    % Maximum mutation value
    Max_mut_factor = .2;
    maxMutationArray = [Max_mut_factor * boundsArray(1:Nparam-2*Nstates,2) / 2; ...
    0.1 * ( boundsArray(Nparam-2*Nstates+1:end,2) - boundsArray(Nparam-2*Nstates+1:end,1)) ];
%     sigma_A1 = 0.032;
%     sigma_A2 = 0.060;
%     sigma_A3 = 0.027;
%     sigma_A4 = 0.060;
%     sigma_A = [sigma_A1; sigma_A2; sigma_A3; sigma_A4];

    param_strings = {'t12','t13','t14','t21','t23','t24','t31','t34','t43','A1','A2','A3','A4','sig1','sig2','sig3','sig4'};
    % Note: the order of the strings needs to be the same as the variable
    % list so we can extract the parameters properly.
    % tij_strings = param_strings(1:(Nparam - Nstates));
    % A_strings = param_strings((Nparam - Nstates + 1):Nparam);
    %----------------------------------------------------------------------
    % Update model
    %----------------------------------------------------------------------
    if length(tijs) > (Nparam - 2*Nstates)
        tijs = tijs(1:(Nparam - 2*Nstates));
    end
    
    rates = 1 ./ tijs;
    
    [k12, k13, k14, k21, k23, k24, k31, k34, k43] = feval(@(x) x{:}, num2cell(reshape(rates,[1, length(rates)])));
    
    k32 = (k12 * k23 * k31)/(k13 * k21);
    k41 = (k31 * k14 * k43)/ (k13 * k34);
    k42 = (k24 * k43 * k31 * k12)/(k34 * k13 * k21);
    % [A1, A2, A3, A4] = feval(@(x) x{:}, num2cell(As));
    
    K = [-(k12 + k13 + k14), k21,               k31,               k41;...
        k12,              -(k21 + k23 + k24), k32,               k42;...
        k13,               k23,              -(k31 + k32 + k34), k43;...
        k14,               k24,                k34,             -(k41 + k42 + k43)];
   
    
% UPDATE: 03-16-2022 - change over this entry in modelBuilder for the 10ms
% case with NEW maxMutate array and NEW sigma fitting 
    
    elseif strcmp(model_name, 'outer_loop3124') == 1
     Nstates = 4;
%     t12_bounds = [10e-6,500e-3];  %Paramater #1 is high--> med    
%     t13_bounds = [1e-3,3000e-3];    %Paramater #2 is high --> low
%     t14_bounds = [100e-6,900e-3];
    t21_bounds = [10e-6,500e-3];%Paramater #3 is med --> high
    t13_bounds = [10e-6,500e-3];%Paramater #4 is med --> low
    t31_bounds = [10e-6,500e-3];
%     t31_bounds = [50e-6,1000e-3];  %Paramater #5 is low --> Medium = [1e-3,10e-3];%Paramater #5 %t32 is low --> Medium
    t24_bounds = [10e-6,500e-3];
    t42_bounds = [10e-6,500e-3];
    t34_bounds = [10e-6,500e-3];
    t43_bounds = [10e-6,500e-3];
    % *t32 wll be determined by the other rates
    
% IN THIS SCHEME, 2 AND 3 ARE THE LOW AMP MOSTLY SYMM/ANTI WHILE 1 AND 4
% ARE THE HIGH AP STATES
    % A1_bounds = [0.65,0.85];%Paramater #6 % HIGH fret State
    A1_bounds = [0.04,0.055];%low amp
    A2_bounds = [0.08,0.13];%2nd lowest
    A3_bounds = [0.15,0.2];%2nd high amp
    A4_bounds = [0.25,0.32]; %1st high amp  
    
    %Bounds for standard deviation
     sigma_A1_bounds = [0.028,0.042]; %Paramater #8 %Low Amp state standard deviation
     sigma_A2_bounds = [0.05,0.07];%Paramater #9 %Med Amp state standard deviation
     sigma_A3_bounds = [0.065,0.09];%Paramater #10 %Med2 Amp state standard deviation
     sigma_A4_bounds = [0.09,0.14];%Paramater #10 %High Amp state standard deviation
    
    boundsArray = [t21_bounds; t13_bounds; t31_bounds; t24_bounds; t42_bounds; t34_bounds; t43_bounds;...
        A1_bounds; A2_bounds; A3_bounds; A4_bounds; ...
        sigma_A1_bounds; sigma_A2_bounds; sigma_A3_bounds ; sigma_A4_bounds];
    
    Nparam = length(boundsArray);
    boundsArray_noAdj = boundsArray;
    
    [boundsArrayAdj,atEdges]= boundsAdjuster_v3(boundsArray,windowMag,paramAdjs,absLow,absHigh,Nparam,Nstates,fixedDecades);
    boundsArray=boundsArrayAdj;
    
    maxMutationCountsArray = zeros(1,Nparam);
    minMutationCountsArray = zeros(1,Nparam);
    
    
    % Maximum mutation value
    Max_mut_factor = .2;
    
    maxMutationArray = [Max_mut_factor * boundsArray(1:Nparam-2*Nstates,2) / 2; ...
        0.1 * ( boundsArray(Nparam-2*Nstates+1:end,2) - boundsArray(Nparam-2*Nstates+1:end,1)) ];
    
%     
    
    param_strings = {'t21','t13','t31','t24','t42','t34','t43','A1','A2','A3','A4','sig1','sig2','sig3','sig4'};
    % Note: the order of the strings needs to be the same as the variable
    % list so we can extract the parameters properly.
    
    % tij_strings = param_strings(1:(Nparam - Nstates));
    % A_strings = param_strings((Nparam - Nstates + 1):Nparam);
    
    %----------------------------------------------------------------------
    % Update model
    %----------------------------------------------------------------------
    if length(tijs) > (Nparam - 2*Nstates)
        tijs = tijs(1:(Nparam - 2*Nstates));
    end
    rates = 1 ./ tijs;
    
    [k21, k13, k31, k24, k42, k34, k43] = feval(@(x) x{:}, num2cell(reshape(rates,[1, length(rates)])));
    k12 = (k21 * k13 * k34 * k42)/(k24 * k43 * k31);
%     k41 = (k31 * k14 * k43)/ (k13 * k34);
%     k42 = (k24 * k43 * k31 * k12)/(k34 * k13 * k21);
    
    % [A1, A2, A3, A4] = feval(@(x) x{:}, num2cell(As));
    
    K = [-(k12 + k13),     k21,               k31,                  0;...
        k12,              -(k21 + k24),      0,                     k42;...
        k13,               0,              -(k31 + k34 ),            k43;...
        0,               k24,                k34,             -(k42 +k43)];
    
 elseif strcmp(model_name, 'outer_loop1324') == 1 % this is the fully connected model
    Nstates = 4;
    t13_bounds = [200e-6,50000e-3];  %Paramater #1 is high--> med
    t14_bounds = [200e-6,50000e-3];
    t23_bounds = [200e-6,50000e-3];%Paramater #3 is med --> high
    t24_bounds = [200e-6,50000e-3];%Paramater #4 is med --> low
    t31_bounds = [200e-6,50000e-3];
    t32_bounds = [200e-6,50000e-3];
    t42_bounds = [200e-6,50000e-3];
    % *t32 wll be determined by the other rates
    
   % IN THIS SCHEME, 2 AND 3 ARE THE LOW AMP MOSTLY SYMM/ANTI WHILE 1 AND 4
% ARE THE HIGH AP STATES
    % A1_bounds = [0.65,0.85];%Paramater #6 % HIGH fret State
    A1_bounds = [0.001,0.055];%low amp
    A2_bounds = [0.08,0.13];%2nd lowest
    A3_bounds = [0.13,0.2];%2nd high amp
    A4_bounds = [0.25,0.4]; %1st high amp  
    
    %Bounds for standard deviation
     sigma_A1_bounds = [0.028,0.042]; %Paramater #8 %Low Amp state standard deviation
     sigma_A2_bounds = [0.05,0.07];%Paramater #9 %Med Amp state standard deviation
     sigma_A3_bounds = [0.065,0.09];%Paramater #10 %Med2 Amp state standard deviation
     sigma_A4_bounds = [0.09,0.14];%Paramater #10 %High Amp state standard deviation
    
    boundsArray = [t13_bounds; t14_bounds; t23_bounds; t24_bounds;  t31_bounds; t32_bounds; t42_bounds; ...
    A1_bounds; A2_bounds; A3_bounds; A4_bounds; ...
    sigma_A1_bounds; sigma_A2_bounds; sigma_A3_bounds ; sigma_A4_bounds];
    
    
    Nparam = length(boundsArray);
    boundsArray_noAdj = boundsArray;
    
    [boundsArrayAdj,atEdges]= boundsAdjuster_v3(boundsArray,windowMag,paramAdjs,absLow,absHigh,Nparam,Nstates,fixedDecades);
    boundsArray=boundsArrayAdj;
    
    maxMutationCountsArray = zeros(1,Nparam);
    minMutationCountsArray = zeros(1,Nparam);
    
    
    % Maximum mutation value
    Max_mut_factor = .2;
    
    maxMutationArray = [Max_mut_factor * boundsArray(1:Nparam-2*Nstates,2) / 2; ...
        0.1 * ( boundsArray(Nparam-2*Nstates+1:end,2) - boundsArray(Nparam-2*Nstates+1:end,1)) ];
    
    
    param_strings = {'t13','t14','t23','t24','t31','t32','t42','A1','A2','A3','A4','sig1','sig2','sig3','sig4'};
    % Note: the order of the strings needs to be the same as the variable
    % list so we can extract the parameters properly.
    
    % tij_strings = param_strings(1:(Nparam - Nstates));
    % A_strings = param_strings((Nparam - Nstates + 1):Nparam);
    
    %----------------------------------------------------------------------
    % Update model
    %----------------------------------------------------------------------
    if length(tijs) > (Nparam - 2*Nstates)
        tijs = tijs(1:(Nparam - 2*Nstates));
    end
    rates = 1 ./ tijs;

    [k13, k14, k23, k24, k31, k32, k42] = feval(@(x) x{:}, num2cell(reshape(rates,[1, length(rates)])));
    k41 = (k14 * k42 * k23 * k31)/ (k13 * k32 * k24);
    
    
    K = [-(k13+k14),      0,                k31,                  k41;...
        0,              -(k23 +k24 ),       k32,                   k42;...
        k13,               k23,              -( k32 + k31),       0;...
        k14,               k24,                 0,                -(k41  + k42)];
    
    
   elseif strcmp(model_name, 'center1') == 1
     Nstates = 4;
    t12_bounds = [10e-6,1000e-3];  %Paramater #1 is high--> med
%     t13_bounds = [1e-3,3000e-3];    %Paramater #2 is high --> low
%     t14_bounds = [100e-6,900e-3];
    t21_bounds = [10e-6,1000e-3];%Paramater #3 is med --> high
    t13_bounds = [10e-6,1000e-3];%Paramater #4 is med --> low
%     t32_bounds = [100e-6,4000e-3];
    t31_bounds = [10e-6,1000e-3];  %Paramater #5 is low --> Medium = [1e-3,10e-3];%Paramater #5 %t32 is low --> Medium
    t14_bounds = [10e-6,1000e-3];
    t41_bounds = [10e-6,100e-3];
    % *t32 wll be determined by the other rates
    
% IN THIS SCHEME, 2 AND 3 ARE THE LOW AMP MOSTLY SYMM/ANTI WHILE 1 AND 4
% ARE THE HIGH AP STATES
    % A1_bounds = [0.65,0.85];%Paramater #6 % HIGH fret State
    A1_bounds = [0.04,0.075];%low amp
    A2_bounds = [0.1,0.15];%2nd lowest
    A3_bounds = [0.15,0.24];%2nd high amp
    A4_bounds = [0.3,0.4]; %1st high amp 
    
    boundsArray = [t12_bounds; t21_bounds; t13_bounds; t31_bounds; t14_bounds; t41_bounds; A1_bounds; A2_bounds; A3_bounds; A4_bounds];
    Nparam = length(boundsArray);
    boundsArray_noAdj = boundsArray;
    
    maxMutationCountsArray = zeros(1,Nparam);
    minMutationCountsArray = zeros(1,Nparam);
    
    % Maximum mutation value
    Max_mut_factor = .2;
    t12_mutate = Max_mut_factor*t12_bounds(2)/2;%100e-6;
%     t13_mutate = Max_mut_factor*t13_bounds(2)/2;%100e-6;
%     t14_mutate = Max_mut_factor*t14_bounds(2)/2;%100e-6;
    t21_mutate = Max_mut_factor*t21_bounds(2)/2;%500e-6;
    t13_mutate = Max_mut_factor*t13_bounds(2)/2;%10e-6;
    t31_mutate = Max_mut_factor*t31_bounds(2)/2;%10e-6;
%     t31_mutate = Max_mut_factor*t31_bounds(2)/2;%100e-6;
    t14_mutate = Max_mut_factor*t14_bounds(2)/2;%100e-6;
    t41_mutate = Max_mut_factor*t41_bounds(2)/2;%100e-6;
    
    A1_mutate = (A1_bounds(2) - A1_bounds(1))*0.1;%0.05;
    A2_mutate = (A2_bounds(2) - A2_bounds(1))*0.1;%0.05;0.05;
    A3_mutate = (A3_bounds(2) - A3_bounds(1))*0.1;%0.05;0.05;
    A4_mutate = (A4_bounds(2) - A4_bounds(1))*0.1;%0.05;0.05;
    maxMutationArray = [t12_mutate;t21_mutate;t13_mutate;t31_mutate;t14_mutate;t41_mutate;A1_mutate;A2_mutate;A3_mutate;A4_mutate];
    
    sigma_A1 = 0.032;
    sigma_A2 = 0.060;
    sigma_A3 = 0.027;
    sigma_A4 = 0.060;
    sigma_A = [sigma_A1; sigma_A2; sigma_A3; sigma_A4];
    
    
    param_strings = {'t12','t21','t13','t31','t14','t41','A1','A2','A3','A4'};
    % Note: the order of the strings needs to be the same as the variable
    % list so we can extract the parameters properly.
    
    % tij_strings = param_strings(1:(Nparam - Nstates));
    % A_strings = param_strings((Nparam - Nstates + 1):Nparam);
    
    %----------------------------------------------------------------------
    % Update model
    %----------------------------------------------------------------------
    if length(tijs) > (Nparam - 2*Nstates)
        tijs = tijs(1:(Nparam - 2*Nstates));
    end
    rates = 1 ./ tijs;
    
    [k12, k21, k13, k31, k14, k41] = feval(@(x) x{:}, num2cell(reshape(rates,[1, length(rates)])));
%     k32 = (k12 * k23 * k31)/(k13 * k21);
%     k41 = (k31 * k14 * k43)/ (k13 * k34);
%     k42 = (k24 * k43 * k31 * k12)/(k34 * k13 * k21);
    
    % [A1, A2, A3, A4] = feval(@(x) x{:}, num2cell(As));
    
    
     K = [-(k12 + k13 + k14), k21,               k31,               k41;...
        k12,              -(k21),               0,                   0;...
        k13,               0,                    -(k31),             0;...
        k14,               0,                     0,               -(k41)];
   

elseif strcmp(model_name, 'outer_loop') == 1 % this is the fully connected model
    Nstates = 4;
    t12_bounds = [10e-6,1000e-3];  %Paramater #1 is high--> med
    t14_bounds = [10e-6,1000e-3];
    t21_bounds = [10e-6,1000e-3];%Paramater #3 is med --> high
    t23_bounds = [10e-6,1000e-3];%Paramater #4 is med --> low
    t32_bounds = [10e-6,1000e-3];
    t34_bounds = [10e-6,1000e-3];
    t43_bounds = [10e-6,1000e-3];
    % *t32 wll be determined by the other rates
    
    % A1_bounds = [0.65,0.85];%Paramater #6 % HIGH fret State
    A1_bounds = [0.04,0.075];%low amp
    A2_bounds = [0.1,0.15];%2nd lowest
    A3_bounds = [0.15,0.24];%2nd high amp
    A4_bounds = [0.3,0.4]; %1st high amp 
    
    boundsArray = [t12_bounds; t14_bounds; t21_bounds; t23_bounds;  t32_bounds; t34_bounds; t43_bounds; A1_bounds; A2_bounds; A3_bounds; A4_bounds];
    Nparam = length(boundsArray);
    boundsArray_noAdj = boundsArray;
    
    [boundsArrayAdj,atEdges]= boundsAdjuster(boundsArray,windowMag,paramAdjs,absLow,absHigh,Nparam,Nstates);
    boundsArray=boundsArrayAdj;
    
    maxMutationCountsArray = zeros(1,Nparam);
    minMutationCountsArray = zeros(1,Nparam);
    
    % Maximum mutation value
    Max_mut_factor = .2;
    t12_mutate = Max_mut_factor*boundsArray(1,2)/2;%100e-6;
    t14_mutate = Max_mut_factor*boundsArray(2,2)/2;%100e-6;
    t21_mutate = Max_mut_factor*boundsArray(3,2)/2;%500e-6;
    t23_mutate = Max_mut_factor*boundsArray(4,2)/2;%10e-6;
    t32_mutate = Max_mut_factor*boundsArray(5,2)/2;%100e-6;
    t34_mutate = Max_mut_factor*boundsArray(6,2)/2;%100e-6;
    t43_mutate = Max_mut_factor*boundsArray(7,2)/2;%100e-6;
    
    A1_mutate = (A1_bounds(2) - A1_bounds(1))*0.1;%0.05;
    A2_mutate = (A2_bounds(2) - A2_bounds(1))*0.1;%0.05;0.05;
    A3_mutate = (A3_bounds(2) - A3_bounds(1))*0.1;%0.05;0.05;
    A4_mutate = (A4_bounds(2) - A4_bounds(1))*0.1;%0.05;0.05;
    maxMutationArray = [t12_mutate;t14_mutate;t21_mutate;t23_mutate;t32_mutate;t34_mutate;t43_mutate;A1_mutate;A2_mutate;A3_mutate;A4_mutate];
    
    sigma_A1 = 0.032;
    sigma_A2 = 0.060;
    sigma_A3 = 0.027;
    sigma_A4 = 0.060;
    sigma_A = [sigma_A1; sigma_A2; sigma_A3; sigma_A4];
    
    
    param_strings = {'t12','t14','t21','t23','t32','t34','t43','A1','A2','A3','A4'};
    % Note: the order of the strings needs to be the same as the variable
    % list so we can extract the parameters properly.
    
    % tij_strings = param_strings(1:(Nparam - Nstates));
    % A_strings = param_strings((Nparam - Nstates + 1):Nparam);
    
    %----------------------------------------------------------------------
    % Update model
    %----------------------------------------------------------------------
    if length(tijs) > (Nparam - 2*Nstates)
        tijs = tijs(1:(Nparam - 2*Nstates));
    end
    rates = 1 ./ tijs;
    
    [k12, k14, k21, k23, k32, k34, k43] = feval(@(x) x{:}, num2cell(reshape(rates,[1, length(rates)])));
    k41 = (k14 * k43 * k32 * k21)/ (k12 * k23 * k34);
    
    
    K = [-(k12+k14),    k21,                0,                  k41;...
        k12,              -(k23 +k21 ),       k32,                0;...
        0,               k23,              -( k32 + k34),       k43;...
        k14,               0,                 k34,                -(k41  + k43)];
    

elseif strcmp(model_name, 'fully_connected_minus41') == 1 % this is the fully connected model
    Nstates = 4;
    t12_bounds = [200e-6,50000e-3];  %Paramater #1 is high--> med
    t13_bounds = [200e-6,50000e-3];    %Paramater #2 is high --> low
    t21_bounds = [200e-6,50000e-3];%Paramater #3 is med --> high
    t23_bounds = [200e-6,50000e-3];%Paramater #4 is med --> low
    t24_bounds = [200e-6,50000e-3];
    t31_bounds = [200e-6,50000e-3];  %Paramater #5 is low --> Medium = [1e-3,10e-3];%Paramater #5 %t32 is low --> Medium
    t34_bounds = [200e-6,50000e-3];
    t43_bounds = [200e-6,50000e-3];
    % *t32 wll be determined by the other rates
    
     % A1_bounds = [0.65,0.85];%Paramater #6 % HIGH fret State
    A1_bounds = [0.001,0.055];%low amp
    A2_bounds = [0.08,0.13];%2nd lowest
    A3_bounds = [0.13,0.2];%2nd high amp
    A4_bounds = [0.25,0.4]; %1st high amp  
    
    %Bounds for standard deviation
     sigma_A1_bounds = [0.028,0.042]; %Paramater #8 %Low Amp state standard deviation
     sigma_A2_bounds = [0.05,0.07];%Paramater #9 %Med Amp state standard deviation
     sigma_A3_bounds = [0.065,0.09];%Paramater #10 %Med2 Amp state standard deviation
     sigma_A4_bounds = [0.09,0.14];%Paramater #10 %High Amp state standard deviation
     
    boundsArray = [t12_bounds; t13_bounds;  t21_bounds; t23_bounds; t24_bounds; t31_bounds; t34_bounds; t43_bounds; ...];
    A1_bounds; A2_bounds; A3_bounds; A4_bounds; ...
    sigma_A1_bounds; sigma_A2_bounds; sigma_A3_bounds ; sigma_A4_bounds];
    
    Nparam = length(boundsArray);
    boundsArray_noAdj = boundsArray;
    
    [boundsArrayAdj,atEdges]= boundsAdjuster_v3(boundsArray,windowMag,paramAdjs,absLow,absHigh,Nparam,Nstates,fixedDecades);
    boundsArray=boundsArrayAdj;
    
    maxMutationCountsArray = zeros(1,Nparam);
    minMutationCountsArray = zeros(1,Nparam);
    
    
    % Maximum mutation value
    Max_mut_factor = .2;
    
    maxMutationArray = [Max_mut_factor * boundsArray(1:Nparam-2*Nstates,2) / 2; ...
        0.1 * ( boundsArray(Nparam-2*Nstates+1:end,2) - boundsArray(Nparam-2*Nstates+1:end,1)) ]; 
    
    
    param_strings = {'t12','t13','t21','t23','t24','t31','t34','t43','A1','A2','A3','A4','sig1','sig2','sig3','sig4'};
    % Note: the order of the strings needs to be the same as the variable
    % list so we can extract the parameters properly.
    
    % tij_strings = param_strings(1:(Nparam - Nstates));
    % A_strings = param_strings((Nparam - Nstates + 1):Nparam);
    
    %----------------------------------------------------------------------
    % Update model
    %----------------------------------------------------------------------
    if length(tijs) > (Nparam - 2*Nstates)
        tijs = tijs(1:(Nparam - 2*Nstates));
    end
    rates = 1 ./ tijs;
    
    [k12, k13, k21, k23, k24, k31, k34, k43] = feval(@(x) x{:}, num2cell(reshape(rates,[1, length(rates)])));
    k32 = (k12 * k23 * k31)/(k13 * k21);
    k42 = (k24 * k43 * k31 * k12)/(k34 * k13 * k21);
    
    % [A1, A2, A3, A4] = feval(@(x) x{:}, num2cell(As));
    
    
    K = [-(k12 + k13),      k21,               k31,                0;...
        k12,                -(k21 + k23 + k24), k32,                k42;...
        k13,                k23,              -(k31 + k32 + k34),   k43;...
        0,                  k24,                k34,                -( k42 + k43)];

    
    elseif strcmp(model_name, 'outer_loop_plus13') == 1 % this is the fully connected model
    Nstates = 4;
    t12_bounds = [10e-6,1000e-3];  %Paramater #1 is high--> med
    t14_bounds = [10e-6,1000e-3];
    t13_bounds = [10e-6,1000e-3];
    t21_bounds = [10e-6,1000e-3];%Paramater #3 is med --> high
    t23_bounds = [10e-6,1000e-3];%Paramater #4 is med --> low
    t32_bounds = [10e-6,1000e-3];
    t34_bounds = [10e-6,1000e-3];
    t43_bounds = [10e-6,1000e-3];
    % *t32 wll be determined by the other rates
    
    % A1_bounds = [0.65,0.85];%Paramater #6 % HIGH fret State
    A1_bounds = [0.04,0.075];%low amp
    A2_bounds = [0.1,0.15];%2nd lowest
    A3_bounds = [0.15,0.24];%2nd high amp
    A4_bounds = [0.3,0.4]; %1st high amp 
    
    boundsArray = [t12_bounds; t14_bounds; t13_bounds; t21_bounds; t23_bounds;  t32_bounds; t34_bounds; t43_bounds; A1_bounds; A2_bounds; A3_bounds; A4_bounds];
    Nparam = length(boundsArray);
    boundsArray_noAdj = boundsArray;
    
    [boundsArrayAdj,atEdges]= boundsAdjuster(boundsArray,windowMag,paramAdjs,absLow,absHigh,Nparam,Nstates);
    boundsArray=boundsArrayAdj;
    
    maxMutationCountsArray = zeros(1,Nparam);
    minMutationCountsArray = zeros(1,Nparam);
    
    % Maximum mutation value
    Max_mut_factor = .2;
    t12_mutate = Max_mut_factor*boundsArray(1,2)/2;%100e-6;
    t14_mutate = Max_mut_factor*boundsArray(2,2)/2;%100e-6;
    t13_mutate = Max_mut_factor*boundsArray(3,2)/2;
    t21_mutate = Max_mut_factor*boundsArray(4,2)/2;%500e-6;
    t23_mutate = Max_mut_factor*boundsArray(5,2)/2;%10e-6;
    t32_mutate = Max_mut_factor*boundsArray(7,2)/2;%100e-6;
    t34_mutate = Max_mut_factor*boundsArray(8,2)/2;%100e-6;
    t43_mutate = Max_mut_factor*boundsArray(9,2)/2;%100e-6;
    
    A1_mutate = (A1_bounds(2) - A1_bounds(1))*0.1;%0.05;
    A2_mutate = (A2_bounds(2) - A2_bounds(1))*0.1;%0.05;0.05;
    A3_mutate = (A3_bounds(2) - A3_bounds(1))*0.1;%0.05;0.05;
    A4_mutate = (A4_bounds(2) - A4_bounds(1))*0.1;%0.05;0.05;
    maxMutationArray = [t12_mutate;t14_mutate;t13_mutate;t21_mutate;t23_mutate;t32_mutate;t34_mutate;t43_mutate;A1_mutate;A2_mutate;A3_mutate;A4_mutate];
    
    sigma_A1 = 0.032;
    sigma_A2 = 0.060;
    sigma_A3 = 0.027;
    sigma_A4 = 0.060;
    sigma_A = [sigma_A1; sigma_A2; sigma_A3; sigma_A4];
    
    
    param_strings = {'t12','t14','t13','t21','t23','t32','t34','t43','A1','A2','A3','A4'};
    % Note: the order of the strings needs to be the same as the variable
    % list so we can extract the parameters properly.
    
    % tij_strings = param_strings(1:(Nparam - Nstates));
    % A_strings = param_strings((Nparam - Nstates + 1):Nparam);
    
    %----------------------------------------------------------------------
    % Update model
    %----------------------------------------------------------------------
    if length(tijs) > (Nparam - 2*Nstates)
        tijs = tijs(1:(Nparam - 2*Nstates));
    end
    rates = 1 ./ tijs;
    
    [k12, k14, k21, k23, k32, k34, k43] = feval(@(x) x{:}, num2cell(reshape(rates,[1, length(rates)])));
    k41 = (k14 * k43 * k32 * k21)/ (k12 * k23 * k34);
    k31 = (k13 * k21 * k32) / (k12 * k23);
    
    K = [-(k12+k14+k13),    k21,               k31,                    k41;...
        k12,              -(k23 +k21 ),       k32,                      0;...
        k13,               k23,              -(k31+ k32 + k34),         k43;...
        k14,               0,                 k34,                      -(k41  + k43)];
    
elseif strcmp(model_name, 'fully_connected_minus43') == 1 % this is the fully connected model
    Nstates = 4;
    t12_bounds = [200e-6,50000e-3];  %Paramater #1 is high--> med
    t21_bounds = [200e-6,50000e-3];    %Paramater #2 is high --> low
    t13_bounds = [200e-6,50000e-3];
    t31_bounds = [200e-6,50000e-3];%Paramater #3 is med --> high
    t14_bounds = [200e-6,50000e-3];%Paramater #4 is med --> low
%     t41_bounds = [100e-6,100e-3];
    t23_bounds = [200e-6,50000e-3];  %Paramater #5 is low --> Medium = [1e-3,10e-3];%Paramater #5 %t32 is low --> Medium
%     t32_bounds = [100e-6,1000e-3];
    t24_bounds = [200e-6,50000e-3];
    t42_bounds = [200e-6,50000e-3];
%     t34_bounds = [100e-6,1000e-3];
%     t43_bounds = [100e-6,1000e-3];
    

    % *t32 wll be determined by the other rates
    
    % A1_bounds = [0.65,0.85];%Paramater #6 % HIGH fret State
    A1_bounds = [0.001,0.055];%low amp
    A2_bounds = [0.08,0.13];%2nd lowest
    A3_bounds = [0.13,0.2];%2nd high amp
    A4_bounds = [0.25,0.4]; %1st high amp  
    
    %Bounds for standard deviation
     sigma_A1_bounds = [0.028,0.042]; %Paramater #8 %Low Amp state standard deviation
     sigma_A2_bounds = [0.05,0.07];%Paramater #9 %Med Amp state standard deviation
     sigma_A3_bounds = [0.065,0.09];%Paramater #10 %Med2 Amp state standard deviation
     sigma_A4_bounds = [0.09,0.14];%Paramater #10 %High Amp state standard deviation
    
    boundsArray = [t12_bounds; t21_bounds; t13_bounds; t31_bounds; t14_bounds; t23_bounds; ...
        t24_bounds; t42_bounds; ...
        A1_bounds; A2_bounds; A3_bounds; A4_bounds; ...
        sigma_A1_bounds; sigma_A2_bounds; sigma_A3_bounds ; sigma_A4_bounds];
    
    Nparam = length(boundsArray);
    boundsArray_noAdj = boundsArray;
    
    [boundsArrayAdj,atEdges]= boundsAdjuster_v3(boundsArray,windowMag,paramAdjs,absLow,absHigh,Nparam,Nstates,fixedDecades);
    boundsArray=boundsArrayAdj;
    
    maxMutationCountsArray = zeros(1,Nparam);
    minMutationCountsArray = zeros(1,Nparam);
    
    
    % Maximum mutation value
    Max_mut_factor = .2;
    
    maxMutationArray = [Max_mut_factor * boundsArray(1:Nparam-2*Nstates,2) / 2; ...
        0.1 * ( boundsArray(Nparam-2*Nstates+1:end,2) - boundsArray(Nparam-2*Nstates+1:end,1)) ]; 
    
    param_strings = {'t12','t21','t13','t31','t14','t23','t24','t42','A1','A2','A3','A4','sig1','sig2','sig3','sig4'};
    % Note: the order of the strings needs to be the same as the variable
    % list so we can extract the parameters properly.
    
    % tij_strings = param_strings(1:(Nparam - Nstates));
    % A_strings = param_strings((Nparam - Nstates + 1):Nparam);
    
    %----------------------------------------------------------------------
    % Update model
    %----------------------------------------------------------------------
    if length(tijs) > (Nparam - 2*Nstates)
        tijs = tijs(1:(Nparam - 2*Nstates));
    end
    rates = 1 ./ tijs;
    
    [k12, k21, k13, k31, k14, k23, k24, k42] = feval(@(x) x{:}, num2cell(reshape(rates,[1, length(rates)])));
    k32 = (k12 * k23 * k31)/(k13 * k21);
    k41 = (k21 * k42 * k14)/ (k12 * k24);
%     k42 = (k24 * k43 * k31 * k12)/(k34 * k13 * k21);
    
    % [A1, A2, A3, A4] = feval(@(x) x{:}, num2cell(As));
    
    
    K = [-(k12 + k13 + k14), k21,               k31,             k41;...
        k12,              -(k21 + k23 + k24), k32,               k42;...
        k13,               k23,              -(k31 + k32 ),      0;...
        k14,               k24,                0,              -(k41 + k42 )];

    elseif strcmp(model_name, 'fully_connected_minus23') == 1 % this is the fully connected model
    Nstates = 4;
    t12_bounds = [200e-6,50000e-3];  %Paramater #1 is high--> med
    t21_bounds = [200e-6,50000e-3];    %Paramater #2 is high --> low
%     t13_bounds = [5e-6,700e-3];
    t31_bounds = [200e-6,50000e-3];%Paramater #3 is med --> high
    t14_bounds = [200e-6,50000e-3];%Paramater #4 is med --> low
    t41_bounds = [200e-6,50000e-3];
%     t23_bounds = [80e-6,500e-3];  %Paramater #5 is low --> Medium = [1e-3,10e-3];%Paramater #5 %t32 is low --> Medium
%     t32_bounds = [100e-6,1000e-3];
    t24_bounds = [200e-6,50000e-3];
%     t42_bounds = [100e-6,1000e-3];
    t34_bounds = [200e-6,50000e-3];
    t43_bounds = [200e-6,50000e-3];
    

    % *t32 wll be determined by the other rates
    
    % A1_bounds = [0.65,0.85];%Paramater #6 % HIGH fret State
    A1_bounds = [0.001,0.055];%low amp
    A2_bounds = [0.08,0.13];%2nd lowest
    A3_bounds = [0.13,0.2];%2nd high amp
    A4_bounds = [0.25,0.4]; %1st high amp  
    
    %Bounds for standard deviation
     sigma_A1_bounds = [0.028,0.042]; %Paramater #8 %Low Amp state standard deviation
     sigma_A2_bounds = [0.05,0.07];%Paramater #9 %Med Amp state standard deviation
     sigma_A3_bounds = [0.065,0.09];%Paramater #10 %Med2 Amp state standard deviation
     sigma_A4_bounds = [0.09,0.14];%Paramater #10 %High Amp state standard deviation
    
    boundsArray = [t12_bounds; t21_bounds; t31_bounds; t14_bounds; t41_bounds; ...
        t24_bounds; t34_bounds; t43_bounds; ...
        A1_bounds; A2_bounds; A3_bounds; A4_bounds; ...
        sigma_A1_bounds; sigma_A2_bounds; sigma_A3_bounds ; sigma_A4_bounds];
    
    Nparam = length(boundsArray);
    boundsArray_noAdj = boundsArray;
    
    [boundsArrayAdj,atEdges]= boundsAdjuster_v3(boundsArray,windowMag,paramAdjs,absLow,absHigh,Nparam,Nstates,fixedDecades);
    boundsArray=boundsArrayAdj;
    
    maxMutationCountsArray = zeros(1,Nparam);
    minMutationCountsArray = zeros(1,Nparam);
    
    
    % Maximum mutation value
    Max_mut_factor = .2;
    
    maxMutationArray = [Max_mut_factor * boundsArray(1:Nparam-2*Nstates,2) / 2; ...
        0.1 * ( boundsArray(Nparam-2*Nstates+1:end,2) - boundsArray(Nparam-2*Nstates+1:end,1)) ]; 
    
    param_strings = {'t12','t21','t31','t14','t41','t24','t34','t43','A1','A2','A3','A4','sig1','sig2','sig3','sig4'};
    % Note: the order of the strings needs to be the same as the variable
    % list so we can extract the parameters properly.
    
    % tij_strings = param_strings(1:(Nparam - Nstates));
    % A_strings = param_strings((Nparam - Nstates + 1):Nparam);
    
    %----------------------------------------------------------------------
    % Update model
    %----------------------------------------------------------------------
    if length(tijs) > (Nparam - 2*Nstates)
        tijs = tijs(1:(Nparam - 2*Nstates));
    end
    rates = 1 ./ tijs;
    
    [k12, k21, k31, k14, k41, k24, k34, k43] = feval(@(x) x{:}, num2cell(reshape(rates,[1, length(rates)])));
    k42 = (k12 * k24 * k41)/(k14 * k21);
    k13 = (k31 * k14 * k43)/ (k41 * k34);
%     k42 = (k24 * k43 * k31 * k12)/(k34 * k13 * k21);
    
    % [A1, A2, A3, A4] = feval(@(x) x{:}, num2cell(As));
    
    
    K = [-(k12 + k13 + k14), k21,               k31,               k41;...
        k12,              -(k21 + k24),          0,               k42;...
        k13,               0,                -(k31 + k34),        k43;...
        k14,               k24,                k34,             -(k41 + k42 + k43)];
    
     elseif strcmp(model_name, 'fully_connected_minus24') == 1 % this is the fully connected model
    Nstates = 4;
    t12_bounds = [200e-6,50000e-3];  %Paramater #1 is high--> med
    t21_bounds = [200e-6,50000e-3];    %Paramater #2 is high --> low
    t13_bounds = [200e-6,50000e-3];
    t31_bounds = [200e-6,50000e-3];%Paramater #3 is med --> high
    t14_bounds = [200e-6,50000e-3];%Paramater #4 is med --> low
%     t41_bounds = [100e-6,100e-3];
    t23_bounds = [200e-6,50000e-3];  %Paramater #5 is low --> Medium = [1e-3,10e-3];%Paramater #5 %t32 is low --> Medium
%     t32_bounds = [100e-6,1000e-3];
%     t24_bounds = [100e-6,1000e-3];
%     t42_bounds = [100e-6,1000e-3];
    t34_bounds = [200e-6,50000e-3];
    t43_bounds = [200e-6,50000e-3];
    

    % *t32 wll be determined by the other rates
    
    % A1_bounds = [0.65,0.85];%Paramater #6 % HIGH fret State
    A1_bounds = [0.001,0.055];%low amp
    A2_bounds = [0.08,0.13];%2nd lowest
    A3_bounds = [0.13,0.2];%2nd high amp
    A4_bounds = [0.25,0.4]; %1st high amp  
    
    %Bounds for standard deviation
     sigma_A1_bounds = [0.028,0.042]; %Paramater #8 %Low Amp state standard deviation
     sigma_A2_bounds = [0.05,0.07];%Paramater #9 %Med Amp state standard deviation
     sigma_A3_bounds = [0.065,0.09];%Paramater #10 %Med2 Amp state standard deviation
     sigma_A4_bounds = [0.09,0.14];%Paramater #10 %High Amp state standard deviation
    
    boundsArray = [t12_bounds; t21_bounds; t13_bounds; t31_bounds; t14_bounds; t23_bounds; ...
         t34_bounds; t43_bounds; ...
        A1_bounds; A2_bounds; A3_bounds; A4_bounds; ...
        sigma_A1_bounds; sigma_A2_bounds; sigma_A3_bounds ; sigma_A4_bounds];
    
    Nparam = length(boundsArray);
    boundsArray_noAdj = boundsArray;
    
    [boundsArrayAdj,atEdges]= boundsAdjuster_v3(boundsArray,windowMag,paramAdjs,absLow,absHigh,Nparam,Nstates,fixedDecades);
    boundsArray=boundsArrayAdj;
    
    maxMutationCountsArray = zeros(1,Nparam);
    minMutationCountsArray = zeros(1,Nparam);
    
    
    % Maximum mutation value
    Max_mut_factor = .2;
    
    maxMutationArray = [Max_mut_factor * boundsArray(1:Nparam-2*Nstates,2) / 2; ...
        0.1 * ( boundsArray(Nparam-2*Nstates+1:end,2) - boundsArray(Nparam-2*Nstates+1:end,1)) ]; 
    
    param_strings = {'t12','t21','t13','t31','t14','t23','t34','t43','A1','A2','A3','A4','sig1','sig2','sig3','sig4'};
    % Note: the order of the strings needs to be the same as the variable
    % list so we can extract the parameters properly.
    
    % tij_strings = param_strings(1:(Nparam - Nstates));
    % A_strings = param_strings((Nparam - Nstates + 1):Nparam);
    
    %----------------------------------------------------------------------
    % Update model
    %----------------------------------------------------------------------
    if length(tijs) > (Nparam - 2*Nstates)
        tijs = tijs(1:(Nparam - 2*Nstates));
    end
    rates = 1 ./ tijs;
    
    [k12, k21, k13, k31, k14, k23, k34, k43] = feval(@(x) x{:}, num2cell(reshape(rates,[1, length(rates)])));
    k32 = (k12 * k23 * k31)/(k13 * k21);
    k41 = (k31 * k14 * k43)/ (k13 * k34);
    
    % [A1, A2, A3, A4] = feval(@(x) x{:}, num2cell(As));
    
    
    K = [-(k12 + k13 + k14), k21,               k31,               k41;...
        k12,              -(k21 + k23),         k32,               0;...
        k13,               k23,              -(k31 + k32 + k34),  k43;...
        k14,               0,                   k34,              -(k41 + k43)];
    
    
 elseif strcmp(model_name, 'fully_connected_TEMPLATE') == 1 % this is the fully connected model
    Nstates = 4;
    t12_bounds = [50e-6,1000e-3];  %Paramater #1 is high--> med
    t21_bounds = [10e-3,1000e-3];    %Paramater #2 is high --> low
    t13_bounds = [5e-6,700e-3];
    t31_bounds = [10e-6,100e-3];%Paramater #3 is med --> high
    t14_bounds = [10e-6,100e-3];%Paramater #4 is med --> low
    t41_bounds = [100e-6,100e-3];
    t23_bounds = [80e-6,500e-3];  %Paramater #5 is low --> Medium = [1e-3,10e-3];%Paramater #5 %t32 is low --> Medium
    t32_bounds = [100e-6,1000e-3];
    t24_bounds = [100e-6,1000e-3];
    t42_bounds = [100e-6,1000e-3];
    t34_bounds = [100e-6,1000e-3];
    t43_bounds = [100e-6,1000e-3];
    

    % *t32 wll be determined by the other rates
    
    % A1_bounds = [0.65,0.85];%Paramater #6 % HIGH fret State
    A1_bounds = [0.04,0.055];%low amp
    A2_bounds = [0.08,0.13];%2nd lowest
    A3_bounds = [0.15,0.2];%2nd high amp
    A4_bounds = [0.25,0.32]; %1st high amp  
    
    %Bounds for standard deviation
     sigma_A1_bounds = [0.028,0.042]; %Paramater #8 %Low Amp state standard deviation
     sigma_A2_bounds = [0.05,0.07];%Paramater #9 %Med Amp state standard deviation
     sigma_A3_bounds = [0.065,0.09];%Paramater #10 %Med2 Amp state standard deviation
     sigma_A4_bounds = [0.09,0.14];%Paramater #10 %High Amp state standard deviation
    
    boundsArray = [t12_bounds; t21_bounds; t13_bounds; t31_bounds; t14_bounds; t41_bounds; t23_bounds; ...
         t32_bounds; t24_bounds; t42_bounds; t34_bounds; t43_bounds; ...
        A1_bounds; A2_bounds; A3_bounds; A4_bounds; ...
        sigma_A1_bounds; sigma_A2_bounds; sigma_A3_bounds ; sigma_A4_bounds];
    
    Nparam = length(boundsArray);
    boundsArray_noAdj = boundsArray;
    
    [boundsArrayAdj,atEdges]= boundsAdjuster_v3(boundsArray,windowMag,paramAdjs,absLow,absHigh,Nparam,Nstates,fixedDecades);
    boundsArray=boundsArrayAdj;
    
    maxMutationCountsArray = zeros(1,Nparam);
    minMutationCountsArray = zeros(1,Nparam);
    
    
    % Maximum mutation value
    Max_mut_factor = .2;
    
    maxMutationArray = [Max_mut_factor * boundsArray(1:Nparam-2*Nstates,2) / 2; ...
        0.1 * ( boundsArray(Nparam-2*Nstates+1:end,2) - boundsArray(Nparam-2*Nstates+1:end,1)) ]; 
    
    param_strings = {'t12','t21','t13','t31','t14','t41','t23','t32','t24','t42','t34','t43','A1','A2','A3','A4','sig1','sig2','sig3','sig4'};
    % Note: the order of the strings needs to be the same as the variable
    % list so we can extract the parameters properly.
    
    % tij_strings = param_strings(1:(Nparam - Nstates));
    % A_strings = param_strings((Nparam - Nstates + 1):Nparam);
    
    %----------------------------------------------------------------------
    % Update model
    %----------------------------------------------------------------------
    if length(tijs) > (Nparam - 2*Nstates)
        tijs = tijs(1:(Nparam - 2*Nstates));
    end
    rates = 1 ./ tijs;
    
    [k12, k21, k13, k31, k14, k41, k23, k32, k24, k42, k34, k43] = feval(@(x) x{:}, num2cell(reshape(rates,[1, length(rates)])));
%     k32 = (k12 * k23 * k31)/(k13 * k21);
%     k41 = (k31 * k14 * k43)/ (k13 * k34);
%     k42 = (k24 * k43 * k31 * k12)/(k34 * k13 * k21);
    
    % [A1, A2, A3, A4] = feval(@(x) x{:}, num2cell(As));
    
    
    K = [-(k12 + k13 + k14), k21,               k31,               k41;...
        k12,              -(k21 + k23 + k24), k32,               k42;...
        k13,               k23,              -(k31 + k32 + k34), k43;...
        k14,               k24,                k34,             -(k41 + k42 + k43)];



%     ************* ANABEL LINEAR CHAINS START HERE *****************
  elseif strcmp(model_name, 'linear1324') == 1
     Nstates = 4;
    t13_bounds = [200e-6,50000e-3];  %Paramater #1 is high--> med    
%     t13_bounds = [1e-3,3000e-3];    %Paramater #2 is high --> low
%     t14_bounds = [100e-6,900e-3];
    t31_bounds = [200e-6,50000e-3];%Paramater #3 is med --> high
    t23_bounds = [200e-6,50000e-3];%Paramater #4 is med --> low
    t32_bounds = [200e-6,50000e-3];
%     t31_bounds = [50e-6,1000e-3];  %Paramater #5 is low --> Medium = [1e-3,10e-3];%Paramater #5 %t32 is low --> Medium
    t24_bounds = [200e-6,50000e-3];
    t42_bounds = [200e-6,50000e-3];
    % *t32 wll be determined by the other rates
    
% IN THIS SCHEME, 2 AND 3 ARE THE LOW AMP MOSTLY SYMM/ANTI WHILE 1 AND 4
% ARE THE HIGH AP STATES
    % A1_bounds = [0.65,0.85];%Paramater #6 % HIGH fret State
    A1_bounds = [0.04,0.055];%low amp
    A2_bounds = [0.08,0.13];%2nd lowest
    A3_bounds = [0.15,0.2];%2nd high amp
    A4_bounds = [0.25,0.32]; %1st high amp  
    
    %Bounds for standard deviation
     sigma_A1_bounds = [0.028,0.042]; %Paramater #8 %Low Amp state standard deviation
     sigma_A2_bounds = [0.05,0.07];%Paramater #9 %Med Amp state standard deviation
     sigma_A3_bounds = [0.065,0.09];%Paramater #10 %Med2 Amp state standard deviation
     sigma_A4_bounds = [0.09,0.14];%Paramater #10 %High Amp state standard deviation
    
    boundsArray = [t13_bounds; t31_bounds; t23_bounds; t32_bounds; t24_bounds; t42_bounds; ...
        A1_bounds; A2_bounds; A3_bounds; A4_bounds; ...
        sigma_A1_bounds; sigma_A2_bounds; sigma_A3_bounds ; sigma_A4_bounds];
    
    Nparam = length(boundsArray);
    boundsArray_noAdj = boundsArray;
    
    [boundsArrayAdj,atEdges]= boundsAdjuster_v3(boundsArray,windowMag,paramAdjs,absLow,absHigh,Nparam,Nstates,fixedDecades);
    boundsArray=boundsArrayAdj;
    
    maxMutationCountsArray = zeros(1,Nparam);
    minMutationCountsArray = zeros(1,Nparam);
    
    
    % Maximum mutation value
    Max_mut_factor = .2;
    
    maxMutationArray = [Max_mut_factor * boundsArray(1:Nparam-2*Nstates,2) / 2; ...
        0.1 * ( boundsArray(Nparam-2*Nstates+1:end,2) - boundsArray(Nparam-2*Nstates+1:end,1)) ];    
    
    param_strings = {'t13','t31','t23','t32','t24','t42','A1','A2','A3','A4','sig1','sig2','sig3','sig4'};
    % Note: the order of the strings needs to be the same as the variable
    % list so we can extract the parameters properly.
    
    % tij_strings = param_strings(1:(Nparam - Nstates));
    % A_strings = param_strings((Nparam - Nstates + 1):Nparam);
    
    %----------------------------------------------------------------------
    % Update model
    %----------------------------------------------------------------------
    if length(tijs) > (Nparam - 2*Nstates)
        tijs = tijs(1:(Nparam - 2*Nstates));
    end
    rates = 1 ./ tijs;
    
    [k13, k31, k23, k32, k24, k42] = feval(@(x) x{:}, num2cell(reshape(rates,[1, length(rates)])));
%     k32 = (k12 * k23 * k31)/(k13 * k21);
%     k41 = (k31 * k14 * k43)/ (k13 * k34);
%     k42 = (k24 * k43 * k31 * k12)/(k34 * k13 * k21);
    
    % [A1, A2, A3, A4] = feval(@(x) x{:}, num2cell(As));
    
    K = [-(k13),          0,                k31,               0;...
        0,              -(k23 + k24),       k32,               k42;...
        k13,             k23,              -(k31 + k32 ),       0;...
        0,               k24,                0,             -(k42)];
    
    elseif strcmp(model_name, 'linear1342') == 1
     Nstates = 4;
    t13_bounds = [200e-6,50000e-3];  %Paramater #1 is high--> med    
%     t13_bounds = [1e-3,3000e-3];    %Paramater #2 is high --> low
%     t14_bounds = [100e-6,900e-3];
    t31_bounds = [200e-6,50000e-3];%Paramater #3 is med --> high
    t34_bounds = [200e-6,50000e-3];%Paramater #4 is med --> low
    t43_bounds = [200e-6,50000e-3];
%     t31_bounds = [50e-6,1000e-3];  %Paramater #5 is low --> Medium = [1e-3,10e-3];%Paramater #5 %t32 is low --> Medium
    t42_bounds = [200e-6,50000e-3];
    t24_bounds = [200e-6,50000e-3];
    % *t32 wll be determined by the other rates
    
% IN THIS SCHEME, 2 AND 3 ARE THE LOW AMP MOSTLY SYMM/ANTI WHILE 1 AND 4
% ARE THE HIGH AP STATES
    % A1_bounds = [0.65,0.85];%Paramater #6 % HIGH fret State
    A1_bounds = [0.04,0.055];%low amp
    A2_bounds = [0.08,0.13];%2nd lowest
    A3_bounds = [0.15,0.2];%2nd high amp
    A4_bounds = [0.25,0.32]; %1st high amp  
    
    %Bounds for standard deviation
     sigma_A1_bounds = [0.028,0.042]; %Paramater #8 %Low Amp state standard deviation
     sigma_A2_bounds = [0.05,0.07];%Paramater #9 %Med Amp state standard deviation
     sigma_A3_bounds = [0.065,0.09];%Paramater #10 %Med2 Amp state standard deviation
     sigma_A4_bounds = [0.09,0.14];%Paramater #10 %High Amp state standard deviation
    
    boundsArray = [t13_bounds; t31_bounds; t34_bounds; t43_bounds; t42_bounds; t24_bounds; ...
        A1_bounds; A2_bounds; A3_bounds; A4_bounds; ...
        sigma_A1_bounds; sigma_A2_bounds; sigma_A3_bounds ; sigma_A4_bounds];
    
    Nparam = length(boundsArray);
    boundsArray_noAdj = boundsArray;
    
    [boundsArrayAdj,atEdges]= boundsAdjuster_v3(boundsArray,windowMag,paramAdjs,absLow,absHigh,Nparam,Nstates,fixedDecades);
    boundsArray=boundsArrayAdj;
    
    maxMutationCountsArray = zeros(1,Nparam);
    minMutationCountsArray = zeros(1,Nparam);
    
    
    % Maximum mutation value
    Max_mut_factor = .2;
    
    maxMutationArray = [Max_mut_factor * boundsArray(1:Nparam-2*Nstates,2) / 2; ...
        0.1 * ( boundsArray(Nparam-2*Nstates+1:end,2) - boundsArray(Nparam-2*Nstates+1:end,1)) ];    
    
    param_strings = {'t13','t31','t34','t43','t42','t24','A1','A2','A3','A4','sig1','sig2','sig3','sig4'};
    % Note: the order of the strings needs to be the same as the variable
    % list so we can extract the parameters properly.
    
    % tij_strings = param_strings(1:(Nparam - Nstates));
    % A_strings = param_strings((Nparam - Nstates + 1):Nparam);
    
    %----------------------------------------------------------------------
    % Update model
    %----------------------------------------------------------------------
    if length(tijs) > (Nparam - 2*Nstates)
        tijs = tijs(1:(Nparam - 2*Nstates));
    end
    rates = 1 ./ tijs;
    
    [k13, k31, k34, k43, k42, k24] = feval(@(x) x{:}, num2cell(reshape(rates,[1, length(rates)])));
%     k32 = (k12 * k23 * k31)/(k13 * k21);
%     k41 = (k31 * k14 * k43)/ (k13 * k34);
%     k42 = (k24 * k43 * k31 * k12)/(k34 * k13 * k21);
    
    % [A1, A2, A3, A4] = feval(@(x) x{:}, num2cell(As));
    
    K = [-(k13),          0,                k31,               0;...
        0,              -(k24),              0,               k42;...
        k13,             0,              -(k31 + k34),         k43;...
        0,               k24,                k34,           -(k42 + k43)];
    
    
     elseif strcmp(model_name, 'linear1432') == 1
     Nstates = 4;
    t14_bounds = [200e-6,50000e-3];  %Paramater #1 is high--> med    
%     t13_bounds = [1e-3,3000e-3];    %Paramater #2 is high --> low
%     t14_bounds = [100e-6,900e-3];
    t41_bounds = [200e-6,50000e-3];%Paramater #3 is med --> high
    t43_bounds = [200e-6,50000e-3];%Paramater #4 is med --> low
    t34_bounds = [200e-6,50000e-3];
%     t31_bounds = [50e-6,1000e-3];  %Paramater #5 is low --> Medium = [1e-3,10e-3];%Paramater #5 %t32 is low --> Medium
    t32_bounds = [200e-6,50000e-3];
    t23_bounds = [200e-6,50000e-3];
    % *t32 wll be determined by the other rates
    
% IN THIS SCHEME, 2 AND 3 ARE THE LOW AMP MOSTLY SYMM/ANTI WHILE 1 AND 4
% ARE THE HIGH AP STATES
    % A1_bounds = [0.65,0.85];%Paramater #6 % HIGH fret State
    A1_bounds = [0.04,0.055];%low amp
    A2_bounds = [0.08,0.13];%2nd lowest
    A3_bounds = [0.15,0.2];%2nd high amp
    A4_bounds = [0.25,0.32]; %1st high amp  
    
    %Bounds for standard deviation
     sigma_A1_bounds = [0.028,0.042]; %Paramater #8 %Low Amp state standard deviation
     sigma_A2_bounds = [0.05,0.07];%Paramater #9 %Med Amp state standard deviation
     sigma_A3_bounds = [0.065,0.09];%Paramater #10 %Med2 Amp state standard deviation
     sigma_A4_bounds = [0.09,0.14];%Paramater #10 %High Amp state standard deviation
    
    boundsArray = [t14_bounds; t41_bounds; t43_bounds; t34_bounds; t32_bounds; t23_bounds; ...
        A1_bounds; A2_bounds; A3_bounds; A4_bounds; ...
        sigma_A1_bounds; sigma_A2_bounds; sigma_A3_bounds ; sigma_A4_bounds];
    
    Nparam = length(boundsArray);
    boundsArray_noAdj = boundsArray;
    
    [boundsArrayAdj,atEdges]= boundsAdjuster_v3(boundsArray,windowMag,paramAdjs,absLow,absHigh,Nparam,Nstates,fixedDecades);
    boundsArray=boundsArrayAdj;
    
    maxMutationCountsArray = zeros(1,Nparam);
    minMutationCountsArray = zeros(1,Nparam);
    
    
    % Maximum mutation value
    Max_mut_factor = .2;
    
    maxMutationArray = [Max_mut_factor * boundsArray(1:Nparam-2*Nstates,2) / 2; ...
        0.1 * ( boundsArray(Nparam-2*Nstates+1:end,2) - boundsArray(Nparam-2*Nstates+1:end,1)) ];    
    
    param_strings = {'t14','t41','t43','t34','t32','t23','A1','A2','A3','A4','sig1','sig2','sig3','sig4'};
    % Note: the order of the strings needs to be the same as the variable
    % list so we can extract the parameters properly.
    
    % tij_strings = param_strings(1:(Nparam - Nstates));
    % A_strings = param_strings((Nparam - Nstates + 1):Nparam);
    
    %----------------------------------------------------------------------
    % Update model
    %----------------------------------------------------------------------
    if length(tijs) > (Nparam - 2*Nstates)
        tijs = tijs(1:(Nparam - 2*Nstates));
    end
    rates = 1 ./ tijs;
    
    [k14, k41, k43, k34, k32, k23] = feval(@(x) x{:}, num2cell(reshape(rates,[1, length(rates)])));
%     k32 = (k12 * k23 * k31)/(k13 * k21);
%     k41 = (k31 * k14 * k43)/ (k13 * k34);
%     k42 = (k24 * k43 * k31 * k12)/(k34 * k13 * k21);
    
    % [A1, A2, A3, A4] = feval(@(x) x{:}, num2cell(As));
    
    K = [-(k14),          0,                 0,                k41;...
        0,              -(k23),              k32,               0;...
        0,                k23,              -(k32 + k34),       k43;...
        k14,               0,                k34,           -(k41 + k43)];
    
  elseif strcmp(model_name, 'linear1234') == 1
     Nstates = 4;
    t12_bounds = [200e-6,50000e-3];  %Paramater #1 is high--> med    
%     t13_bounds = [1e-3,3000e-3];    %Paramater #2 is high --> low
%     t14_bounds = [100e-6,900e-3];
    t21_bounds = [200e-6,50000e-3];%Paramater #3 is med --> high
    t23_bounds = [200e-6,50000e-3];%Paramater #4 is med --> low
    t32_bounds = [200e-6,50000e-3];
%     t31_bounds = [50e-6,1000e-3];  %Paramater #5 is low --> Medium = [1e-3,10e-3];%Paramater #5 %t32 is low --> Medium
    t34_bounds = [200e-6,50000e-3];
    t43_bounds = [200e-6,50000e-3];
    % *t32 wll be determined by the other rates
    
% IN THIS SCHEME, 2 AND 3 ARE THE LOW AMP MOSTLY SYMM/ANTI WHILE 1 AND 4
% ARE THE HIGH AP STATES
    % A1_bounds = [0.65,0.85];%Paramater #6 % HIGH fret State
    A1_bounds = [0.04,0.055];%low amp
    A2_bounds = [0.08,0.13];%2nd lowest
    A3_bounds = [0.15,0.2];%2nd high amp
    A4_bounds = [0.25,0.32]; %1st high amp  
    
    %Bounds for standard deviation
     sigma_A1_bounds = [0.028,0.042]; %Paramater #8 %Low Amp state standard deviation
     sigma_A2_bounds = [0.05,0.07];%Paramater #9 %Med Amp state standard deviation
     sigma_A3_bounds = [0.065,0.09];%Paramater #10 %Med2 Amp state standard deviation
     sigma_A4_bounds = [0.09,0.14];%Paramater #10 %High Amp state standard deviation
    
    boundsArray = [t12_bounds; t21_bounds; t23_bounds; t32_bounds; t34_bounds; t43_bounds; ...
        A1_bounds; A2_bounds; A3_bounds; A4_bounds; ...
        sigma_A1_bounds; sigma_A2_bounds; sigma_A3_bounds ; sigma_A4_bounds];
    
    Nparam = length(boundsArray);
    boundsArray_noAdj = boundsArray;
    
    [boundsArrayAdj,atEdges]= boundsAdjuster_v3(boundsArray,windowMag,paramAdjs,absLow,absHigh,Nparam,Nstates,fixedDecades);
    boundsArray=boundsArrayAdj;
    
    maxMutationCountsArray = zeros(1,Nparam);
    minMutationCountsArray = zeros(1,Nparam);
    
    
    % Maximum mutation value
    Max_mut_factor = .2;
    
    maxMutationArray = [Max_mut_factor * boundsArray(1:Nparam-2*Nstates,2) / 2; ...
        0.1 * ( boundsArray(Nparam-2*Nstates+1:end,2) - boundsArray(Nparam-2*Nstates+1:end,1)) ];    
    
    param_strings = {'t12','t21','t23','t32','t34','t43','A1','A2','A3','A4','sig1','sig2','sig3','sig4'};
    % Note: the order of the strings needs to be the same as the variable
    % list so we can extract the parameters properly.
    
    % tij_strings = param_strings(1:(Nparam - Nstates));
    % A_strings = param_strings((Nparam - Nstates + 1):Nparam);
    
    %----------------------------------------------------------------------
    % Update model
    %----------------------------------------------------------------------
    if length(tijs) > (Nparam - 2*Nstates)
        tijs = tijs(1:(Nparam - 2*Nstates));
    end
    rates = 1 ./ tijs;
    
    [k12, k21, k23, k32, k34, k43] = feval(@(x) x{:}, num2cell(reshape(rates,[1, length(rates)])));
%     k32 = (k12 * k23 * k31)/(k13 * k21);
%     k41 = (k31 * k14 * k43)/ (k13 * k34);
%     k42 = (k24 * k43 * k31 * k12)/(k34 * k13 * k21);
    
    % [A1, A2, A3, A4] = feval(@(x) x{:}, num2cell(As));
    
    K = [-(k12),          k21,                 0,                0;...
        k12,              -(k21 + k23),        k32,               0;...
        0,                k23,              -(k32 + k34),       k43;...
        0,                 0,                 k34,             -(k43)];
    
      elseif strcmp(model_name, 'linear1423') == 1
     Nstates = 4;
    t14_bounds = [200e-6,50000e-3];  %Paramater #1 is high--> med    
%     t13_bounds = [1e-3,3000e-3];    %Paramater #2 is high --> low
%     t14_bounds = [100e-6,900e-3];
    t41_bounds = [200e-6,50000e-3];%Paramater #3 is med --> high
    t42_bounds = [200e-6,50000e-3];%Paramater #4 is med --> low
    t24_bounds = [200e-6,50000e-3];
%     t31_bounds = [50e-6,1000e-3];  %Paramater #5 is low --> Medium = [1e-3,10e-3];%Paramater #5 %t32 is low --> Medium
    t23_bounds = [200e-6,50000e-3];
    t32_bounds = [200e-6,50000e-3];
    % *t32 wll be determined by the other rates
    
% IN THIS SCHEME, 2 AND 3 ARE THE LOW AMP MOSTLY SYMM/ANTI WHILE 1 AND 4
% ARE THE HIGH AP STATES
    % A1_bounds = [0.65,0.85];%Paramater #6 % HIGH fret State
    A1_bounds = [0.04,0.055];%low amp
    A2_bounds = [0.08,0.13];%2nd lowest
    A3_bounds = [0.15,0.2];%2nd high amp
    A4_bounds = [0.25,0.32]; %1st high amp  
    
    %Bounds for standard deviation
     sigma_A1_bounds = [0.028,0.042]; %Paramater #8 %Low Amp state standard deviation
     sigma_A2_bounds = [0.05,0.07];%Paramater #9 %Med Amp state standard deviation
     sigma_A3_bounds = [0.065,0.09];%Paramater #10 %Med2 Amp state standard deviation
     sigma_A4_bounds = [0.09,0.14];%Paramater #10 %High Amp state standard deviation
    
    boundsArray = [t14_bounds; t41_bounds; t42_bounds; t24_bounds; t23_bounds; t32_bounds; ...
        A1_bounds; A2_bounds; A3_bounds; A4_bounds; ...
        sigma_A1_bounds; sigma_A2_bounds; sigma_A3_bounds ; sigma_A4_bounds];
    
    Nparam = length(boundsArray);
    boundsArray_noAdj = boundsArray;
    
    [boundsArrayAdj,atEdges]= boundsAdjuster_v3(boundsArray,windowMag,paramAdjs,absLow,absHigh,Nparam,Nstates,fixedDecades);
    boundsArray=boundsArrayAdj;
    
    maxMutationCountsArray = zeros(1,Nparam);
    minMutationCountsArray = zeros(1,Nparam);
    
    
    % Maximum mutation value
    Max_mut_factor = .2;
    
    maxMutationArray = [Max_mut_factor * boundsArray(1:Nparam-2*Nstates,2) / 2; ...
        0.1 * ( boundsArray(Nparam-2*Nstates+1:end,2) - boundsArray(Nparam-2*Nstates+1:end,1)) ];    
    
    param_strings = {'t14','t41','t42','t24','t23','t32','A1','A2','A3','A4','sig1','sig2','sig3','sig4'};
    % Note: the order of the strings needs to be the same as the variable
    % list so we can extract the parameters properly.
    
    % tij_strings = param_strings(1:(Nparam - Nstates));
    % A_strings = param_strings((Nparam - Nstates + 1):Nparam);
    
    %----------------------------------------------------------------------
    % Update model
    %----------------------------------------------------------------------
    if length(tijs) > (Nparam - 2*Nstates)
        tijs = tijs(1:(Nparam - 2*Nstates));
    end
    rates = 1 ./ tijs;
    
    [k14, k41, k42, k24, k23, k32] = feval(@(x) x{:}, num2cell(reshape(rates,[1, length(rates)])));
%     k32 = (k12 * k23 * k31)/(k13 * k21);
%     k41 = (k31 * k14 * k43)/ (k13 * k34);
%     k42 = (k24 * k43 * k31 * k12)/(k34 * k13 * k21);
    
    % [A1, A2, A3, A4] = feval(@(x) x{:}, num2cell(As));
    
    K = [-(k14),          0,                 0,                k41;...
        0,              -(k23 + k24),        k32,               k42;...
        0,                k23,              -(k32),             0;...
        k14,              k24,               0,             -(k41 + k42)];

    elseif strcmp(model_name, 'linear1243') == 1
     Nstates = 4;
    t12_bounds = [200e-6,50000e-3];  %Paramater #1 is high--> med    
%     t13_bounds = [1e-3,3000e-3];    %Paramater #2 is high --> low
%     t14_bounds = [100e-6,900e-3];
    t21_bounds = [200e-6,50000e-3];%Paramater #3 is med --> high
    t24_bounds = [200e-6,50000e-3];%Paramater #4 is med --> low
    t42_bounds = [200e-6,50000e-3];
%     t31_bounds = [50e-6,1000e-3];  %Paramater #5 is low --> Medium = [1e-3,10e-3];%Paramater #5 %t32 is low --> Medium
    t43_bounds = [200e-6,50000e-3];
    t34_bounds = [200e-6,50000e-3];
    % *t32 wll be determined by the other rates
    
% IN THIS SCHEME, 2 AND 3 ARE THE LOW AMP MOSTLY SYMM/ANTI WHILE 1 AND 4
% ARE THE HIGH AP STATES
    % A1_bounds = [0.65,0.85];%Paramater #6 % HIGH fret State
    A1_bounds = [0.04,0.055];%low amp
    A2_bounds = [0.08,0.13];%2nd lowest
    A3_bounds = [0.15,0.2];%2nd high amp
    A4_bounds = [0.25,0.32]; %1st high amp  
    
    %Bounds for standard deviation
     sigma_A1_bounds = [0.028,0.042]; %Paramater #8 %Low Amp state standard deviation
     sigma_A2_bounds = [0.05,0.07];%Paramater #9 %Med Amp state standard deviation
     sigma_A3_bounds = [0.065,0.09];%Paramater #10 %Med2 Amp state standard deviation
     sigma_A4_bounds = [0.09,0.14];%Paramater #10 %High Amp state standard deviation
    
    boundsArray = [t12_bounds; t21_bounds; t24_bounds; t42_bounds; t43_bounds; t34_bounds; ...
        A1_bounds; A2_bounds; A3_bounds; A4_bounds; ...
        sigma_A1_bounds; sigma_A2_bounds; sigma_A3_bounds ; sigma_A4_bounds];
    
    Nparam = length(boundsArray);
    boundsArray_noAdj = boundsArray;
    
    [boundsArrayAdj,atEdges]= boundsAdjuster_v3(boundsArray,windowMag,paramAdjs,absLow,absHigh,Nparam,Nstates,fixedDecades);
    boundsArray=boundsArrayAdj;
    
    maxMutationCountsArray = zeros(1,Nparam);
    minMutationCountsArray = zeros(1,Nparam);
    
    
    % Maximum mutation value
    Max_mut_factor = .2;
    
    maxMutationArray = [Max_mut_factor * boundsArray(1:Nparam-2*Nstates,2) / 2; ...
        0.1 * ( boundsArray(Nparam-2*Nstates+1:end,2) - boundsArray(Nparam-2*Nstates+1:end,1)) ];    
    
    param_strings = {'t12','t21','t24','t42','t43','t34','A1','A2','A3','A4','sig1','sig2','sig3','sig4'};
    % Note: the order of the strings needs to be the same as the variable
    % list so we can extract the parameters properly.
    
    % tij_strings = param_strings(1:(Nparam - Nstates));
    % A_strings = param_strings((Nparam - Nstates + 1):Nparam);
    
    %----------------------------------------------------------------------
    % Update model
    %----------------------------------------------------------------------
    if length(tijs) > (Nparam - 2*Nstates)
        tijs = tijs(1:(Nparam - 2*Nstates));
    end
    rates = 1 ./ tijs;
    
    [k12, k21, k24, k42, k43, k34] = feval(@(x) x{:}, num2cell(reshape(rates,[1, length(rates)])));
%     k32 = (k12 * k23 * k31)/(k13 * k21);
%     k41 = (k31 * k14 * k43)/ (k13 * k34);
%     k42 = (k24 * k43 * k31 * k12)/(k34 * k13 * k21);
    
    % [A1, A2, A3, A4] = feval(@(x) x{:}, num2cell(As));
    
    K = [-(k12),          k21,                 0,                0;...
        k12,              -(k21 + k24),        0,                 k42;...
        0,                0,                -(k34),               k43;...
        0,                k24,                k34,             -(k42 + k43)];
    
    
     elseif strcmp(model_name, 'linear2134') == 1
     Nstates = 4;
    t21_bounds = [200e-6,50000e-3];  %Paramater #1 is high--> med    
%     t13_bounds = [1e-3,3000e-3];    %Paramater #2 is high --> low
%     t14_bounds = [100e-6,900e-3];
    t12_bounds = [200e-6,50000e-3];%Paramater #3 is med --> high
    t13_bounds = [200e-6,50000e-3];%Paramater #4 is med --> low
    t31_bounds = [200e-6,50000e-3];
%     t31_bounds = [50e-6,1000e-3];  %Paramater #5 is low --> Medium = [1e-3,10e-3];%Paramater #5 %t32 is low --> Medium
    t34_bounds = [200e-6,50000e-3];
    t43_bounds = [200e-6,50000e-3];
    % *t32 wll be determined by the other rates
    
% IN THIS SCHEME, 2 AND 3 ARE THE LOW AMP MOSTLY SYMM/ANTI WHILE 1 AND 4
% ARE THE HIGH AP STATES
    % A1_bounds = [0.65,0.85];%Paramater #6 % HIGH fret State
    A1_bounds = [0.04,0.055];%low amp
    A2_bounds = [0.08,0.13];%2nd lowest
    A3_bounds = [0.15,0.2];%2nd high amp
    A4_bounds = [0.25,0.32]; %1st high amp  
    
    %Bounds for standard deviation
     sigma_A1_bounds = [0.028,0.042]; %Paramater #8 %Low Amp state standard deviation
     sigma_A2_bounds = [0.05,0.07];%Paramater #9 %Med Amp state standard deviation
     sigma_A3_bounds = [0.065,0.09];%Paramater #10 %Med2 Amp state standard deviation
     sigma_A4_bounds = [0.09,0.14];%Paramater #10 %High Amp state standard deviation
    
    boundsArray = [t21_bounds; t12_bounds; t13_bounds; t31_bounds; t34_bounds; t43_bounds; ...
        A1_bounds; A2_bounds; A3_bounds; A4_bounds; ...
        sigma_A1_bounds; sigma_A2_bounds; sigma_A3_bounds ; sigma_A4_bounds];
    
    Nparam = length(boundsArray);
    boundsArray_noAdj = boundsArray;
    
    [boundsArrayAdj,atEdges]= boundsAdjuster_v3(boundsArray,windowMag,paramAdjs,absLow,absHigh,Nparam,Nstates,fixedDecades);
    boundsArray=boundsArrayAdj;
    
    maxMutationCountsArray = zeros(1,Nparam);
    minMutationCountsArray = zeros(1,Nparam);
    
    
    % Maximum mutation value
    Max_mut_factor = .2;
    
    maxMutationArray = [Max_mut_factor * boundsArray(1:Nparam-2*Nstates,2) / 2; ...
        0.1 * ( boundsArray(Nparam-2*Nstates+1:end,2) - boundsArray(Nparam-2*Nstates+1:end,1)) ];    
    
    param_strings = {'t21','t12','t13','t31','t34','t43','A1','A2','A3','A4','sig1','sig2','sig3','sig4'};
    % Note: the order of the strings needs to be the same as the variable
    % list so we can extract the parameters properly.
    
    % tij_strings = param_strings(1:(Nparam - Nstates));
    % A_strings = param_strings((Nparam - Nstates + 1):Nparam);
    
    %----------------------------------------------------------------------
    % Update model
    %----------------------------------------------------------------------
    if length(tijs) > (Nparam - 2*Nstates)
        tijs = tijs(1:(Nparam - 2*Nstates));
    end
    rates = 1 ./ tijs;
    
    [k21, k12, k13, k31, k34, k43] = feval(@(x) x{:}, num2cell(reshape(rates,[1, length(rates)])));
%     k32 = (k12 * k23 * k31)/(k13 * k21);
%     k41 = (k31 * k14 * k43)/ (k13 * k34);
%     k42 = (k24 * k43 * k31 * k12)/(k34 * k13 * k21);
    
    % [A1, A2, A3, A4] = feval(@(x) x{:}, num2cell(As));
    
    K = [-(k12 + k13),       k21,               k31,                0;...
        k12,              -(k21),                0,                 0;...
        k13,                0,                -(k31 + k34),          k43;...
        0,                   0,                k34,               -(k43)];
    
    
    elseif strcmp(model_name, 'linear2143') == 1
     Nstates = 4;
    t21_bounds = [200e-6,50000e-3];  %Paramater #1 is high--> med    
%     t13_bounds = [1e-3,3000e-3];    %Paramater #2 is high --> low
%     t14_bounds = [100e-6,900e-3];
    t12_bounds = [200e-6,50000e-3];%Paramater #3 is med --> high
    t14_bounds = [200e-6,50000e-3];%Paramater #4 is med --> low
    t41_bounds = [200e-6,50000e-3];
%     t31_bounds = [50e-6,1000e-3];  %Paramater #5 is low --> Medium = [1e-3,10e-3];%Paramater #5 %t32 is low --> Medium
    t43_bounds = [200e-6,50000e-3];
    t34_bounds = [200e-6,50000e-3];
    % *t32 wll be determined by the other rates
    
% IN THIS SCHEME, 2 AND 3 ARE THE LOW AMP MOSTLY SYMM/ANTI WHILE 1 AND 4
% ARE THE HIGH AP STATES
    % A1_bounds = [0.65,0.85];%Paramater #6 % HIGH fret State
    A1_bounds = [0.04,0.055];%low amp
    A2_bounds = [0.08,0.13];%2nd lowest
    A3_bounds = [0.15,0.2];%2nd high amp
    A4_bounds = [0.25,0.32]; %1st high amp  
    
    %Bounds for standard deviation
     sigma_A1_bounds = [0.028,0.042]; %Paramater #8 %Low Amp state standard deviation
     sigma_A2_bounds = [0.05,0.07];%Paramater #9 %Med Amp state standard deviation
     sigma_A3_bounds = [0.065,0.09];%Paramater #10 %Med2 Amp state standard deviation
     sigma_A4_bounds = [0.09,0.14];%Paramater #10 %High Amp state standard deviation
    
    boundsArray = [t21_bounds; t12_bounds; t14_bounds; t41_bounds; t43_bounds; t34_bounds; ...
        A1_bounds; A2_bounds; A3_bounds; A4_bounds; ...
        sigma_A1_bounds; sigma_A2_bounds; sigma_A3_bounds ; sigma_A4_bounds];
    
    Nparam = length(boundsArray);
    boundsArray_noAdj = boundsArray;
    
    [boundsArrayAdj,atEdges]= boundsAdjuster_v3(boundsArray,windowMag,paramAdjs,absLow,absHigh,Nparam,Nstates,fixedDecades);
    boundsArray=boundsArrayAdj;
    
    maxMutationCountsArray = zeros(1,Nparam);
    minMutationCountsArray = zeros(1,Nparam);
    
    
    % Maximum mutation value
    Max_mut_factor = .2;
    
    maxMutationArray = [Max_mut_factor * boundsArray(1:Nparam-2*Nstates,2) / 2; ...
        0.1 * ( boundsArray(Nparam-2*Nstates+1:end,2) - boundsArray(Nparam-2*Nstates+1:end,1)) ];    
    
    param_strings = {'t21','t12','t14','t41','t43','t34','A1','A2','A3','A4','sig1','sig2','sig3','sig4'};
    % Note: the order of the strings needs to be the same as the variable
    % list so we can extract the parameters properly.
    
    % tij_strings = param_strings(1:(Nparam - Nstates));
    % A_strings = param_strings((Nparam - Nstates + 1):Nparam);
    
    %----------------------------------------------------------------------
    % Update model
    %----------------------------------------------------------------------
    if length(tijs) > (Nparam - 2*Nstates)
        tijs = tijs(1:(Nparam - 2*Nstates));
    end
    rates = 1 ./ tijs;
    
    [k21, k12, k14, k41, k43, k34] = feval(@(x) x{:}, num2cell(reshape(rates,[1, length(rates)])));
%     k32 = (k12 * k23 * k31)/(k13 * k21);
%     k41 = (k31 * k14 * k43)/ (k13 * k34);
%     k42 = (k24 * k43 * k31 * k12)/(k34 * k13 * k21);
    
    % [A1, A2, A3, A4] = feval(@(x) x{:}, num2cell(As));
    
    K = [-(k12 + k14),       k21,               0,                k41;...
        k12,              -(k21),                0,                 0;...
        0,                   0,                -(k34),          k43;...
        k14,                 0,                k34,               -(k41 + k43)];
    
    elseif strcmp(model_name, 'linear2314') == 1
     Nstates = 4;
    t23_bounds = [200e-6,50000e-3];  %Paramater #1 is high--> med    
%     t13_bounds = [1e-3,3000e-3];    %Paramater #2 is high --> low
%     t14_bounds = [100e-6,900e-3];
    t32_bounds = [200e-6,50000e-3];%Paramater #3 is med --> high
    t31_bounds = [200e-6,50000e-3];%Paramater #4 is med --> low
    t13_bounds = [200e-6,50000e-3];
%     t31_bounds = [50e-6,1000e-3];  %Paramater #5 is low --> Medium = [1e-3,10e-3];%Paramater #5 %t32 is low --> Medium
    t14_bounds = [200e-6,50000e-3];
    t41_bounds = [200e-6,50000e-3];
    % *t32 wll be determined by the other rates
    
% IN THIS SCHEME, 2 AND 3 ARE THE LOW AMP MOSTLY SYMM/ANTI WHILE 1 AND 4
% ARE THE HIGH AP STATES
    % A1_bounds = [0.65,0.85];%Paramater #6 % HIGH fret State
    A1_bounds = [0.04,0.055];%low amp
    A2_bounds = [0.08,0.13];%2nd lowest
    A3_bounds = [0.15,0.2];%2nd high amp
    A4_bounds = [0.25,0.32]; %1st high amp  
    
    %Bounds for standard deviation
     sigma_A1_bounds = [0.028,0.042]; %Paramater #8 %Low Amp state standard deviation
     sigma_A2_bounds = [0.05,0.07];%Paramater #9 %Med Amp state standard deviation
     sigma_A3_bounds = [0.065,0.09];%Paramater #10 %Med2 Amp state standard deviation
     sigma_A4_bounds = [0.09,0.14];%Paramater #10 %High Amp state standard deviation
    
    boundsArray = [t23_bounds; t32_bounds; t31_bounds; t13_bounds; t14_bounds; t41_bounds; ...
        A1_bounds; A2_bounds; A3_bounds; A4_bounds; ...
        sigma_A1_bounds; sigma_A2_bounds; sigma_A3_bounds ; sigma_A4_bounds];
    
    Nparam = length(boundsArray);
    boundsArray_noAdj = boundsArray;
    
    [boundsArrayAdj,atEdges]= boundsAdjuster_v3(boundsArray,windowMag,paramAdjs,absLow,absHigh,Nparam,Nstates,fixedDecades);
    boundsArray=boundsArrayAdj;
    
    maxMutationCountsArray = zeros(1,Nparam);
    minMutationCountsArray = zeros(1,Nparam);
    
    
    % Maximum mutation value
    Max_mut_factor = .2;
    
    maxMutationArray = [Max_mut_factor * boundsArray(1:Nparam-2*Nstates,2) / 2; ...
        0.1 * ( boundsArray(Nparam-2*Nstates+1:end,2) - boundsArray(Nparam-2*Nstates+1:end,1)) ];    
    
    param_strings = {'t23','t32','t31','t13','t14','t41','A1','A2','A3','A4','sig1','sig2','sig3','sig4'};
    % Note: the order of the strings needs to be the same as the variable
    % list so we can extract the parameters properly.
    
    % tij_strings = param_strings(1:(Nparam - Nstates));
    % A_strings = param_strings((Nparam - Nstates + 1):Nparam);
    
    %----------------------------------------------------------------------
    % Update model
    %----------------------------------------------------------------------
    if length(tijs) > (Nparam - 2*Nstates)
        tijs = tijs(1:(Nparam - 2*Nstates));
    end
    rates = 1 ./ tijs;
    
    [k23, k32, k31, k13, k14, k41] = feval(@(x) x{:}, num2cell(reshape(rates,[1, length(rates)])));
%     k32 = (k12 * k23 * k31)/(k13 * k21);
%     k41 = (k31 * k14 * k43)/ (k13 * k34);
%     k42 = (k24 * k43 * k31 * k12)/(k34 * k13 * k21);
    
    % [A1, A2, A3, A4] = feval(@(x) x{:}, num2cell(As));
    
    K = [-(k13 + k14),       0,               k31,                k41;...
        0,              -(k23),                k32,                 0;...
        k13,                  k23,           -(k31 + k32),          0;...
        k14,                 0,                0,               -(k41)];
    
    
     elseif strcmp(model_name, 'linear2413') == 1
     Nstates = 4;
    t24_bounds = [200e-6,50000e-3];  %Paramater #1 is high--> med    
%     t13_bounds = [1e-3,3000e-3];    %Paramater #2 is high --> low
%     t14_bounds = [100e-6,900e-3];
    t42_bounds = [200e-6,50000e-3];%Paramater #3 is med --> high
    t41_bounds = [200e-6,50000e-3];%Paramater #4 is med --> low
    t14_bounds = [200e-6,50000e-3];
%     t31_bounds = [50e-6,1000e-3];  %Paramater #5 is low --> Medium = [1e-3,10e-3];%Paramater #5 %t32 is low --> Medium
    t13_bounds = [200e-6,50000e-3];
    t31_bounds = [200e-6,50000e-3];
    % *t32 wll be determined by the other rates
    
% IN THIS SCHEME, 2 AND 3 ARE THE LOW AMP MOSTLY SYMM/ANTI WHILE 1 AND 4
% ARE THE HIGH AP STATES
    % A1_bounds = [0.65,0.85];%Paramater #6 % HIGH fret State
    A1_bounds = [0.04,0.055];%low amp
    A2_bounds = [0.08,0.13];%2nd lowest
    A3_bounds = [0.15,0.2];%2nd high amp
    A4_bounds = [0.25,0.32]; %1st high amp  
    
    %Bounds for standard deviation
     sigma_A1_bounds = [0.028,0.042]; %Paramater #8 %Low Amp state standard deviation
     sigma_A2_bounds = [0.05,0.07];%Paramater #9 %Med Amp state standard deviation
     sigma_A3_bounds = [0.065,0.09];%Paramater #10 %Med2 Amp state standard deviation
     sigma_A4_bounds = [0.09,0.14];%Paramater #10 %High Amp state standard deviation
    
    boundsArray = [t24_bounds; t42_bounds; t41_bounds; t14_bounds; t13_bounds; t31_bounds; ...
        A1_bounds; A2_bounds; A3_bounds; A4_bounds; ...
        sigma_A1_bounds; sigma_A2_bounds; sigma_A3_bounds ; sigma_A4_bounds];
    
    Nparam = length(boundsArray);
    boundsArray_noAdj = boundsArray;
    
    [boundsArrayAdj,atEdges]= boundsAdjuster_v3(boundsArray,windowMag,paramAdjs,absLow,absHigh,Nparam,Nstates,fixedDecades);
    boundsArray=boundsArrayAdj;
    
    maxMutationCountsArray = zeros(1,Nparam);
    minMutationCountsArray = zeros(1,Nparam);
    
    
    % Maximum mutation value
    Max_mut_factor = .2;
    
    maxMutationArray = [Max_mut_factor * boundsArray(1:Nparam-2*Nstates,2) / 2; ...
        0.1 * ( boundsArray(Nparam-2*Nstates+1:end,2) - boundsArray(Nparam-2*Nstates+1:end,1)) ];    
    
    param_strings = {'t24','t42','t41','t14','t13','t31','A1','A2','A3','A4','sig1','sig2','sig3','sig4'};
    % Note: the order of the strings needs to be the same as the variable
    % list so we can extract the parameters properly.
    
    % tij_strings = param_strings(1:(Nparam - Nstates));
    % A_strings = param_strings((Nparam - Nstates + 1):Nparam);
    
    %----------------------------------------------------------------------
    % Update model
    %----------------------------------------------------------------------
    if length(tijs) > (Nparam - 2*Nstates)
        tijs = tijs(1:(Nparam - 2*Nstates));
    end
    rates = 1 ./ tijs;
    
    [k24, k42, k41, k14, k13, k31] = feval(@(x) x{:}, num2cell(reshape(rates,[1, length(rates)])));
%     k32 = (k12 * k23 * k31)/(k13 * k21);
%     k41 = (k31 * k14 * k43)/ (k13 * k34);
%     k42 = (k24 * k43 * k31 * k12)/(k34 * k13 * k21);
    
    % [A1, A2, A3, A4] = feval(@(x) x{:}, num2cell(As));
    
    K = [-(k13 + k14),       0,               k31,                k41;...
        0,              -(k24),                0,                 k42;...
        k13,                  0,           -(k31),                0;...
        k14,                 k24,                0,               -(k41 + k42)];

    
     elseif strcmp(model_name, 'linear3124') == 1
     Nstates = 4;
    t31_bounds = [200e-6,50000e-3];  %Paramater #1 is high--> med    
%     t13_bounds = [1e-3,3000e-3];    %Paramater #2 is high --> low
%     t14_bounds = [100e-6,900e-3];
    t13_bounds = [200e-6,50000e-3];%Paramater #3 is med --> high
    t12_bounds = [200e-6,50000e-3];%Paramater #4 is med --> low
    t21_bounds = [200e-6,50000e-3];
%     t31_bounds = [50e-6,1000e-3];  %Paramater #5 is low --> Medium = [1e-3,10e-3];%Paramater #5 %t32 is low --> Medium
    t24_bounds = [200e-6,50000e-3];
    t42_bounds = [200e-6,50000e-3];
    % *t32 wll be determined by the other rates
    
% IN THIS SCHEME, 2 AND 3 ARE THE LOW AMP MOSTLY SYMM/ANTI WHILE 1 AND 4
% ARE THE HIGH AP STATES
    % A1_bounds = [0.65,0.85];%Paramater #6 % HIGH fret State
    A1_bounds = [0.04,0.055];%low amp
    A2_bounds = [0.08,0.13];%2nd lowest
    A3_bounds = [0.15,0.2];%2nd high amp
    A4_bounds = [0.25,0.32]; %1st high amp  
    
    %Bounds for standard deviation
     sigma_A1_bounds = [0.028,0.042]; %Paramater #8 %Low Amp state standard deviation
     sigma_A2_bounds = [0.05,0.07];%Paramater #9 %Med Amp state standard deviation
     sigma_A3_bounds = [0.065,0.09];%Paramater #10 %Med2 Amp state standard deviation
     sigma_A4_bounds = [0.09,0.14];%Paramater #10 %High Amp state standard deviation
    
    boundsArray = [t31_bounds; t13_bounds; t12_bounds; t21_bounds; t24_bounds; t42_bounds; ...
        A1_bounds; A2_bounds; A3_bounds; A4_bounds; ...
        sigma_A1_bounds; sigma_A2_bounds; sigma_A3_bounds ; sigma_A4_bounds];
    
    Nparam = length(boundsArray);
    boundsArray_noAdj = boundsArray;
    
    [boundsArrayAdj,atEdges]= boundsAdjuster_v3(boundsArray,windowMag,paramAdjs,absLow,absHigh,Nparam,Nstates,fixedDecades);
    boundsArray=boundsArrayAdj;
    
    maxMutationCountsArray = zeros(1,Nparam);
    minMutationCountsArray = zeros(1,Nparam);
    
    
    % Maximum mutation value
    Max_mut_factor = .2;
    
    maxMutationArray = [Max_mut_factor * boundsArray(1:Nparam-2*Nstates,2) / 2; ...
        0.1 * ( boundsArray(Nparam-2*Nstates+1:end,2) - boundsArray(Nparam-2*Nstates+1:end,1)) ];    
    
    param_strings = {'t31','t13','t12','t21','t24','t42','A1','A2','A3','A4','sig1','sig2','sig3','sig4'};
    % Note: the order of the strings needs to be the same as the variable
    % list so we can extract the parameters properly.
    
    % tij_strings = param_strings(1:(Nparam - Nstates));
    % A_strings = param_strings((Nparam - Nstates + 1):Nparam);
    
    %----------------------------------------------------------------------
    % Update model
    %----------------------------------------------------------------------
    if length(tijs) > (Nparam - 2*Nstates)
        tijs = tijs(1:(Nparam - 2*Nstates));
    end
    rates = 1 ./ tijs;
    
    [k31, k13, k12, k21, k24, k42] = feval(@(x) x{:}, num2cell(reshape(rates,[1, length(rates)])));
%     k32 = (k12 * k23 * k31)/(k13 * k21);
%     k41 = (k31 * k14 * k43)/ (k13 * k34);
%     k42 = (k24 * k43 * k31 * k12)/(k34 * k13 * k21);
    
    % [A1, A2, A3, A4] = feval(@(x) x{:}, num2cell(As));
    
    K = [-(k12 + k13),       k21,               k31,                0;...
        k12,              -(k21 + k24),           0,               k42;...
        k13,                  0,           -(k31),                   0;...
        0,                 k24,                0,               -(k42)];
    
    elseif strcmp(model_name, 'linear3214') == 1
     Nstates = 4;
    t32_bounds = [200e-6,50000e-3];  %Paramater #1 is high--> med    
%     t13_bounds = [1e-3,3000e-3];    %Paramater #2 is high --> low
%     t14_bounds = [100e-6,900e-3];
    t23_bounds = [200e-6,50000e-3];%Paramater #3 is med --> high
    t21_bounds = [200e-6,50000e-3];%Paramater #4 is med --> low
    t12_bounds = [200e-6,50000e-3];
%     t31_bounds = [50e-6,1000e-3];  %Paramater #5 is low --> Medium = [1e-3,10e-3];%Paramater #5 %t32 is low --> Medium
    t14_bounds = [200e-6,50000e-3];
    t41_bounds = [200e-6,50000e-3];
    % *t32 wll be determined by the other rates
    
% IN THIS SCHEME, 2 AND 3 ARE THE LOW AMP MOSTLY SYMM/ANTI WHILE 1 AND 4
% ARE THE HIGH AP STATES
    % A1_bounds = [0.65,0.85];%Paramater #6 % HIGH fret State
    A1_bounds = [0.04,0.055];%low amp
    A2_bounds = [0.08,0.13];%2nd lowest
    A3_bounds = [0.15,0.2];%2nd high amp
    A4_bounds = [0.25,0.32]; %1st high amp  
    
    %Bounds for standard deviation
     sigma_A1_bounds = [0.028,0.042]; %Paramater #8 %Low Amp state standard deviation
     sigma_A2_bounds = [0.05,0.07];%Paramater #9 %Med Amp state standard deviation
     sigma_A3_bounds = [0.065,0.09];%Paramater #10 %Med2 Amp state standard deviation
     sigma_A4_bounds = [0.09,0.14];%Paramater #10 %High Amp state standard deviation
    
    boundsArray = [t32_bounds; t23_bounds; t21_bounds; t12_bounds; t14_bounds; t41_bounds; ...
        A1_bounds; A2_bounds; A3_bounds; A4_bounds; ...
        sigma_A1_bounds; sigma_A2_bounds; sigma_A3_bounds ; sigma_A4_bounds];
    
    Nparam = length(boundsArray);
    boundsArray_noAdj = boundsArray;
    
    [boundsArrayAdj,atEdges]= boundsAdjuster_v3(boundsArray,windowMag,paramAdjs,absLow,absHigh,Nparam,Nstates,fixedDecades);
    boundsArray=boundsArrayAdj;
    
    maxMutationCountsArray = zeros(1,Nparam);
    minMutationCountsArray = zeros(1,Nparam);
    
    
    % Maximum mutation value
    Max_mut_factor = .2;
    
    maxMutationArray = [Max_mut_factor * boundsArray(1:Nparam-2*Nstates,2) / 2; ...
        0.1 * ( boundsArray(Nparam-2*Nstates+1:end,2) - boundsArray(Nparam-2*Nstates+1:end,1)) ];    
    
    param_strings = {'t32','t23','t21','t12','t14','t41','A1','A2','A3','A4','sig1','sig2','sig3','sig4'};
    % Note: the order of the strings needs to be the same as the variable
    % list so we can extract the parameters properly.
    
    % tij_strings = param_strings(1:(Nparam - Nstates));
    % A_strings = param_strings((Nparam - Nstates + 1):Nparam);
    
    %----------------------------------------------------------------------
    % Update model
    %----------------------------------------------------------------------
    if length(tijs) > (Nparam - 2*Nstates)
        tijs = tijs(1:(Nparam - 2*Nstates));
    end
    rates = 1 ./ tijs;
    
    [k32, k23, k21, k12, k14, k41] = feval(@(x) x{:}, num2cell(reshape(rates,[1, length(rates)])));
%     k32 = (k12 * k23 * k31)/(k13 * k21);
%     k41 = (k31 * k14 * k43)/ (k13 * k34);
%     k42 = (k24 * k43 * k31 * k12)/(k34 * k13 * k21);
    
    % [A1, A2, A3, A4] = feval(@(x) x{:}, num2cell(As));
    
    K = [-(k12 + k14),       k21,               0,                 k41;...
        k12,              -(k21 + k23),         k32,               0;...
        0,                   k23,           -(k32),                  0;...
        k14,                 0,                0,               -(k41)];

    
 elseif strcmp(model_name, 'pyramid13') == 1
     Nstates = 4;
    t13_bounds = [200e-6,50000e-3];  %Paramater #1 is high--> med    
%     t13_bounds = [1e-3,3000e-3];    %Paramater #2 is high --> low
%     t14_bounds = [100e-6,900e-3];
    t31_bounds = [200e-6,50000e-3];%Paramater #3 is med --> high
    t23_bounds = [200e-6,50000e-3];%Paramater #4 is med --> low
    t32_bounds = [200e-6,50000e-3];
%     t31_bounds = [50e-6,1000e-3];  %Paramater #5 is low --> Medium = [1e-3,10e-3];%Paramater #5 %t32 is low --> Medium
    t24_bounds = [200e-6,50000e-3];
    t42_bounds = [200e-6,50000e-3];
%     t34_bounds = [200e-6,50000e-3];
    t43_bounds = [200e-6,50000e-3];
    
% IN THIS SCHEME, 2 AND 3 ARE THE LOW AMP MOSTLY SYMM/ANTI WHILE 1 AND 4
% ARE THE HIGH AP STATES
    % A1_bounds = [0.65,0.85];%Paramater #6 % HIGH fret State
    A1_bounds = [0.04,0.055];%low amp
    A2_bounds = [0.08,0.13];%2nd lowest
    A3_bounds = [0.15,0.2];%2nd high amp
    A4_bounds = [0.25,0.32]; %1st high amp  
    
    %Bounds for standard deviation
     sigma_A1_bounds = [0.028,0.042]; %Paramater #8 %Low Amp state standard deviation
     sigma_A2_bounds = [0.05,0.07];%Paramater #9 %Med Amp state standard deviation
     sigma_A3_bounds = [0.065,0.09];%Paramater #10 %Med2 Amp state standard deviation
     sigma_A4_bounds = [0.09,0.14];%Paramater #10 %High Amp state standard deviation
    
    boundsArray = [t13_bounds; t31_bounds; t23_bounds; t32_bounds; t24_bounds; t42_bounds; t43_bounds;...
        A1_bounds; A2_bounds; A3_bounds; A4_bounds; ...
        sigma_A1_bounds; sigma_A2_bounds; sigma_A3_bounds ; sigma_A4_bounds];
    
    Nparam = length(boundsArray);
    boundsArray_noAdj = boundsArray;
    
    [boundsArrayAdj,atEdges]= boundsAdjuster_v3(boundsArray,windowMag,paramAdjs,absLow,absHigh,Nparam,Nstates,fixedDecades);
    boundsArray=boundsArrayAdj;
    
    maxMutationCountsArray = zeros(1,Nparam);
    minMutationCountsArray = zeros(1,Nparam);
    
    
    % Maximum mutation value
    Max_mut_factor = .2;
    
    maxMutationArray = [Max_mut_factor * boundsArray(1:Nparam-2*Nstates,2) / 2; ...
        0.1 * ( boundsArray(Nparam-2*Nstates+1:end,2) - boundsArray(Nparam-2*Nstates+1:end,1)) ];    
    
    param_strings = {'t13','t31','t23','t32','t24','t42','t43','A1','A2','A3','A4','sig1','sig2','sig3','sig4'};
    % Note: the order of the strings needs to be the same as the variable
    % list so we can extract the parameters properly.
    
    % tij_strings = param_strings(1:(Nparam - Nstates));
    % A_strings = param_strings((Nparam - Nstates + 1):Nparam);
    
    %----------------------------------------------------------------------
    % Update model
    %----------------------------------------------------------------------
    if length(tijs) > (Nparam - 2*Nstates)
        tijs = tijs(1:(Nparam - 2*Nstates));
    end
    rates = 1 ./ tijs;
    
    [k13, k31, k23, k32, k24, k42, k43] = feval(@(x) x{:}, num2cell(reshape(rates,[1, length(rates)])));
    k34 = (k32 * k24 * k43)/(k42 * k23);
%     k41 = (k31 * k14 * k43)/ (k13 * k34);
%     k42 = (k24 * k43 * k31 * k12)/(k34 * k13 * k21);
    
    % [A1, A2, A3, A4] = feval(@(x) x{:}, num2cell(As));
    
    K = [-(k13),          0,                k31,               0;...
        0,              -(k23 + k24),       k32,               k42;...
        k13,             k23,         -(k31 + k32 + k34),       k43;...
        0,               k24,                k34,             -(k42 + k43)];  
    
    
    
      elseif strcmp(model_name, 'pyramid12') == 1
     Nstates = 4;
    t12_bounds = [200e-6,50000e-3];  %Paramater #1 is high--> med    
%     t13_bounds = [1e-3,3000e-3];    %Paramater #2 is high --> low
%     t14_bounds = [100e-6,900e-3];
    t21_bounds = [200e-6,50000e-3];%Paramater #3 is med --> high
    t23_bounds = [200e-6,50000e-3];%Paramater #4 is med --> low
    t32_bounds = [200e-6,50000e-3];
%     t31_bounds = [50e-6,1000e-3];  %Paramater #5 is low --> Medium = [1e-3,10e-3];%Paramater #5 %t32 is low --> Medium
    t34_bounds = [200e-6,50000e-3];
    t43_bounds = [200e-6,50000e-3];
%     t24_bounds = [200e-6,50000e-3];
    t42_bounds = [200e-6,50000e-3];
    % *t32 wll be determined by the other rates
    
% IN THIS SCHEME, 2 AND 3 ARE THE LOW AMP MOSTLY SYMM/ANTI WHILE 1 AND 4
% ARE THE HIGH AP STATES
    % A1_bounds = [0.65,0.85];%Paramater #6 % HIGH fret State
    A1_bounds = [0.04,0.055];%low amp
    A2_bounds = [0.08,0.13];%2nd lowest
    A3_bounds = [0.15,0.2];%2nd high amp
    A4_bounds = [0.25,0.32]; %1st high amp  
    
    %Bounds for standard deviation
     sigma_A1_bounds = [0.028,0.042]; %Paramater #8 %Low Amp state standard deviation
     sigma_A2_bounds = [0.05,0.07];%Paramater #9 %Med Amp state standard deviation
     sigma_A3_bounds = [0.065,0.09];%Paramater #10 %Med2 Amp state standard deviation
     sigma_A4_bounds = [0.09,0.14];%Paramater #10 %High Amp state standard deviation
    
    boundsArray = [t12_bounds; t21_bounds; t23_bounds; t32_bounds; t34_bounds; t43_bounds; t42_bounds;...
        A1_bounds; A2_bounds; A3_bounds; A4_bounds; ...
        sigma_A1_bounds; sigma_A2_bounds; sigma_A3_bounds ; sigma_A4_bounds];
    
    Nparam = length(boundsArray);
    boundsArray_noAdj = boundsArray;
    
    [boundsArrayAdj,atEdges]= boundsAdjuster_v3(boundsArray,windowMag,paramAdjs,absLow,absHigh,Nparam,Nstates,fixedDecades);
    boundsArray=boundsArrayAdj;
    
    maxMutationCountsArray = zeros(1,Nparam);
    minMutationCountsArray = zeros(1,Nparam);
    
    
    % Maximum mutation value
    Max_mut_factor = .2;
    
    maxMutationArray = [Max_mut_factor * boundsArray(1:Nparam-2*Nstates,2) / 2; ...
        0.1 * ( boundsArray(Nparam-2*Nstates+1:end,2) - boundsArray(Nparam-2*Nstates+1:end,1)) ];    
    
    param_strings = {'t12','t21','t23','t32','t34','t43','t42','A1','A2','A3','A4','sig1','sig2','sig3','sig4'};
    % Note: the order of the strings needs to be the same as the variable
    % list so we can extract the parameters properly.
    
    % tij_strings = param_strings(1:(Nparam - Nstates));
    % A_strings = param_strings((Nparam - Nstates + 1):Nparam);
    
    %----------------------------------------------------------------------
    % Update model
    %----------------------------------------------------------------------
    if length(tijs) > (Nparam - 2*Nstates)
        tijs = tijs(1:(Nparam - 2*Nstates));
    end
    rates = 1 ./ tijs;
    
    [k12, k21, k23, k32, k34, k43, k42] = feval(@(x) x{:}, num2cell(reshape(rates,[1, length(rates)])));
    k24 = (k23 * k34 * k42)/(k43 * k32);
%     k41 = (k31 * k14 * k43)/ (k13 * k34);
%     k42 = (k24 * k43 * k31 * k12)/(k34 * k13 * k21);
    
    % [A1, A2, A3, A4] = feval(@(x) x{:}, num2cell(As));
    
    K = [-(k12),          k21,                 0,                0;...
        k12,       -(k21 + k23 + k24),        k32,               k42;...
        0,                k23,              -(k32 + k34),       k43;...
        0,                 k24,                 k34,             -(k43 + k42)];
    
    
    
     elseif strcmp(model_name, 'pyramid14') == 1
     Nstates = 4;
    t14_bounds = [200e-6,50000e-3];  %Paramater #1 is high--> med    
%     t13_bounds = [1e-3,3000e-3];    %Paramater #2 is high --> low
%     t14_bounds = [100e-6,900e-3];
    t41_bounds = [200e-6,50000e-3];%Paramater #3 is med --> high
    t42_bounds = [200e-6,50000e-3];%Paramater #4 is med --> low
    t24_bounds = [200e-6,50000e-3];
%     t31_bounds = [50e-6,1000e-3];  %Paramater #5 is low --> Medium = [1e-3,10e-3];%Paramater #5 %t32 is low --> Medium
    t23_bounds = [200e-6,50000e-3];
    t32_bounds = [200e-6,50000e-3];
%     t43_bounds = [200e-6,50000e-3];
    t34_bounds = [200e-6,50000e-3];
    % *t32 wll be determined by the other rates
    
% IN THIS SCHEME, 2 AND 3 ARE THE LOW AMP MOSTLY SYMM/ANTI WHILE 1 AND 4
% ARE THE HIGH AP STATES
    % A1_bounds = [0.65,0.85];%Paramater #6 % HIGH fret State
    A1_bounds = [0.04,0.055];%low amp
    A2_bounds = [0.08,0.13];%2nd lowest
    A3_bounds = [0.15,0.2];%2nd high amp
    A4_bounds = [0.25,0.32]; %1st high amp  
    
    %Bounds for standard deviation
     sigma_A1_bounds = [0.028,0.042]; %Paramater #8 %Low Amp state standard deviation
     sigma_A2_bounds = [0.05,0.07];%Paramater #9 %Med Amp state standard deviation
     sigma_A3_bounds = [0.065,0.09];%Paramater #10 %Med2 Amp state standard deviation
     sigma_A4_bounds = [0.09,0.14];%Paramater #10 %High Amp state standard deviation
    
    boundsArray = [t14_bounds; t41_bounds; t42_bounds; t24_bounds; t23_bounds; t32_bounds; t34_bounds;...
        A1_bounds; A2_bounds; A3_bounds; A4_bounds; ...
        sigma_A1_bounds; sigma_A2_bounds; sigma_A3_bounds ; sigma_A4_bounds];
    
    Nparam = length(boundsArray);
    boundsArray_noAdj = boundsArray;
    
    [boundsArrayAdj,atEdges]= boundsAdjuster_v3(boundsArray,windowMag,paramAdjs,absLow,absHigh,Nparam,Nstates,fixedDecades);
    boundsArray=boundsArrayAdj;
    
    maxMutationCountsArray = zeros(1,Nparam);
    minMutationCountsArray = zeros(1,Nparam);
    
    
    % Maximum mutation value
    Max_mut_factor = .2;
    
    maxMutationArray = [Max_mut_factor * boundsArray(1:Nparam-2*Nstates,2) / 2; ...
        0.1 * ( boundsArray(Nparam-2*Nstates+1:end,2) - boundsArray(Nparam-2*Nstates+1:end,1)) ];    
    
    param_strings = {'t14','t41','t42','t24','t23','t32','t34','A1','A2','A3','A4','sig1','sig2','sig3','sig4'};
    % Note: the order of the strings needs to be the same as the variable
    % list so we can extract the parameters properly.
    
    % tij_strings = param_strings(1:(Nparam - Nstates));
    % A_strings = param_strings((Nparam - Nstates + 1):Nparam);
    
    %----------------------------------------------------------------------
    % Update model
    %----------------------------------------------------------------------
    if length(tijs) > (Nparam - 2*Nstates)
        tijs = tijs(1:(Nparam - 2*Nstates));
    end
    rates = 1 ./ tijs;
    
    [k14, k41, k42, k24, k23, k32, k34] = feval(@(x) x{:}, num2cell(reshape(rates,[1, length(rates)])));
    k43 = (k42 * k23 * k34)/(k32 * k24);
%     k41 = (k31 * k14 * k43)/ (k13 * k34);
%     k42 = (k24 * k43 * k31 * k12)/(k34 * k13 * k21);
    
    % [A1, A2, A3, A4] = feval(@(x) x{:}, num2cell(As));
    
    K = [-(k14),          0,                 0,                k41; ...
        0,              -(k23 + k24),        k32,               k42; ...
        0,                k23,          -(k32 + k34),             k43; ...
        k14,              k24,               k34,         -(k41 + k42 + k43)];
    
    
    
     elseif strcmp(model_name, 'pyramid21') == 1
     Nstates = 4;
    t21_bounds = [200e-6,50000e-3];  %Paramater #1 is high--> med    
%     t13_bounds = [1e-3,3000e-3];    %Paramater #2 is high --> low
%     t14_bounds = [100e-6,900e-3];
    t12_bounds = [200e-6,50000e-3];%Paramater #3 is med --> high
    t13_bounds = [200e-6,50000e-3];%Paramater #4 is med --> low
    t31_bounds = [200e-6,50000e-3];
%     t31_bounds = [50e-6,1000e-3];  %Paramater #5 is low --> Medium = [1e-3,10e-3];%Paramater #5 %t32 is low --> Medium
    t34_bounds = [200e-6,50000e-3];
    t43_bounds = [200e-6,50000e-3];
%     t14_bounds = [200e-6,50000e-3];
    t41_bounds = [200e-6,50000e-3];
    % *t32 wll be determined by the other rates
    
% IN THIS SCHEME, 2 AND 3 ARE THE LOW AMP MOSTLY SYMM/ANTI WHILE 1 AND 4
% ARE THE HIGH AP STATES
    % A1_bounds = [0.65,0.85];%Paramater #6 % HIGH fret State
    A1_bounds = [0.04,0.055];%low amp
    A2_bounds = [0.08,0.13];%2nd lowest
    A3_bounds = [0.15,0.2];%2nd high amp
    A4_bounds = [0.25,0.32]; %1st high amp  
    
    %Bounds for standard deviation
     sigma_A1_bounds = [0.028,0.042]; %Paramater #8 %Low Amp state standard deviation
     sigma_A2_bounds = [0.05,0.07];%Paramater #9 %Med Amp state standard deviation
     sigma_A3_bounds = [0.065,0.09];%Paramater #10 %Med2 Amp state standard deviation
     sigma_A4_bounds = [0.09,0.14];%Paramater #10 %High Amp state standard deviation
    
    boundsArray = [t21_bounds; t12_bounds; t13_bounds; t31_bounds; t34_bounds; t43_bounds; t41_bounds; ...
        A1_bounds; A2_bounds; A3_bounds; A4_bounds; ...
        sigma_A1_bounds; sigma_A2_bounds; sigma_A3_bounds ; sigma_A4_bounds];
    
    Nparam = length(boundsArray);
    boundsArray_noAdj = boundsArray;
    
    [boundsArrayAdj,atEdges]= boundsAdjuster_v3(boundsArray,windowMag,paramAdjs,absLow,absHigh,Nparam,Nstates,fixedDecades);
    boundsArray=boundsArrayAdj;
    
    maxMutationCountsArray = zeros(1,Nparam);
    minMutationCountsArray = zeros(1,Nparam);
    
    
    % Maximum mutation value
    Max_mut_factor = .2;
    
    maxMutationArray = [Max_mut_factor * boundsArray(1:Nparam-2*Nstates,2) / 2; ...
        0.1 * ( boundsArray(Nparam-2*Nstates+1:end,2) - boundsArray(Nparam-2*Nstates+1:end,1)) ];    
    
    param_strings = {'t21','t12','t13','t31','t34','t43','t41','A1','A2','A3','A4','sig1','sig2','sig3','sig4'};
    % Note: the order of the strings needs to be the same as the variable
    % list so we can extract the parameters properly.
    
    % tij_strings = param_strings(1:(Nparam - Nstates));
    % A_strings = param_strings((Nparam - Nstates + 1):Nparam);
    
    %----------------------------------------------------------------------
    % Update model
    %----------------------------------------------------------------------
    if length(tijs) > (Nparam - 2*Nstates)
        tijs = tijs(1:(Nparam - 2*Nstates));
    end
    rates = 1 ./ tijs;
    
    [k21, k12, k13, k31, k34, k43, k41] = feval(@(x) x{:}, num2cell(reshape(rates,[1, length(rates)])));
    k14 = (k13 * k34 * k41)/(k43 * k31);
%     k41 = (k31 * k14 * k43)/ (k13 * k34);
%     k42 = (k24 * k43 * k31 * k12)/(k34 * k13 * k21);
    
    % [A1, A2, A3, A4] = feval(@(x) x{:}, num2cell(As));
    
    K = [-(k12 + k13 + k14),   k21,               k31,                k41;...
        k12,              -(k21),                0,                 0;...
        k13,                0,                -(k31 + k34),          k43;...
        k14,                   0,                k34,               -(k43 + k41)];
   
    
     elseif strcmp(model_name, 'pyramid23') == 1
     Nstates = 4;
    t23_bounds = [200e-6,50000e-3];  %Paramater #1 is high--> med    
%     t13_bounds = [1e-3,3000e-3];    %Paramater #2 is high --> low
%     t14_bounds = [100e-6,900e-3];
    t32_bounds = [200e-6,50000e-3];%Paramater #3 is med --> high
    t31_bounds = [200e-6,50000e-3];%Paramater #4 is med --> low
    t13_bounds = [200e-6,50000e-3];
%     t31_bounds = [50e-6,1000e-3];  %Paramater #5 is low --> Medium = [1e-3,10e-3];%Paramater #5 %t32 is low --> Medium
    t14_bounds = [200e-6,50000e-3];
    t41_bounds = [200e-6,50000e-3];
%     t34_bounds = [200e-6,50000e-3];
    t43_bounds = [200e-6,50000e-3];
    % *t32 wll be determined by the other rates
    
% IN THIS SCHEME, 2 AND 3 ARE THE LOW AMP MOSTLY SYMM/ANTI WHILE 1 AND 4
% ARE THE HIGH AP STATES
    % A1_bounds = [0.65,0.85];%Paramater #6 % HIGH fret State
    A1_bounds = [0.04,0.055];%low amp
    A2_bounds = [0.08,0.13];%2nd lowest
    A3_bounds = [0.15,0.2];%2nd high amp
    A4_bounds = [0.25,0.32]; %1st high amp  
    
    %Bounds for standard deviation
     sigma_A1_bounds = [0.028,0.042]; %Paramater #8 %Low Amp state standard deviation
     sigma_A2_bounds = [0.05,0.07];%Paramater #9 %Med Amp state standard deviation
     sigma_A3_bounds = [0.065,0.09];%Paramater #10 %Med2 Amp state standard deviation
     sigma_A4_bounds = [0.09,0.14];%Paramater #10 %High Amp state standard deviation
    
    boundsArray = [t23_bounds; t32_bounds; t31_bounds; t13_bounds; t14_bounds; t41_bounds; t43_bounds;...
        A1_bounds; A2_bounds; A3_bounds; A4_bounds; ...
        sigma_A1_bounds; sigma_A2_bounds; sigma_A3_bounds ; sigma_A4_bounds];
    
    Nparam = length(boundsArray);
    boundsArray_noAdj = boundsArray;
    
    [boundsArrayAdj,atEdges]= boundsAdjuster_v3(boundsArray,windowMag,paramAdjs,absLow,absHigh,Nparam,Nstates,fixedDecades);
    boundsArray=boundsArrayAdj;
    
    maxMutationCountsArray = zeros(1,Nparam);
    minMutationCountsArray = zeros(1,Nparam);
    
    
    % Maximum mutation value
    Max_mut_factor = .2;
    
    maxMutationArray = [Max_mut_factor * boundsArray(1:Nparam-2*Nstates,2) / 2; ...
        0.1 * ( boundsArray(Nparam-2*Nstates+1:end,2) - boundsArray(Nparam-2*Nstates+1:end,1)) ];    
    
    param_strings = {'t23','t32','t31','t13','t14','t41','t43','A1','A2','A3','A4','sig1','sig2','sig3','sig4'};
    % Note: the order of the strings needs to be the same as the variable
    % list so we can extract the parameters properly.
    
    % tij_strings = param_strings(1:(Nparam - Nstates));
    % A_strings = param_strings((Nparam - Nstates + 1):Nparam);
    
    %----------------------------------------------------------------------
    % Update model
    %----------------------------------------------------------------------
    if length(tijs) > (Nparam - 2*Nstates)
        tijs = tijs(1:(Nparam - 2*Nstates));
    end
    rates = 1 ./ tijs;
    
    [k23, k32, k31, k13, k14, k41, k43] = feval(@(x) x{:}, num2cell(reshape(rates,[1, length(rates)])));
    k34 = (k31 * k14 * k43)/(k41 * k13);
%     k41 = (k31 * k14 * k43)/ (k13 * k34);
%     k42 = (k24 * k43 * k31 * k12)/(k34 * k13 * k21);
    
    % [A1, A2, A3, A4] = feval(@(x) x{:}, num2cell(As));
    
    K = [-(k13 + k14),       0,               k31,                k41;...
        0,              -(k23),                k32,                 0;...
        k13,                  k23,     -(k31 + k32 + k34),          k43;...
        k14,                 0,                k34,            -(k41 + k43)];
    
    
    
    
    elseif strcmp(model_name, 'pyramid24') == 1
     Nstates = 4;
    t24_bounds = [200e-6,50000e-3];  %Paramater #1 is high--> med    
%     t13_bounds = [1e-3,3000e-3];    %Paramater #2 is high --> low
%     t14_bounds = [100e-6,900e-3];
    t42_bounds = [200e-6,50000e-3];%Paramater #3 is med --> high
    t41_bounds = [200e-6,50000e-3];%Paramater #4 is med --> low
    t14_bounds = [200e-6,50000e-3];
%     t31_bounds = [50e-6,1000e-3];  %Paramater #5 is low --> Medium = [1e-3,10e-3];%Paramater #5 %t32 is low --> Medium
    t13_bounds = [200e-6,50000e-3];
    t31_bounds = [200e-6,50000e-3];
%     t43_bounds = [200e-6,50000e-3];
    t34_bounds = [200e-6,50000e-3];
    % *t32 wll be determined by the other rates
    
% IN THIS SCHEME, 2 AND 3 ARE THE LOW AMP MOSTLY SYMM/ANTI WHILE 1 AND 4
% ARE THE HIGH AP STATES
    % A1_bounds = [0.65,0.85];%Paramater #6 % HIGH fret State
    A1_bounds = [0.04,0.055];%low amp
    A2_bounds = [0.08,0.13];%2nd lowest
    A3_bounds = [0.15,0.2];%2nd high amp
    A4_bounds = [0.25,0.32]; %1st high amp  
    
    %Bounds for standard deviation
     sigma_A1_bounds = [0.028,0.042]; %Paramater #8 %Low Amp state standard deviation
     sigma_A2_bounds = [0.05,0.07];%Paramater #9 %Med Amp state standard deviation
     sigma_A3_bounds = [0.065,0.09];%Paramater #10 %Med2 Amp state standard deviation
     sigma_A4_bounds = [0.09,0.14];%Paramater #10 %High Amp state standard deviation
    
    boundsArray = [t24_bounds; t42_bounds; t41_bounds; t14_bounds; t13_bounds; t31_bounds; t34_bounds;...
        A1_bounds; A2_bounds; A3_bounds; A4_bounds; ...
        sigma_A1_bounds; sigma_A2_bounds; sigma_A3_bounds ; sigma_A4_bounds];
    
    Nparam = length(boundsArray);
    boundsArray_noAdj = boundsArray;
    
    [boundsArrayAdj,atEdges]= boundsAdjuster_v3(boundsArray,windowMag,paramAdjs,absLow,absHigh,Nparam,Nstates,fixedDecades);
    boundsArray=boundsArrayAdj;
    
    maxMutationCountsArray = zeros(1,Nparam);
    minMutationCountsArray = zeros(1,Nparam);
    
    
    % Maximum mutation value
    Max_mut_factor = .2;
    
    maxMutationArray = [Max_mut_factor * boundsArray(1:Nparam-2*Nstates,2) / 2; ...
        0.1 * ( boundsArray(Nparam-2*Nstates+1:end,2) - boundsArray(Nparam-2*Nstates+1:end,1)) ];    
    
    param_strings = {'t24','t42','t41','t14','t13','t31','t34','A1','A2','A3','A4','sig1','sig2','sig3','sig4'};
    % Note: the order of the strings needs to be the same as the variable
    % list so we can extract the parameters properly.
    
    % tij_strings = param_strings(1:(Nparam - Nstates));
    % A_strings = param_strings((Nparam - Nstates + 1):Nparam);
    
    %----------------------------------------------------------------------
    % Update model
    %----------------------------------------------------------------------
    if length(tijs) > (Nparam - 2*Nstates)
        tijs = tijs(1:(Nparam - 2*Nstates));
    end
    rates = 1 ./ tijs;
    
    [k24, k42, k41, k14, k13, k31, k34] = feval(@(x) x{:}, num2cell(reshape(rates,[1, length(rates)])));
    k43 = (k41 * k13 * k34)/(k31 * k14);
%     k41 = (k31 * k14 * k43)/ (k13 * k34);
%     k42 = (k24 * k43 * k31 * k12)/(k34 * k13 * k21);
    
    % [A1, A2, A3, A4] = feval(@(x) x{:}, num2cell(As));
    
    K = [-(k13 + k14),       0,               k31,                k41;...
        0,              -(k24),                0,                 k42;...
        k13,                  0,           -(k31 + k34),          k43;...
        k14,                 k24,         k34,               -(k41 + k42 + k43)];

    
    elseif strcmp(model_name, 'pyramid31') == 1
     Nstates = 4;
    t31_bounds = [200e-6,50000e-3];  %Paramater #1 is high--> med    
%     t13_bounds = [1e-3,3000e-3];    %Paramater #2 is high --> low
%     t14_bounds = [100e-6,900e-3];
    t13_bounds = [200e-6,50000e-3];%Paramater #3 is med --> high
    t12_bounds = [200e-6,50000e-3];%Paramater #4 is med --> low
    t21_bounds = [200e-6,50000e-3];
%     t31_bounds = [50e-6,1000e-3];  %Paramater #5 is low --> Medium = [1e-3,10e-3];%Paramater #5 %t32 is low --> Medium
    t24_bounds = [200e-6,50000e-3];
    t42_bounds = [200e-6,50000e-3];
%     t14_bounds = [200e-6,50000e-3];
    t41_bounds = [200e-6,50000e-3];
    % *t32 wll be determined by the other rates
    
% IN THIS SCHEME, 2 AND 3 ARE THE LOW AMP MOSTLY SYMM/ANTI WHILE 1 AND 4
% ARE THE HIGH AP STATES
    % A1_bounds = [0.65,0.85];%Paramater #6 % HIGH fret State
    A1_bounds = [0.04,0.055];%low amp
    A2_bounds = [0.08,0.13];%2nd lowest
    A3_bounds = [0.15,0.2];%2nd high amp
    A4_bounds = [0.25,0.32]; %1st high amp  
    
    %Bounds for standard deviation
     sigma_A1_bounds = [0.028,0.042]; %Paramater #8 %Low Amp state standard deviation
     sigma_A2_bounds = [0.05,0.07];%Paramater #9 %Med Amp state standard deviation
     sigma_A3_bounds = [0.065,0.09];%Paramater #10 %Med2 Amp state standard deviation
     sigma_A4_bounds = [0.09,0.14];%Paramater #10 %High Amp state standard deviation
    
    boundsArray = [t31_bounds; t13_bounds; t12_bounds; t21_bounds; t24_bounds; t42_bounds; t41_bounds; ...
        A1_bounds; A2_bounds; A3_bounds; A4_bounds; ...
        sigma_A1_bounds; sigma_A2_bounds; sigma_A3_bounds ; sigma_A4_bounds];
    
    Nparam = length(boundsArray);
    boundsArray_noAdj = boundsArray;
    
    [boundsArrayAdj,atEdges]= boundsAdjuster_v3(boundsArray,windowMag,paramAdjs,absLow,absHigh,Nparam,Nstates,fixedDecades);
    boundsArray=boundsArrayAdj;
    
    maxMutationCountsArray = zeros(1,Nparam);
    minMutationCountsArray = zeros(1,Nparam);
    
    
    % Maximum mutation value
    Max_mut_factor = .2;
    
    maxMutationArray = [Max_mut_factor * boundsArray(1:Nparam-2*Nstates,2) / 2; ...
        0.1 * ( boundsArray(Nparam-2*Nstates+1:end,2) - boundsArray(Nparam-2*Nstates+1:end,1)) ];    
    
    param_strings = {'t31','t13','t12','t21','t24','t42','t41','A1','A2','A3','A4','sig1','sig2','sig3','sig4'};
    % Note: the order of the strings needs to be the same as the variable
    % list so we can extract the parameters properly.
    
    % tij_strings = param_strings(1:(Nparam - Nstates));
    % A_strings = param_strings((Nparam - Nstates + 1):Nparam);
    
    %----------------------------------------------------------------------
    % Update model
    %----------------------------------------------------------------------
    if length(tijs) > (Nparam - 2*Nstates)
        tijs = tijs(1:(Nparam - 2*Nstates));
    end
    rates = 1 ./ tijs;
    
    [k31, k13, k12, k21, k24, k42, k41] = feval(@(x) x{:}, num2cell(reshape(rates,[1, length(rates)])));
    k14 = (k12 * k24 * k41)/(k42 * k21);
%     k41 = (k31 * k14 * k43)/ (k13 * k34);
%     k42 = (k24 * k43 * k31 * k12)/(k34 * k13 * k21);
    
    % [A1, A2, A3, A4] = feval(@(x) x{:}, num2cell(As));
    
    K = [-(k12 + k13 + k14),  k21,               k31,                k41;...
        k12,              -(k21 + k24),           0,               k42;...
        k13,                  0,           -(k31),                   0;...
        k14,                 k24,                0,             -(k42 + k41)];
    
    elseif strcmp(model_name, 'pyramid32') == 1
     Nstates = 4;
    t32_bounds = [200e-6,50000e-3];  %Paramater #1 is high--> med    
%     t13_bounds = [1e-3,3000e-3];    %Paramater #2 is high --> low
%     t14_bounds = [100e-6,900e-3];
    t23_bounds = [200e-6,50000e-3];%Paramater #3 is med --> high
    t21_bounds = [200e-6,50000e-3];%Paramater #4 is med --> low
    t12_bounds = [200e-6,50000e-3];
%     t31_bounds = [50e-6,1000e-3];  %Paramater #5 is low --> Medium = [1e-3,10e-3];%Paramater #5 %t32 is low --> Medium
    t14_bounds = [200e-6,50000e-3];
    t41_bounds = [200e-6,50000e-3];
%     t24_bounds = [200e-6,50000e-3];
    t42_bounds = [200e-6,50000e-3];
    % *t32 wll be determined by the other rates
    
% IN THIS SCHEME, 2 AND 3 ARE THE LOW AMP MOSTLY SYMM/ANTI WHILE 1 AND 4
% ARE THE HIGH AP STATES
    % A1_bounds = [0.65,0.85];%Paramater #6 % HIGH fret State
    A1_bounds = [0.04,0.055];%low amp
    A2_bounds = [0.08,0.13];%2nd lowest
    A3_bounds = [0.15,0.2];%2nd high amp
    A4_bounds = [0.25,0.32]; %1st high amp  
    
    %Bounds for standard deviation
     sigma_A1_bounds = [0.028,0.042]; %Paramater #8 %Low Amp state standard deviation
     sigma_A2_bounds = [0.05,0.07];%Paramater #9 %Med Amp state standard deviation
     sigma_A3_bounds = [0.065,0.09];%Paramater #10 %Med2 Amp state standard deviation
     sigma_A4_bounds = [0.09,0.14];%Paramater #10 %High Amp state standard deviation
    
    boundsArray = [t32_bounds; t23_bounds; t21_bounds; t12_bounds; t14_bounds; t41_bounds; t42_bounds;...
        A1_bounds; A2_bounds; A3_bounds; A4_bounds; ...
        sigma_A1_bounds; sigma_A2_bounds; sigma_A3_bounds ; sigma_A4_bounds];
    
    Nparam = length(boundsArray);
    boundsArray_noAdj = boundsArray;
    
    [boundsArrayAdj,atEdges]= boundsAdjuster_v3(boundsArray,windowMag,paramAdjs,absLow,absHigh,Nparam,Nstates,fixedDecades);
    boundsArray=boundsArrayAdj;
    
    maxMutationCountsArray = zeros(1,Nparam);
    minMutationCountsArray = zeros(1,Nparam);
    
    
    % Maximum mutation value
    Max_mut_factor = .2;
    
    maxMutationArray = [Max_mut_factor * boundsArray(1:Nparam-2*Nstates,2) / 2; ...
        0.1 * ( boundsArray(Nparam-2*Nstates+1:end,2) - boundsArray(Nparam-2*Nstates+1:end,1)) ];    
    
    param_strings = {'t32','t23','t21','t12','t14','t41','t42','A1','A2','A3','A4','sig1','sig2','sig3','sig4'};
    % Note: the order of the strings needs to be the same as the variable
    % list so we can extract the parameters properly.
    
    % tij_strings = param_strings(1:(Nparam - Nstates));
    % A_strings = param_strings((Nparam - Nstates + 1):Nparam);
    
    %----------------------------------------------------------------------
    % Update model
    %----------------------------------------------------------------------
    if length(tijs) > (Nparam - 2*Nstates)
        tijs = tijs(1:(Nparam - 2*Nstates));
    end
    rates = 1 ./ tijs;
    
    [k32, k23, k21, k12, k14, k41, k42] = feval(@(x) x{:}, num2cell(reshape(rates,[1, length(rates)])));
    k24 = (k21 * k14 * k42)/(k41 * k12);
%     k41 = (k31 * k14 * k43)/ (k13 * k34);
%     k42 = (k24 * k43 * k31 * k12)/(k34 * k13 * k21);
    
    % [A1, A2, A3, A4] = feval(@(x) x{:}, num2cell(As));
    
    K = [-(k12 + k14),       k21,               0,                 k41;...
        k12,         -(k21 + k23 + k24),         k32,               k42;...
        0,                   k23,           -(k32),                  0;...
        k14,                 k24,                0,             -(k41 + k42)];
    
    
    elseif strcmp(model_name, 'pyramid34') == 1
     Nstates = 4;
    t12_bounds = [200e-6,50000e-3];  %Paramater #1 is high--> med    
%     t13_bounds = [1e-3,3000e-3];    %Paramater #2 is high --> low
%     t14_bounds = [100e-6,900e-3];
    t21_bounds = [200e-6,50000e-3];%Paramater #3 is med --> high
    t24_bounds = [200e-6,50000e-3];%Paramater #4 is med --> low
    t42_bounds = [200e-6,50000e-3];
%     t31_bounds = [50e-6,1000e-3];  %Paramater #5 is low --> Medium = [1e-3,10e-3];%Paramater #5 %t32 is low --> Medium
    t43_bounds = [200e-6,50000e-3];
    t34_bounds = [200e-6,50000e-3];
%   t14_bounds = [200e-6,50000e-3];
    t41_bounds = [200e-6,50000e-3];
    % *t32 wll be determined by the other rates
    
% IN THIS SCHEME, 2 AND 3 ARE THE LOW AMP MOSTLY SYMM/ANTI WHILE 1 AND 4
% ARE THE HIGH AP STATES
    % A1_bounds = [0.65,0.85];%Paramater #6 % HIGH fret State
    A1_bounds = [0.04,0.055];%low amp
    A2_bounds = [0.08,0.13];%2nd lowest
    A3_bounds = [0.15,0.2];%2nd high amp
    A4_bounds = [0.25,0.32]; %1st high amp  
    
    %Bounds for standard deviation
     sigma_A1_bounds = [0.028,0.042]; %Paramater #8 %Low Amp state standard deviation
     sigma_A2_bounds = [0.05,0.07];%Paramater #9 %Med Amp state standard deviation
     sigma_A3_bounds = [0.065,0.09];%Paramater #10 %Med2 Amp state standard deviation
     sigma_A4_bounds = [0.09,0.14];%Paramater #10 %High Amp state standard deviation
    
    boundsArray = [t12_bounds; t21_bounds; t24_bounds; t42_bounds; t43_bounds; t34_bounds; t41_bounds;...
        A1_bounds; A2_bounds; A3_bounds; A4_bounds; ...
        sigma_A1_bounds; sigma_A2_bounds; sigma_A3_bounds ; sigma_A4_bounds];
    
    Nparam = length(boundsArray);
    boundsArray_noAdj = boundsArray;
    
    [boundsArrayAdj,atEdges]= boundsAdjuster_v3(boundsArray,windowMag,paramAdjs,absLow,absHigh,Nparam,Nstates,fixedDecades);
    boundsArray=boundsArrayAdj;
    
    maxMutationCountsArray = zeros(1,Nparam);
    minMutationCountsArray = zeros(1,Nparam);
    
    
    % Maximum mutation value
    Max_mut_factor = .2;
    
    maxMutationArray = [Max_mut_factor * boundsArray(1:Nparam-2*Nstates,2) / 2; ...
        0.1 * ( boundsArray(Nparam-2*Nstates+1:end,2) - boundsArray(Nparam-2*Nstates+1:end,1)) ];    
    
    param_strings = {'t12','t21','t24','t42','t43','t34','t41','A1','A2','A3','A4','sig1','sig2','sig3','sig4'};
    % Note: the order of the strings needs to be the same as the variable
    % list so we can extract the parameters properly.
    
    % tij_strings = param_strings(1:(Nparam - Nstates));
    % A_strings = param_strings((Nparam - Nstates + 1):Nparam);
    
    %----------------------------------------------------------------------
    % Update model
    %----------------------------------------------------------------------
    if length(tijs) > (Nparam - 2*Nstates)
        tijs = tijs(1:(Nparam - 2*Nstates));
    end
    rates = 1 ./ tijs;
    
    [k12, k21, k24, k42, k43, k34, k41] = feval(@(x) x{:}, num2cell(reshape(rates,[1, length(rates)])));
    k14 = (k41 * k12 * k24)/(k42 * k21);
%     k41 = (k31 * k14 * k43)/ (k13 * k34);
%     k42 = (k24 * k43 * k31 * k12)/(k34 * k13 * k21);
    
    % [A1, A2, A3, A4] = feval(@(x) x{:}, num2cell(As));
    
    K = [-(k12 + k14),          k21,                 0,                k41;...
        k12,              -(k21 + k24),        0,                 k42;...
        0,                0,                -(k34),               k43;...
        k14,                k24,                k34,             -(k42 + k43 + k41)];
    
    
    
     elseif strcmp(model_name, 'pyramid41') == 1
     Nstates = 4;
    t23_bounds = [200e-6,50000e-3];  %Paramater #1 is high--> med    
%     t13_bounds = [1e-3,3000e-3];    %Paramater #2 is high --> low
%     t14_bounds = [100e-6,900e-3];
    t32_bounds = [200e-6,50000e-3];%Paramater #3 is med --> high
    t31_bounds = [200e-6,50000e-3];%Paramater #4 is med --> low
    t13_bounds = [200e-6,50000e-3];
%     t31_bounds = [50e-6,1000e-3];  %Paramater #5 is low --> Medium = [1e-3,10e-3];%Paramater #5 %t32 is low --> Medium
    t14_bounds = [200e-6,50000e-3];
    t41_bounds = [200e-6,50000e-3];
    t12_bounds = [200e-6,50000e-3];
%     t21_bounds = [200e-6,50000e-3];
    % *t32 wll be determined by the other rates
    
% IN THIS SCHEME, 2 AND 3 ARE THE LOW AMP MOSTLY SYMM/ANTI WHILE 1 AND 4
% ARE THE HIGH AP STATES
    % A1_bounds = [0.65,0.85];%Paramater #6 % HIGH fret State
    A1_bounds = [0.04,0.055];%low amp
    A2_bounds = [0.08,0.13];%2nd lowest
    A3_bounds = [0.15,0.2];%2nd high amp
    A4_bounds = [0.25,0.32]; %1st high amp  
    
    %Bounds for standard deviation
     sigma_A1_bounds = [0.028,0.042]; %Paramater #8 %Low Amp state standard deviation
     sigma_A2_bounds = [0.05,0.07];%Paramater #9 %Med Amp state standard deviation
     sigma_A3_bounds = [0.065,0.09];%Paramater #10 %Med2 Amp state standard deviation
     sigma_A4_bounds = [0.09,0.14];%Paramater #10 %High Amp state standard deviation
    
    boundsArray = [t23_bounds; t32_bounds; t31_bounds; t13_bounds; t14_bounds; t41_bounds; t12_bounds; ...
        A1_bounds; A2_bounds; A3_bounds; A4_bounds; ...
        sigma_A1_bounds; sigma_A2_bounds; sigma_A3_bounds ; sigma_A4_bounds];
    
    Nparam = length(boundsArray);
    boundsArray_noAdj = boundsArray;
    
    [boundsArrayAdj,atEdges]= boundsAdjuster_v3(boundsArray,windowMag,paramAdjs,absLow,absHigh,Nparam,Nstates,fixedDecades);
    boundsArray=boundsArrayAdj;
    
    maxMutationCountsArray = zeros(1,Nparam);
    minMutationCountsArray = zeros(1,Nparam);
    
    
    % Maximum mutation value
    Max_mut_factor = .2;
    
    maxMutationArray = [Max_mut_factor * boundsArray(1:Nparam-2*Nstates,2) / 2; ...
        0.1 * ( boundsArray(Nparam-2*Nstates+1:end,2) - boundsArray(Nparam-2*Nstates+1:end,1)) ];    
    
    param_strings = {'t23','t32','t31','t13','t14','t41','t12','A1','A2','A3','A4','sig1','sig2','sig3','sig4'};
    % Note: the order of the strings needs to be the same as the variable
    % list so we can extract the parameters properly.
    
    % tij_strings = param_strings(1:(Nparam - Nstates));
    % A_strings = param_strings((Nparam - Nstates + 1):Nparam);
    
    %----------------------------------------------------------------------
    % Update model
    %----------------------------------------------------------------------
    if length(tijs) > (Nparam - 2*Nstates)
        tijs = tijs(1:(Nparam - 2*Nstates));
    end
    rates = 1 ./ tijs;
    
    [k23, k32, k31, k13, k14, k41, k12] = feval(@(x) x{:}, num2cell(reshape(rates,[1, length(rates)])));
    k21 = (k12 * k23 * k31)/(k13 * k32);
%     k41 = (k31 * k14 * k43)/ (k13 * k34);
%     k42 = (k24 * k43 * k31 * k12)/(k34 * k13 * k21);
    
    % [A1, A2, A3, A4] = feval(@(x) x{:}, num2cell(As));
    
    K = [-(k13 + k14 + k12),  k21,               k31,                k41;...
        k12,           -(k23 + k21),             k32,              0;...
        k13,               k23,           -(k31 + k32),          0;...
        k14,                 0,                0,               -(k41)];
    
    
    elseif strcmp(model_name, 'pyramid42') == 1
     Nstates = 4;
    t13_bounds = [200e-6,50000e-3];  %Paramater #1 is high--> med    
%     t13_bounds = [1e-3,3000e-3];    %Paramater #2 is high --> low
%     t14_bounds = [100e-6,900e-3];
    t31_bounds = [200e-6,50000e-3];%Paramater #3 is med --> high
    t23_bounds = [200e-6,50000e-3];%Paramater #4 is med --> low
    t32_bounds = [200e-6,50000e-3];
%     t31_bounds = [50e-6,1000e-3];  %Paramater #5 is low --> Medium = [1e-3,10e-3];%Paramater #5 %t32 is low --> Medium
    t24_bounds = [200e-6,50000e-3];
    t42_bounds = [200e-6,50000e-3];
%     t12_bounds = [200e-6,50000e-3];
    t21_bounds = [200e-6,50000e-3];
    % *t32 wll be determined by the other rates
    
% IN THIS SCHEME, 2 AND 3 ARE THE LOW AMP MOSTLY SYMM/ANTI WHILE 1 AND 4
% ARE THE HIGH AP STATES
    % A1_bounds = [0.65,0.85];%Paramater #6 % HIGH fret State
    A1_bounds = [0.04,0.055];%low amp
    A2_bounds = [0.08,0.13];%2nd lowest
    A3_bounds = [0.15,0.2];%2nd high amp
    A4_bounds = [0.25,0.32]; %1st high amp  
    
    %Bounds for standard deviation
     sigma_A1_bounds = [0.028,0.042]; %Paramater #8 %Low Amp state standard deviation
     sigma_A2_bounds = [0.05,0.07];%Paramater #9 %Med Amp state standard deviation
     sigma_A3_bounds = [0.065,0.09];%Paramater #10 %Med2 Amp state standard deviation
     sigma_A4_bounds = [0.09,0.14];%Paramater #10 %High Amp state standard deviation
    
    boundsArray = [t13_bounds; t31_bounds; t23_bounds; t32_bounds; t24_bounds; t42_bounds; t21_bounds; ...
        A1_bounds; A2_bounds; A3_bounds; A4_bounds; ...
        sigma_A1_bounds; sigma_A2_bounds; sigma_A3_bounds ; sigma_A4_bounds];
    
    Nparam = length(boundsArray);
    boundsArray_noAdj = boundsArray;
    
    [boundsArrayAdj,atEdges]= boundsAdjuster_v3(boundsArray,windowMag,paramAdjs,absLow,absHigh,Nparam,Nstates,fixedDecades);
    boundsArray=boundsArrayAdj;
    
    maxMutationCountsArray = zeros(1,Nparam);
    minMutationCountsArray = zeros(1,Nparam);
    
    
    % Maximum mutation value
    Max_mut_factor = .2;
    
    maxMutationArray = [Max_mut_factor * boundsArray(1:Nparam-2*Nstates,2) / 2; ...
        0.1 * ( boundsArray(Nparam-2*Nstates+1:end,2) - boundsArray(Nparam-2*Nstates+1:end,1)) ];    
    
    param_strings = {'t13','t31','t23','t32','t24','t21','t42','A1','A2','A3','A4','sig1','sig2','sig3','sig4'};
    % Note: the order of the strings needs to be the same as the variable
    % list so we can extract the parameters properly.
    
    % tij_strings = param_strings(1:(Nparam - Nstates));
    % A_strings = param_strings((Nparam - Nstates + 1):Nparam);
    
    %----------------------------------------------------------------------
    % Update model
    %----------------------------------------------------------------------
    if length(tijs) > (Nparam - 2*Nstates)
        tijs = tijs(1:(Nparam - 2*Nstates));
    end
    rates = 1 ./ tijs;
    
    [k13, k31, k23, k32, k24, k42, k21] = feval(@(x) x{:}, num2cell(reshape(rates,[1, length(rates)])));
      k12 = (k21 * k13 * k32)/(k23 * k31);
%     k41 = (k31 * k14 * k43)/ (k13 * k34);
%     k42 = (k24 * k43 * k31 * k12)/(k34 * k13 * k21);
    
    % [A1, A2, A3, A4] = feval(@(x) x{:}, num2cell(As));
    
    K = [-(k13 + k12),     k21,                k31,               0;...
        k12,          -(k23 + k24 + k21),       k32,               k42;...
        k13,             k23,              -(k31 + k32 ),       0;...
        0,               k24,                0,             -(k42)];
    
    
    elseif strcmp(model_name, 'pyramid43') == 1
     Nstates = 4;
    t12_bounds = [200e-6,50000e-3];  %Paramater #1 is high--> med    
%     t13_bounds = [1e-3,3000e-3];    %Paramater #2 is high --> low
%     t14_bounds = [100e-6,900e-3];
    t21_bounds = [200e-6,50000e-3];%Paramater #3 is med --> high
    t23_bounds = [200e-6,50000e-3];%Paramater #4 is med --> low
    t32_bounds = [200e-6,50000e-3];
%     t31_bounds = [50e-6,1000e-3];  %Paramater #5 is low --> Medium = [1e-3,10e-3];%Paramater #5 %t32 is low --> Medium
    t34_bounds = [200e-6,50000e-3];
    t43_bounds = [200e-6,50000e-3];
%     t13_bounds = [200e-6,50000e-3];
    t31_bounds = [200e-6,50000e-3];
    % *t32 wll be determined by the other rates
    
% IN THIS SCHEME, 2 AND 3 ARE THE LOW AMP MOSTLY SYMM/ANTI WHILE 1 AND 4
% ARE THE HIGH AP STATES
    % A1_bounds = [0.65,0.85];%Paramater #6 % HIGH fret State
    A1_bounds = [0.04,0.055];%low amp
    A2_bounds = [0.08,0.13];%2nd lowest
    A3_bounds = [0.15,0.2];%2nd high amp
    A4_bounds = [0.25,0.32]; %1st high amp  
    
    %Bounds for standard deviation
     sigma_A1_bounds = [0.028,0.042]; %Paramater #8 %Low Amp state standard deviation
     sigma_A2_bounds = [0.05,0.07];%Paramater #9 %Med Amp state standard deviation
     sigma_A3_bounds = [0.065,0.09];%Paramater #10 %Med2 Amp state standard deviation
     sigma_A4_bounds = [0.09,0.14];%Paramater #10 %High Amp state standard deviation
    
    boundsArray = [t12_bounds; t21_bounds; t23_bounds; t32_bounds; t34_bounds; t43_bounds; t31_bounds;...
        A1_bounds; A2_bounds; A3_bounds; A4_bounds; ...
        sigma_A1_bounds; sigma_A2_bounds; sigma_A3_bounds ; sigma_A4_bounds];
    
    Nparam = length(boundsArray);
    boundsArray_noAdj = boundsArray;
    
    [boundsArrayAdj,atEdges]= boundsAdjuster_v3(boundsArray,windowMag,paramAdjs,absLow,absHigh,Nparam,Nstates,fixedDecades);
    boundsArray=boundsArrayAdj;
    
    maxMutationCountsArray = zeros(1,Nparam);
    minMutationCountsArray = zeros(1,Nparam);
    
    
    % Maximum mutation value
    Max_mut_factor = .2;
    
    maxMutationArray = [Max_mut_factor * boundsArray(1:Nparam-2*Nstates,2) / 2; ...
        0.1 * ( boundsArray(Nparam-2*Nstates+1:end,2) - boundsArray(Nparam-2*Nstates+1:end,1)) ];    
    
    param_strings = {'t12','t21','t23','t32','t34','t43','t31','A1','A2','A3','A4','sig1','sig2','sig3','sig4'};
    % Note: the order of the strings needs to be the same as the variable
    % list so we can extract the parameters properly.
    
    % tij_strings = param_strings(1:(Nparam - Nstates));
    % A_strings = param_strings((Nparam - Nstates + 1):Nparam);
    
    %----------------------------------------------------------------------
    % Update model
    %----------------------------------------------------------------------
    if length(tijs) > (Nparam - 2*Nstates)            
        tijs = tijs(1:(Nparam - 2*Nstates));
    end
    rates = 1 ./ tijs;
    
    [k12, k21, k23, k32, k34, k43, k31] = feval(@(x) x{:}, num2cell(reshape(rates,[1, length(rates)])));
    k13 = (k31 * k12 * k23)/(k32 * k21);
%     k41 = (k31 * k14 * k43)/ (k13 * k34);
%     k42 = (k24 * k43 * k31 * k12)/(k34 * k13 * k21);
    
    % [A1, A2, A3, A4] = feval(@(x) x{:}, num2cell(As));
    
    K = [-(k12 + k13),          k21,                 k31,                0;...
        k12,              -(k21 + k23),        k32,               0;...
        k13,                k23,              -(k32 + k34 + k31),       k43;...
        0,                 0,                 k34,             -(k43)];
 

elseif strcmp(model_name, 'outer_loop1234') == 1
     Nstates = 4;
    t12_bounds = [200e-6,50000e-3];  %Paramater #1 is high--> med
%     t13_bounds = [200e-6,50000e-3];    %Paramater #2 is high --> low
%     t14_bounds = [200e-6,50000e-3];
    t21_bounds = [200e-6,50000e-3];%Paramater #3 is med --> high
    t23_bounds = [200e-6,50000e-3];%Paramater #4 is med --> low
    t32_bounds = [200e-6,50000e-3];
%     t31_bounds = [200e-6,50000e-3];  %Paramater #5 is low --> Medium = [1e-3,10e-3];%Paramater #5 %t32 is low --> Medium
    t34_bounds = [200e-6,50000e-3];
    t43_bounds = [200e-6,50000e-3];
    t41_bounds = [200e-6,50000e-3];

    % *t32 wll be determined by the other rates

% IN THIS SCHEME, 2 AND 3 ARE THE LOW AMP MOSTLY SYMM/ANTI WHILE 1 AND 4
% ARE THE HIGH AP STATES
    % A1_bounds = [0.65,0.85];%Paramater #6 % HIGH fret State
    A1_bounds = [0.001,0.055];%low amp
    A2_bounds = [0.08,0.13];%2nd lowest
    A3_bounds = [0.13,0.2];%2nd high amp
    A4_bounds = [0.25,0.4]; %1st high amp

    %Bounds for standard deviation
     sigma_A1_bounds = [0.028,0.042]; %Paramater #8 %Low Amp state standard deviation
     sigma_A2_bounds = [0.05,0.07];%Paramater #9 %Med Amp state standard deviation
     sigma_A3_bounds = [0.065,0.09];%Paramater #10 %Med2 Amp state standard deviation
     sigma_A4_bounds = [0.09,0.14];%Paramater #10 %High Amp state standard deviation

    boundsArray = [t12_bounds; t21_bounds; t23_bounds; t32_bounds; t34_bounds; t43_bounds; t41_bounds;...
        A1_bounds; A2_bounds; A3_bounds; A4_bounds; ...
        sigma_A1_bounds; sigma_A2_bounds; sigma_A3_bounds ; sigma_A4_bounds];

    Nparam = length(boundsArray);
    boundsArray_noAdj = boundsArray;

    [boundsArrayAdj,atEdges]= boundsAdjuster_v3(boundsArray,windowMag,paramAdjs,absLow,absHigh,Nparam,Nstates,fixedDecades);
    boundsArray=boundsArrayAdj;

    maxMutationCountsArray = zeros(1,Nparam);
    minMutationCountsArray = zeros(1,Nparam);


    % Maximum mutation value
    Max_mut_factor = .2;

    maxMutationArray = [Max_mut_factor * boundsArray(1:Nparam-2*Nstates,2) / 2; ...
        0.1 * ( boundsArray(Nparam-2*Nstates+1:end,2) - boundsArray(Nparam-2*Nstates+1:end,1)) ];

    param_strings = {'t12','t21','t23','t32','t34','t43','t41','A1','A2','A3','A4','sig1','sig2','sig3','sig4'};
    % Note: the order of the strings needs to be the same as the variable
    % list so we can extract the parameters properly.

    % tij_strings = param_strings(1:(Nparam - Nstates));
    % A_strings = param_strings((Nparam - Nstates + 1):Nparam);

    %----------------------------------------------------------------------
    % Update model
    %----------------------------------------------------------------------
    if length(tijs) > (Nparam - 2*Nstates)
        tijs = tijs(1:(Nparam - 2*Nstates));
    end
    rates = 1 ./ tijs;

    [k12, k21, k23, k32, k34, k43, k41] = feval(@(x) x{:}, num2cell(reshape(rates,[1, length(rates)])));
%     k32 = (k12 * k23 * k31)/(k13 * k21);
    k14 = (k12 * k23 * k34 * k41)/ (k43 * k32 * k21);
%     k42 = (k24 * k43 * k31 * k12)/(k34 * k13 * k21);

    % [A1, A2, A3, A4] = feval(@(x) x{:}, num2cell(As));

    K = [-(k12 + k14),          k21,                 0,                k41;...
        k12,              -(k21 + k23),        k32,               0;...
        0,                k23,              -(k32 + k34),       k43;...
        k14,                 0,                 k34,             -(k43 + k41)];
    
    elseif strcmp(model_name, 'outer_loop1243') == 1
     Nstates = 4;
    t12_bounds = [200e-6,50000e-3];  %Paramater #1 is high--> med    
%     t13_bounds = [1e-3,3000e-3];    %Paramater #2 is high --> low
%     t14_bounds = [100e-6,900e-3];
    t21_bounds = [200e-6,50000e-3];%Paramater #3 is med --> high
    t24_bounds = [200e-6,50000e-3];%Paramater #4 is med --> low
    t42_bounds = [200e-6,50000e-3];
%     t31_bounds = [50e-6,1000e-3];  %Paramater #5 is low --> Medium = [1e-3,10e-3];%Paramater #5 %t32 is low --> Medium
    t43_bounds = [200e-6,50000e-3];
    t34_bounds = [200e-6,50000e-3];
    t13_bounds = [200e-6,50000e-3]; 
    % *t32 wll be determined by the other rates
    
% IN THIS SCHEME, 2 AND 3 ARE THE LOW AMP MOSTLY SYMM/ANTI WHILE 1 AND 4
% ARE THE HIGH AP STATES
    % A1_bounds = [0.65,0.85];%Paramater #6 % HIGH fret State
    A1_bounds = [0.001,0.055];%low amp
    A2_bounds = [0.08,0.13];%2nd lowest
    A3_bounds = [0.13,0.2];%2nd high amp
    A4_bounds = [0.25,0.4]; %1st high amp  
    
    %Bounds for standard deviation
     sigma_A1_bounds = [0.028,0.042]; %Paramater #8 %Low Amp state standard deviation
     sigma_A2_bounds = [0.05,0.07];%Paramater #9 %Med Amp state standard deviation
     sigma_A3_bounds = [0.065,0.09];%Paramater #10 %Med2 Amp state standard deviation
     sigma_A4_bounds = [0.09,0.14];%Paramater #10 %High Amp state standard deviation
    
    boundsArray = [t12_bounds; t21_bounds; t24_bounds; t42_bounds; t43_bounds; t34_bounds; t13_bounds; ...
        A1_bounds; A2_bounds; A3_bounds; A4_bounds; ...
        sigma_A1_bounds; sigma_A2_bounds; sigma_A3_bounds ; sigma_A4_bounds];
    
    Nparam = length(boundsArray);
    boundsArray_noAdj = boundsArray;
    
    [boundsArrayAdj,atEdges]= boundsAdjuster_v3(boundsArray,windowMag,paramAdjs,absLow,absHigh,Nparam,Nstates,fixedDecades);
    boundsArray=boundsArrayAdj;
    
    maxMutationCountsArray = zeros(1,Nparam);
    minMutationCountsArray = zeros(1,Nparam);
    
    
    % Maximum mutation value
    Max_mut_factor = .2;
    
    maxMutationArray = [Max_mut_factor * boundsArray(1:Nparam-2*Nstates,2) / 2; ...
        0.1 * ( boundsArray(Nparam-2*Nstates+1:end,2) - boundsArray(Nparam-2*Nstates+1:end,1)) ];    
    
    param_strings = {'t12','t21','t24','t42','t43','t34','t13','A1','A2','A3','A4','sig1','sig2','sig3','sig4'};
    % Note: the order of the strings needs to be the same as the variable
    % list so we can extract the parameters properly.
    
    % tij_strings = param_strings(1:(Nparam - Nstates));
    % A_strings = param_strings((Nparam - Nstates + 1):Nparam);
    
    %----------------------------------------------------------------------
    % Update model
    %----------------------------------------------------------------------
    if length(tijs) > (Nparam - 2*Nstates)
        tijs = tijs(1:(Nparam - 2*Nstates));
    end
    rates = 1 ./ tijs;
    
    [k12, k21, k24, k42, k43, k34, k13] = feval(@(x) x{:}, num2cell(reshape(rates,[1, length(rates)])));
    k31 = (k13 * k34 * k42 *k21)/(k12 * k24 *k43);
%     k41 = (k31 * k14 * k43)/ (k13 * k34);
%     k42 = (k24 * k43 * k31 * k12)/(k34 * k13 * k21);
    
    % [A1, A2, A3, A4] = feval(@(x) x{:}, num2cell(As));
    
    K = [-(k12 + k13),          k21,           k31,                0;...
        k12,              -(k21 + k24),        0,                 k42;...
        k13,                0,                -(k34 +k31),               k43;...
        0,                k24,                k34,             -(k42 + k43)];
    
    elseif strcmp(model_name, 'fully_connected_minus31') == 1 % this is the fully connected model
    Nstates = 4;
    t12_bounds = [200e-6,50000e-3];  %Paramater #1 is high--> med
    t21_bounds = [200e-6,50000e-3];    %Paramater #2 is high --> low
%     t13_bounds = [5e-6,700e-3];
%     t31_bounds = [10e-6,100e-3];%Paramater #3 is med --> high
%     t14_bounds = [10e-6,100e-3];%Paramater #4 is med --> low
    t41_bounds = [200e-6,50000e-3];
    t23_bounds = [200e-6,50000e-3];  %Paramater #5 is low --> Medium = [1e-3,10e-3];%Paramater #5 %t32 is low --> Medium
%     t32_bounds = [100e-6,1000e-3];
    t24_bounds = [200e-6,50000e-3];
    t42_bounds = [200e-6,50000e-3];
    t34_bounds = [200e-6,50000e-3];
    t43_bounds = [200e-6,50000e-3];
   

    % *t32 wll be determined by the other rates
   
    % A1_bounds = [0.65,0.85];%Paramater #6 % HIGH fret State
    A1_bounds = [0.001,0.055];%low amp
    A2_bounds = [0.08,0.13];%2nd lowest
    A3_bounds = [0.13,0.2];%2nd high amp
    A4_bounds = [0.25,0.4]; %1st high amp  
   
    %Bounds for standard deviation
     sigma_A1_bounds = [0.028,0.042]; %Paramater #8 %Low Amp state standard deviation
     sigma_A2_bounds = [0.05,0.07];%Paramater #9 %Med Amp state standard deviation
     sigma_A3_bounds = [0.065,0.09];%Paramater #10 %Med2 Amp state standard deviation
     sigma_A4_bounds = [0.09,0.14];%Paramater #10 %High Amp state standard deviation
   
    boundsArray = [t12_bounds; t21_bounds; t41_bounds; t23_bounds; ...
         t24_bounds; t42_bounds; t34_bounds; t43_bounds; ...
        A1_bounds; A2_bounds; A3_bounds; A4_bounds; ...
        sigma_A1_bounds; sigma_A2_bounds; sigma_A3_bounds ; sigma_A4_bounds];
   
    Nparam = length(boundsArray);
    boundsArray_noAdj = boundsArray;
   
    [boundsArrayAdj,atEdges]= boundsAdjuster_v3(boundsArray,windowMag,paramAdjs,absLow,absHigh,Nparam,Nstates,fixedDecades);
    boundsArray=boundsArrayAdj;
   
    maxMutationCountsArray = zeros(1,Nparam);
    minMutationCountsArray = zeros(1,Nparam);
   
   
    % Maximum mutation value
    Max_mut_factor = .2;
   
    maxMutationArray = [Max_mut_factor * boundsArray(1:Nparam-2*Nstates,2) / 2; ...
        0.1 * ( boundsArray(Nparam-2*Nstates+1:end,2) - boundsArray(Nparam-2*Nstates+1:end,1)) ];
   
    param_strings = {'t12','t21','t41','t23','t24','t42','t34','t43','A1','A2','A3','A4','sig1','sig2','sig3','sig4'};
    % Note: the order of the strings needs to be the same as the variable
    % list so we can extract the parameters properly.
   
    % tij_strings = param_strings(1:(Nparam - Nstates));
    % A_strings = param_strings((Nparam - Nstates + 1):Nparam);
   
    %----------------------------------------------------------------------
    % Update model
    %----------------------------------------------------------------------
    if length(tijs) > (Nparam - 2*Nstates)
        tijs = tijs(1:(Nparam - 2*Nstates));
    end
    rates = 1 ./ tijs;
   
    [k12, k21,  k41, k23, k24, k42, k34, k43] = feval(@(x) x{:}, num2cell(reshape(rates,[1, length(rates)])));
    k32 = (k34 * k42 * k23)/(k24 * k43);
    k14 = (k41 * k12 * k24)/ (k21 * k42);
%     k42 = (k24 * k43 * k31 * k12)/(k34 * k13 * k21);
   
    % [A1, A2, A3, A4] = feval(@(x) x{:}, num2cell(As));
   
   
    K = [-(k12 + k14),     k21,               0,               k41;...
        k12,              -(k21 + k23 + k24), k32,               k42;...
        0,               k23,              -( k32 + k34),       k43;...
        k14,               k24,                k34,             -(k41 + k42 + k43)];
    
    elseif strcmp(model_name, 'fully_connected_minus21') == 1 % this is the fully connected model
    Nstates = 4;
    %t12_bounds = [50e-6,1000e-3];  %Paramater #1 is high--> med
    %t21_bounds = [10e-3,1000e-3];    %Paramater #2 is high --> low
    t13_bounds = [200e-6,50000e-3];
    t31_bounds = [200e-6,50000e-3];%Paramater #3 is med --> high
    t14_bounds = [200e-6,50000e-3];%Paramater #4 is med --> low
    t41_bounds = [200e-6,50000e-3];
    %t23_bounds = [80e-6,500e-3];  %Paramater #5 is low --> Medium = [1e-3,10e-3];%Paramater #5 %t32 is low --> Medium
    t32_bounds = [200e-6,50000e-3];
    t24_bounds = [200e-6,50000e-3];
    t42_bounds = [200e-6,50000e-3];
    %t34_bounds = [200e-6,50000e-3];
    t43_bounds = [200e-6,50000e-3];
    

    % *t32 wll be determined by the other rates
    
    % A1_bounds = [0.65,0.85];%Paramater #6 % HIGH fret State
    A1_bounds = [0.001,0.055];%low amp
    A2_bounds = [0.08,0.13];%2nd lowest
    A3_bounds = [0.13,0.2];%2nd high amp
    A4_bounds = [0.25,0.4]; %1st high amp  
    
    %Bounds for standard deviation
     sigma_A1_bounds = [0.028,0.042]; %Paramater #8 %Low Amp state standard deviation
     sigma_A2_bounds = [0.05,0.07];%Paramater #9 %Med Amp state standard deviation
     sigma_A3_bounds = [0.065,0.09];%Paramater #10 %Med2 Amp state standard deviation
     sigma_A4_bounds = [0.09,0.14];%Paramater #10 %High Amp state standard deviation
    
    boundsArray = [t13_bounds; t31_bounds; t14_bounds; t41_bounds; ...
         t32_bounds; t24_bounds; t42_bounds; t43_bounds; ...
        A1_bounds; A2_bounds; A3_bounds; A4_bounds; ...
        sigma_A1_bounds; sigma_A2_bounds; sigma_A3_bounds ; sigma_A4_bounds];
    
    Nparam = length(boundsArray);
    boundsArray_noAdj = boundsArray;
    
    [boundsArrayAdj,atEdges]= boundsAdjuster_v3(boundsArray,windowMag,paramAdjs,absLow,absHigh,Nparam,Nstates,fixedDecades);
    boundsArray=boundsArrayAdj;
    
    maxMutationCountsArray = zeros(1,Nparam);
    minMutationCountsArray = zeros(1,Nparam);
    
    
    % Maximum mutation value
    Max_mut_factor = .2;
    
    maxMutationArray = [Max_mut_factor * boundsArray(1:Nparam-2*Nstates,2) / 2; ...
        0.1 * ( boundsArray(Nparam-2*Nstates+1:end,2) - boundsArray(Nparam-2*Nstates+1:end,1)) ]; 
    
    param_strings = {'t13','t31','t14','t41','t32','t24','t42','t43','A1','A2','A3','A4','sig1','sig2','sig3','sig4'};
    % Note: the order of the strings needs to be the same as the variable
    % list so we can extract the parameters properly.
    
    % tij_strings = param_strings(1:(Nparam - Nstates));
    % A_strings = param_strings((Nparam - Nstates + 1):Nparam);
    
    %----------------------------------------------------------------------
    % Update model
    %----------------------------------------------------------------------
    if length(tijs) > (Nparam - 2*Nstates)
        tijs = tijs(1:(Nparam - 2*Nstates));
    end
    rates = 1 ./ tijs;
    
    [k13, k31, k14, k41, k32, k24, k42, k43] = feval(@(x) x{:}, num2cell(reshape(rates,[1, length(rates)])));
    k34 = (k14 * k43 * k31)/(k13 * k41);
    k23 = (k43 * k32 * k24)/ (k34 * k42);
%     k42 = (k24 * k43 * k31 * k12)/(k34 * k13 * k21);
    
    % [A1, A2, A3, A4] = feval(@(x) x{:}, num2cell(As));
    
    
    K = [-( k13 + k14), 0,                  k31,                  k41;...
        0,              -( k23 + k24),      k32,                  k42;...
        k13,               k23,           -(k31 + k32 + k34),     k43;...
        k14,               k24,             k34,                -(k41 + k42 + k43)];
 
    strcmp(model_name, 'outer_loop123') == 1
    %----------------------------------------------------------------------
    % Initialzing parameters
    %----------------------------------------------------------------------
    Nstates = 3;
    t12_bounds = [20e-6,70000e-3];  %Paramater #1 is high--> med
    t13_bounds = [20e-6,70000e-3];    %Paramater #2 is high --> low
    t21_bounds = [20e-6,70000e-3];%Paramater #3 is med --> high
    t23_bounds = [20e-6,70000e-3];%Paramater #4 is med --> low
    t31_bounds = [20e-6,70000e-3];  %Paramater #5 is low --> Medium = [1e-3,10e-3];%Paramater #5 %t32 is low --> Medium
    % *t32 wll be determined by the other rates
    
%     % A1_bounds = [0.65,0.85];%Paramater #6 % HIGH fret State
%     A1_bounds = [0.04,0.08];%Paramater #6 % HIGH fret State
%     A2_bounds = [0.08,0.25];%Paramater #7 % Med FRET state
%     A3_bounds = [0.25,0.4];%Paramater #8 %Low FRET state
    
    A1_bounds = [0.001,0.055];%low amp
    A2_bounds = [0.06,0.15];%2nd lowest
    A3_bounds = [0.15,0.3];%2nd high amp
%     A4_bounds = [0.25,0.4]; %1st high amp  
    
    %Bounds for standard deviation
     sigma_A1_bounds = [0.028,0.07]; %Paramater #8 %Low Amp state standard deviation
     sigma_A2_bounds = [0.03,0.07];%Paramater #9 %Med Amp state standard deviation
     sigma_A3_bounds = [0.03,0.09];%Paramater #10 %Med2 Amp state standard deviation
%      sigma_A4_bounds = [0.09,0.14];%Paramater #10 %High Amp state standard deviation
    
      
    boundsArray = [t12_bounds; t13_bounds; t21_bounds; t23_bounds; ...
         t31_bounds; A1_bounds; A2_bounds; A3_bounds;  ...
        sigma_A1_bounds; sigma_A2_bounds; sigma_A3_bounds ];
    
    Nparam = length(boundsArray);
    boundsArray_noAdj = boundsArray;
    
    [boundsArrayAdj,atEdges]= boundsAdjuster_v3(boundsArray,windowMag,paramAdjs,absLow,absHigh,Nparam,Nstates,fixedDecades);
    boundsArray=boundsArrayAdj;
    
    maxMutationCountsArray = zeros(1,Nparam);
    minMutationCountsArray = zeros(1,Nparam);
    
    
    % Maximum mutation value
    Max_mut_factor = .2;
    
    maxMutationArray = [Max_mut_factor * boundsArray(1:Nparam-2*Nstates,2) / 2; ...
        0.1 * ( boundsArray(Nparam-2*Nstates+1:end,2) - boundsArray(Nparam-2*Nstates+1:end,1)) ]; 
    
    param_strings = {'t12','t13','t21','t23','t31','A1','A2','A3','sig1','sig2','sig3'};
    % Note: the order of the strings needs to be the same as the variable
    % list so we can extract the parameters properly.
    
    % tij_strings = param_strings(1:(Nparam - Nstates));
    % A_strings = param_strings((Nparam - Nstates + 1):Nparam);
    
    %----------------------------------------------------------------------
    % Update model
    %----------------------------------------------------------------------
    if length(tijs) > (Nparam - 2*Nstates)
        tijs = tijs(1:(Nparam - 2*Nstates));
    end
    
%     param_strings = {'t12','t13','t21','t23','t31','A1','A2','A3'};
    
    %     disp(['model_builder shape tijs: ',num2str(size(tijs))])
    rates = 1 ./ tijs;
    %     disp(['model_builder shape rates: ',num2str(size(rates))])
    [k12, k13, k21, k23, k31] = feval(@(x) x{:}, num2cell(reshape(rates,[1, length(rates)])));
    k32 = (k12 * k23 * k31)/(k13 * k21);
    
    K = [-(k12 + k13), k21,           k31;...
        k12,          -(k21 + k23),   k32;...
        k13,           k23,          -(k31 + k32)];
  
    
    strcmp(model_name, 'linear123') == 1
    %----------------------------------------------------------------------
    % Initialzing parameters
    %----------------------------------------------------------------------
    Nstates = 3;
    t12_bounds = [20e-6,70000e-3];  %Paramater #1 is high--> med
%     t13_bounds = [20e-6,70000e-3];    %Paramater #2 is high --> low
    t21_bounds = [20e-6,70000e-3];%Paramater #3 is med --> high
    t23_bounds = [20e-6,70000e-3];%Paramater #4 is med --> low
    t32_bounds = [20e-6,70000e-3];  %Paramater #5 is low --> Medium = [1e-3,10e-3];%Paramater #5 %t32 is low --> Medium
    % *t32 wll be determined by the other rates
    
%     % A1_bounds = [0.65,0.85];%Paramater #6 % HIGH fret State
%     A1_bounds = [0.04,0.08];%Paramater #6 % HIGH fret State
%     A2_bounds = [0.08,0.25];%Paramater #7 % Med FRET state
%     A3_bounds = [0.25,0.4];%Paramater #8 %Low FRET state
    
    A1_bounds = [0.001,0.055];%low amp
    A2_bounds = [0.06,0.15];%2nd lowest
    A3_bounds = [0.15,0.3];%2nd high amp
%     A4_bounds = [0.25,0.4]; %1st high amp  
    
    %Bounds for standard deviation
     sigma_A1_bounds = [0.028,0.07]; %Paramater #8 %Low Amp state standard deviation
     sigma_A2_bounds = [0.03,0.07];%Paramater #9 %Med Amp state standard deviation
     sigma_A3_bounds = [0.03,0.09];%Paramater #10 %Med2 Amp state standard deviation
%      sigma_A4_bounds = [0.09,0.14];%Paramater #10 %High Amp state standard deviation
    
      
    boundsArray = [t12_bounds; t21_bounds; t23_bounds; ...
         t32_bounds; A1_bounds; A2_bounds; A3_bounds;  ...
        sigma_A1_bounds; sigma_A2_bounds; sigma_A3_bounds ];
    
    Nparam = length(boundsArray);
    boundsArray_noAdj = boundsArray;
    
    [boundsArrayAdj,atEdges]= boundsAdjuster_v3(boundsArray,windowMag,paramAdjs,absLow,absHigh,Nparam,Nstates,fixedDecades);
    boundsArray=boundsArrayAdj;
    
    maxMutationCountsArray = zeros(1,Nparam);
    minMutationCountsArray = zeros(1,Nparam);
    
    
    % Maximum mutation value
    Max_mut_factor = .2;
    
    maxMutationArray = [Max_mut_factor * boundsArray(1:Nparam-2*Nstates,2) / 2; ...
        0.1 * ( boundsArray(Nparam-2*Nstates+1:end,2) - boundsArray(Nparam-2*Nstates+1:end,1)) ]; 
    
    param_strings = {'t12','t21','t23','t32','A1','A2','A3','sig1','sig2','sig3'};
    % Note: the order of the strings needs to be the same as the variable
    % list so we can extract the parameters properly.
    
    % tij_strings = param_strings(1:(Nparam - Nstates));
    % A_strings = param_strings((Nparam - Nstates + 1):Nparam);
    
    %----------------------------------------------------------------------
    % Update model
    %----------------------------------------------------------------------
    if length(tijs) > (Nparam - 2*Nstates)
        tijs = tijs(1:(Nparam - 2*Nstates));
    end
    
%     param_strings = {'t12','t13','t21','t23','t31','A1','A2','A3'};
    
    %     disp(['model_builder shape tijs: ',num2str(size(tijs))])
    rates = 1 ./ tijs;
    %     disp(['model_builder shape rates: ',num2str(size(rates))])
    [k12, k21, k23, k32] = feval(@(x) x{:}, num2cell(reshape(rates,[1, length(rates)])));
%     k32 = (k12 * k23 * k31)/(k13 * k21);
    
    K = [-(k12 ),       k21,           0;...
        k12,          -(k21 + k23),   k32;...
        0,             k23,           -(k32)];
  
    strcmp(model_name, 'linear231') == 1
    %----------------------------------------------------------------------
    % Initializing parameters
    %----------------------------------------------------------------------
    Nstates = 3;
    t23_bounds = [20e-6,70000e-3];  %Paramater #1 is high--> med
    t32_bounds = [20e-6,70000e-3];    %Paramater #2 is high --> low
    t31_bounds = [20e-6,70000e-3];%Paramater #3 is med --> high
    t13_bounds = [20e-6,70000e-3];%Paramater #4 is med --> low
%     t31_bounds = [20e-6,70000e-3];  %Paramater #5 is low --> Medium = [1e-3,10e-3];%Paramater #5 %t32 is low --> Medium
    % *t32 wll be determined by the other rates
    
%     % A1_bounds = [0.65,0.85];%Paramater #6 % HIGH fret State
%     A1_bounds = [0.04,0.08];%Paramater #6 % HIGH fret State
%     A2_bounds = [0.08,0.25];%Paramater #7 % Med FRET state
%     A3_bounds = [0.25,0.4];%Paramater #8 %Low FRET state
    
    A1_bounds = [0.001,0.055];%low amp
    A2_bounds = [0.06,0.15];%2nd lowest
    A3_bounds = [0.15,0.3];%2nd high amp
%     A4_bounds = [0.25,0.4]; %1st high amp  
    
    %Bounds for standard deviation
     sigma_A1_bounds = [0.028,0.07]; %Paramater #8 %Low Amp state standard deviation
     sigma_A2_bounds = [0.03,0.07];%Paramater #9 %Med Amp state standard deviation
     sigma_A3_bounds = [0.03,0.09];%Paramater #10 %Med2 Amp state standard deviation
%      sigma_A4_bounds = [0.09,0.14];%Paramater #10 %High Amp state standard deviation
    
      
    boundsArray = [t23_bounds; t32_bounds; t31_bounds; t13_bounds; ...
          A1_bounds; A2_bounds; A3_bounds;  ...
        sigma_A1_bounds; sigma_A2_bounds; sigma_A3_bounds ];
    
    Nparam = length(boundsArray);
    boundsArray_noAdj = boundsArray;
    
    [boundsArrayAdj,atEdges]= boundsAdjuster_v3(boundsArray,windowMag,paramAdjs,absLow,absHigh,Nparam,Nstates,fixedDecades);
    boundsArray=boundsArrayAdj;
    
    maxMutationCountsArray = zeros(1,Nparam);
    minMutationCountsArray = zeros(1,Nparam);
    
    
    % Maximum mutation value
    Max_mut_factor = .2;
    
    maxMutationArray = [Max_mut_factor * boundsArray(1:Nparam-2*Nstates,2) / 2; ...
        0.1 * ( boundsArray(Nparam-2*Nstates+1:end,2) - boundsArray(Nparam-2*Nstates+1:end,1)) ]; 
    
    param_strings = {'t23','t32','t31','t13','A1','A2','A3','sig1','sig2','sig3'};
    % Note: the order of the strings needs to be the same as the variable
    % list so we can extract the parameters properly.
    
    % tij_strings = param_strings(1:(Nparam - Nstates));
    % A_strings = param_strings((Nparam - Nstates + 1):Nparam);
    
    %----------------------------------------------------------------------
    % Update model
    %----------------------------------------------------------------------
    if length(tijs) > (Nparam - 2*Nstates)
        tijs = tijs(1:(Nparam - 2*Nstates));
    end
    
%     param_strings = {'t12','t13','t21','t23','t31','A1','A2','A3'};
    
    %     disp(['model_builder shape tijs: ',num2str(size(tijs))])
    rates = 1 ./ tijs;
    %     disp(['model_builder shape rates: ',num2str(size(rates))])
    [k23, k32, k31, k13] = feval(@(x) x{:}, num2cell(reshape(rates,[1, length(rates)])));
%     k32 = (k12 * k23 * k31)/(k13 * k21);
    
    K = [-(k13),        0,           k31;...
        0,          -(k23),          k32;...
        k13,           k23,          -(k31 + k32)];
    
    strcmp(model_name, 'linear312') == 1
    %----------------------------------------------------------------------
    % Initialzing parameters
    %----------------------------------------------------------------------
    Nstates = 3;
    t31_bounds = [20e-6,70000e-3];  %Paramater #1 is high--> med
    t13_bounds = [20e-6,70000e-3];    %Paramater #2 is high --> low
    t12_bounds = [20e-6,70000e-3];%Paramater #3 is med --> high
    t21_bounds = [20e-6,70000e-3];%Paramater #4 is med --> low
%     t31_bounds = [20e-6,70000e-3];  %Paramater #5 is low --> Medium = [1e-3,10e-3];%Paramater #5 %t32 is low --> Medium
    % *t32 wll be determined by the other rates
    
%     % A1_bounds = [0.65,0.85];%Paramater #6 % HIGH fret State
%     A1_bounds = [0.04,0.08];%Paramater #6 % HIGH fret State
%     A2_bounds = [0.08,0.25];%Paramater #7 % Med FRET state
%     A3_bounds = [0.25,0.4];%Paramater #8 %Low FRET state
    
    A1_bounds = [0.001,0.055];%low amp
    A2_bounds = [0.06,0.15];%2nd lowest
    A3_bounds = [0.15,0.3];%2nd high amp
%     A4_bounds = [0.25,0.4]; %1st high amp  
    
    %Bounds for standard deviation
     sigma_A1_bounds = [0.028,0.07]; %Paramater #8 %Low Amp state standard deviation
     sigma_A2_bounds = [0.03,0.07];%Paramater #9 %Med Amp state standard deviation
     sigma_A3_bounds = [0.03,0.09];%Paramater #10 %Med2 Amp state standard deviation
%      sigma_A4_bounds = [0.09,0.14];%Paramater #10 %High Amp state standard deviation
    
      
    boundsArray = [t31_bounds; t13_bounds; t12_bounds; t21_bounds; ...
         A1_bounds; A2_bounds; A3_bounds;  ...
        sigma_A1_bounds; sigma_A2_bounds; sigma_A3_bounds ];
    
    Nparam = length(boundsArray);
    boundsArray_noAdj = boundsArray;
    
    [boundsArrayAdj,atEdges]= boundsAdjuster_v3(boundsArray,windowMag,paramAdjs,absLow,absHigh,Nparam,Nstates,fixedDecades);
    boundsArray=boundsArrayAdj;
    
    maxMutationCountsArray = zeros(1,Nparam);
    minMutationCountsArray = zeros(1,Nparam);
    
    
    % Maximum mutation value
    Max_mut_factor = .2;
    
    maxMutationArray = [Max_mut_factor * boundsArray(1:Nparam-2*Nstates,2) / 2; ...
        0.1 * ( boundsArray(Nparam-2*Nstates+1:end,2) - boundsArray(Nparam-2*Nstates+1:end,1)) ]; 
    
    param_strings = {'t31','t13','t12','t21','A1','A2','A3','sig1','sig2','sig3'};
    % Note: the order of the strings needs to be the same as the variable
    % list so we can extract the parameters properly.
    
    % tij_strings = param_strings(1:(Nparam - Nstates));
    % A_strings = param_strings((Nparam - Nstates + 1):Nparam);
    
    %----------------------------------------------------------------------
    % Update model
    %----------------------------------------------------------------------
    if length(tijs) > (Nparam - 2*Nstates)
        tijs = tijs(1:(Nparam - 2*Nstates));
    end
    
%     param_strings = {'t12','t13','t21','t23','t31','A1','A2','A3'};
    
    %     disp(['model_builder shape tijs: ',num2str(size(tijs))])
    rates = 1 ./ tijs;
    %     disp(['model_builder shape rates: ',num2str(size(rates))])
    [k31, k13, k12, k21] = feval(@(x) x{:}, num2cell(reshape(rates,[1, length(rates)])));
%     k32 = (k12 * k23 * k31)/(k13 * k21);
    
    K = [-(k12 + k13), k21,           k31;...
        k12,          -(k21),           0;...
        k13,           0,          -(k31)];    
    
end

global varyCtrlMode NctrlParams ctrlMutateAmpc2 ctrlMutateTau ctrlMutateAmpc4 C2_exp_y C4_tau2eq0_exp% added to vary the control timescale
    if varyCtrlMode == 1
    %   NctrlParams = 3;
    %   ctrlMutate = 0.1;
    %   c2ctrlAmp = 9.25e-3;
    %   c2ctrlTau = 2.41;
      c2ctrlAmp_bounds = c2ctrlAmp * [0.5, 1.5];
      c2ctrlTau_bounds = c2ctrlTau * [0.75, 1.25];
      c4ctrlAmp_bounds = c4ctrlAmp * [0.25, 1.75];
       boundsArray = [boundsArray; c2ctrlAmp_bounds; c2ctrlTau_bounds; c4ctrlAmp_bounds];
       Nparam = length(boundsArray);
       maxMutationArray = [maxMutationArray ; (ctrlMutateAmpc2*(boundsArray(Nparam-2, 2) - boundsArray(Nparam-2, 1)))];
       maxMutationArray = [maxMutationArray ; (ctrlMutateTau*(boundsArray(Nparam-1, 2) - boundsArray(Nparam-1, 1)))];
       maxMutationArray = [maxMutationArray ; (ctrlMutateAmpc4*(boundsArray(Nparam, 2) - boundsArray(Nparam, 1)))];
       param_strings = {param_strings{:}, 'c2ctrlAmp','c2ctrlTau', 'c4ctrlAmp'};
       
    end

end