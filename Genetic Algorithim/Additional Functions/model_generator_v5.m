%--------------------------------------------------------------------------
% AUTHORS: Claire & Jack
% CREATED: Summer 2022
% PURPOSE: generate cells containing the model information for N = 3,4,5,6
%   - So far we can account for all linear chains, loops, and some
%   branching
%--------------------------------------------------------------------------

function [model_lin,model_loop1,model_loopN,model_brch] = model_generator_v5(N,model_set)
% cc
% N = 6; 
saveMode = 0;

% select the number of states you want
% this currently works for 3,4,5,6
chains_all = perms(linspace(1,N,N));

%% Remove reduntant forward backward chains
chains = chains_all;
first_occur = [];
second_occur = [];
pairs = zeros(length(chains),2);
for i=1:length(chains)
    pairs(i,:) = [i, find((sum(flip(chains(i,:)) == chains,2) == N)==1,1)];
end

pairs_firsthalf = pairs(1:length(chains)/2,:);
pairs_secondhalf = pairs(length(chains)/2+1:length(chains),:);
% can't quite just take the pairs in either the first or second half
% because there are some redundant and some still missing.

% find the pairs with indicies either both > or < the halfway point, those
% are the redundant/missing ones

pairs2rm = find(sum(pairs_firsthalf <= length(chains)/2,2) == 2);
%     disp('pairs2rm: ')
%     disp(pairs(pairs2rm,:))
pairs2rm_unique = pairs2rm(1:length(pairs2rm)/2,:);
%     disp('pairs2rm_unique: ')
%     disp(pairs(pairs2rm_unique,:))

pairs2add = find(sum(pairs_secondhalf >= length(chains)/2,2) == 2);
%     disp('pairs2add: ')
%     disp(pairs(pairs2add+(length(pairs))/2,:))
pairs2add_unique = pairs2add(1:(length(pairs2add)/2),:);
%     disp('pairs2add_unique: ')
%     disp(pairs(pairs2add_unique+(length(pairs)/2),:))

pairs_unique = pairs(1:length(pairs)/2,:);
pairs_unique(pairs2rm_unique,:) = zeros(size(pairs_unique(pairs2rm_unique,:)));
pairs_unique =[pairs_unique; pairs_secondhalf(pairs2add_unique,:)];

pairs_unique(find(pairs_unique == [0,0])) = [];
pairs_unique = reshape(pairs_unique,[length(pairs)/2,2]);

chains = chains(pairs_unique(:,1),:);
% this has the possible linear chain arrays with redundancies removed
% Note: this is in state space (not index space)

%
% this part isn't general yet... not sure if we can make it that way. but
% not too bad so far.
if N == 3
    str_mat = [ "0","k21","k31";...
        "k12","0","k32";...
        "k13","k23","0"];
elseif N == 4
    str_mat = [ "0","k21","k31","k41";...
        "k12","0","k32","k42";...
        "k13","k23","0","k43";...
        "k14","k24","k34","0"];
elseif N == 5
    str_mat = [ "0","k21","k31","k41","k51";...
        "k12","0","k32","k42","k52";...
        "k13","k23","0","k43","k53";...
        "k14","k24","k34","0","k54";...
        "k15","k25","k35","k45","0"];
elseif N == 6
    str_mat = [ "0","k21","k31","k41","k51","k61";...
        "k12","0","k32","k42","k52","k62";...
        "k13","k23","0","k43","k53","k63";...
        "k14","k24","k34","0","k54","k64";...
        "k15","k25","k35","k45","0","k65";...
        "k16","k26","k36","k46","k56","0"];
end
str_list = reshape(str_mat, [N*N,1]);

% Using the linear chains calculate the model specific info
model_lin = cell(length(chains),6);
for j = 1:length(chains)
    bool_mat = zeros(N,N); % bool_mat is going to be the boolean masks that we build a library of to see if the calculated model already exists or not.
    chain = chains(j,:);
    %     disp(['chain: ',num2str(chain)])
    for i = 1:length(chain)-1
        bool_mat(chain(i), chain(i+1)) = 1;  % build the matrix with the rates that are in the chain
        bool_mat(chain(i+1), chain(i)) = 1; % put 1 and the i,j and j,i positions in the matrix for the sucesstive connections in the chain
    end
    %     disp(rate_mat)
    rates_str_woDB = str_list(find(bool_mat == 1));  % select the rate strings based on the occupation of the matrix
    rates_idx_woDB = find(bool_mat == 1);
    
    rate_idx = rates_idx_woDB;
    Nrates = length(rate_idx);
    rates = linspace(1,Nrates,Nrates);
    % logRandUniform...
    
    idxs = reshape(linspace(1,N*N,N*N),[N,N]);
    
    K = bool_mat;
    K(rates_idx_woDB) = rates;
    diagK = -sum(K,1);
    K(diag(idxs)) = diagK;
    if sum(K,1) ~= zeros(1,N)
        error(['ERROR: columns of K do not sum to zero!'])
    end
    
    
    model_name = strrep(['chain' strtrim(num2str(chain(1:end)))],' ','');
    %     model_lin{j,1} = chain;     % the linear connectivity
    %     model_lin{j,2} = bool_mat;  % the matrix which sets the rates that exist
    %     model_lin{j,3} = rates_str_woDB; % the rates that exist in the system
    %     model_lin{j,4} = rates_idx_woDB;
    %     model_lin{j,5} = K;
    %     model_lin{j,6} = model_name;
    
    model_lin{j,1} = model_name;
    model_lin{j,2} = K;
    model_lin{j,3} = rates_idx_woDB;
    model_lin{j,4} = rates_str_woDB;
    model_lin{j,5} = chain;
    model_lin{j,6} = bool_mat;
end


disp(['Found ' num2str(size(model_lin,1)) ' linear connectivities for N = ' num2str(N)])


%% Build K from an array of rates

% currently linear models do not exist in bool_stack... come back later if
% need a comprehensive bool_stack array.

% [chain, bool_mat,  rates_idx_woDB, rates_str_woDB, model_name] = model_lin{1,:};
%
% rate_idx = rates_idx_woDB;

% rate_idx = find(rate_mat == 1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% I moved this part above so the model_lin cell array will contain the K matricies. %%%
% Nrates = length(rate_idx);
% rates = linspace(1,Nrates,Nrates);
%
% idxs = reshape(linspace(1,N*N,N*N),[N,N]);
%
% K = bool_mat;
% K(rate_idx) = rates;
% diagK = -sum(K,1);
% K(diag(idxs)) = diagK;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% graph the connectivity
% figure(1);clf;
% plot(graph(bool_mat)) % way to visualize the extra connections


%% Try adding one loop
plotGraphMode = 1;


% rate_mat_loop = rate_mat;
model_loop1 = {};
bool_stack = [];

idxs = reshape(linspace(1,N*N,N*N),[N,N]); % general array of indices for the chosen chain length

for model_idx = 1:size(model_lin,1)
    %     [chain, bool_mat, rates_str_woDB, rates_idx_woDB] = model_lin{model_idx,:};
%     [model_name, K, rates_idx_woDB, rates_str_woDB, chain,bool_mat] = model_lin{model_idx,:};
%     rate_idx = rates_idx_woDB;
    
    for k = 1:length(chain)
        [model_name, K, rates_idx_woDB, rates_str_woDB, chain, bool_mat] = model_lin{model_idx,:};
        rate_idx = rates_idx_woDB;
        start_idx = k;
        cxns = linspace(1,length(chain),length(chain));
        allowed = cxns < start_idx-1;
        allowed = allowed + (cxns > start_idx+1);
        cxns = cxns .* allowed;
        cxns = cxns(nonzeros(cxns));
%         disp(['possible connections for start ' num2str(chain(start_idx)) ': ', num2str(chain(cxns))])
        cxns_all = chain(cxns);
        
        for i = 1:length(cxns)            
            [model_name, K, rates_idx_woDB, rates_str_woDB, chain, bool_mat] = model_lin{model_idx,:};
            bool_mat_loop = bool_mat;
            stop_idx = cxns(i);
            cxns_cur = [chain(start_idx), chain(stop_idx)];
            %         disp(['adding connection: ',num2str(chain(start_idx)), ' ', num2str(chain(stop_idx))])
            %         disp(['making loop of: ', num2str(chain(min(start_idx,stop_idx):max(start_idx,stop_idx)))])
            
            bool_mat_loop(chain(start_idx), chain(stop_idx)) = 1;
            bool_mat_loop(chain(stop_idx), chain(start_idx)) = 1;
            
            if isempty(bool_stack) == 1
                bool_stack = bool_mat_loop;
            elseif sum(sum(sum(bool_stack - bool_mat_loop ==0))== N*N) == 0
                bool_stack = cat(3,bool_stack, bool_mat_loop);
            else
                continue
            end
            
            rate_idx = [rate_idx; idxs(chain(start_idx),chain(stop_idx));idxs(chain(stop_idx),chain(start_idx))];
            Nrates = length(rate_idx);
            rates = linspace(1,Nrates,Nrates); %rand(1,Nrates);
            % logRandUniform will probably do this step... need to figure
            % out how to set bounds
            
            K = bool_mat_loop;
            
            K(rate_idx) = rates;
            
            
            bool_mat_loop(chain(stop_idx), chain(start_idx)) = 2; % use identify rate for detailed balance?
            rates_str_DB = str_list(find(bool_mat_loop == 2));
            rates_idx_DB = find(bool_mat_loop == 2);
            
            rates_str_woDB = str_list(find(bool_mat_loop == 1));  % select the rate strings based on the occupation of the matrix
            rates_idx_woDB = find(bool_mat_loop == 1);

            
            dbCell = {};
            loop_memb = chain(min(start_idx,stop_idx):max(start_idx,stop_idx));
            for_prod = 1;
            for_num = 0;
            for_prod_idx_arr = [];
            for i=1:length(loop_memb)
                if i == length(loop_memb) && bool_mat_loop(loop_memb(1),loop_memb(length(loop_memb))) ~= 2
                    for_prod = for_prod * K(loop_memb(1), loop_memb(length(loop_memb)));
                    for_num = for_num + 1;
                    for_prod_idx_arr = [for_prod_idx_arr, idxs(loop_memb(1), loop_memb(length(loop_memb)))];
                elseif i < length(loop_memb) && bool_mat_loop(loop_memb(i),loop_memb(i+1)) ~= 2
                    for_prod = for_prod * K(loop_memb(i), loop_memb(i+1));
                    for_num = for_num + 1;
                    for_prod_idx_arr = [for_prod_idx_arr, idxs(loop_memb(i), loop_memb(i+1))];
                end
            end
            rev_prod = 1;
            rev_num = 0;
            rev_prod_idx_arr = [];
            for i=1:length(loop_memb)
                if i == length(loop_memb) && bool_mat_loop(loop_memb(length(loop_memb)),loop_memb(1)) ~= 2
                    rev_prod = rev_prod * K(loop_memb(length(loop_memb)), loop_memb(1));
                    rev_num = rev_num + 1;
                    rev_prod_idx_arr = [rev_prod_idx_arr, idxs(loop_memb(length(loop_memb)),loop_memb(1))];
                elseif i < length(loop_memb) && bool_mat_loop(loop_memb(i+1),loop_memb(i)) ~= 2
                    rev_prod = rev_prod * K(loop_memb(i+1), loop_memb(i));
                    rev_num = rev_num + 1;
                    rev_prod_idx_arr = [rev_prod_idx_arr, idxs(loop_memb(i+1),loop_memb(i))];
                end
            end
            
            if isempty(dbCell) == 0
                dbCell =[dbCell; {for_prod_idx_arr, rev_prod_idx_arr}];
            else
                dbCell = {for_prod_idx_arr, rev_prod_idx_arr};
            end
            
            if for_num < rev_num
                K(rates_idx_DB(end))= rev_prod/for_prod;
            elseif for_num > rev_num
                K(rates_idx_DB(end))= for_prod/rev_prod;
            else
                error(['DB error: forward & backward loop are same size'])
            end

            
            diagK = -sum(K,1);
            K(diag(idxs)) = diagK;
            if sum(K,1) ~= zeros(1,N)
                error(['ERROR: columns of K do not sum to zero!'])
            end
            
            K = zeros(N, N);
            
            model_name = strrep(['chain' strtrim(num2str(chain(1:end))), '_cxn', num2str([chain(start_idx) chain(stop_idx)])],' ','');
            numLoops = 1;
            %             model_loop1 = [model_loop1; {K, rates_str_woDB, rates_str_DB, rates_idx_DB, rates_idx_woDB, bool_mat_loop, cxns_all, cxns_cur, model_name,numLoops, chain}];
            model_loop1 = [model_loop1; {model_name, K, rates_idx_woDB, rates_str_woDB, chain, bool_mat_loop, rates_idx_DB, rates_str_DB, dbCell, cxns_all, cxns_cur, numLoops}];
            
        end
        
    end
end

disp(['Found ' num2str(size(model_loop1,1)) ' single loop connectivities for N = ' num2str(N)])


%%
if strcmp(model_set, 'model_loopN') == 1
model_loopN = {};

for n = 1:size(model_loop1,1)

    state_idx = linspace(1,N,N);
    start_idxs = state_idx(sum(bool_mat_loop > 0,2)< N-1);
    
    for z = 1:length(start_idxs)
        [model_name_loop, K, rates_idx_woDB, rates_str_woDB, chain, bool_mat_loop, rates_idx_DB, rates_str_DB, dbCell, cxns_all, cxns_cur, numLoops] = model_loop1{n,:};

        
        % all in terms of indexes, not state labels
        start_idx = start_idxs(z);
        cxns = linspace(1,length(chain),length(chain));
        allowed = cxns < start_idx-1;
        allowed = allowed + (cxns > start_idx+1);
        cxns = cxns .* allowed;
        cxns = cxns(nonzeros(cxns));
        % convert back to state labels here
%         disp(['possible connections for start ' num2str(chain(start_idx)) ': ', num2str(chain(cxns))])
        cxns_all = chain(cxns);
        
        
        existing_cxns = state_idx(bool_mat_loop(:,chain(start_idx)) >0);
        %     start_idx = find(chain == cxns_cur(1)); %%%%%%%%%%%%%%%% need to choose different start indices %%%%%%%%%%%%%%%%
        %     existing_cxn = cxns_cur(2);

        rem_cxns = setdiff(cxns_all, existing_cxns); %cxns_all(cxns_all ~= existing_cxn);
        rates_idx = [rates_idx_woDB; rates_idx_DB];
%         disp(['remaining connections for start at ' num2str(chain(start_idx)) ': ', num2str(rem_cxns)])

        %while numLoops <= length(cxns_all)
        for j = 1:length(rem_cxns)
            identifyDB = j+2; % identify detailed balance condition by number of loops + 1
            stop_idx = find(chain == rem_cxns(j));
            if j == 1
                bool_mat_loopN = bool_mat_loop;
            end
            
            cxns_cur = [cxns_cur, rem_cxns(j)];
            
%             disp(['start: ' num2str(chain(start_idx))])
%             disp(['stop: ' num2str(chain(stop_idx))])
            bool_mat_loopN(chain(start_idx), chain(stop_idx)) = 1;
            bool_mat_loopN(chain(stop_idx), chain(start_idx)) = 1;
            
            bool_mat_loopN(chain(stop_idx), chain(start_idx)) = identifyDB; % use identify rate for detailed balance
            rates_str_DB = [rates_str_DB ;str_list(find(bool_mat_loopN == identifyDB))];
            rates_idx_DB = [rates_idx_DB; find(bool_mat_loopN == identifyDB)];
            
            %     if isempty(bool_stack) == 1
            %         bool_stack = bool_mat_loop;
            
            
           
            
            rates_idx = [rates_idx; idxs(chain(start_idx),chain(stop_idx));idxs(chain(stop_idx),chain(start_idx))];
            Nrates = length(rates_idx);
            rates = reshape(rand(1,Nrates), size(rates_idx));
            % longRand here ?
            
            K = bool_mat_loopN;
            
            K(rates_idx) = rates;
            

            rates_str_woDB = str_list(find(bool_mat_loopN == 1));  % select the rate strings based on the occupation of the matrix
            rates_idx_woDB = find(bool_mat_loopN == 1);

            loop_memb = chain(min(start_idx,stop_idx):max(start_idx,stop_idx)); % FIND A BETTER WAY TO DO THIS
            for_prod = 1;
            for_num = 0;
            for_prod_idx_arr = [];
            for i=1:length(loop_memb)
                if i == length(loop_memb) && (bool_mat_loopN(loop_memb(1),loop_memb(length(loop_memb))) ==1 || bool_mat_loopN(loop_memb(1),loop_memb(length(loop_memb))) < identifyDB)%~= identifyDB %~= 2
                    for_prod = for_prod * K(loop_memb(1), loop_memb(length(loop_memb)));
                    for_num = for_num + 1;
                    for_prod_idx_arr = [for_prod_idx_arr, idxs(loop_memb(1), loop_memb(length(loop_memb)))];
                elseif i < length(loop_memb) && (bool_mat_loopN(loop_memb(i),loop_memb(i+1)) == 1 || bool_mat_loopN(loop_memb(i),loop_memb(i+1)) < identifyDB)%~= identifyDB %~= 2
                    for_prod = for_prod * K(loop_memb(i), loop_memb(i+1));
                    for_num = for_num + 1;
                    for_prod_idx_arr = [for_prod_idx_arr, idxs(loop_memb(i), loop_memb(i+1))];
                end
            end
            rev_prod = 1;
            rev_num = 0;
            rev_prod_idx_arr = [];
            for i=1:length(loop_memb)
                if i == length(loop_memb) && ( bool_mat_loopN(loop_memb(length(loop_memb)),loop_memb(1)) == 1 || bool_mat_loopN(loop_memb(length(loop_memb)),loop_memb(1)) < identifyDB) %~= identifyDB % ~= 2
                    rev_prod = rev_prod * K(loop_memb(length(loop_memb)), loop_memb(1));
                    rev_num = rev_num + 1;
                    rev_prod_idx_arr = [rev_prod_idx_arr, idxs(loop_memb(length(loop_memb)),loop_memb(1))];
                elseif i < length(loop_memb) && (bool_mat_loopN(loop_memb(i+1),loop_memb(i)) == 1 || bool_mat_loopN(loop_memb(i+1),loop_memb(i)) < identifyDB) %~= identifyDB %~= 2
                    rev_prod = rev_prod * K(loop_memb(i+1), loop_memb(i));
                    rev_num = rev_num + 1;
                    rev_prod_idx_arr = [rev_prod_idx_arr, idxs(loop_memb(i+1),loop_memb(i))];
                end
            end

            if isempty(dbCell) == 0
                dbCell =[dbCell; {for_prod_idx_arr, rev_prod_idx_arr}]; % DOESN'T SEEM TO FILL FOR >1 ADDITIONAL CONNECTION
            else
                dbCell = {for_prod_idx_arr, rev_prod_idx_arr};
            end

%             disp(dbCell{1})
%             disp(dbCell{2})
            
            if for_num < rev_num
                K(rates_idx_DB(end))= rev_prod/for_prod;
            elseif for_num > rev_num
                K(rates_idx_DB(end))= for_prod/rev_prod;
            else
                error(['DB error: forward & backward loop are same size'])
            end
            
            diagK = -sum(K,1);
            K(diag(idxs)) = diagK;
            if sum(K,1) ~= zeros(1,N)
                error(['ERROR: columns of K do not sum to zero!'])
            end
            
            if sum(sum(sum(bool_stack - (bool_mat_loopN > 0) ==0))== N*N) == 0
                bool_stack = cat(3,bool_stack, (bool_mat_loopN>0));
            else
                continue
            end
            

            [row,col] = ind2sub(size(bool_mat_loopN),find(bool_mat_loopN(:,:) > 2));
            for i=1:length(row)
                if i == 1
                    model_name = strrep([model_name_loop  '_cxn', num2str([row(i) col(i)])],' ','');
                else
                  model_name = strrep([model_name  '_cxn', num2str([row(i) col(i)])],' ','');  
                end
            end
            
            K = zeros(N, N);
            
            numLoops = j+1;
            %         model_loopN = [model_loopN; {K, rates_str_woDB, rates_str_DB, rates_idx_DB, rates_idx_woDB, bool_mat_loop, cxns_all, cxns_cur, model_name, numLoops, chain}];
            model_loopN = [model_loopN; {model_name, K, rates_idx_woDB, rates_str_woDB, chain, bool_mat_loopN, rates_idx_DB, rates_str_DB, dbCell, cxns_all, cxns_cur,numLoops }];
            %
        end
    end
    
end
% for i = 1:size(bool_stack,3)
%     figure();clf;
%     plot(graph(bool_stack(:,:,i)))
% end

disp(['Found ' num2str(size(model_loopN,1)) ' loopN connectivities for N = ' num2str(N)])
else
    model_loopN = [];
end

%%

if strcmp(model_set, 'model_brch') == 1

model_brch = {};
idxs = reshape(linspace(1,N*N,N*N),[N,N]);

for q = 1:size(model_lin,1)
    [model_name_lin, K, rates_idx_woDB, rates_str_woDB, chain,bool_mat] = model_lin{q,:};
    % in terms of indices
    brk_idxs = [1, length(chain)];
    for l = 1:length(brk_idxs)
        start_idx = brk_idxs(l); % 1 or end
        brchs = linspace(1,length(chain),length(chain));
        
        if start_idx == 1
            allowed = brchs <= length(brchs)-1;
            allowed = allowed .* (brchs >= 3);
        else
            allowed = brchs <= length(brchs)-2;
            allowed = allowed .* (brchs >= 2);
        end
        brchs = brchs(nonzeros(brchs .* allowed));
        % from here on in terms of state labels
        poss_brchs = chain(brchs);
        
        if start_idx == 1
            brk_cxn = [chain(start_idx) chain(start_idx + 1)];
        else
            brk_cxn = [chain(start_idx-1) chain(start_idx)];
        end
        
        % start loop over possible brchs
        for i = 1:length(poss_brchs)
            
            bool_mat_brch = bool_mat;
            bool_mat_brch(brk_cxn(1), brk_cxn(2)) = 0;
            bool_mat_brch(brk_cxn(2), brk_cxn(1)) = 0;
            
            bool_mat_brch(chain(start_idx), poss_brchs(i)) = 1;
            bool_mat_brch(poss_brchs(i), chain(start_idx)) = 1;
            
            if sum(sum(sum(bool_stack - (bool_mat_brch > 0) ==0))== N*N) == 0
                bool_stack = cat(3,bool_stack, (bool_mat_brch>0));
            else
                continue
            end
            
            rates_str_woDB = str_list(find(bool_mat_brch == 1));  % select the rate strings based on the occupation of the matrix
            rates_idx_woDB = find(bool_mat_brch == 1);
            
            Nrates = length(rates_idx_woDB);
            rates = linspace(1,Nrates,Nrates);
            % logRandUniform here!
            
            K = bool_mat_brch;
            K(rates_idx_woDB) = rates;
            diagK = -sum(K,1);
            K(diag(idxs)) = diagK;
            
            if start_idx == 1
                model_name = [model_name_lin '_del' num2str(chain(start_idx)) num2str(chain(start_idx+1)) '_add' num2str(chain(start_idx)) num2str(poss_brchs(i))];
            else
                model_name = [model_name_lin '_del' num2str(chain(start_idx-1)) num2str(chain(start_idx)) '_add' num2str(chain(start_idx)) num2str(poss_brchs(i))];
            end
            
            %         model_name_brch = [model_name_brch; {K, rates_str_woDB, rates_str_DB, rates_idx_DB, rates_idx_woDB, bool_mat_brch, poss_brchs, [chain(start_idx), poss_brchs(i)], model_name, chain}];
            model_brch = [model_brch; {model_name, K, rates_idx_woDB, rates_str_woDB, chain, bool_mat_brch, poss_brchs, [chain(start_idx), poss_brchs(i)]}];
            % * do we want to come back and make the output array the same for
            % all models? Or have different output arrays for different models
            
            %    model_loopN = [model_loopN; {model_name, K, rates_idx_woDB, rates_str_woDB, chain, bool_mat_loop, rates_idx_DB, rates_str_DB,for_prod_idx_arr, rev_prod_idx_arr, cxns_all, cxns_cur,numLoops }];
            
        end
        
    end
end

disp(['Found ' num2str(size(model_brch,1)) ' branch connectivities for N = ' num2str(N)])
else
    model_brch = [];
end

% for i = 1:length(model_name_brch(:,6))/2
%     figure();clf;
%     plot(graph(model_name_brch{i,6}))
% end


% use product of allowed masks for start_idx and 1 and end to identify
% place for second (or nth?) branch point within an existing branched model



%% Saving libraries
% model_lin
% model_loop1
% model_loopN
% model_brch

if saveMode == 1
    libraryPath = '/Users/calbrecht/Dropbox/MATLAB_programs/claire_programs/from_Jack/updatedgenalgcodes';
    outputFold = [libraryPath filesep() 'modelLibraries'];
if exist(outputFold, 'dir')~= 7
    mkdir(outputFold);
end
    save([outputFold filesep() 'lib_N' num2str(N)], 'model_lin','model_loop1','model_loopN','model_brch')
    disp(['...saving model libraries in: ',outputFold])
end
