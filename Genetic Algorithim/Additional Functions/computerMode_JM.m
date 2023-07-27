% AUTHOR:   Claire Albrecht
% PURPOSE:  Identify the computer running the program and give that
%           informatin to the genetic algorithm codes so it can easily find files
%           properly.
%
% MODIFICATIONS:
%   CA: Added the DropboxLocationPrefix output to accomodate the fmincon
%   code
%----------------------------------------------------------------------------

function [computer_terminal_str, terminalID, DropboxLocationPrefix] = computerMode_JM()
global useFigurePosnMode
switch nargin
    case 0
end
wd = pwd;
%----------------------------------------------------------------------------
% User Options
%----------------------------------------------------------------------------
verboseMode = 0;

% Define identifying strings in each computer's command line
baimi = 'C:\Users\baimi\';
bisraels = '/Users/bisraels/';
claire = '/Users/clairealbrecht/';
calbrecht = '/Users/calbrecht/';
Brett = 'X:\Users\Brett\';
JackLaptop='C:\Users\Ryzen 5\';
labTower='C:\Users\jmaurer3';

% Name the computers and create string for  terminal string
computer_IDs = char(baimi, bisraels, claire, calbrecht, Brett, JackLaptop, labTower);
computer_names = char('baimi', 'bisraels', 'claire','calbrecht','Brett', 'JackLaptop', 'labTower');
% mode_str  = char('computer_baimi_mode', 'computer_bisraels_mode','computer_claire_mode')

numComputers = numel(computer_IDs(:,1));
for compIdx = 1:numComputers
    if contains(wd, deblank(computer_IDs(compIdx,:)))   % Find identifying piece of working directory
        user = deblank(computer_names(compIdx,:));  % Associate the user name
        terminalID = deblank(computer_IDs(compIdx,:));
        if verboseMode == 1
            disp(['function computerMode(): Assigning terminalID as ' terminalID]);
        end
        break;
    else
        user = 'general';           % If not one of the normal computers - general computer
        computer_general_mode = 1;
        useFigurePosnMode = 0;
          if verboseMode == 1
            disp(['function computerMode(): Assigning useFigurePosnMode = 0' ...
                10 'User = general']);
        end
    end
end

computer_terminal_str = ['computer_' user '_mode'];  % Define the output

% if user=='JackLaptop'
% DropboxLocationPrefix = [terminalID 'Dropbox' filesep() 'chosenAPD2mat_output'];   
% else
% DropboxLocationPrefix = [terminalID 'Dropbox' filesep() 'MarcusLab'];  %%% ADDED THIS LINE JANUARY 2020 %%%
% end
end

    