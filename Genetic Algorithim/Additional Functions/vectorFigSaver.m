function vectorFigSaver(foutName,folderNameOut,extension)
%__________________________________________________________________________
% Author: Brett Israels
%
% FUNCTION: Saves figure as a vector graphic. Will be the same size as it
% is on the screen. Small screen size = small file size, vice versa too.
%
% MODIFICATION LOG:
%  BI  9/29/16  Creaton.
%  BI  4/30/21  Documentation
%__________________________________________________________________________

%% Program options
% Folder name to store traces in
folderNameOut_default = '.';
foutName_default = 'vectorPDF';
extension_default = 'pdf';
verboseMode = 1;
saveFigureMode = 1;

%% Check for user input
switch nargin
    case 0
        % foutName is the name of the file you are about to write
        foutName = input('Enter the name of the file to be saved: ','s');
        if isempty(foutName)
            foutName = foutName_default;
            disp(['Setting foutName  = ' foutName]);
        end
        %         folderNameOut = input(['Enter the name of folder to be save the file (Press ENTER for ' folderNameOut_default '): '],'s');
        %         if isempty(folderNameOut),
        %             folderNameOut = folderNameOut_default;
        %             dis            p(['Setting output Folder  = ' folderNameOut]);
        %         end
        folderNameOut = folderNameOut_default;
        extension = extension_default;
    case 1
        folderNameOut = folderNameOut_default;
        %         folderNameOut = input(['Enter the name of folder to be save the file (Press ENTER for ' folderNameOut_default '): '],'s');
        %          if isempty(folderNameOut),
        %             folderNameOut = folderNameOut_default;
        %             disp(['Setting output Folder  = ' folderNameOut]);
        %          end
        extension = extension_default;
        
    case 2
        extension = extension_default;
end


%% Check to see if the folder already exist
if exist(folderNameOut,'dir') ~= 7
    % If it does not exist, make the folder.
    mkdir(folderNameOut);
    disp(['     *Creating a folder named ' folderNameOut ' to store the figure in.']);
end

%% Save the PDF
fig = gcf;
set(gcf,'Units','Inches');
pos = get(gcf,'Position');
set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3),pos(4)]);
% print(fig,[folderNameOut filesep() foutName],['-d' extension],'-r0');
print(fig,[folderNameOut filesep() foutName '.' extension],['-d' extension],'-r0');

if verboseMode == 1
    disp(['      **Saving figure as ' foutName '.' extension ' in ' folderNameOut]);
end
%make the figure units in pixels again
set(gcf,'Units','pixels');

%% Save the MATLAB Figure
if saveFigureMode == 1
   savefig([folderNameOut filesep() foutName '.fig']);
end
