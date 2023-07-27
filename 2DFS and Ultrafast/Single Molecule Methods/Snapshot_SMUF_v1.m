%"SNapshot" viewer of t21 data sets for deciding on an Autocut of all
%trajectories or a manual cut to eliminate errors

fileFolder = 'C:\Users\Ryzen 5\Dropbox\SM_UF_Data\1DTimePhaseTags_2DFPGA\RhodB\02-28-2023';
outFolderGen= 'C:\Users\Ryzen 5\Dropbox\SM_UF_Data\1DTimePhaseTags_2DFPGA\OutputData\ChosenMolecules';
ConstructName='RhodB'; 

cd(fileFolder); 

filesList = dir('*2DFPGA'); 
numFiles = numel(filesList); 
curPath = pwd;
usres=30000;

  for i=1:length(filesList)

    
 
        fileFolder= [curPath filesep() filesList(i).name]; 
        scanID = filesList(i).name(1:15);

        tb1 = dlmread(fullfile(fileFolder, 'timebase1.txt'));
        tb2 = dlmread(fullfile(fileFolder, 'timebase2.txt'));
        xl = length(tb1);
        yl = length(tb2);
        
        if i==1
            OGTimeInc = (tb2(2)-tb2(1))*1000;
            OGTimeInc=round(OGTimeInc); 
            halfBinWidth=OGTimeInc/4;           
            arbTimeAxis=linspace(1,length(tb2),length(tb2));
            arbTimeAxis=[0,arbTimeAxis*OGTimeInc]; 
            lowerBounds=arbTimeAxis-halfBinWidth; 
            upperBounds=arbTimeAxis+halfBinWidth;
            globBounds=[];
            foldTable=[];
            
            for j=1:length(lowerBounds)
                if j==length(lowerBounds)
                globBounds=[globBounds;[lowerBounds(j), upperBounds(j)]]; 
                else
                globBounds=[globBounds;[lowerBounds(j), upperBounds(j)] ; [upperBounds(j), lowerBounds(j+1)]];
                end
                
                if j==length(lowerBounds)
                outFolder1 = ([outFolderGen filesep() ConstructName '_' num2str(usres) 'usec'  filesep() 'Tau=' num2str(lowerBounds(j)) '-' num2str(upperBounds(j)) 'fs']);
                foldTable=[foldTable, struct('f',outFolder1)];
                else
                outFolder1 = ([outFolderGen filesep() ConstructName '_' num2str(usres) 'usec'  filesep() 'Tau=' num2str(lowerBounds(j)) '-' num2str(upperBounds(j)) 'fs']);
                outFolder2 = ([outFolderGen filesep() ConstructName '_' num2str(usres) 'usec'  filesep() 'Tau=' num2str(upperBounds(j)) '-' num2str(lowerBounds(j+1)) 'fs']);
                foldTable=[foldTable, struct('f',outFolder1),struct('f',outFolder2)];

                end
                
            end
            % create a structure to hold the output folder names for use within the
            % loop
%             foldTable= struct('f1',outFolder1,'f2',outFolder2,'f3',outFolder3,'f4',outFolder4,'f5',outFolder5,'f6',outFolder6,'f7',outFolder7,'f8',outFolder8,'f9',outFolder9,'f10',outFolder10);             
        end
        

%       zero the remaining steps to start from step 1=0fs
%       ( need the -1 since the X-values are in the
%       -x direction )
        tb1=-1*(tb1 - tb1(1));
        tb2=(tb2 - tb2(1));
        tb1FS=tb1*1000;
        tb2FS=tb2*1000;        
        
        XmatTemp = zeros(xl, yl);
        YmatTemp = zeros(xl, yl);
%         
%         if j==1
%            Xmat= zeros(xl, yl);
%            Ymat= zeros(xl, yl);
%         end 
        
%%


        for yi = 1:yl
            for xi = 1:xl
                p1file = fullfile(fileFolder, [num2str(xi-1,'%02.0f') num2str(yi-1, '%02.0f') 'p1.bin']);
                p2file = fullfile(fileFolder, [num2str(xi-1,'%02.0f') num2str(yi-1, '%02.0f') 'p2.bin']);
                timefile = fullfile(fileFolder, [num2str(xi-1,'%02.0f') num2str(yi-1, '%02.0f') 'time.bin']);
                p1ID = fopen(p1file);
                p2ID = fopen(p2file);
                timeID = fopen(timefile);
                p1 = fread(p1ID,Inf,'float64=>double',0,'s');
                p2 = fread(p2ID,Inf,'float64=>double',0,'s');
                time = fread(timeID,Inf,'uint64=>uint64',0,'s');
                fclose('all');
                
                times = ((double (time-time(1))./(8e7))+1);
                % take out the first time from all later times, making event 1 occur at
                % t=0sec, then convert the times to usec in another array
                times=times-times(1);
                timesUsec=times*1e6; 
                % scan length is in seconds, scanLengths is just a holder for all scan
                % times at a glance
                scanLength= times(end);
                
                numWindows=floor((scanLength*1e6)/usres); 
                plotTime=linspace(0,scanLength,numWindows); 
                countsTraj = zeros(1,numWindows);
                PhaseFctTrajP1 = zeros(1,numWindows);
                
                
                      XquadP1 = mean(cos(p1.*2*pi));
                      YquadP1 = mean(sin(p1.*2*pi));                   
                      shiftAngleP1=atan2(YquadP1,XquadP1);
                      p1 = (p1.*2*pi)-(shiftAngleP1);
                      
                        for n=1:numWindows
                           Phases1= p1(find(timesUsec<=n*usres & timesUsec>=(n-1)*usres ));                    
                        % create the 30ms trajectories for plotting purposes. 
                           countsTraj(n)= numel(Phases1);
                           PhaseFctTrajP1(n)= mean(exp(1i*(Phases1)));
                        end 
                        
                        %here the individual traj has been processed and is
                        %stored in the arrays above (lines 115-120)
                        %NOTE: try using histcounts() array to fill in the
                        %photon stats for this data set 
                        
                        
            end
        end
        
  end