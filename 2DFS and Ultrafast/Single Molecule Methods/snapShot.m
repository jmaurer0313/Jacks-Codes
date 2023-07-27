function [elim_IDlist] = snapShot(fileFolder,usres)

cd(fileFolder); 

filesList = dir('*2DFPGA'); 
numFiles = numel(filesList); 
curPath = pwd;

min_counts = [];
max_counts = [];
  
  for i=1:length(filesList)

        fileFolder= [curPath filesep() filesList(i).name]; 
        scanID = filesList(i).name(1:15);

        tb1 = dlmread(fullfile(fileFolder, 'timebase1.txt'));
        tb2 = dlmread(fullfile(fileFolder, 'timebase2.txt'));
        xl = length(tb1); %1
        yl = length(tb2); %7

        min_occs = [];
        max_occs =[];
        for yi = 1:yl
            for xi = 1:xl
               
                p1file = fullfile(fileFolder, [num2str(xi-1,'%02.0f') num2str(yi-1, '%02.0f') 'p1.bin']);
                p2file = fullfile(fileFolder, [num2str(xi-1,'%02.0f') num2str(yi-1, '%02.0f') 'p2.bin']);
                timefile = fullfile(fileFolder, [num2str(xi-1,'%02.0f') num2str(yi-1, '%02.0f') 'time.bin']);
                p1ID = fopen(p1file);
                p2ID = fopen(p2file);
                timeID = fopen(timefile);
                p1 = fread(p1ID,Inf,'float64=>double',0,'s'); % Inf - Column vector, with each element containing a value in the file
                p2 = fread(p2ID,Inf,'float64=>double',0,'s');
                time = fread(timeID,Inf,'uint64=>uint64',0,'s');
                fclose('all');
                
                times = ((double (time-time(1))./(8e7))+1);
                times=times-times(1);
                timesUsec=times*1e6; 
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
                        scanIDfull = [scanID '_' num2str(xi*yi)];
                        minocc = min(countsTraj);
                        maxocc = max(countsTraj);

                        minoccs = {minocc; scanIDfull};
                        maxoccs = {maxocc; scanIDfull};
            end
        min_occs = [min_occs minoccs];
        max_occs = [max_occs maxoccs];
        end
   min_counts = [min_counts min_occs];
   max_counts = [max_counts max_occs];
  end

  hmin = cell2mat(min_counts(1,:));
  hmax = cell2mat(max_counts(1,:));

%%% Condition 1: minimum

  pd1 = fitdist(hmin.','Normal');
  x_values1 = min(hmin):1:max(hmin);
  y1 = pdf(pd1,x_values1);
%   plot(x_values1,y1)
  sigma_min = pd1.sigma
  lowb1 = mean(hmin) - 2*sigma_min;
  upb1 = mean(hmin) + 2*sigma_min;
  [V,I]=find(hmin < lowb1 | hmin > upb1);
  scanIDSkip1=min_counts(2,I);

%%% Condition 2: maximum
  pd2 = fitdist(hmax.','Normal');
  x_values2 = min(hmax):1:max(hmax);
  y2 = pdf(pd2,x_values2);
%   plot(x_values2,y2)
  sigma_max = pd2.sigma;
  lowb2 = mean(hmax) - 2*sigma_max;
  upb2 = mean(hmax) + 2*sigma_max;

  [V,I]=find(hmax < lowb2 | hmax > upb2);
  scanIDSkip2=max_counts(2,I);

  finallist = [scanIDSkip1 scanIDSkip2];
  elim_IDlist = unique(finallist);
  Nscans = length(filesList)*xl*yl;
  if length(elim_IDlist) < 0.05*Nscans;
      disp(['Autocut is approved']);
  else
      disp(['Manual cut is recommended']);
      disp(['The number of suspicious scan = ', num2str(length(elim_IDlist)), ' out of total ', num2str(Nscans), ' scans']);
  end
  
end



