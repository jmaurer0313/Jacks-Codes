function [RPmat, NRPmat, tb1, tb2,times] = avg2DFPGA_Lapse(intTime,homeFolder)
startFolder= pwd;
cd(homeFolder); 
filesList= dir('*2DFPGA'); 
numFiles=numel(filesList); 
curPath=pwd;

for j=1:numFiles
    fileFolder= [curPath filesep() filesList(j).name]; 

        tb1 = dlmread(fullfile(fileFolder, 'timebase1.txt'));
        tb2 = dlmread(fullfile(fileFolder, 'timebase2.txt'));
        xl = length(tb1);
        yl = length(tb2);
        
        RPmatTemp = zeros(xl, yl);
        NRPmatTemp = zeros(xl, yl);
        
        if j==1
            RPmat= zeros(xl, yl);
            NRPmat= zeros(xl, yl);
        end 
        

        for yi = 1:yl
            for xi = 1:xl
                timeID = fullfile(fileFolder, [num2str(xi-1,'%02.0f') num2str(yi-1, '%02.0f') 'time.bin']);
                p1file = fullfile(fileFolder, [num2str(xi-1,'%02.0f') num2str(yi-1, '%02.0f') 'p1.bin']);
                p2file = fullfile(fileFolder, [num2str(xi-1,'%02.0f') num2str(yi-1, '%02.0f') 'p2.bin']);
                timeFile = fopen(timeID);
                p1ID = fopen(p1file);
                p2ID = fopen(p2file);
                time = fread(timeFile,Inf,'uint64=>uint64',0,'s');
                times = ((double (time-time(1))./(8e7))+1);
                times=times-times(1);
                p1 = fread(p1ID,Inf,'float64=>double',0,'s');
                p2 = fread(p2ID,Inf,'float64=>double',0,'s');
                [intWind]=find(times<=intTime);
%                At this point isolate the phase and time lists from zero
%                up to the integration time set
                p1=p1(1:intWind(end));
                p2=p2(1:intWind(end));   
                fclose('all');
                RPmatTemp(xi, yi) = mean(exp(2i*pi*(p1-p2)));
                NRPmatTemp(xi, yi) = mean(exp(2i*pi*(p1+p2)));
            end
        end
        
        RPmat= RPmat + RPmatTemp;
        NRPmat= NRPmat + NRPmatTemp;
            
           

end 
    RPmat= RPmat./numFiles;
    NRPmat= NRPmat./numFiles;
% Return time as a femtosecond column vector
tb1 = tb1.'*1e3;
tb2 = tb2.'*1e3;
cd(startFolder); 

end 