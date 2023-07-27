% avg2DFPGA_IP (individual phasing)

function [RPmat, NRPmat, tb1, tb2] = avg2DFPGA_IP(homeFolder)
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
                p1file = fullfile(fileFolder, [num2str(xi-1,'%02.0f') num2str(yi-1, '%02.0f') 'p1.bin']);
                p2file = fullfile(fileFolder, [num2str(xi-1,'%02.0f') num2str(yi-1, '%02.0f') 'p2.bin']);
                p1ID = fopen(p1file);
                p2ID = fopen(p2file);
                p1 = fread(p1ID,Inf,'float64=>double',0,'s');
                p2 = fread(p2ID,Inf,'float64=>double',0,'s');
                fclose('all');
                RPmatTemp(xi, yi) = mean(exp(2i*pi*(p1-p2)));
                NRPmatTemp(xi, yi) = mean(exp(2i*pi*(p1+p2)));
            end
        end
        angle(RPmatTemp(2,1))
        angle(NRPmatTemp(2,1))
        RPmatTemp =RPmatTemp*exp(1i*(angle(RPmatTemp(2,1))));
        NRPmatTemp = NRPmatTemp*exp(1i*(angle(NRPmatTemp(2,1))));
        
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