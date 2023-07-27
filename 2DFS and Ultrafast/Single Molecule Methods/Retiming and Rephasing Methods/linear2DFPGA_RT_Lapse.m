% linear2DFPGA_RT

%11-17-2022: this code simply processes the linear demodualtion channels
%from the 2D data grid, phases the time domain data, transforms them, and gives back the phased time domain data
% transposes are done to both the X and Y dimension - early pre-2021 scans
% will be 10x9 with 10 points in t21 - so that can be used to identify the
% axis in time corresponding to Xmat/Ymat and tb1/tb2
%

function [XmatTemp, YmatTemp, tb1, tb2, photGrid] = linear2DFPGA_RT_Lapse(filePath, intTime)

zeroPad=10; 

%     
    fileFolder= filePath; 

        tb1 = dlmread(fullfile(fileFolder, 'timebase1.txt'));
        tb2 = dlmread(fullfile(fileFolder, 'timebase2.txt'));
        xl = length(tb1);
        yl = length(tb2);
        
        XmatTemp = zeros(xl, yl);
        YmatTemp = zeros(xl, yl);
        photGrid= zeros(xl, yl);
%         
%         if j==1
%            Xmat= zeros(xl, yl);
%            Ymat= zeros(xl, yl);
%         end 
        

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
                photGrid(xi,yi)=numel(p1); 
                fclose('all');
%                 This will form a matrix of the 5 (x) and 8 (y) Khz components at
%                 every time step. The linear Abs Spectra can be obtained
%                 by transformation of the rows/columns. The 5Khz will only
%                 be defined along the rows and the 8Khz only defined
%                 alogn the columns, as the fringe pattern of the 5Khz is
%                 static along columns as is the 8Khz fringe pattern along
%                 the rows. 
%         

                %Xmat and yMat here are the linear signal componenets at 5
                %and 8kHz respectively 
                XmatTemp(xi, yi) = mean(exp(2i*pi*(p1)));
                YmatTemp(xi, yi) = mean(exp(2i*pi*(p2)));
            end
        end
%         
            

%           Throw out the first row of each Matrix which is the extra X
%           step for timing purposes 

%        NOTE: for the purposes of retiming each scan individually, the first row
%        should not be thrown out, since it is possible we missed the true
%        time zero toward the negative time domain.
%             XmatTemp(1,:)= [];
%             YmatTemp(1,:)= [];
%         

%         Transpose the Xmat so that it is the rows in both cases being
%         transformed (For instance in the case of 9 steps in Y and 10
%         steps in X (1 extra), you have 9 instances of 10 X data points (10 X-points for each Y-point),
%         such that you will have 9 lin abs spectra in the X channel and 10
%         in the y channel. 
            XmatTemp =XmatTemp.'; 
            
            

%        For each row of the Xmat and each row of the Ymat, obtain the
%        phase (angle) of the first time point and phase the entire vector(
%       multiply by conj(phasefactor at time zero) 

        for i=1:size(XmatTemp,1)
            XmatTemp(i,:)=XmatTemp(i,:)*conj(XmatTemp(i,1));
        end 
        
        for m=1:size(YmatTemp,1)
            YmatTemp(m,:)=YmatTemp(m,:)*conj(YmatTemp(m,1));
        end 
%         
        
        
% %      fft(X,n,dim) returns the Fourier transform along the dimension dim. 
% %      For example, if X is a matrix, then fft(X,n,2) returns the n-point 
% %      Fourier transform of each row.

%         Transform for each file to yield the resulting lin abs spectra 
        
          Xabsp = fftshift(fft(XmatTemp,size(XmatTemp,1)+zeroPad, 2),2);       
          Yabsp = fftshift(fft(YmatTemp,size(YmatTemp,1)+zeroPad, 2),2);   
          
          XabspAVG = fftshift(fft(mean(XmatTemp,1),size(XmatTemp,1)+zeroPad));       
          YabspAVG = fftshift(fft(mean(YmatTemp,1),size(YmatTemp,1)+zeroPad));          
   
          YmatTemp =YmatTemp.'; 
          
%           figure(1);            
%           for n=1:size(Xabsp,1)
%              plot(abs(Xabsp(n,:)));
%              hold on          
%           end
%            title('X Lin Abs');
%             xlabel('Arb. Freq');
%             ylabel('Amp');
%            
%           
%           figure(2);            
%           for n=1:size(Yabsp,1)
%              plot(abs(Yabsp(n,:)));
%              hold on          
%           end
%           title('Y Lin Abs');
%             xlabel('Arb. Freq');
%             ylabel('Amp');
          
            
%           figure(3);            
%           plot(abs(XabspAVG));
%           hold on
%           plot(abs(YabspAVG));
%           title('X and Y Avg Lin Abs');
%             xlabel('Arb. Freq');
%             ylabel('Amp');       
%           
%           
%           saveas(figure(1),[pwd '/XLinearSpectra/' filesList(j).name '_Xabs.fig']);
%           saveas(figure(2),[pwd '/YLinearSpectra/' filesList(j).name '_Yabs.fig']);
%           saveas(figure(3),[pwd '/AvgLinearSpectra/' filesList(j).name '_XYavgAbs.fig']);
            
%           close all

end 
