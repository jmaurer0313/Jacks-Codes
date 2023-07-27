% linear2DFPGA_RT

function [XmatTemp, YmatTemp, tb1, tb2] = linear2DFPGA_RT(filePath)

zeroPad=10; 

%     
    fileFolder= filePath; 

        tb1 = dlmread(fullfile(fileFolder, 'timebase1.txt'));
        tb2 = dlmread(fullfile(fileFolder, 'timebase2.txt'));
        xl = length(tb1);
        yl = length(tb2);
        
        XmatTemp = zeros(xl, yl);
        YmatTemp = zeros(xl, yl);
%         
%         if j==1
%            Xmat= zeros(xl, yl);
%            Ymat= zeros(xl, yl);
%         end 
        

        for yi = 1:yl
            for xi = 1:xl
                p1file = fullfile(fileFolder, [num2str(xi-1,'%02.0f') num2str(yi-1, '%02.0f') 'p1.bin']);
                p2file = fullfile(fileFolder, [num2str(xi-1,'%02.0f') num2str(yi-1, '%02.0f') 'p2.bin']);
                p1ID = fopen(p1file);
                p2ID = fopen(p2file);
                p1 = fread(p1ID,Inf,'float64=>double',0,'s');
                p2 = fread(p2ID,Inf,'float64=>double',0,'s');
                fclose('all');
%                 This will form a matrix of the 5 (x) and 8 (y) Khz components at
%                 every time step. The linear Abs Spectra can be obtained
%                 by transformation of the rows/columns. The 5Khz will only
%                 be defined along the rows and the 8Khz only defined
%                 alogn the columns, as the fringe pattern of the 5Khz is
%                 static along columns as is the 8Khz fringe pattern along
%                 the rows. 
%         
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
% %      For mexaple, if X is a matrix, then fft(X,n,2) returns the n-point 
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
