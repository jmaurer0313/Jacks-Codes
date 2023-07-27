% GOAL: obtain the columns and rows of the 2D-FPGA data as linear abs
% spectra of the molecule. Given a 10x10 grid there should be a total of 20
% lin abs spectra possible, 10 for 5Khz along axis 1 10 for 8Khz along axis
% 2. The dimension along which the x-stage is
% stepped should be considered first as these steps are made sequentially
% rather than being seperated by 10 waiting periods like the Y-stage steps.
% 


function [XmatTemp, YmatTemp, Xabsp, Yabsp, XabspAVG, YabspAVG, tb1, tb2] = linear2DFPGA(homeFolder)
startFolder = pwd;
cd(homeFolder); 

outFolder1 = (['XLinearSpectra']);
if exist(outFolder1, 'dir')~=7
    mkdir(outFolder1)
end

outFolder2 = (['YLinearSpectra']);
if exist(outFolder2, 'dir')~=7
    mkdir(outFolder2)
end

outFolder3 = (['AvgLinearSpectra']);
if exist(outFolder3, 'dir')~=7
    mkdir(outFolder3)
end
    
filesList = dir('*2DFPGA'); 
numFiles = numel(filesList); 
curPath = pwd;
zeroPad=32; 
c0 = 0.000299792458; % speed of light mm/fs


for j=1:numFiles
%     
    name= filesList(j).name(1:end-7);
    fileFolder= [curPath filesep() filesList(j).name]; 

        tb1 = dlmread(fullfile(fileFolder, 'timebase1.txt'));
        tb2 = dlmread(fullfile(fileFolder, 'timebase2.txt'));
        xl = length(tb1);
        yl = length(tb2);
        
        tx= tb1*(-1)*(1e3);
        ty= tb2*1e3; 
        Ts=abs(tx(1)-tx(2));
        
        windowX = delayedGaussian(tx, 15, 25);
        windowY = delayedGaussian(ty, 0, 45);
        
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
            XmatPlot= XmatTemp;
            XmatPlot(1,:)= [];
            
            YmatPlot= YmatTemp;
            YmatPlot(1,:)= [];
%             YmatTemp(1,:)= [];
%         

%         Transpose the Xmat so that it is the rows in both cases being
%         transformed (For instance in the case of 9 steps in Y and 10
%         steps in X (1 extra), you have 9 instances of 10 X data points (10 X-points for each Y-point),
%         such that you will have 9 lin abs spectra in the X channel and 10
%         in the y channel. 
            XmatTemp =XmatTemp.';
            XmatPlot =XmatPlot.';
%             angle(XmatTemp(1,2))
%             angle(YmatTemp(2,1))
            figure(4)
            plot(tx, abs(XmatTemp.*windowX))
            title('X time domain Windowed');
            
            figure(5)            
            plot(ty, abs(YmatTemp.*windowY))
            title('Y time domain Windowed');
            
            
%             multiply the first entry in the time domain by 0.5 to handle
%             the single sidedness
        XmatPlot(:,1)=XmatPlot(:,1)*0.5;
        XmatPlot(:,end)=XmatPlot(:,end)*0.5;
        
        YmatPlot(:,1)=YmatPlot(:,1)*0.5;
        YmatPlot(:,end)=YmatPlot(:,end)*0.5;
%         XmatTemp(:,1)=XmatTemp(:,1)*0.5;
%         YmatTemp(:,1)=YmatTemp(:,1)*0.5;
        
        NFFT=zeroPad+size(XmatPlot,1); 
        f = ceil((-NFFT/2:NFFT/2-1))/(Ts*NFFT); % freq axis in PHz
        
        lam_mono = 550; % monochromater wavelength in nanometers
        Fmono = 1e7 / lam_mono;
        % T = 7/6/1e3*2/c0;
        Fs = 10/(c0*Ts);
        Faxis = Fs/2*linspace(-1,1,NFFT);

%        For each row of the Xmat and each row of the Ymat, obtain the
%        phase (angle) of the first time point and phase the entire vector(
%       multiply by conj(phasefactor at time zero) 
        for i=1:size(XmatPlot,1)
            XmatPlot(i,:)=XmatPlot(i,:).*exp(-1i*angle(XmatPlot(i,1)));
        end 
        
        for i=1:size(YmatPlot,1)
            YmatPlot(i,:)=YmatPlot(i,:).*exp(-1i*angle(YmatPlot(i,1)));
        end 
        
        for i=1:size(XmatTemp,1)
            XmatTemp(i,:)=XmatTemp(i,:).*exp(-1i*angle(XmatTemp(i,1)));
        end 
        
        for m=1:size(YmatTemp,1)
            YmatTemp(m,:)=YmatTemp(m,:).*exp(-1i*angle(YmatTemp(m,1)));
        end 
%         
        figure(6)            
            plotAbs = plot(Faxis+Fmono,abs(fftshift(fft((XmatPlot),NFFT, 2),2)).','color',[1 1 0]);
            hold on
%           Faxis+Fmono,real(fftshift(fft((XmatPlot),NFFT, 2),2)).',Faxis+Fmono,imag(fftshift(fft((XmatPlot),NFFT, 2),2)).',
            plotReal = plot(Faxis+Fmono,real(fftshift(fft((XmatPlot),NFFT, 2),2)).','color',[0.8500 0.3250 0.0980]);
            plotImag = plot(Faxis+Fmono,imag(fftshift(fft((XmatPlot),NFFT, 2),2)).','color',[0 0.4470 0.7410]);
            ax = gca;
            ax.XAxis.Exponent = 3;
            title(['X Scans Frequency Domain - from 2D ' name],'FontSize',15);
            xlim([17000 20500]);
            xlabel('Wavenumbers (cm^{-1})');
            ylabel('Amplitude (AU)');
            legSet= [plotAbs(1) plotReal(1) plotImag(1)];
            legend(legSet, 'Abs', 'Real','Imag');
            hold off
            
            figure(7)
            plotAbs = plot(Faxis+Fmono,abs(fftshift(fft((YmatPlot),NFFT, 2),2)).','color',[1 1 0]);
            hold on
%           Faxis+Fmono,real(fftshift(fft((XmatPlot),NFFT, 2),2)).',Faxis+Fmono,imag(fftshift(fft((XmatPlot),NFFT, 2),2)).',
            plotReal = plot(Faxis+Fmono,real(fftshift(fft((YmatPlot),NFFT, 2),2)).','color',[0.8500 0.3250 0.0980]);
            plotImag = plot(Faxis+Fmono,imag(fftshift(fft((YmatPlot),NFFT, 2),2)).','color',[0 0.4470 0.7410]);
            ax = gca;
            ax.XAxis.Exponent = 3;
            title(['Y Scans Frequency Domain - from 2D ' name],'FontSize',15);
            xlim([17000 20500]);
            xlabel('Wavenumbers (cm^{-1})');
            ylabel('Amplitude (AU)');
            legSet= [plotAbs(1) plotReal(1) plotImag(1)];
            legend(legSet, 'Abs', 'Real','Imag');
            hold off
            
        
% %      fft(X,n,dim) returns the Fourier transform along the dimension dim. 
% %      For mexaple, if X is a matrix, then fft(X,n,2) returns the n-point 
% %      Fourier transform of each row.

%         Transform for each file to yield the resulting lin abs spectra 
        
          Xabsp = fftshift(fft(XmatTemp,size(XmatTemp,1)+zeroPad, 2),2);       
          Yabsp = fftshift(fft(YmatTemp,size(YmatTemp,1)+zeroPad, 2),2);   
          
%           XabspAVG = fftshift(fft(mean(XmatTemp,1),size(XmatTemp,1)+zeroPad));       
%           YabspAVG = fftshift(fft(mean(YmatTemp,1),size(YmatTemp,1)+zeroPad));
          
          XabspAVG = fftshift(fft(mean(XmatPlot,1),size(XmatPlot,1)+zeroPad));       
          YabspAVG = fftshift(fft(mean(YmatPlot,1),size(YmatPlot,1)+zeroPad));
   
          YmatTemp =YmatTemp.'; 
          
          figure(1);            
          for n=1:size(Xabsp,1)
             plot(abs(Xabsp(n,:)));
             hold on          
          end
           title('X Lin Abs');
            xlabel('Arb. Freq');
            ylabel('Amp');
           
          
          figure(2);            
          for n=1:size(Yabsp,1)
             plot(abs(Yabsp(n,:)));
             hold on          
          end
          title('Y Lin Abs');
            xlabel('Arb. Freq');
            ylabel('Amp');
          
            
          figure(3);            
          plot(real(XabspAVG));
          hold on
          plot(real(YabspAVG));
          title('X and Y Avg Lin Abs');
            xlabel('Arb. Freq');
            ylabel('Amp');       
          
          
          saveas(figure(6),[pwd '/XLinearSpectra/' filesList(j).name '_Xabs.fig']);
          saveas(figure(7),[pwd '/YLinearSpectra/' filesList(j).name '_Yabs.fig']);
          saveas(figure(3),[pwd '/AvgLinearSpectra/' filesList(j).name '_XYavgAbs.fig']);
            
%           close all

end 

end 
