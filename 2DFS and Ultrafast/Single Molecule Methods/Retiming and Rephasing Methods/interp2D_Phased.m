% interp2D_Phased

% function which will take: 1) time domain data 2) Exp step size 3) desired
% step size (post interpolation) 4) Dim along which the interpolation will
% be performed. Function should interpolate both the 2D sideband and Linear
% data, so that the sideband data can be retimed and shifted according to
% the maximum of the abs(Linear data) .


% NOTE: the use of a Ts could be problematic given roundign done on the
% side of the LAbView to get an actual step size which is near to, but not
% exactly equal to the Ts of say 7fs (real is 6.9381fs for example)

function[ NRPtimed , RPtimed , tnX2 , tnY2 , window2Dfinal ] = interp2D_Phased( plotOpt, NRPmat, RPmat, tb1, tb2, NIFFT , Xs, Ys, Xzero, retiming, rephase ) 
% Conditional checks to ensure inputs are proper
% if size(S2D,Dim)==1
%     disp('Dimension specified is singular');
%        return;
% end
% 
% if Tdes<1e-3
%     disp('Desired time resolution is too small');
%     return; 
% end 
% 
% if abs(Tdes)>abs(Ts)
%     disp('Desired time resolution must be less than input time resolution');
%     return; 
% end 
    
%Constants to be defined for use througout the script
c0 = 299.792458; % speed of light in nm/fs
% as long as the NRP/RP mats do not get signifncantly larger this will
% work. since nextpow2 of 9 and 10 are both 4
NFFT = 2^(nextpow2(size(NRPmat,2)));
% NFFT=32; 
Npad = NIFFT-NFFT; % Number of zeros to pad, based on the difference of the 
%  size of the intial 2D time domain compared to the incoming interpolated
%  linear time domain.

% Set of time axes to handle a) the general plotting b) the windows for lin
% and 2D. Takes care of the extra step in X to allow for the possibity of
% retiming by +/- nfs

% if tb1(2)-tb1(1)<0.1
%  tx=(-1*(tb1))*1e3;
% ty=tb2*1e3; 
% Ts=abs(ty(2)-ty(1));
% t = (0:9)*Ts;   
% else

tx=(-1*(tb1));
ty=tb2; 
Ts=abs(tx(2)-tx(1));
t = (0:9)*Ts;
% end


% Introduce a windowing function to be applied to both the linear and 2D
% data in order to eliminate any ringing in the interploated time domain as
% a result of nonzero signal at the edge of the original raw data in time.
% IMPORTANTLY: the window must be applied to the data prior to the first
% FFT, otherwise the reuslt of padding by a small number of zeros
% NFFT-length(S) will give a shapr edge in the intial time domain, which
% rings in the frequency domain. THE VALUES CHOSEN HERE ARE MEANT TO
% INTRODUCE A GENTLE WINDOW WHICH STARTS VERY EARLY BUT DECAYS QUIT SLOWLY.
% GIVEN OUR COARSELY SAMPLE DATA THIS APPEARS TO THE ONLY SCHEME WHICH
% WORKS EFFECTIVELY
wX= repmat(tx',1,size(NRPmat,2)) - repmat(Xs(:,1)',size(NRPmat,1),1);
wY= repmat(ty,size(NRPmat,1),1) - repmat(Ys(:,1),1,size(NRPmat,2));
% [wX, wY] = meshgrid(tx2D,ty);
% standard window used for data print outs is 25,40
% window2D = delayedGaussian(sqrt(wX.^2 + wY.^2), 25, 40);
window2D = delayedGaussian(sqrt(wX.^2 + wY.^2), 40, 32);

if plotOpt
figure(15)
contourf(window2D);
title('Initial Window');
end

% reduce time domain values at zero along the first dimension being transformed by 0.5 to
% handle the cutoff from negative freqs to zero. Window the NRP and RP
% individually 
SRP= RPmat;
SNRP= NRPmat;

% SRP=SRP*exp(1i*(-angle(SRP(2,1))));
% SNRP=SNRP*exp(1i*(-angle(SNRP(2,1))));

SRP= SRP.*window2D;
SNRP= SNRP.*window2D;

SRP=SRP.';
SNRP=SNRP.';

% % % % TEST KIT: multiply the last point of each row and each column by 0.5
% % % % (last point in the X and Y runs) to see if the appearance of side band
% % % % crosspeaks is an artifact of sharp edged
% SRP(:,size(SRP,2))=SRP(:,size(SRP,2))*0.5;
% SRP(size(SRP,1),:)=SRP(size(SRP,1),:)*0.5;
% 
% SNRP(:,size(SNRP,2))=SNRP(:,size(SNRP,2))*0.5;
% SNRP(size(SNRP,1),:)=SNRP(size(SNRP,1),:)*0.5;


% make a plot of the X and Y nonlinear slices after windowing to examine
% the extent to which they are damped by the end of the array
if plotOpt
figure(1) % X slices
for j=1:size(SRP,1)
    plot(tx,real(SRP(j,:)));
    hold on
    plot(tx, real(SNRP(j,:)));
    title('Nonlinear T21 Signal - Pre-interpolation - Windowed')
    xlabel('Time (fs)');
    ylabel('Magnitude (arb.)'); 
end 
Xdamp = (mean(abs(SNRP(:,size(SNRP,2)))./abs(SNRP(:,1)))+ mean(abs(SRP(:,size(SRP,2)))./abs(SRP(:,1))))/2;

figure(20) % Y slices
for j=1:size(SRP,2)
    plot(ty',real(SRP(:,j)));
    hold on
    plot(ty', real(SNRP(:,j)));
    title('Nonlinear T43 Signal - Pre-interpolation - Windowed')
    xlabel('Time (fs)');
    ylabel('Magnitude (arb.)'); 
end 
end 

NRPrawplot=SNRP;
RPrawplot=SRP;
f = ceil((-NFFT/2:NFFT/2-1))/(Ts*NFFT); % freq axis in PHz

% if Dim==2
%     this is the case of rows for FFT

% try phasing the individual rows here and storing their values in a vector
% for reverse correction after interpolation
xRPphases=zeros(1,size(SRP,1));
xNRPphases=zeros(1,size(SNRP,1)); 

for m=1:length(xRPphases)
    xRPphases(m)=angle(SRP(m,1)); 
    xNRPphases(m)=angle(SNRP(m,1)); 
end 

SRP=SRP.*(exp(1i*(-1*xRPphases')));
SNRP=SNRP.*(exp(1i*(-1*xNRPphases')));
% at this point the NRP and RP mat are phased to the first point along the
% x dimension, such that what enters the first FFT should be purely real at
% the first time point
RPrawPhased=SRP;
NRPrawPhased=SNRP;

SNRP(:,1) = SNRP(:,1)*0.5;
SRP(:,1) = SRP(:,1)*0.5;

% FIRST TRANSFORM FROM TIME to OMEGA
SFRP = fftshift(fft(SRP, NFFT, 2),2); % Signal in the frequency domain, FFT shift applied 
SFNRP = fftshift(fft(SNRP, NFFT, 2),2); % Signal in the frequency domain, FFT shift applied 


if plotOpt
% first transform of the RP time domain along the rows (x runs)
figure(2)
contourf(real(SFRP)); 
title('RP 2D Freq Domain Data(Real)');

% first transform of the NRP time domain along the rows (x runs)
figure(3)
contourf(real(SFNRP)); 
title('NRP 2D Freq Domain Data(Real)');

figure(4)
contourf(real(SNRP));
title('NRP 2D Time Domain Data(Real)');

figure(5)
contourf(real(SRP));
title('RP 2D Time Domain Data(Real)');
end 

% try dropping the zeros at the edge of the freq domain after having windowed, OR the wrap around
% point where freq goes + to - 
SFRPzeros=[real(SFRP) , zeros(size(SFRP,1),Npad)];
SFNRPzeros=[real(SFNRP) , zeros(size(SFNRP,1),Npad)];

% sizeRPF= size(SFRPzeros)
% determine the new number of points and their values in the time domain for
% the interpolated data. Ensure that the new number of points along the X
% dimension matches the number of points from the interplated linear data
NIFFTX = length(SFRPzeros);
if NIFFT~=NIFFTX
    disp('Interpolated Linear and Nonlinear dimensions do not match in size');
    return;
end 


% Shift the newly zero-padded array by the amount which it was off from
% zero to start, then calculate the proper complex vector for each row/col
% (depending on Dim) from the real valued vectors in the freq domain.
% hilbert operates on columns when a matirx is passed.
unshiftedRP= circshift(SFRPzeros, -floor(NFFT/2),2);
unshiftedNRP= circshift(SFNRPzeros, -floor(NFFT/2),2);

% LAST KNOWN CONTACT HERE


% real and complex unshifted matrices, transposed here to handle the hilbert
% operation properly
RunshRP = unshiftedRP;
RunshNRP = unshiftedNRP;

SnRP = (ifft(RunshRP, NIFFT, 2));
SnNRP = (ifft(RunshNRP, NIFFT, 2));

SnRP(:,((ceil(NIFFT/2)+1):end))=0;
SnNRP(:,((ceil(NIFFT/2)+1):end))=0;

% MIGHT NEED A FACTOR OF 2 TWICE HERE, ONCE FOR THE 0.5 ENTERING THE FFT
% AND ONCE FOR CHOPPING THE NEGATIVE FREQS OFF THE IFFT RESULT
SnRP=SnRP.*2;
SnNRP=SnNRP.*2;

RPinterpPhased=SnRP;
NRPinterpPhased= SnNRP;

% reapply the original phases to the x dimension
SnRP=SnRP.*(exp(1i*(xRPphases')));
SnNRP=SnNRP.*(exp(1i*(xNRPphases')));

% sizeSt=size(SnRP)



% CunshRP = complex(RunshRP, imag(hilbert(-RunshRP)));
% CunshNRP = complex(RunshNRP, imag(hilbert(-RunshNRP)));
% 
% CunshRP= CunshRP.';
% CunshNRP= CunshNRP.';

% these are the new time axes for the interpolated 2D data BUT TAKEN FROM
% ZERO. Until the X-axis is time corrected, plotting against these values
% will not be indivative of the true value of the signal, since X starts
% from -Ts and goes to 9*Ts, rather than from 0 to 10*Ts.
Tsn = (Ts*NFFT)/NIFFT;
tn = ((0:NIFFT-1))*Tsn;
tnX=tn;
tnY=tn;

if plotOpt
figure(13)
contourf((real(RunshRP)));

figure(17)
contourf((real(RunshNRP)));
end

% 
% SnRP = (ifft(CunshRP, NIFFT, 2));
% SnNRP = (ifft(CunshNRP, NIFFT, 2));

% SnRP(:,1) = SnRP(:,1)*2;
% SnNRP(:,1) = SnNRP(:,1)*2;

if plotOpt
figure(8)
contourf((real(SnRP)));

figure(9)
contourf((real(SnNRP)));

end 

if plotOpt
    figure(18)
   
    for k=1:3
        % plot((t),imag(NRPrawplot(k,:)),'o',(t),real(NRPrawplot(k,:)),'o')
% hold on
% plot(tnX, imag(SnNRP(k,:))*(NIFFT/NFFT),tnX, real(SnNRP(k,:))*(NIFFT/NFFT))

plot((t),real(NRPrawplot(k,:)),'o');
hold on
plot(tnX, real(SnNRP(k,:))*(NIFFT/NFFT));
    end
    title('Nonlinear T21 Signal - Unretimed - Real Part - Unwindowed')
    xlabel('Time (fs)')
    ylabel('Magnitude (arb.)') 

    figure(21)
    title('Nonlinear T21 Signal - Interpolated - Phased to Pure Real');
    xlabel('Time (fs)');
    ylabel('Magnitude (arb.)');
    for k=1:3
plot((t),imag(NRPrawPhased(k,:)),'o',(t),real(NRPrawPhased(k,:)),'o')
hold on
plot(tnX, imag(NRPinterpPhased(k,:))*(NIFFT/NFFT),tnX, real(NRPinterpPhased(k,:))*(NIFFT/NFFT))
end 
end 


if retiming
% Now the time shifts can be applied to the resulting interplated time
% domain signals of the 2D data, for the X case apply the shifts
% individually by first zeroing out all data from start to Max-1, then
% circshift to set the maximum at the first index of each row
% sizeRPpreShifts=size(SnRP)
for i=1:size(SnRP,1)
    if (Xs(i,2))~=0
    SnRP(i,(1:Xs(i,2)-1))=0; 
    SnRP(i,:)=circshift((SnRP(i,:)),-(Xs(i,2)-1),2);  
    end
end

for n=1:size(SnNRP,1)
    if (Xs(n,2))~=0
    SnNRP(n,(1:Xs(n,2)-1))=0; 
    SnNRP(n,:)=circshift((SnNRP(n,:)),-(Xs(n,2)-1),2);
    end
end

end

if plotOpt
    
    figure(22)
    title('Nonlinear T21 Signal - Retimed - Real and Abs - Unwindowed')
    xlabel('Time (fs)')
    ylabel('Magnitude (arb.)')
    for k=1:3
%  plot(tnX, abs(SnNRP(k,:))*(NIFFT/NFFT),'o');
% plot((t-Xs(k,1)),imag(NRPrawplot(k,:)),'o',(t-Xs(k,1)),real(NRPrawplot(k,:)),'o')
hold on
% plot(tnX, imag(SnNRP(k,:))*(NIFFT/NFFT),tnX, real(SnNRP(k,:))*(NIFFT/NFFT))
plot(tnX, real(SnNRP(k,:))*(NIFFT/NFFT));

    end 
    
end 


% sizeRPpostShifts=size(SnRP)

% ******AT THIS POINT THE RP AND NRP HAVE BEEN PROPERLY INTERPOLATED*****
% AND ZEROD FOR X, NOW THEY MUST BE INTERPOLATED AND SHIFTED FOR Y BY THE
% AVERAGE Y SHIFT

% pre-phase the y dimension
yRPphases=zeros(1,size(SnRP,2));
yNRPphases=zeros(1,size(SnNRP,2)); 

for s=1:length(yRPphases)
    yRPphases(s)=angle(SnRP(1,s)); 
    yNRPphases(s)=angle(SnNRP(1,s)); 
end 
% Get the "raw" (pre second interp) data for comparison to post interp
NRP2raw=SnNRP;
RP2raw= SnRP;

% phase the y runs, 
SnRP=SnRP.*(exp(1i*(-1*yRPphases)));
SnNRP=SnNRP.*(exp(1i*(-1*yNRPphases)));

% START THE PROCESS OVER AGAIN FOR THE Y DIMENSION
SnNRP(1,:) = SnNRP(1,:)*0.5;
SnRP(1,:) = SnRP(1,:)*0.5;
% 
% NFFT
% sizeRPpreFFT=size(SnRP)
SFnRP = fftshift(fft(SnRP, NFFT, 1),1); % Signal in the frequency domain, FFT shift applied
SFnNRP = fftshift(fft(SnNRP, NFFT, 1),1); % Signal in the frequency domain, FFT shift applied

% apply the zeros to the new matrices, the correction dimensions should be
% Npad-1 (because Y runs are one element longer than X runs) by the size of 
%  the interpolated dimension along X
% Npad
% sizeRPpostFFT=size(SFnRP)

SFnRPzeros=[real(SFnRP) ; zeros(Npad,size(SFnRP,2))];
SFnNRPzeros=[real(SFnNRP) ; zeros(Npad,size(SFnNRP,2))];

% unshift the zero padded Mats by the right amount... MAY BE THAT EITHER X
% OR Y IS OFF BY ONE HERE due to the 10x9 size of the RP/NRP mats
unshiftednRP= circshift(SFnRPzeros, -floor(NFFT/2),1);
unshiftednNRP= circshift(SFnNRPzeros, -floor(NFFT/2),1);

RunshRPn = unshiftednRP;
RunshNRPn = unshiftednNRP;

Sn2RP = (ifft(RunshRPn, NIFFT, 1));
Sn2NRP = (ifft(RunshNRPn, NIFFT, 1));

Sn2RP(((ceil(NIFFT/2)+1):end),:)=0;
Sn2NRP(((ceil(NIFFT/2)+1):end),:)=0;

Sn2RP=Sn2RP*2;
Sn2NRP=Sn2NRP*2;

% reapply the y phases to get back to the original data phases
Sn2RP=Sn2RP.*(exp(1i*(yRPphases)));
Sn2NRP=Sn2NRP.*(exp(1i*(yNRPphases)));

% CunshRPn = complex(RunshRPn, imag(hilbert(-RunshRPn)));
% CunshNRPn = complex(RunshNRPn, imag(hilbert(-RunshNRPn)));
% size(RunshRPn,1)
% size(RunshRPn,2)
% NIFFTX
if size(RunshRPn,1)==size(RunshRPn,2) && size(RunshRPn,1)==NIFFTX
    NIFFT=size(RunshRPn,1);
else
    disp('2D data post interpolation has asymmetric number of indices'); 
    return;
end

if plotOpt
    figure(19)
    title('Y dimension - Interpolated - Imag Part'); 
    xlabel('Time (fs)'); 
    ylabel('Magnitude (arb)'); 
    for w=1:6
plot((ty)',imag(NRP2raw(:,w)),'o');
hold on
plot(tnY', imag(Sn2NRP(:,w))*(NIFFT/NFFT));
    end 

    figure(23)
    title('Y dimension - Raw - Imag Part'); 
    xlabel('Time (fs)'); 
    ylabel('Magnitude (arb)'); 
    for w=1:6
        plot((ty)',imag(NRPrawplot(:,w)));
        hold on
    end
end 

if plotOpt
figure(13)
contourf((real(RunshRPn)));

figure(17)
contourf((real(RunshNRPn)));
end

% Sn2RP = (ifft(CunshRPn, NIFFT, 1));
% Sn2NRP = (ifft(CunshNRPn, NIFFT, 1));
% 
% Sn2RP(1,:) = Sn2RP(1,:)*2;
% Sn2NRP(1,:) = Sn2NRP(1,:)*2;


% define an average Yshift to be applied

yAvg= round(mean(Ys(:,2)));

if retiming
% % zero out the matrices until the Yavg shift
if yAvg~=0
    Sn2RP((1:yAvg-1),:)=0;
    Sn2NRP((1:yAvg-1),:)=0;
    
    % % Now the Y shift can be applied as the average Y shift along all the
% % columns of the RP/NRP mats
Sn2RP=circshift(Sn2RP,-(yAvg-1),1);
Sn2NRP=circshift(Sn2NRP,-(yAvg-1),1);
    
end

end 

if plotOpt
    figure(24)
    title('Y dimension - Interpolated - Retimed - Imag Part'); 
    xlabel('Time (fs)'); 
    ylabel('Magnitude (arb)'); 
    for w=1:6
        plot(tnY', imag(Sn2NRP(:,w))*(NIFFT/NFFT));
        hold on;
    end 
    
    
end 

 


% ****AT THIS POINT BOTH AXES ARE INTERPOLATED AND RETIMED******
% NOW THE RP AND NRP MATS MUST BE PHASED.

% Y-data is the first index here
RPangle=angle(Sn2RP(1,1))
NRPangle=angle(Sn2NRP(1,1))
avgTheta = (angle(Sn2RP(1,1))+angle(Sn2NRP(1,1)))/2;

if rephase
Sn2RP = exp(1i*(-angle(Sn2RP(1,1))))*Sn2RP;
Sn2NRP = exp(1i*(-angle(Sn2NRP(1,1))))*Sn2NRP;
% Sn2RP = exp(1i*(pi+pi/10))*Sn2RP;
% Sn2NRP = exp(1i*(avgTheta))*Sn2NRP;
end

% scale back up by the right factor
Sn2RP= Sn2RP.*(NIFFT/NFFT);
Sn2NRP= Sn2NRP.*(NIFFT/NFFT);

% DECIMATE THE INTERPOLATED TIME DOMAIN DATA TO REDUCE THE SIZE OF THE
% MATRIX BY A FACTOR OF 2 ALONG EACH AXIS
% Sn2RPdec= decimate(Sn2RP,2);  
% Sn2NRPdec= decimate(Sn2NRP,2);  
% 
% figure(12)
% contourf(real(fftshift(fft2(Sn2RPdec))));
% axis equal
% title('Real part of RP Freq Decimated');
% 
% figure(13)
% contourf(real(fftshift(fft2(Sn2NRPdec))));
% axis equal
% title('Real part of NRP Freq Decimated');



% WINDOW THE DECIMATED DATA AND THE UNDECIMATED DATA TO COMPARE THE OUTPUTS
wYn= repmat(tnY.',1,size(Sn2NRP,1));
wXn= repmat(tnX,size(Sn2NRP,2),1);
% [wX, wY] = meshgrid(tx2D,ty);
window2Dec = delayedGaussian(sqrt(wXn.^2 + wYn.^2), 45, 10);


% % zero out the window beyond a certain point
% for k=1:size(window2Dec,1)
%     for w=1:size(window2Dec,2)
%         if abs(window2Dec(k,w))<(0.005)
%             window2Dec(k,w)=0;
%         end
%         
%     end 
%     
% end 


if plotOpt
figure(16)
contourf(real(Sn2RP.*window2Dec));
axis equal
title('Real part of RP - windowed');

figure(14)
contourf(real(Sn2NRP.*window2Dec));
axis equal
title('Real part of NRP - windowed');
end 
% freq axis for the transform of the NON zero padded interpolated time domain data to
% freq domain 
f2 = ceil((-NIFFT/2:NIFFT/2-1))/(Tsn*NIFFT); 

% ZERO PAD THE WINDOWED DECIMATED AND WINDOWED UNDECIMATED DATA PRIOR TO
% FINAL FFT 

Sn2RPzeros = [Sn2RP ; zeros(2^nextpow2(size(Sn2RP,1)),size(Sn2RP,2))];
Sn2RPzeros = [Sn2RPzeros , zeros(size(Sn2RPzeros,1),2^nextpow2(size(Sn2RP,1)))];

Sn2NRPzeros = [Sn2NRP ; zeros(2^nextpow2(size(Sn2NRP,1)),size(Sn2NRP,2))];
Sn2NRPzeros = [Sn2NRPzeros , zeros(size(Sn2NRPzeros,1),2^nextpow2(size(Sn2NRP,1)))];

NFFT2=size(Sn2RPzeros,1);
% new length of the time domain post-zero padding, increment is still the
% same
tn2 = ((0:NFFT2-1))*Tsn;
tnX2=tn2;
tnY2=tn2;

% a third window for the zero padded data
wXn2= repmat(tnX2.',1,size(Sn2NRPzeros,2));
wYn2= repmat(tnY2,size(Sn2NRPzeros,1),1);
% [wX, wY] = meshgrid(tx2D,ty);
window2Dfinal = delayedGaussian(sqrt(wXn2.^2 + wYn2.^2), 75, 45);

if plotOpt
figure(12)
contourf(window2Dfinal);
axis equal tight

figure(6)
contourf(real(Sn2RPzeros.*window2Dfinal));
axis equal tight
title('Real part of RP');

figure(7)
contourf(real(Sn2NRPzeros.*window2Dfinal));
axis equal
title('Real part of NRP');

figure(10)
contourf(real(fftshift(fft2(Sn2RPzeros.*window2Dfinal))));
axis equal
title('Real part of RP Freq');

figure(11)
contourf(real(fftshift(fft2(Sn2NRPzeros.*window2Dfinal))));
axis equal
title('Real part of NRP Freq');
end 

% % Second set of time axes for the newly zero padded interpolated data ONLY 
% % IF THE BACK TRANFORM OF THE FREQ DOMAIN IS TAKEN. 
% Tsn2 = (Tsn*NIFFT)/NIFFT2;
% tn2 = ((0:NIFFT2-1))*Tsn2;
% tnX2=tn2;
% tnY2=tn2;
Sn2NRPzeros=Sn2NRPzeros.';
Sn2RPzeros=Sn2RPzeros.';

NRPtimed= Sn2NRPzeros;
RPtimed= Sn2RPzeros;


end 