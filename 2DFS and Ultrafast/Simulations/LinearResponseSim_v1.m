% Linear Response and Spectral Analysis of rotating polarized pulses
% incident on a strongly coupled dimer. 

% Step wise Build up:
% 1. Develop the equations for the electric field of two pulses and
% establish their interference and phase effects 

%******* equation for the j-th pulse in the pulse train (starting with just 2 for now)*******
% IN TIME
% sigma = 50; % target is 30 fs decay
% a = 1/(2*(sigma.^2)); % input width to gaussian function

Polz=1; 
QWPmat=0;
simp=1;
useHandles=1;

%%
if ~Polz && ~simp
syms w
syms E1(t,t1,T,AO1,sig1,wL);
% syms E1test(t,T,AO1,sig1,wL);
syms E2(t,t2,T,AO2,sig2,wL);
% gauss = @(x, a) sqrt(a/pi)*exp(-a*x.^2); % general gaussian of std = 1/sqrt(2a);

% in this description T is the "lab time" or "absolute time" since it
% drives the evolution of the phase cycle. While t and t1/t2 are the
% arrival time coordinates of the pulses, which are not neccesarily
% evolving in the frame of the lab. For this reason, the time variable that
% should be allowed to evolve during a particular integration period - i.e.
% a fixed set of t1/t2 - is the time variable T
E1(t,t1,T,AO1,sig1,wL)= ((1/(sig1*sqrt(2*pi)))*exp(-((t-t1).^2)/(2*sig1^2))).*(exp(1i*(wL*(t-t1)+ 2*pi*AO1*T))); 
% E1test(t,T,AO1,sig1,wL)= ((1/(sig1*sqrt(2*pi)))*exp(-((t).^2)/(2*sig1^2))).*(exp(1i*(wL*(t)+ AO1*T))); 
E2(t,t2,T,AO2,sig2,wL)= ((1/(sig2*sqrt(2*pi)))*exp(-((t-t2).^2)/(2*sig2^2))).*(exp(1i*(wL*(t-t2)+ 2*pi*AO2*T))); 

taxis=linspace(-200e-15,200e-15,1000);
t1a= 0e-15; 
t2a= 20e-15;
sig1a= 40e-15;
sig2a= 40e-15;
% 532nm is the same as 563THz
wLa=563e12;
% wLa=530e12;

% AO1 and T is going to play the role of the AOM phase, Where AO1 is the
% freq of imparted by the AOM and T is the time within the modulaiton
% cycle. Notably, the t-space for the pulse arrivals and the T-space for
% the AOM cycle are not interdependent at all. This will manifest itself in
% the transform. 
AO1a=1.005e6; 
AO2a=1.0e6;
Ta=0e-4;
Taxis=linspace(0,2e-4,15);


figure(1)
plot(taxis,real(E1(taxis,t1a,Ta,AO1a,sig1a,wLa)+E2(taxis,t2a,Ta,AO2a,sig2a,wLa)),taxis,imag(E1(taxis,t1a,Ta,AO1a,sig1a,wLa)+E2(taxis,t2a,Ta,AO2a,sig2a,wLa)));
title('Abs E1 + E2 - Time domain');
xlabel('Time (seconds)');
ylabel('Amp (arb)');

% attempting to compute the FT of the pulse time profile
assume([t1,T,AO1,sig1,wL],'real');
E1w=fourier(E1,t,w);
E1w_simp=simplify(E1w); 
assume([t1,T,AO1,sig1,wL],'clear');
% syms E1w_final(w,t1,T,AO1,sig1,wL);
E1w_final(w,t1,T,AO1,sig1,wL)=E1w; 
E1w_finalS(w,t1,T,AO1,sig1,wL)=E1w_simp; 

assume([t2,T,AO2,sig2,wL],'real');
E2w=fourier(E2,t,w);
E2w_simp=simplify(E2w); 
assume([t2,T,AO2,sig2,wL],'clear');
% syms E2w_final(w,t2,T,AO2,sig2,wL);
E2w_final(w,t2,T,AO2,sig2,wL)=E2w; 
E2w_finalS(w,t2,T,AO2,sig2,wL)=E2w_simp; 

wAxis=linspace(430e12,750e12,1000); 
figure(2)
plot(wAxis,real(E1w_finalS(wAxis,t1a,Ta,AO1a,sig1a,wLa)),wAxis,imag(E1w_finalS(wAxis,t1a,Ta,AO1a,sig1a,wLa)));
hold on
plot(wAxis,abs(E1w_finalS(wAxis,t1a,Ta,AO1a,sig1a,wLa)));
title('E1 - Freq domain');
xlabel('Frequency (Hz)');
ylabel('Amp (arb)');
legend('real part','imag part','Abs'); 

figure(3)
plot(wAxis,real(E2w_finalS(wAxis,t2a,Ta,AO2a,sig2a,wLa)),wAxis,imag(E2w_finalS(wAxis,t2a,Ta,AO2a,sig2a,wLa)));
hold on
plot(wAxis,abs(E2w_finalS(wAxis,t2a,Ta,AO2a,sig2a,wLa)));
title('E2 - Freq domain');
xlabel('Frequency (Hz)');
ylabel('Amp (arb)');
legend('real part','imag part','Abs');


Etot=E1(t,t1,T,AO1,sig1,wL)+E2(t,t2,T,AO2,sig2,wL);
% compute the spectral profile of the combined 2 pulse train
% *******This moethod works thus far. However, it seems that in order to
% capture the phase sweep, the product of these two fields must be taken
% (to combine their complex phase factors). Despite the transform accepting
% this sum of E1 and E2, it does not appear to combine the fields in the
% same way as the reference articles. Something to note....
assume([t1,T,AO1,sig1,wL,t2,AO2,sig2],'real');
EwT=fourier(Etot,t,w);
Ew_simpT=simplify(EwT); 
assume([t1,T,AO1,sig1,wL,t2,AO2,sig2],'clear');
% syms E1w_final(w,T,AO1,sig1,wL,t2,AO2,sig2)
Ew_finalTot(w,t1,T,AO1,sig1,wL,t2,AO2,sig2)=EwT; 
Ew_finalSTot(w,t1,T,AO1,sig1,wL,t2,AO2,sig2)=Ew_simpT; 

wAxis=linspace(430e12,750e12,1000); 
figure(4)
plot(wAxis,real(Ew_finalSTot(wAxis,t1a,Ta,AO1a,sig1a,wLa,t2a,AO2a,sig2a)),wAxis,imag(Ew_finalSTot(wAxis,t1a,Ta,AO1a,sig1a,wLa,t2a,AO2a,sig2a)));
hold on
plot(wAxis,abs(Ew_finalSTot(wAxis,t1a,Ta,AO1a,sig1a,wLa,t2a,AO2a,sig2a)));
title('Etotal = E1 + E2 - Freq domain');
xlabel('Frequency (Hz)');
ylabel('Amp (arb)');
legend('real part','imag part','Abs');
%%

% Now to create a manifold for the sym/anti transistions of the dimer.
% Overlap it with the two fields, and then consider their inner products.
% See Appendix Eq.9,10 in WPI 2006 artcile 
syms symm(w,w00,w01A,w01B,C00,C01A,C01B,gamma);
syms anti(w,w00,w01A,w01B,C00,C01A,C01B,gamma);

% % Sum of Gaussian curves
symm(w,w00,w01A,w01B,C00,C01A,C01B,gamma)= ...
    C00*(1/(gamma*sqrt(2*pi)))*exp(-((w-w00).^2)/(2*gamma^2)) + ...
    C01A*(1/(gamma*sqrt(2*pi)))*exp(-((w-w01A).^2)/(2*gamma^2)) + ...
    C01B*(1/(gamma*sqrt(2*pi)))*exp(-((w-w01B).^2)/(2*gamma^2));

anti(w,w00,w01A,w01B,C00,C01A,C01B,gamma)= ...
    C00*(1/(gamma*sqrt(2*pi)))*exp(-((w-w00).^2)/(2*gamma^2)) + ...
    C01A*(1/(gamma*sqrt(2*pi)))*exp(-((w-w01A).^2)/(2*gamma^2)) + ...
    C01B*(1/(gamma*sqrt(2*pi)))*exp(-((w-w01B).^2)/(2*gamma^2));

% THESE VALUES REPRESENT A DUPLEX DIMER AT 15DEG C (OBTAINED FROM
% LIN/NONLIN FIT OF THE DATA)
% assign the center frequncies of the lorenztians 
w00s=536.63e12;
w01As=570.12e12;
w01Bs=579.74e12;

w00a=552.34e12;
w01Aa=584.36e12;
w01Ba=585.97e12;
% assign the frequncy widths of the lorenztians (should be in THz)
gammaA=10.07e12;
% gammaAa=10e12;
% assign the peak heights
C00s= 0.0934;
C01As= 0.0252;
C01Bs= 0.0063;

C00a= 0.0556;
C01Aa= 0.0771;
C01Ba= 0.0194;


figure(4)
plot(wAxis,symm(wAxis,w00s,w01As,w01Bs,C00s,C01As,C01Bs,gammaA)./max(symm(wAxis,w00s,w01As,w01Bs,C00s,C01As,C01Bs,gammaA)),wAxis,anti(wAxis,w00a,w01Aa,w01Ba,C00a,C01Aa,C01Ba,gammaA)./max(anti(wAxis,w00a,w01Aa,w01Ba,C00a,C01Aa,C01Ba,gammaA)));
hold on;
plot(wAxis,abs(Ew_finalSTot(wAxis,t1a,Ta,AO1a,sig1a,wLa,t2a,AO2a,sig2a)));
plot(wAxis,abs(E1w_finalS(wAxis,t1a,Ta,AO1a,sig1a,wLa)));
hold off;
title('Pure Molecular Spectra - Freq domain');
xlabel('Frequency (Hz)');
ylabel('Amp (arb)');
legend('Symm Abs Spec','Anti Abs Spec','E1+E2 Laser Spectrum','E1 Spectrum'); 

% test an overlap calucation between the field and either the symm/anti
% manifold
% E1wSymm(w,t1,T,AO1,sig1,wL,w0,gamma)= (int((E1w_finalS(w,t1,T,AO1,sig1,wL).*symm(w,w0,gamma)),w)); 
% E1wSymm_Simp=simplify(E1wSymm);

% try overlapping the combined fields with the molecualr spectrum to get
% the fianl sigal out
% EtotSymm(w,t1,T,AO1,sig1,wL,t2,AO2,sig2,w0,gamma)= int(conj(Ew_finalSTot(w,t1,T,AO1,sig1,wL,t2,AO2,sig2).*symm(w,w0,gamma)).*(Ew_finalSTot(w,t1,T,AO1,sig1,wL,t2,AO2,sig2).*symm(w,w0,gamma)),w);
% EtotSymm_Simp=simplify(EtotSymm);
symmSpec= symm(wAxis,w00s,w01As,w01Bs,C00s,C01As,C01Bs,gammaA);
antiSpec= anti(wAxis,w00a,w01Aa,w01Ba,C00a,C01Aa,C01Ba,gammaA);
% Overlap of the two fields with the molecualr spectra.
% See Appendix Eq.9,10 in WPI 2006 article - should be equivalent to Eqn 1b
% in Tamimi 2020 as well. 
sigTot= (conj(E1w_finalS(wAxis,t1a,Ta,AO1a,sig1a,wLa)).*E1w_finalS(wAxis,t1a,Ta,AO1a,sig1a,wLa).*((symmSpec+antiSpec).^2)) + ...
        (conj(E2w_finalS(wAxis,t2a,Ta,AO2a,sig2a,wLa)).*E2w_finalS(wAxis,t2a,Ta,AO2a,sig2a,wLa).*((symmSpec+antiSpec).^2)) + ...
        2*real((conj(E2w_finalS(wAxis,t2a,Ta,AO2a,sig2a,wLa)).*E1w_finalS(wAxis,t1a,Ta,AO1a,sig1a,wLa).*((symmSpec+antiSpec).^2))); 

interfComp=2*real((conj(E2w_finalS(wAxis,t2a,Ta,AO2a,sig2a,wLa)).*E1w_finalS(wAxis,t1a,Ta,AO1a,sig1a,wLa).*((symmSpec+antiSpec).^2)));
figure(5)
plot(wAxis,sigTot,wAxis,abs(interfComp));
title('Total Signal - Freq domain');
xlabel('Frequency (Hz)');
ylabel('Amp (arb)');
legend('Total Signal','Interference Component');

numPlots=length(Taxis);
    DimRow= floor(sqrt(numPlots)); 
    DimCol= floor(sqrt(numPlots)); 
    
    if ((DimRow*DimCol)< numPlots)
        DimRow=DimRow+1;
        if(DimRow*DimCol)< numPlots
            DimCol=DimCol+1; 
        end
    end
%     create a vector to hold the magnitudes of the polz components over
%     the mod cycle
    innerProdComp = zeros(2,numel(Taxis)); 
    
   
    for j=1:numPlots
         demodCompTot=((conj(E2w_finalS(wAxis,t2a,Taxis(j),AO2a,sig2a,wLa)).*E1w_finalS(wAxis,t1a,Taxis(j),AO1a,sig1a,wLa).*((symmSpec+antiSpec).^2)));
         demodCompSymm=((conj(E2w_finalS(wAxis,t2a,Taxis(j),AO2a,sig2a,wLa)).*E1w_finalS(wAxis,t1a,Taxis(j),AO1a,sig1a,wLa).*((symmSpec).^2)));
         demodCompAnti=((conj(E2w_finalS(wAxis,t2a,Taxis(j),AO2a,sig2a,wLa)).*E1w_finalS(wAxis,t1a,Taxis(j),AO1a,sig1a,wLa).*((antiSpec).^2)));

         innerProdSym= ((E2w_finalS(wAxis,t2a,Taxis(j),AO2a,sig2a,wLa)).*(symmSpec))*((E1w_finalS(wAxis,t1a,Taxis(j),AO1a,sig1a,wLa)').*(symmSpec'));
         innerProdAnti= ((E2w_finalS(wAxis,t2a,Taxis(j),AO2a,sig2a,wLa)).*(antiSpec))*((E1w_finalS(wAxis,t1a,Taxis(j),AO1a,sig1a,wLa)').*(antiSpec'));

        
         figure(6)
        sgtitle(['Real Overlap Spectra - Freq domain']); 
        subplot(DimRow,DimCol,(j));
        plot(wAxis,real(demodCompSymm),wAxis,real(demodCompAnti));
        title(['AOM Phase=' num2str((Taxis(j)./(Taxis(end)))*2*pi) 'radians']);
        xlabel('Frequency (Hz)');
        ylabel('Amp (arb)');
        ylim([-8e-30,8e-30]);
        legend('Symm Component','Anti Component');
        
        
        figure(7)
        sgtitle(['Imag Overlap Spectra - Freq domain']); 
        subplot(DimRow,DimCol,(j));
        plot(wAxis,imag(demodCompSymm),wAxis,imag(demodCompAnti));
        title(['AOM Phase=' num2str((Taxis(j)./(Taxis(end)))*2*pi) 'radians']);
        xlabel('Frequency (Hz)');
        ylabel('Amp (arb)');
        ylim([-8e-30,8e-30]);
        legend('Symm Component','Anti Component');
        
        innerProdComp(1,j)=innerProdSym;
        innerProdComp(2,j)=innerProdAnti;
        
    end
    
    figure(8)
    plot((Taxis),real(innerProdComp(1,:)),(Taxis),imag(innerProdComp(1,:))); 
    hold on
    plot((Taxis),real(innerProdComp(2,:)),(Taxis),imag(innerProdComp(2,:)));
    title('Inner Prod components over AOM Phase Cycle');
    xlabel('AOM Clock Cycle (sec)');
    ylabel('Mag (arb)');
    legend('Symm Real','Symm Imag','Anti Real','Anti Imag');
    

end

%%
% Same calculation as before, but now performed with Cross polarized pulses
% setting up a 2x1 vector to overlap with the cross polarized Symm/anti
% manifolds

if Polz && ~QWPmat && ~simp
 syms w
syms E1(t,t1,T,AO1,AO2,sig1,wL);
% syms E1test(t,T,AO1,sig1,wL);
syms E2(t,t2,T,AO1,AO2,sig2,wL);
% gauss = @(x, a) sqrt(a/pi)*exp(-a*x.^2); % general gaussian of std = 1/sqrt(2a);

% in this description T is the "lab time" or "absolute time" since it
% drives the evolution of the phase cycle. While t and t1/t2 are the
% arrival time coordinates of the pulses, which are not neccesarily
% evolving in the frame of the lab. For this reason, the time variable that
% should be allowed to evolve during a particular integration period - i.e.
% a fixed set of t1/t2 - is the time variable T
% WE WILL ASSUME E1 IS POLARIZED ALONG V AND E2 ALOING H 
E1(t,t1,T,AO1,AO2,sig1,wL)= ((1/(sig1*sqrt(2*pi)))*exp(-((t-t1).^2)/(2*sig1^2))).*(exp(1i*(wL*(t-t1)))).*exp(((1i*2*pi*(AO1-AO2)*T)/2) - pi/4).*cos((2*pi*((AO1-AO2)*T)/2) - pi/4); 
% E1test(t,T,AO1,sig1,wL)= ((1/(sig1*sqrt(2*pi)))*exp(-((t).^2)/(2*sig1^2))).*(exp(1i*(wL*(t)+ AO1*T))); 
E2(t,t2,T,AO1,AO2,sig2,wL)= ((1/(sig2*sqrt(2*pi)))*exp(-((t-t2).^2)/(2*sig2^2))).*(exp(1i*(wL*(t-t2)))).*exp(((1i*2*pi*(AO1-AO2)*T)/2) - pi/4).*(-sin((2*pi*((AO1-AO2)*T)/2) - pi/4)); 

taxis=linspace(-200e-15,200e-15,1000);
t1a= 0e-15; 
t2a= 20e-15;
sig1a= 40e-15;
sig2a= 40e-15;
% 532nm is the same as 563THz
wLa=563e12;
% wLa=530e12;

% AO1 and T is going to play the role of the AOM phase, Where AO1 is the
% freq of imparted by the AOM and T is the time within the modulaiton
% cycle. Notably, the t-space for the pulse arrivals and the T-space for
% the AOM cycle are not interdependent at all. This will manifest itself in
% the transform. 
AO1a=1.005e6; 
AO2a=1.0e6;
Ta=0e-4;
% another time axis for the modulation cycle evolution
Taxis=linspace(0,2e-4,6);


figure(1)
plot(taxis,real(E1(taxis,t1a,Ta,AO1a,AO2a,sig1a,wLa)+E2(taxis,t2a,Ta,AO1a,AO2a,sig2a,wLa)));
hold on
plot(taxis,imag(E1(taxis,t1a,Ta,AO1a,AO2a,sig1a,wLa)+E2(taxis,t2a,Ta,AO1a,AO2a,sig2a,wLa)));
title('E1 + E2 - Time domain - assumes proj onto 45 deg axis');
xlabel('Time (seconds)');
ylabel('Amp (arb)');
legend('real part', 'imag part'); 

% attempting to compute the FT of the pulse time profile for both E1 and E2
% - The cross polarization can be handled afterwards, with E1w and E2w
% being stored in orthogonal vectors 
assume([t1,T,AO1,AO2,sig1,wL],'real');
E1w=fourier(E1,t,w);
E1w_simp=simplify(E1w); 
assume([t1,T,AO1,AO2,sig1,wL],'clear');
% syms E1w_final(w,t1,T,AO1,sig1,wL);
E1w_final(w,t1,T,AO1,AO2,sig1,wL)=E1w; 
E1w_finalS(w,t1,T,AO1,AO2,sig1,wL)=E1w_simp; 

assume([t2,T,AO1,AO2,sig2,wL],'real');
E2w=fourier(E2,t,w);
E2w_simp=simplify(E2w); 
assume([t2,T,AO1,AO2,sig2,wL],'clear');
% syms E2w_final(w,t2,T,AO2,sig2,wL);
E2w_final(w,t2,T,AO1,AO2,sig2,wL)=E2w; 
E2w_finalS(w,t2,T,AO1,AO2,sig2,wL)=E2w_simp; 

wAxis=linspace(430e12,750e12,1000); 
figure(2)
plot(wAxis,real(E1w_finalS(wAxis,t1a,Ta,AO1a,AO2a,sig1a,wLa)),wAxis,imag(E1w_finalS(wAxis,t1a,Ta,AO1a,AO2a,sig1a,wLa)));
hold on
plot(wAxis,abs(E1w_finalS(wAxis,t1a,Ta,AO1a,AO2a,sig1a,wLa)));
title('E1 - Freq domain');
xlabel('Frequency (Hz)');
ylabel('Amp (arb)');
legend('real part','imag part','Abs'); 

figure(3)
plot(wAxis,real(E2w_finalS(wAxis,t2a,Ta,AO1a,AO2a,sig2a,wLa)),wAxis,imag(E2w_finalS(wAxis,t2a,Ta,AO1a,AO2a,sig2a,wLa)));
hold on
plot(wAxis,abs(E2w_finalS(wAxis,t2a,Ta,AO1a,AO2a,sig2a,wLa)));
title('E2 - Freq domain');
xlabel('Frequency (Hz)');
ylabel('Amp (arb)');
legend('real part','imag part','Abs'); 

% FOR THIS CROSS POLARIZED CASE NO ETOTAL SHOULD BE COSNDIERED, AS THE TWO
% FIELDS ARE ORTHOGONAL AND WILL NOT PRODUCE THE SAME OUTCOME AS TWO
% MUTUALLY POLARIZED FIELDS 


%%

% Now to create a manifold for the sym/anti transistions of the dimer.
% Overlap it with the two fields, and then consider their inner products.
% See Appendix Eq.9,10 in WPI 2006 artcile 
% syms symm(w,w0,gamma);
% syms anti(w,w0,gamma);

syms symm(w,w00,w01A,w01B,C00,C01A,C01B,gamma);
syms anti(w,w00,w01A,w01B,C00,C01A,C01B,gamma);

% % Lorentzian curves
% symm(w,w0,gamma)= 1./((pi*gamma)*(1+((w-w0)/gamma).^2));  
% anti(w,w0,gamma)= (1/(pi*gamma))*((gamma^2)/((w-w0).^2 + gamma^2));

% % Gaussian curves
% symm(w,w0,gamma)=(1/(gamma*sqrt(2*pi)))*exp(-((w-w0).^2)/(2*gamma^2));
% anti(w,w0,gamma)=(1/(gamma*sqrt(2*pi)))*exp(-((w-w0).^2)/(2*gamma^2));

% % Sum of Gaussian curves
symm(w,w00,w01A,w01B,C00,C01A,C01B,gamma)= ...
    C00*(1/(gamma*sqrt(2*pi)))*exp(-((w-w00).^2)/(2*gamma^2)) + ...
    C01A*(1/(gamma*sqrt(2*pi)))*exp(-((w-w01A).^2)/(2*gamma^2)) + ...
    C01B*(1/(gamma*sqrt(2*pi)))*exp(-((w-w01B).^2)/(2*gamma^2));

anti(w,w00,w01A,w01B,C00,C01A,C01B,gamma)= ...
    C00*(1/(gamma*sqrt(2*pi)))*exp(-((w-w00).^2)/(2*gamma^2)) + ...
    C01A*(1/(gamma*sqrt(2*pi)))*exp(-((w-w01A).^2)/(2*gamma^2)) + ...
    C01B*(1/(gamma*sqrt(2*pi)))*exp(-((w-w01B).^2)/(2*gamma^2));


% % assign the center frequncies of the lorenztians 
% w0sa=530e12;
% w0aa=580e12;
% % assign the frequncy widths of the lorenztians (should be in THz)
% gammaSa=10e12;
% gammaAa=10e12;

% THESE VALUES REPRESENT A DUPLEX DIMER AT 15DEG C (OBTAINED FROM
% LIN/NONLIN FIT OF THE DATA)
% assign the center frequncies of the lorenztians 
w00s=536.63e12;
w01As=570.12e12;
w01Bs=579.74e12;

w00a=552.34e12;
w01Aa=584.36e12;
w01Ba=585.97e12;
% assign the frequncy widths of the lorenztians (should be in THz)
gammaA=10.07e12;
% gammaAa=10e12;
% assign the peak heights
C00s= 0.0934;
C01As= 0.0252;
C01Bs= 0.0063;

C00a= 0.0556;
C01Aa= 0.0771;
C01Ba= 0.0194;

% WE WILL ASSUME THAT THE SYMM MANIFOLD IS POLARIZED ALONG V AND THE ANTI
% ALONG H


figure(4)
plot(wAxis,symm(wAxis,w00s,w01As,w01Bs,C00s,C01As,C01Bs,gammaA)./max(symm(wAxis,w00s,w01As,w01Bs,C00s,C01As,C01Bs,gammaA)),wAxis,anti(wAxis,w00a,w01Aa,w01Ba,C00a,C01Aa,C01Ba,gammaA)./max(anti(wAxis,w00a,w01Aa,w01Ba,C00a,C01Aa,C01Ba,gammaA)));
hold on;
plot(wAxis,abs(E1w_finalS(wAxis,t1a,Ta,AO1a,AO2a,sig1a,wLa)));
hold off;
title('Pure Molecular Spectra - Freq domain');
xlabel('Frequency (Hz)');
ylabel('Amp (arb)');
legend('Symm Abs Spec','Anti Abs Spec','Laser Spectrum'); 

% % Cast the symm/anti as a 2x1 vector, as well as the field. Determine the
% overlap of the two.
moleSpec=[symm(wAxis,w00s,w01As,w01Bs,C00s,C01As,C01Bs,gammaA) ; anti(wAxis,w00a,w01Aa,w01Ba,C00a,C01Aa,C01Ba,gammaA)];
 
% Overlap of the two fields with the molecualr spectra.
% See Appendix Eq.9,10 in WPI 2006 article - should be equivalent to Eqn 1b
% in Tamimi 2020 as well. 
% sigTot= (conj(E1w_finalS(wAxis,t1a,Ta,AO1a,sig1a,wLa)).*E1w_finalS(wAxis,t1a,Ta,AO1a,sig1a,wLa).*((symm(wAxis,w0sa,gammaSa)+anti(wAxis,w0aa,gammaAa)).^2)) + ...
%         (conj(E2w_finalS(wAxis,t2a,Ta,AO2a,sig2a,wLa)).*E2w_finalS(wAxis,t2a,Ta,AO2a,sig2a,wLa).*((symm(wAxis,w0sa,gammaSa)+anti(wAxis,w0aa,gammaAa)).^2)) + ...
%         2*real((conj(E2w_finalS(wAxis,t2a,Ta,AO2a,sig2a,wLa)).*E1w_finalS(wAxis,t1a,Ta,AO1a,sig1a,wLa).*((symm(wAxis,w0sa,gammaSa)+anti(wAxis,w0aa,gammaAa)).^2))); 

% interfComp=2*real((conj(E2w_finalS(wAxis,t2a,Ta,AO2a,sig2a,wLa)).*E1w_finalS(wAxis,t1a,Ta,AO1a,sig1a,wLa).*((symm(wAxis,w0sa,gammaSa)+anti(wAxis,w0aa,gammaAa)).^2)));

% figure(5)
% plot(wAxis,overlapSpec(1,:),wAxis,overlapSpec(2,:));
% title('Overlap Spectra - Freq domain');
% xlabel('Frequency (Hz)');
% ylabel('Amp (arb)');
% legend('Symm Component','Anti Component');   
    
numPlots=length(Taxis);
    DimRow= floor(sqrt(numPlots)); 
    DimCol= floor(sqrt(numPlots)); 
    
    if ((DimRow*DimCol)< numPlots)
        DimRow=DimRow+1;
        if(DimRow*DimCol)< numPlots
            DimCol=DimCol+1; 
        end
    end
%     create a vector to hold the magnitudes of the polz components over
%     the mod cycle
    polzComp = zeros(2,numel(Taxis)); 
    innerProdComp = zeros(2,numel(Taxis));
    
   
    for j=1:numPlots
%         *******pulses polarized at 45 degrees sich that their projections
%         are equal... This scheme assumes the product of the two fields
%         onto the dipole captures the interference. 
%         fieldComp=[conj(E2w_finalS(wAxis,t2a,Taxis(j),AO1a,AO2a,sig2a,wLa)).*(E1w_finalS(wAxis,t1a,Taxis(j),AO1a,AO2a,sig1a,wLa)) ; ...
%             conj(E2w_finalS(wAxis,t2a,Taxis(j),AO1a,AO2a,sig2a,wLa)).*(E1w_finalS(wAxis,t1a,Taxis(j),AO1a,AO2a,sig1a,wLa))];

%        ******pulses polarized at an angle of 45 degrees relative to the dimer
%        axes. Such that one pulse (E1) has postive x and y component, but
%        E2 has a postive x and negative y component. Tryign here to modle
%        the total field as the sum, lettign the product be taken in the
%        overlap spectra
            fieldComp=[(1/sqrt(2))*(E1w_finalS(wAxis,t1a,Taxis(j),AO1a,AO2a,sig1a,wLa) - E2w_finalS(wAxis,t2a,Taxis(j),AO1a,AO2a,sig2a,wLa)) ...
            ; (1/sqrt(2))*(E1w_finalS(wAxis,t1a,Taxis(j),AO1a,AO2a,sig1a,wLa) + E2w_finalS(wAxis,t2a,Taxis(j),AO1a,AO2a,sig2a,wLa))];

        % for now the overlap spectra should contain everything - would be ideal to
        % ioslate the interference components as well 
        overlapSpec= (fieldComp.*moleSpec).*conj(fieldComp.*moleSpec);
%         overlapSpec= (fieldComp.*moleSpec).*conj(moleSpec);
        

%          overlapSpec2= ((fieldComp')*(fieldComp))*((moleSpec')*(moleSpec))';
%          innerProdSym= ((E2w_finalS(wAxis,t2a,Taxis(j),AO1a,AO2a,sig2a,wLa)).*(moleSpec(1,:)))*((E1w_finalS(wAxis,t1a,Taxis(j),AO1a,AO2a,sig1a,wLa)').*(moleSpec(1,:)'));
%          innerProdAnti= ((E2w_finalS(wAxis,t2a,Taxis(j),AO1a,AO2a,sig2a,wLa)).*(moleSpec(2,:)))*((E1w_finalS(wAxis,t1a,Taxis(j),AO1a,AO2a,sig1a,wLa)').*(moleSpec(2,:)'));
%          
         innerProdSym= ((fieldComp(1,:)).*(moleSpec(1,:)))*((fieldComp(1,:)').*(moleSpec(1,:)'));
         innerProdAnti= ((fieldComp(2,:)).*(moleSpec(2,:)))*((fieldComp(2,:)').*(moleSpec(2,:)'));
         
         innerProdComp(1,j)=innerProdSym;
         innerProdComp(2,j)=innerProdAnti;
        
         figure(6)
        sgtitle(['Real Overlap Spectra - Freq domain']); 
        subplot(DimRow,DimCol,(j));
        plot(wAxis,real(overlapSpec(1,:)),wAxis,real(overlapSpec(2,:)));
        title(['AOM Phase=' num2str((Taxis(j)./(Taxis(end)))*2*pi) 'radians']);
        xlabel('Frequency (Hz)');
        ylabel('Amp (arb)');
        legend('Symm Component','Anti Component');
        
        
        figure(7)
        sgtitle(['Imag Overlap Spectra - Freq domain']); 
        subplot(DimRow,DimCol,(j));
        plot(wAxis,imag(overlapSpec(1,:)),wAxis,imag(overlapSpec(2,:)));
        title(['AOM Phase=' num2str((Taxis(j)./(Taxis(end)))*2*pi) 'radians']);
        xlabel('Frequency (Hz)');
        ylabel('Amp (arb)');
        legend('Symm Component','Anti Component');
        
        polzComp(1,j)=sum(abs(fieldComp(1,:)));
        polzComp(2,j)=sum(abs(fieldComp(2,:)));
        
    end
    
    figure(8)
        plot(Taxis,(polzComp(1,:)),Taxis,(polzComp(2,:)));
        title(['Polarization Component Magnitudes']);
        xlabel('AOM Phase');
        ylabel('Mag (arb)');
        legend('Vert Component','Horz Component');
    
         figure(9)
        plot(Taxis,real((innerProdComp(1,:))),Taxis,imag((innerProdComp(1,:))));
        hold on
        plot(Taxis,real((innerProdComp(2,:))),Taxis,imag((innerProdComp(2,:))));
        title(['Inner Prod Components']);
        xlabel('AOM Clock Cycle (sec)');
        ylabel('Mag (arb)');
        legend('Symm Real','Symm Imag','Anti Real','Anti Imag');
        
        
%          figure(9)
%         plot(Taxis,sqrt(real(polzComp(1,:)).^2 + imag(polzComp(1,:)).^2),Taxis,sqrt(real(polzComp(2,:)).^2 + imag(polzComp(2,:)).^2));
%         title(['Polarization Component Magnitudes']);
%         xlabel('AOM Phase');
%         ylabel('Mag (arb)');
%         legend('Vert Component','Horz Component');
    
end

if Polz && QWPmat && ~simp   

% ALTERNATE APPROACH TO CROSS POLARIZED EFFECTS - KEEP THE FIELDS THE SAME
% AS BEFORE (IN THE PARALLEL CASE) THEN SIMPLY APPLY THE MATRIX OPERATION
% OF THE QWP
syms w
syms E1(t,t1,T,AO1,sig1,wL);
% syms E1test(t,T,AO1,sig1,wL);
syms E2(t,t2,T,AO2,sig2,wL);
% gauss = @(x, a) sqrt(a/pi)*exp(-a*x.^2); % general gaussian of std = 1/sqrt(2a);

% in this description T is the "lab time" or "absolute time" since it
% drives the evolution of the phase cycle. While t and t1/t2 are the
% arrival time coordinates of the pulses, which are not neccesarily
% evolving in the frame of the lab. For this reason, the time variable that
% should be allowed to evolve during a particular integration period - i.e.
% a fixed set of t1/t2 - is the time variable T
E1(t,t1,T,AO1,sig1,wL)= ((1/(sig1*sqrt(2*pi)))*exp(-((t-t1).^2)/(2*sig1^2))).*(exp(1i*(wL*(t-t1)+ 2*pi*AO1*T)))*( (1/2)*exp(-1i*pi/4)*((1+1i) + (1-1i)) ); 
% E1test(t,T,AO1,sig1,wL)= ((1/(sig1*sqrt(2*pi)))*exp(-((t).^2)/(2*sig1^2))).*(exp(1i*(wL*(t)+ AO1*T))); 
E2(t,t2,T,AO2,sig2,wL)= ((1/(sig2*sqrt(2*pi)))*exp(-((t-t2).^2)/(2*sig2^2))).*(exp(1i*(wL*(t-t2)+ 2*pi*AO2*T)))*( (1/2)*exp(-1i*pi/4)*((1+1i) + (1-1i)) ); 

% cant use syms to do matrix operation, so simpyl going to apply the
% results directly above^^


taxis=linspace(-200e-15,200e-15,1000);
t1a= 0e-15; 
t2a= 20e-15;
sig1a= 40e-15;
sig2a= 40e-15;
% 532nm is the same as 563THz
wLa=563e12;
% wLa=530e12;

% AO1 and T is going to play the role of the AOM phase, Where AO1 is the
% freq of imparted by the AOM and T is the time within the modulaiton
% cycle. Notably, the t-space for the pulse arrivals and the T-space for
% the AOM cycle are not interdependent at all. This will manifest itself in
% the transform. 
AO1a=1.005e6; 
AO2a=1.0e6;
Ta=0e-4;
Taxis=linspace(0,2e-4,6);


figure(1)
plot(taxis,real(E1(taxis,t1a,Ta,AO1a,sig1a,wLa)+E2(taxis,t2a,Ta,AO2a,sig2a,wLa)),taxis,imag(E1(taxis,t1a,Ta,AO1a,sig1a,wLa)+E2(taxis,t2a,Ta,AO2a,sig2a,wLa)));
title('Abs E1 + E2 - Time domain');
xlabel('Time (seconds)');
ylabel('Amp (arb)');

% attempting to compute the FT of the pulse time profile
assume([t1,T,AO1,sig1,wL],'real');
E1w=fourier(E1,t,w);
E1w_simp=simplify(E1w); 
assume([t1,T,AO1,sig1,wL],'clear');
% syms E1w_final(w,t1,T,AO1,sig1,wL);
E1w_final(w,t1,T,AO1,sig1,wL)=E1w; 
E1w_finalS(w,t1,T,AO1,sig1,wL)=E1w_simp; 

assume([t2,T,AO2,sig2,wL],'real');
E2w=fourier(E2,t,w);
E2w_simp=simplify(E2w); 
assume([t2,T,AO2,sig2,wL],'clear');
% syms E2w_final(w,t2,T,AO2,sig2,wL);
E2w_final(w,t2,T,AO2,sig2,wL)=E2w; 
E2w_finalS(w,t2,T,AO2,sig2,wL)=E2w_simp; 

wAxis=linspace(430e12,750e12,1000); 
figure(2)
plot(wAxis,real(E1w_finalS(wAxis,t1a,Ta,AO1a,sig1a,wLa)),wAxis,imag(E1w_finalS(wAxis,t1a,Ta,AO1a,sig1a,wLa)));
hold on
plot(wAxis,abs(E1w_finalS(wAxis,t1a,Ta,AO1a,sig1a,wLa)));
title('E1 - Freq domain');
xlabel('Frequency (Hz)');
ylabel('Amp (arb)');
legend('real part','imag part','Abs'); 

figure(3)
plot(wAxis,real(E2w_finalS(wAxis,t2a,Ta,AO2a,sig2a,wLa)),wAxis,imag(E2w_finalS(wAxis,t2a,Ta,AO2a,sig2a,wLa)));
hold on
plot(wAxis,abs(E2w_finalS(wAxis,t2a,Ta,AO2a,sig2a,wLa)));
title('E2 - Freq domain');
xlabel('Frequency (Hz)');
ylabel('Amp (arb)');
legend('real part','imag part','Abs');

%%

% Now to create a manifold for the sym/anti transistions of the dimer.
% Overlap it with the two fields, and then consider their inner products.
% See Appendix Eq.9,10 in WPI 2006 artcile 
syms symm(w,w00,w01A,w01B,C00,C01A,C01B,gamma);
syms anti(w,w00,w01A,w01B,C00,C01A,C01B,gamma);

% % Sum of Gaussian curves
symm(w,w00,w01A,w01B,C00,C01A,C01B,gamma)= ...
    C00*(1/(gamma*sqrt(2*pi)))*exp(-((w-w00).^2)/(2*gamma^2)) + ...
    C01A*(1/(gamma*sqrt(2*pi)))*exp(-((w-w01A).^2)/(2*gamma^2)) + ...
    C01B*(1/(gamma*sqrt(2*pi)))*exp(-((w-w01B).^2)/(2*gamma^2));

anti(w,w00,w01A,w01B,C00,C01A,C01B,gamma)= ...
    C00*(1/(gamma*sqrt(2*pi)))*exp(-((w-w00).^2)/(2*gamma^2)) + ...
    C01A*(1/(gamma*sqrt(2*pi)))*exp(-((w-w01A).^2)/(2*gamma^2)) + ...
    C01B*(1/(gamma*sqrt(2*pi)))*exp(-((w-w01B).^2)/(2*gamma^2));

% THESE VALUES REPRESENT A DUPLEX DIMER AT 15DEG C (OBTAINED FROM
% LIN/NONLIN FIT OF THE DATA)
% assign the center frequncies of the lorenztians 
w00s=536.63e12;
w01As=570.12e12;
w01Bs=579.74e12;

w00a=552.34e12;
w01Aa=584.36e12;
w01Ba=585.97e12;
% assign the frequncy widths of the lorenztians (should be in THz)
gammaA=10.07e12;
% gammaAa=10e12;
% assign the peak heights
C00s= 0.0934;
C01As= 0.0252;
C01Bs= 0.0063;

C00a= 0.0556;
C01Aa= 0.0771;
C01Ba= 0.0194;


figure(4)
plot(wAxis,symm(wAxis,w00s,w01As,w01Bs,C00s,C01As,C01Bs,gammaA)./max(symm(wAxis,w00s,w01As,w01Bs,C00s,C01As,C01Bs,gammaA)),wAxis,anti(wAxis,w00a,w01Aa,w01Ba,C00a,C01Aa,C01Ba,gammaA)./max(anti(wAxis,w00a,w01Aa,w01Ba,C00a,C01Aa,C01Ba,gammaA)));
hold on;
% plot(wAxis,abs(Ew_finalSTot(wAxis,t1a,Ta,AO1a,sig1a,wLa,t2a,AO2a,sig2a)));
plot(wAxis,abs(E1w_finalS(wAxis,t1a,Ta,AO1a,sig1a,wLa)));
hold off;
title('Pure Molecular Spectra - Freq domain');
xlabel('Frequency (Hz)');
ylabel('Amp (arb)');
legend('Symm Abs Spec','Anti Abs Spec','E1 Spectrum'); 


% % Cast the symm/anti as a 2x1 vector, as well as the field. Determine the
% overlap of the two.
moleSpec=[symm(wAxis,w00s,w01As,w01Bs,C00s,C01As,C01Bs,gammaA) ; anti(wAxis,w00a,w01Aa,w01Ba,C00a,C01Aa,C01Ba,gammaA)];
 
% Overlap of the two fields with the molecualr spectra.
% See Appendix Eq.9,10 in WPI 2006 article - should be equivalent to Eqn 1b
% in Tamimi 2020 as well. 
% sigTot= (conj(E1w_finalS(wAxis,t1a,Ta,AO1a,sig1a,wLa)).*E1w_finalS(wAxis,t1a,Ta,AO1a,sig1a,wLa).*((symm(wAxis,w0sa,gammaSa)+anti(wAxis,w0aa,gammaAa)).^2)) + ...
%         (conj(E2w_finalS(wAxis,t2a,Ta,AO2a,sig2a,wLa)).*E2w_finalS(wAxis,t2a,Ta,AO2a,sig2a,wLa).*((symm(wAxis,w0sa,gammaSa)+anti(wAxis,w0aa,gammaAa)).^2)) + ...
%         2*real((conj(E2w_finalS(wAxis,t2a,Ta,AO2a,sig2a,wLa)).*E1w_finalS(wAxis,t1a,Ta,AO1a,sig1a,wLa).*((symm(wAxis,w0sa,gammaSa)+anti(wAxis,w0aa,gammaAa)).^2))); 

% interfComp=2*real((conj(E2w_finalS(wAxis,t2a,Ta,AO2a,sig2a,wLa)).*E1w_finalS(wAxis,t1a,Ta,AO1a,sig1a,wLa).*((symm(wAxis,w0sa,gammaSa)+anti(wAxis,w0aa,gammaAa)).^2)));

% figure(5)
% plot(wAxis,overlapSpec(1,:),wAxis,overlapSpec(2,:));
% title('Overlap Spectra - Freq domain');
% xlabel('Frequency (Hz)');
% ylabel('Amp (arb)');
% legend('Symm Component','Anti Component');   
    
numPlots=length(Taxis);
    DimRow= floor(sqrt(numPlots)); 
    DimCol= floor(sqrt(numPlots)); 
    
    if ((DimRow*DimCol)< numPlots)
        DimRow=DimRow+1;
        if(DimRow*DimCol)< numPlots
            DimCol=DimCol+1; 
        end
    end
%     create a vector to hold the magnitudes of the polz components over
%     the mod cycle
    polzComp = zeros(2,numel(Taxis)); 
    innerProdComp = zeros(2,numel(Taxis));
    
   
    for j=1:numPlots
%         *******SCHEME 1: pulses polarized at 45 degrees sich that their projections
%         are equal... This scheme assumes the product of the two fields
%         onto the dipole captures the interference. 
%         fieldComp=[conj(E2w_finalS(wAxis,t2a,Taxis(j),AO2a,sig2a,wLa)).*(E1w_finalS(wAxis,t1a,Taxis(j),AO1a,sig1a,wLa)) ; ...
%             conj(E2w_finalS(wAxis,t2a,Taxis(j),AO2a,sig2a,wLa)).*(E1w_finalS(wAxis,t1a,Taxis(j),AO1a,sig1a,wLa))];

%        ******SCHEME 2: pulses polarized at an angle of 45 degrees relative to the dimer
%        axes. Such that one pulse (E1) has postive x and y component, but
%        E2 has a postive x and negative y component. Tryign here to modle
%        the total field as the sum, lettign the product be taken in the
%        overlap spectra
            fieldComp=[(1/sqrt(2))*(E1w_finalS(wAxis,t1a,Taxis(j),AO1a,sig1a,wLa) - E2w_finalS(wAxis,t2a,Taxis(j),AO2a,sig2a,wLa)) ...
            ; (1/sqrt(2))*(E1w_finalS(wAxis,t1a,Taxis(j),AO1a,sig1a,wLa) + E2w_finalS(wAxis,t2a,Taxis(j),AO2a,sig2a,wLa))];

        % for now the overlap spectra should contain everything - would be ideal to
        % ioslate the interference components as well 
%           SCHEME 1
%         overlapSpec= (fieldComp.*moleSpec).*conj(moleSpec);

%           SCHEME 2
            overlapSpec= (fieldComp.*moleSpec).*conj(fieldComp.*moleSpec);
        

%          overlapSpec2= ((fieldComp')*(fieldComp))*((moleSpec')*(moleSpec))';
%           SCHEME 1
%          innerProdSym= ((E2w_finalS(wAxis,t2a,Taxis(j),AO2a,sig2a,wLa)).*(moleSpec(1,:)))*((E1w_finalS(wAxis,t1a,Taxis(j),AO1a,sig1a,wLa)').*(moleSpec(1,:)'));
%          innerProdAnti= ((E2w_finalS(wAxis,t2a,Taxis(j),AO2a,sig2a,wLa)).*(moleSpec(2,:)))*((E1w_finalS(wAxis,t1a,Taxis(j),AO1a,sig1a,wLa)').*(moleSpec(2,:)'));

%           SCHEME 2
         innerProdSym= ((fieldComp(1,:)).*(moleSpec(1,:)))*((fieldComp(1,:)').*(moleSpec(1,:)'));
         innerProdAnti= ((fieldComp(2,:)).*(moleSpec(2,:)))*((fieldComp(2,:)').*(moleSpec(2,:)'));
         
         innerProdComp(1,j)=innerProdSym;
         innerProdComp(2,j)=innerProdAnti;
        
        figure(6)
        sgtitle(['Real Overlap Spectra - Freq domain']); 
        subplot(DimRow,DimCol,(j));
        plot(wAxis,real(overlapSpec(1,:)),wAxis,real(overlapSpec(2,:)));
        title(['AOM Phase=' num2str((Taxis(j)./(Taxis(end)))*2*pi) 'radians']);
        xlabel('Frequency (Hz)');
        ylabel('Amp (arb)');
%         ylim([-7e-30,7e-30]);
        legend('Symm Component','Anti Component');
        
        
        figure(7)
        sgtitle(['Imag Overlap Spectra - Freq domain']); 
        subplot(DimRow,DimCol,(j));
        plot(wAxis,imag(overlapSpec(1,:)),wAxis,imag(overlapSpec(2,:)));
        title(['AOM Phase=' num2str((Taxis(j)./(Taxis(end)))*2*pi) 'radians']);
        xlabel('Frequency (Hz)');
        ylabel('Amp (arb)');
%         ylim([-7e-30,7e-30]);
        legend('Symm Component','Anti Component');
        
        polzComp(1,j)=sum(abs(fieldComp(1,:)));
        polzComp(2,j)=sum(abs(fieldComp(2,:)));
        
    end
    
    figure(8)
        plot(Taxis,(polzComp(1,:)),Taxis,(polzComp(2,:)));
        title(['Polarization Component Magnitudes']);
        xlabel('AOM Phase');
        ylabel('Mag (arb)');
        legend('Vert Component','Horz Component');
    
         figure(9)
        plot(Taxis,real((innerProdComp(1,:))),Taxis,imag((innerProdComp(1,:))));
        hold on
        plot(Taxis,real((innerProdComp(2,:))),Taxis,imag((innerProdComp(2,:))));
        title(['Inner Prod Components']);
        xlabel('AOM Clock Cycle (sec)');
        ylabel('Mag (arb)');
        legend('Symm Real','Symm Imag','Anti Real','Anti Imag');
        

end


if Polz && ~QWPmat && simp && ~useHandles
%%
% METHOD 2: start with the spectral dependence, using eqn 1b from Tamimi
% 2020. Simulate the crosspolarized components of the strongly coupled
% dimer as they interact with the rotating broadband field. 

% Now to create a manifold for the sym/anti transistions of the dimer.
% Overlap it with the two fields, and then consider their inner products.
% See Appendix Eq.9,10 in WPI 2006 artcile 
syms Ew(w,w0,sigma0);
syms symm00(w,w00,C00,gamma);
syms symm01a(w,w01A,C01A,gamma);
syms symm01b(w,w01B,C01B,gamma);

syms anti00(w,w00,C00,gamma);
syms anti01a(w,w01A,C01A,gamma);
syms anti01b(w,w01B,C01B,gamma);


 % % Sum of Gaussian curves
% symm(w,w00,w01A,w01B,C00,C01A,C01B,gamma)= ...
%     C00*(1/(gamma*sqrt(2*pi)))*exp(-((w-w00).^2)/(2*gamma^2)) + ...
%     C01A*(1/(gamma*sqrt(2*pi)))*exp(-((w-w01A).^2)/(2*gamma^2)) + ...
%     C01B*(1/(gamma*sqrt(2*pi)))*exp(-((w-w01B).^2)/(2*gamma^2));

% anti(w,w00,w01A,w01B,C00,C01A,C01B,gamma)= ...
%     C00*(1/(gamma*sqrt(2*pi)))*exp(-((w-w00).^2)/(2*gamma^2)) + ...
%     C01A*(1/(gamma*sqrt(2*pi)))*exp(-((w-w01A).^2)/(2*gamma^2)) + ...
%     C01B*(1/(gamma*sqrt(2*pi)))*exp(-((w-w01B).^2)/(2*gamma^2));

% laser spectrum
Ew(w,w0,sigma0)=(1/(sigma0*sqrt(2*pi)))*exp(-((w-w0).^2)/(2*sigma0^2));

% symmetric states
symm00(w,w00,C00,gamma)= C00*(1/(gamma*sqrt(2*pi)))*exp(-((w-w00).^2)/(2*gamma^2)) ; 
symm01a(w,w01A,C01A,gamma)= C01A*(1/(gamma*sqrt(2*pi)))*exp(-((w-w01A).^2)/(2*gamma^2));
symm01b(w,w01B,C01B,gamma) = C01B*(1/(gamma*sqrt(2*pi)))*exp(-((w-w01B).^2)/(2*gamma^2));

% antisymmetric states
anti00(w,w00,C00,gamma)= C00*(1/(gamma*sqrt(2*pi)))*exp(-((w-w00).^2)/(2*gamma^2));
anti01a(w,w01A,C01A,gamma)= C01A*(1/(gamma*sqrt(2*pi)))*exp(-((w-w01A).^2)/(2*gamma^2));
anti01b(w,w01B,C01B,gamma)= C01B*(1/(gamma*sqrt(2*pi)))*exp(-((w-w01B).^2)/(2*gamma^2));

% THESE VALUES REPRESENT A DUPLEX DIMER AT 15DEG C (OBTAINED FROM
% LIN/NONLIN FIT OF THE DATA)
% assign the center frequncies of the lorenztians 
w00s=536.63e12;
w01As=570.12e12;
w01Bs=579.74e12;

w00a=552.34e12;
w01Aa=584.36e12;
w01Ba=585.97e12;

% assign the frequncy widths of the Gaussians (should be in THz)
gammaA=10.07e12;
% gammaAa=10e12;
% assign the peak heights
C00s= 0.0934;
C01As= 0.0252;
C01Bs= 0.0063;

C00a= 0.0556;
C01Aa= 0.0771;
C01Ba= 0.0194;

% 532nm is the same as 563THz
wLa=563e12;
% this is the standard deviation corresponding to a 35nm FWHM
sigmaL=14.7e12;

% AO1 and T is going to play the role of the AOM phase, Where AO1 is the
% freq of imparted by the AOM and T is the time within the modulaiton
% cycle. Notably, the t-space for the pulse arrivals and the T-space for
% the AOM cycle are not interdependent at all. This will manifest itself in
% the transform. 
AO1a=1.005e6; 
AO2a=1.0e6;
Ta=0e-4;
Taxis=linspace(0,2e-4,20);
t21Axis=linspace(0,60e-15,3); 
wAxis=linspace(430e12,750e12,1000); 

% *****interpulse delay in fs*****
t21=0e-15;

% symmetric manifold background term (we'll assume along the V axis, so
% assign the cos(A01-AO2) to it 
symmBkgd = 2*(cos(((2*pi*(AO1a-AO2a)*Ta)/2) - pi/4)*Ew(wAxis,wLa,sigmaL).*(symm00(wAxis,w00s,C00s,gammaA) + symm01a(wAxis,w01As,C01As,gammaA) + symm01b(wAxis,w01Bs,C01Bs,gammaA))).^2 ;

antiBkgd =  2*(-sin(((2*pi*(AO1a-AO2a)*Ta)/2) - pi/4)*Ew(wAxis,wLa,sigmaL).*(anti00(wAxis,w00a,C00a,gammaA) + anti01a(wAxis,w01Aa,C01Aa,gammaA) + anti01b(wAxis,w01Ba,C01Ba,gammaA))).^2 ;


figure(1)
plot(wAxis,symmBkgd./(max(symmBkgd)))
hold on
plot(wAxis,antiBkgd./(max(antiBkgd)))
plot(wAxis,Ew(wAxis,wLa,sigmaL)./(2*max(Ew(wAxis,wLa,sigmaL))));
title('Background Terms');
xlabel('Frequency (Hz)'); 
ylabel('Amp (arb)');
legend('Symm Bkgd','Anti Bkgd','Laser Spectrum'); 
hold off

% symmetric and antisymmetric modulated terms 
symmMod = 2*(cos(((2*pi*(AO1a-AO2a)*Ta)/2) - pi/4)*Ew(wAxis,wLa,sigmaL).* ...
   (symm00(wAxis,w00s,C00s,gammaA)*cos(2*pi*(AO1a-AO2a)*Ta - t21*w00s) + ...
    symm01a(wAxis,w01As,C01As,gammaA)*cos(2*pi*(AO1a-AO2a)*Ta - t21*w01As) + ...
    symm01b(wAxis,w01Bs,C01Bs,gammaA)*cos(2*pi*(AO1a-AO2a)*Ta - t21*w01Bs))).^2 ;

antiMod = 2*(-sin(((2*pi*(AO1a-AO2a)*Ta)/2) - pi/4)*Ew(wAxis,wLa,sigmaL).* ...
    (anti00(wAxis,w00a,C00a,gammaA)*cos(2*pi*(AO1a-AO2a)*Ta - t21*w00a) + ...
    anti01a(wAxis,w01Aa,C01Aa,gammaA)*cos(2*pi*(AO1a-AO2a)*Ta - t21*w01Aa) + ...
    anti01b(wAxis,w01Ba,C01Ba,gammaA)*cos(2*pi*(AO1a-AO2a)*Ta - t21*w01Ba))).^2 ;

figure(2)
plot(wAxis,symmMod./(max(symmMod)))
hold on
plot(wAxis,antiMod./(max(antiMod)))
plot(wAxis,Ew(wAxis,wLa,sigmaL)./(2*max(Ew(wAxis,wLa,sigmaL))));
title('Modulated Terms');
xlabel('Frequency (Hz)'); 
ylabel('Amp (arb)');
legend('Symm Mod','Anti Mod','Laser Spectrum'); 
hold off

symmModAmp=zeros(numel(t21Axis),numel(Taxis));
antiModAmp=zeros(numel(t21Axis),numel(Taxis));

symmBkgdAmp=zeros(numel(t21Axis),numel(Taxis));
antiBkgdAmp=zeros(numel(t21Axis),numel(Taxis));

symmTotAmp=zeros(numel(t21Axis),numel(Taxis));
antiTotAmp=zeros(numel(t21Axis),numel(Taxis));

for m=1:numel(t21Axis)
for j=1:numel(Taxis)
    % symmetric manifold background term (we'll assume along the V axis, so
    % assign the cos(A01-AO2) to it and the sin() to the H axis
    symmBkgd = 2*(cos(((2*pi*(AO1a-AO2a)*Taxis(j))/2) - pi/4)*Ew(wAxis,wLa,sigmaL).*(symm00(wAxis,w00s,C00s,gammaA) + symm01a(wAxis,w01As,C01As,gammaA) + symm01b(wAxis,w01Bs,C01Bs,gammaA))).^2 ;
    antiBkgd =  2*(-sin(((2*pi*(AO1a-AO2a)*Taxis(j))/2) - pi/4)*Ew(wAxis,wLa,sigmaL).*(anti00(wAxis,w00a,C00a,gammaA) + anti01a(wAxis,w01Aa,C01Aa,gammaA) + anti01b(wAxis,w01Ba,C01Ba,gammaA))).^2 ;
    
    % symmetric and antisymmetric modulated terms 
    symmMod = 2*(cos(((2*pi*(AO1a-AO2a)*Taxis(j))/2) - pi/4)*Ew(wAxis,wLa,sigmaL).* ...
    (symm00(wAxis,w00s,C00s,gammaA)*cos(2*pi*(AO1a-AO2a)*Taxis(j) - t21Axis(m)*w00s) + ...
    symm01a(wAxis,w01As,C01As,gammaA)*cos(2*pi*(AO1a-AO2a)*Taxis(j) - t21Axis(m)*w01As) + ...
    symm01b(wAxis,w01Bs,C01Bs,gammaA)*cos(2*pi*(AO1a-AO2a)*Taxis(j) - t21Axis(m)*w01Bs))).^2 ;

    antiMod = 2*(-sin(((2*pi*(AO1a-AO2a)*Taxis(j))/2) - pi/4)*Ew(wAxis,wLa,sigmaL).* ...
    (anti00(wAxis,w00a,C00a,gammaA)*cos(2*pi*(AO1a-AO2a)*Taxis(j) - t21Axis(m)*w00a) + ...
    anti01a(wAxis,w01Aa,C01Aa,gammaA)*cos(2*pi*(AO1a-AO2a)*Taxis(j) - t21Axis(m)*w01Aa) + ...
    anti01b(wAxis,w01Ba,C01Ba,gammaA)*cos(2*pi*(AO1a-AO2a)*Taxis(j) - t21Axis(m)*w01Ba))).^2 ;

    symmModAmp(m,j)= sum(symmMod);
    antiModAmp(m,j)= sum(antiMod);

    symmBkgdAmp(m,j)= sum(symmBkgd);
    antiBkgdAmp(m,j)= sum(antiBkgd);

    symmTotAmp(m,j)= (sum(symmMod)+sum(symmBkgd));
    antiTotAmp(m,j)=(sum(antiMod)+sum(antiBkgd));
    
    
end

  

end

figure(3)
plot(Taxis,symmModAmp,Taxis,symmBkgdAmp)
title('Symmetric - Modulated and Background terms');
xlabel('AOM clock cycle');
ylabel('Mag (arb)');
legend('Modulated','Background');

figure(4)
plot(Taxis,antiModAmp,Taxis,antiBkgdAmp)
title('Antisymmetric - Modulated and Background terms');
xlabel('AOM clock cycle');
ylabel('Mag (arb)');
legend('Modulated','Background');

figure(5)
for n=1:numel(t21Axis)
plot(Taxis,symmTotAmp(n,:),Taxis,antiTotAmp(n,:))
hold on
end
title('Total Signal - Modulated plus Background terms');
xlabel('AOM clock cycle');
ylabel('Mag (arb)');
legend('Symmetric','Antisymmetric');



end

if Polz && simp && useHandles

%%
%%
% METHOD 2: start with the spectral dependence, using eqn 1b from Tamimi
% 2020. Simulate the crosspolarized components of the strongly coupled
% dimer as they interact with the rotating broadband field. 

% Now to create a manifold for the sym/anti transistions of the dimer.
% Overlap it with the two fields, and then consider their inner products.
% See Appendix Eq.9,10 in WPI 2006 artcile 
% syms Ew(w,w0,sigma0);
% syms symm00(w,w00,C00,gamma);
% syms symm01a(w,w01A,C01A,gamma);
% syms symm01b(w,w01B,C01B,gamma);
% 
% syms anti00(w,w00,C00,gamma);
% syms anti01a(w,w01A,C01A,gamma);
% syms anti01b(w,w01B,C01B,gamma);


 % % Sum of Gaussian curves
% symm(w,w00,w01A,w01B,C00,C01A,C01B,gamma)= ...
%     C00*(1/(gamma*sqrt(2*pi)))*exp(-((w-w00).^2)/(2*gamma^2)) + ...
%     C01A*(1/(gamma*sqrt(2*pi)))*exp(-((w-w01A).^2)/(2*gamma^2)) + ...
%     C01B*(1/(gamma*sqrt(2*pi)))*exp(-((w-w01B).^2)/(2*gamma^2));

% anti(w,w00,w01A,w01B,C00,C01A,C01B,gamma)= ...
%     C00*(1/(gamma*sqrt(2*pi)))*exp(-((w-w00).^2)/(2*gamma^2)) + ...
%     C01A*(1/(gamma*sqrt(2*pi)))*exp(-((w-w01A).^2)/(2*gamma^2)) + ...
%     C01B*(1/(gamma*sqrt(2*pi)))*exp(-((w-w01B).^2)/(2*gamma^2));

% laser spectrum
Ew = @(w,w0,sigma0) (1/(sigma0*sqrt(2*pi)))*exp(-((w-w0).^2)/(2*sigma0^2));

% symmetric states
symm00 = @(w,w00,C00,gamma) C00*(1/(gamma*sqrt(2*pi)))*exp(-((w-w00).^2)/(2*gamma^2)) ; 
symm01a = @(w,w01A,C01A,gamma) C01A*(1/(gamma*sqrt(2*pi)))*exp(-((w-w01A).^2)/(2*gamma^2));
symm01b = @(w,w01B,C01B,gamma) C01B*(1/(gamma*sqrt(2*pi)))*exp(-((w-w01B).^2)/(2*gamma^2));

% antisymmetric states
anti00 = @(w,w00,C00,gamma) C00*(1/(gamma*sqrt(2*pi)))*exp(-((w-w00).^2)/(2*gamma^2));
anti01a = @(w,w01A,C01A,gamma) C01A*(1/(gamma*sqrt(2*pi)))*exp(-((w-w01A).^2)/(2*gamma^2));
anti01b = @(w,w01B,C01B,gamma) C01B*(1/(gamma*sqrt(2*pi)))*exp(-((w-w01B).^2)/(2*gamma^2));

% THESE VALUES REPRESENT A DUPLEX DIMER AT 15DEG C (OBTAINED FROM
% LIN/NONLIN FIT OF THE DATA)
% assign the center frequncies of the lorenztians 
w00a=536.63e12;
w01Aa=570.12e12;
w01Ba=579.74e12;

w00s=552.34e12;
w01As=584.36e12;
w01Bs=585.97e12;

% assign the frequncy widths of the Gaussians (should be in THz)
gammaA=10.07e12;
% gammaAa=10e12;
% assign the peak heights
C00a= 0.0934;
C01Aa= 0.0252;
C01Ba= 0.0063;

C00s= 0.0556;
C01As= 0.0771;
C01Bs= 0.0194;

% 532nm is the same as 563THz
% wLa=563e12;
wLa=555e12;

% this is the standard deviation corresponding to a 35nm FWHM
sigmaL=14.7e12;

% AO1 and T is going to play the role of the AOM phase, Where AO1 is the
% freq of imparted by the AOM and T is the time within the modulaiton
% cycle. Notably, the t-space for the pulse arrivals and the T-space for
% the AOM cycle are not interdependent at all. This will manifest itself in
% the transform. 
AO1a=1.005e6; 
AO2a=1.0e6;
Ta=0e-4;
Taxis=linspace(0,2e-4,20);
t21Axis=linspace(0,90e-15,10);
% t21Axis=0;
wAxis=linspace(430e12,750e12,1000); % t21Axis=linspace(0,60e-15,5); 

% *****interpulse delay in fs*****
t21=0e-15;

% symmetric manifold background term (we'll assume along the V axis, so
% assign the cos(A01-AO2) to it 
symmBkgd = 2*(cos(((2*pi*(AO1a-AO2a)*Ta)/2) - pi/4)*Ew(wAxis,wLa,sigmaL).*(symm00(wAxis,w00s,C00s,gammaA) + symm01a(wAxis,w01As,C01As,gammaA) + symm01b(wAxis,w01Bs,C01Bs,gammaA))).^2 ;

antiBkgd =  2*(-sin(((2*pi*(AO1a-AO2a)*Ta)/2) - pi/4)*Ew(wAxis,wLa,sigmaL).*(anti00(wAxis,w00a,C00a,gammaA) + anti01a(wAxis,w01Aa,C01Aa,gammaA) + anti01b(wAxis,w01Ba,C01Ba,gammaA))).^2 ;


figure(1)
plot(wAxis,symmBkgd./Ew(wAxis,wLa,sigmaL))
hold on
plot(wAxis,antiBkgd./Ew(wAxis,wLa,sigmaL))
plot(wAxis,2e-43*Ew(wAxis,wLa,sigmaL)./(2*max(Ew(wAxis,wLa,sigmaL))));
title('Background Terms');
xlabel('Frequency (Hz)'); 
xlim([500e12,625e12]); 
ylabel('Amp (arb)');
legend('Symm Bkgd','Anti Bkgd','Laser Spectrum'); 
hold off

figure(5)
plot(wAxis,(symm00(wAxis,w00s,C00s,gammaA) + symm01a(wAxis,w01As,C01As,gammaA) + symm01b(wAxis,w01Bs,C01Bs,gammaA)))
hold on
plot(wAxis,(anti00(wAxis,w00a,C00a,gammaA) + anti01a(wAxis,w01Aa,C01Aa,gammaA) + anti01b(wAxis,w01Ba,C01Ba,gammaA)))
plot(wAxis,6e-15*Ew(wAxis,wLa,sigmaL)./(2*max(Ew(wAxis,wLa,sigmaL))));
plot([w00s w00s],[0 max(symm00(wAxis,w00s,C00s,gammaA))],'--b', ...
    [w01As w01As],[0 max(symm01a(wAxis,w01As,C01As,gammaA))],'--b', ...
    [w01Bs w01Bs],[0 max(symm01b(wAxis,w01Bs,C01Bs,gammaA))], '--b',...
    [w00a w00a],[0 max(anti00(wAxis,w00a,C00a,gammaA))],'--r', ... 
    [w01Aa w01Aa],[0 max(anti01a(wAxis,w01Aa,C01Aa,gammaA))],'--r', ... 
    [w01Ba w01Ba],[0 max(anti01b(wAxis,w01Ba,C01Ba,gammaA))],'--r'); 
title('Background Terms');
xlabel('Frequency (Hz)'); 
xlim([500e12,625e12]); 
ylabel('Amp (arb)');
legend('Symm Bkgd','Anti Bkgd','Laser Spectrum'); 
hold off

% symmetric and antisymmetric modulated terms 
symmMod = 2*(cos(((2*pi*(AO1a-AO2a)*Ta)/2) - pi/4)*Ew(wAxis,wLa,sigmaL).* ...
   (symm00(wAxis,w00s,C00s,gammaA)*cos(2*pi*(AO1a-AO2a)*Ta - t21*w00s) + ...
    symm01a(wAxis,w01As,C01As,gammaA)*cos(2*pi*(AO1a-AO2a)*Ta - t21*w01As) + ...
    symm01b(wAxis,w01Bs,C01Bs,gammaA)*cos(2*pi*(AO1a-AO2a)*Ta - t21*w01Bs))).^2 ;

antiMod = 2*(-sin(((2*pi*(AO1a-AO2a)*Ta)/2) - pi/4)*Ew(wAxis,wLa,sigmaL).* ...
    (anti00(wAxis,w00a,C00a,gammaA)*cos(2*pi*(AO1a-AO2a)*Ta - t21*w00a) + ...
    anti01a(wAxis,w01Aa,C01Aa,gammaA)*cos(2*pi*(AO1a-AO2a)*Ta - t21*w01Aa) + ...
    anti01b(wAxis,w01Ba,C01Ba,gammaA)*cos(2*pi*(AO1a-AO2a)*Ta - t21*w01Ba))).^2 ;

figure(2)
plot(wAxis,symmMod./(max(symmMod)))
hold on
plot(wAxis,antiMod./(max(antiMod)))
plot(wAxis,Ew(wAxis,wLa,sigmaL)./(2*max(Ew(wAxis,wLa,sigmaL))));
title('Modulated Terms');
xlabel('Frequency (Hz)'); 
ylabel('Amp (arb)');
legend('Symm Mod','Anti Mod','Laser Spectrum'); 
hold off

numPlots=length(t21Axis);
    DimRow= floor(sqrt(numPlots)); 
    DimCol= floor(sqrt(numPlots)); 
    
    if ((DimRow*DimCol)< numPlots)
        DimRow=DimRow+1;
        if(DimRow*DimCol)< numPlots
            DimCol=DimCol+1; 
        end
    end

symmModAmp=zeros(numel(t21Axis),numel(Taxis));
antiModAmp=zeros(numel(t21Axis),numel(Taxis));

symmBkgdAmp=zeros(numel(t21Axis),numel(Taxis));
antiBkgdAmp=zeros(numel(t21Axis),numel(Taxis));

symmTotAmp=zeros(numel(t21Axis),numel(Taxis));
antiTotAmp=zeros(numel(t21Axis),numel(Taxis));

for m=1:numel(t21Axis)
for j=1:numel(Taxis)
    % symmetric manifold background term (we'll assume along the V axis, so
    % assign the cos(A01-AO2) to it and the sin() to the H axis
    symmBkgd = 2*(cos(((2*pi*(AO1a-AO2a)*Taxis(j))/2) - pi/4)*Ew(wAxis,wLa,sigmaL).*(symm00(wAxis,w00s,C00s,gammaA) + symm01a(wAxis,w01As,C01As,gammaA) + symm01b(wAxis,w01Bs,C01Bs,gammaA))).^2 ;
    antiBkgd =  2*(-sin(((2*pi*(AO1a-AO2a)*Taxis(j))/2) - pi/4)*Ew(wAxis,wLa,sigmaL).*(anti00(wAxis,w00a,C00a,gammaA) + anti01a(wAxis,w01Aa,C01Aa,gammaA) + anti01b(wAxis,w01Ba,C01Ba,gammaA))).^2 ;
    
    % symmetric and antisymmetric modulated terms 
    symmMod = 2*(cos(((2*pi*(AO1a-AO2a)*Taxis(j))/2) - pi/4)*Ew(wAxis,wLa,sigmaL).* ...
    (symm00(wAxis,w00s,C00s,gammaA)*cos(2*pi*(AO1a-AO2a)*Taxis(j) - t21Axis(m)*w00s) + ...
    symm01a(wAxis,w01As,C01As,gammaA)*cos(2*pi*(AO1a-AO2a)*Taxis(j) - t21Axis(m)*w01As) + ...
    symm01b(wAxis,w01Bs,C01Bs,gammaA)*cos(2*pi*(AO1a-AO2a)*Taxis(j) - t21Axis(m)*w01Bs))).^2 ;

    antiMod = 2*(-sin(((2*pi*(AO1a-AO2a)*Taxis(j))/2) - pi/4)*Ew(wAxis,wLa,sigmaL).* ...
    (anti00(wAxis,w00a,C00a,gammaA)*cos(2*pi*(AO1a-AO2a)*Taxis(j) - t21Axis(m)*w00a) + ...
    anti01a(wAxis,w01Aa,C01Aa,gammaA)*cos(2*pi*(AO1a-AO2a)*Taxis(j) - t21Axis(m)*w01Aa) + ...
    anti01b(wAxis,w01Ba,C01Ba,gammaA)*cos(2*pi*(AO1a-AO2a)*Taxis(j) - t21Axis(m)*w01Ba))).^2 ;

    symmModAmp(m,j)= sum(symmMod);
    antiModAmp(m,j)= sum(antiMod);

    symmBkgdAmp(m,j)= sum(symmBkgd);
    antiBkgdAmp(m,j)= sum(antiBkgd);

    symmTotAmp(m,j)= (sum(symmMod)+sum(symmBkgd));
    antiTotAmp(m,j)=(sum(antiMod)+sum(antiBkgd));
    
    
end

  

end

figure(3)
plot(Taxis,symmModAmp)
title('Symmetric - Modulated terms');
xlabel('AOM clock cycle');
ylabel('Mag (arb)');
% legend('Modulated','Background');

figure(4)
plot(Taxis,antiModAmp)
title('Antisymmetric - Modulated terms');
xlabel('AOM clock cycle');
ylabel('Mag (arb)');
% legend('Modulated','Background');

figure(6)
plot(Taxis,symmBkgdAmp)
title('Symmetric - Background terms');
xlabel('AOM clock cycle');
ylabel('Mag (arb)');
% legend('Modulated','Background');

figure(7)
plot(Taxis,antiBkgdAmp)
title('Antisymmetric - Background terms');
xlabel('AOM clock cycle');
ylabel('Mag (arb)');
% legend('Modulated','Background');

figure(8)
for n=1:numel(t21Axis)
    subplot(DimRow,DimCol,n)
    plot(Taxis,symmTotAmp(n,:),'--')
    hold on
    plot(Taxis,symmModAmp(n,:))
    plot(Taxis,symmBkgdAmp(n,:))
    title(['Symmetric Linear Signal - t21=' num2str(t21Axis(n)*1e15) 'fs']);
    xlabel('AOM clock cycle');
    ylabel('Mag (arb)');
    legend('Total Signal','Modulated','Background');

end

figure(9)
for n=1:numel(t21Axis)
    subplot(DimRow,DimCol,n)
    plot(Taxis,antiTotAmp(n,:),'--')
    hold on
    plot(Taxis,antiModAmp(n,:))
    plot(Taxis,antiBkgdAmp(n,:))
    title(['Antisymmetric Linear Signal - t21=' num2str(t21Axis(n)*1e15) 'fs']);
    xlabel('AOM clock cycle');
    ylabel('Mag (arb)');
    legend('Total Signal','Modulated','Background');

end

figure(10)
for n=1:numel(t21Axis)
    subplot(DimRow,DimCol,n)
    plot(Taxis,(antiModAmp(n,:)+symmModAmp(n,:)),'--')
    hold on
    plot(Taxis,symmModAmp(n,:))
    plot(Taxis,antiModAmp(n,:))
    title(['Modulated Linear Signal - t21=' num2str(t21Axis(n)*1e15) 'fs']);
    xlabel('AOM clock cycle');
    ylabel('Mag (arb)');
    legend('Total Mod Signal','Symm Mod','Anti Mod');

end

figure(11)
for n=1:numel(t21Axis)
    subplot(DimRow,DimCol,n)
    plot(Taxis,(antiModAmp(n,:)+symmModAmp(n,:)+antiBkgdAmp(n,:)+symmBkgdAmp(n,:)),'--')
    hold on
    plot(Taxis,symmModAmp(n,:) + symmBkgdAmp(n,:))
    plot(Taxis,antiModAmp(n,:) + antiBkgdAmp(n,:))
    title(['Modulated Linear Signal - t21=' num2str(t21Axis(n)*1e15) 'fs']);
    xlabel('AOM clock cycle');
    ylabel('Mag (arb)');
    legend('Total Signal','Symm Mod + Bkgd','Anti Mod + Bkgd');

end

figure(12)
for n=1:numel(t21Axis)
    subplot(DimRow,DimCol,n)
    plot(Taxis,(antiBkgdAmp(n,:)+symmBkgdAmp(n,:)),'--')
    hold on
    plot(Taxis,symmBkgdAmp(n,:))
    plot(Taxis,antiBkgdAmp(n,:))
    title(['Background Linear Signal - t21=' num2str(t21Axis(n)*1e15) 'fs']);
    xlabel('AOM clock cycle');
    ylabel('Mag (arb)');
    legend('Total Background','Symm Bkgd','Anti Bkgd');

end
% title('Total Signal - Modulated plus Background terms');
% xlabel('AOM clock cycle');
% ylabel('Mag (arb)');
% % legend('Symmetric','Antisymmetric');

end