%%
% METHOD 2: start with the spectral dependence, using eqn 1b from Tamimi
% 2020. Simulate the crosspolarized components of the strongly coupled
% dimer as they interact with the rotating broadband field. 

% UPDATE: 8-26-2021: start from the same form as v3 of the linResp, but now
% pull the function definitions outside the loop and incorporate a dot
% product between the components of the field and symm/anti manifolds (as
% outlined in notebook)

% booleans to control output
showSubplots=0; 
showTotalPlots=1;

sumComp=1;
integrateComp=0;

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

% this angle will be orientation of the symmetric/anti with respect to the
% vertical axis (V). Given a value of 0 radians will mean the Symm axis is
% lined up perfectly with the vertical polarization component - which in
% this code goes with the cos(deltaAOM/2 -pi/4) term - will need to be
% changed in order to reflect Jones vector formalism... 
phi0=0;

% 532nm is the same as 563THz
% wLa=563e12;

% 537nm
% wLa=558.27e12; 

% wLa=550e12;

% 540 nm is 555e12 Hz
% wLa=555e12;

% 528nm
%  wLa=567.78e12;
 
%  522nm
 wLa=574.31e12; 

% frequnecies for the reference to use (going to use 560nm = 535Thz) 
% wR=535.34e12;

% 550 nm ref
% wR=545.07e12;

% 532nm ref
% wR=563e12;

% 515nm ref
wR=582.12e12;

% 520nm ref
% wR=576.52e12;

% 500nm ref
% wR=599.48e12; 

% 570nm ref
% wR=525.95e12;
% integration time for the lock-in - will use a fairly standard 100ms for
% now
intTime= 100e-3;

% Dephasing constant based on 15deg C duplex dimer measuremtns 
TauDephase=106e-15;

% this is the standard deviation corresponding to a 35nm FWHM (14.7e12 THz)
sigmaL=14.7e12;
% sigmaL=94.7e12;


% sigmaL=24.7e12;

% AO1 and T is going to play the role of the AOM phase, Where AO1 is the
% freq of imparted by the AOM and T is the time within the modulaiton
% cycle. Notably, the t-space for the pulse arrivals and the T-space for
% the AOM cycle are not interdependent at all. This will manifest itself in
% the transform. 
AO1a=1.005e6; 
AO2a=1.0e6;
Ta=0e-4;

if sumComp
%     demodulation is sensitive to the spacing and length here - as it
%     seems possible for the AOM cycle to be insuffiently sampled in order
%     to capture the effects of the rotation
Taxis=linspace(0,20e-4,200);
else
Taxis=linspace(0,20e-4,12);
end

% interpulse delay axis in femtoseconds
t21Axis=linspace(0,200e-15,201);


% t21Axis=linspace(0e-15,100e-15,12);

% t21Axis=0;
% frequency axis to evaluate the laser and molecular spectra over 
if sumComp
wAxis=linspace(430e12,750e12,200); 

else
wAxis=linspace(430e12,750e12,1000); 
end
% *****interpulse delay in fs*****
% t21=0e-15;

% symmetric manifold background term (we'll assume along the V axis, so
% assign the cos(A01-AO2) to it - this will need to be altered for other
% configrations of the dimer - might want to implement a proper dot product
% 
% symmBkgd = 2*(cos(((2*pi*(AO1a-AO2a)*Ta)/2) - pi/4)*Ew(wAxis,wLa,sigmaL).*(symm00(wAxis,w00s,C00s,gammaA) + symm01a(wAxis,w01As,C01As,gammaA) + symm01b(wAxis,w01Bs,C01Bs,gammaA))).^2 ;
% symmBkgd =  2*((cos(((2*pi*(AO1a-AO2a)*Ta)/2) - pi/4).*Ew(wAxis,wLa,sigmaL).*symm00(wAxis,w00s,C00s,gammaA)).^2 ...
%             + (cos(((2*pi*(AO1a-AO2a)*Ta)/2) - pi/4).*Ew(wAxis,wLa,sigmaL).*symm01a(wAxis,w01As,C01As,gammaA)).^2 ...
%             + (cos(((2*pi*(AO1a-AO2a)*Ta)/2) - pi/4).*Ew(wAxis,wLa,sigmaL).*symm01b(wAxis,w01Bs,C01Bs,gammaA)).^2) ;
%         
% antiBkgd =  2*((-sin(((2*pi*(AO1a-AO2a)*Ta)/2) - pi/4).*Ew(wAxis,wLa,sigmaL).*anti00(wAxis,w00a,C00a,gammaA)).^2 ...
%             + (-sin(((2*pi*(AO1a-AO2a)*Ta)/2) - pi/4).*Ew(wAxis,wLa,sigmaL).*anti01a(wAxis,w01Aa,C01Aa,gammaA)).^2 ...
%             + (-sin(((2*pi*(AO1a-AO2a)*Ta)/2) - pi/4).*Ew(wAxis,wLa,sigmaL).*anti01b(wAxis,w01Ba,C01Ba,gammaA)).^2) ;

% background signals for demodulation - with the added terms to capture
% equations 3b/a in Tamimi 2020
% symmBkgdDemodX = @(w,tAOM) 2.*cos(2*pi*(AO1a-AO2a)*tAOM - wR*t21).*exp(-tAOM/intTime).*(2/intTime).* (cos(((2*pi*(AO1a-AO2a)*tAOM)/2) - pi/4).*Ew(w,wLa,sigmaL).*(symm00(w,w00s,C00s,gammaA) + symm01a(w,w01As,C01As,gammaA) + symm01b(w,w01Bs,C01Bs,gammaA))).^2 ;
% antiBkgdDemodX =  @(w,tAOM) 2.*cos(2*pi*(AO1a-AO2a)*tAOM - wR*t21).*exp(-tAOM/intTime).*(2/intTime).* (-sin(((2*pi*(AO1a-AO2a)*tAOM)/2) - pi/4).*Ew(w,wLa,sigmaL).*(anti00(w,w00a,C00a,gammaA) + anti01a(w,w01Aa,C01Aa,gammaA) + anti01b(w,w01Ba,C01Ba,gammaA))).^2 ;

% symmBkgdDemodX = @(w,tAOM) 2.*cos(2*pi*(AO1a-AO2a)*tAOM - wR*t21).*exp(-tAOM/intTime).*(2/intTime).* ...
%             ((cos(((2*pi*(AO1a-AO2a)*tAOM)/2) - pi/4).*Ew(w,wLa,sigmaL).*symm00(w,w00s,C00s,gammaA)).^2 ...
%             + (cos(((2*pi*(AO1a-AO2a)*tAOM)/2) - pi/4).*Ew(w,wLa,sigmaL).*symm01a(w,w01As,C01As,gammaA)).^2 ...
%             + (cos(((2*pi*(AO1a-AO2a)*tAOM)/2) - pi/4).*Ew(w,wLa,sigmaL).*symm01b(w,w01Bs,C01Bs,gammaA)).^2) ;
%         
% antiBkgdDemodX =  @(w,tAOM) 2.*cos(2*pi*(AO1a-AO2a)*tAOM - wR*t21).*exp(-tAOM/intTime).*(2/intTime).* ...
%              ((-sin(((2*pi*(AO1a-AO2a)*tAOM)/2) - pi/4).*Ew(w,wLa,sigmaL).*anti00(w,w00a,C00a,gammaA)).^2 ...
%             + (-sin(((2*pi*(AO1a-AO2a)*tAOM)/2) - pi/4).*Ew(w,wLa,sigmaL).*anti01a(w,w01Aa,C01Aa,gammaA)).^2 ...
%             + (-sin(((2*pi*(AO1a-AO2a)*tAOM)/2) - pi/4).*Ew(w,wLa,sigmaL).*anti01b(w,w01Ba,C01Ba,gammaA)).^2) ;
%         
% symmBkgdDemodY = @(w,tAOM) -2.*sin(2*pi*(AO1a-AO2a)*tAOM - wR*t21).*exp(-tAOM/intTime).*(2/intTime).* ...
%             ((cos(((2*pi*(AO1a-AO2a)*tAOM)/2) - pi/4).*Ew(w,wLa,sigmaL).*symm00(w,w00s,C00s,gammaA)).^2 ...
%             + (cos(((2*pi*(AO1a-AO2a)*tAOM)/2) - pi/4).*Ew(w,wLa,sigmaL).*symm01a(w,w01As,C01As,gammaA)).^2 ...
%             + (cos(((2*pi*(AO1a-AO2a)*tAOM)/2) - pi/4).*Ew(w,wLa,sigmaL).*symm01b(w,w01Bs,C01Bs,gammaA)).^2) ;
%         
% antiBkgdDemodY =  @(w,tAOM) -2.*sin(2*pi*(AO1a-AO2a)*tAOM - wR*t21).*exp(-tAOM/intTime).*(2/intTime).* ...
%              ((-sin(((2*pi*(AO1a-AO2a)*tAOM)/2) - pi/4).*Ew(w,wLa,sigmaL).*anti00(w,w00a,C00a,gammaA)).^2 ...
%             + (-sin(((2*pi*(AO1a-AO2a)*tAOM)/2) - pi/4).*Ew(w,wLa,sigmaL).*anti01a(w,w01Aa,C01Aa,gammaA)).^2 ...
%             + (-sin(((2*pi*(AO1a-AO2a)*tAOM)/2) - pi/4).*Ew(w,wLa,sigmaL).*anti01b(w,w01Ba,C01Ba,gammaA)).^2) ;

% symmBkgdDemodY = @(w,tAOM) -2.*sin(2*pi*(AO1a-AO2a)*tAOM - wR*t21).*exp(-tAOM/intTime).*(2/intTime).* (cos(((2*pi*(AO1a-AO2a)*tAOM)/2) - pi/4).*Ew(w,wLa,sigmaL).*(symm00(w,w00s,C00s,gammaA) + symm01a(w,w01As,C01As,gammaA) + symm01b(w,w01Bs,C01Bs,gammaA))).^2 ;
% antiBkgdDemodY =  @(w,tAOM) -2.*sin(2*pi*(AO1a-AO2a)*tAOM - wR*t21).*exp(-tAOM/intTime).*(2/intTime).* (-sin(((2*pi*(AO1a-AO2a)*tAOM)/2) - pi/4).*Ew(w,wLa,sigmaL).*(anti00(w,w00a,C00a,gammaA) + anti01a(w,w01Aa,C01Aa,gammaA) + anti01b(w,w01Ba,C01Ba,gammaA))).^2 ;

figure(1)
plot(wAxis,Ew(wAxis,wLa,sigmaL).*(symm00(wAxis,w00s,C00s,gammaA) + symm01a(wAxis,w01As,C01As,gammaA) + symm01b(wAxis,w01Bs,C01Bs,gammaA)))
hold on
plot(wAxis,Ew(wAxis,wLa,sigmaL).*(anti00(wAxis,w00a,C00a,gammaA) + anti01a(wAxis,w01Aa,C01Aa,gammaA) + anti01b(wAxis,w01Ba,C01Ba,gammaA)))
plot(wAxis,Ew(wAxis,wLa,sigmaL).*(anti00(wAxis,w00a,C00a,gammaA) + anti01a(wAxis,w01Aa,C01Aa,gammaA) + anti01b(wAxis,w01Ba,C01Ba,gammaA) +...
    symm00(wAxis,w00s,C00s,gammaA) + symm01a(wAxis,w01As,C01As,gammaA) + symm01b(wAxis,w01Bs,C01Bs,gammaA)),'--');
plot(wAxis,10e-29*Ew(wAxis,wLa,sigmaL)./(2*max(Ew(wAxis,wLa,sigmaL))));
title('Overlap Spectra - laser and manifold product');
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
title('Laser and Manifolds overlay- no product');
xlabel('Frequency (Hz)'); 
xlim([500e12,625e12]); 
ylabel('Amp (arb)');
legend('Symm Bkgd','Anti Bkgd','Laser Spectrum'); 
hold off

tic
% DEFINE ALL THE RELEVANT RESPONSE FUNCTIONS HERE AS FUINCTIONS OF W, TAOM
% AND T21
% 8-26-2021: move all function defintions outside the loop and parameteried
% by t21 instead. This seems to be only marginally faster than redefining
% the function with each iteration. The outcome of the sim is the same i
% both cases.
% 8-26-2021: going to incoporate the dot product here as outlined in the
% notes
% 8-18-2021: eliminated the minus signs from the Y-quad terms (both for inf
% and bkgd) since it didnt make sense for the background magnitude ot be
% negative 

% the cos(phi+pi/2) terms with the cos(deltaAOM) will be the horz component
% of the polarization while the corresponding sin() terms are the vertical
% component. By convention the symm orientation in geometric space will
% always be equal to the anti (phi0) plus pi/2
symmModX = @(w,tAOM,t21,phi) 2.*exp(-abs(t21)/TauDephase).*( ...
    ((cos(phi+pi/2).*cos(((2*pi*(AO1a-AO2a)*tAOM)/2) - pi/4).*Ew(w,wLa,sigmaL).*symm00(w,w00s,C00s,gammaA)).^2 ).* cos(t21*(w00s-wR)) + ...
    ((cos(phi+pi/2).*cos(((2*pi*(AO1a-AO2a)*tAOM)/2) - pi/4).*Ew(w,wLa,sigmaL).*symm01a(w,w01As,C01As,gammaA)).^2) .* cos(t21*(w01As-wR)) + ...
    ((cos(phi+pi/2).*cos(((2*pi*(AO1a-AO2a)*tAOM)/2) - pi/4).*Ew(w,wLa,sigmaL).*symm01b(w,w01Bs,C01Bs,gammaA)).^2) .* cos(t21*(w01Bs-wR)) + ...
    ((-sin(phi+pi/2).*sin(((2*pi*(AO1a-AO2a)*tAOM)/2) - pi/4).*Ew(w,wLa,sigmaL).*symm00(w,w00s,C00s,gammaA)).^2 ).* cos(t21*(w00s-wR)) + ...
    ((-sin(phi+pi/2).*sin(((2*pi*(AO1a-AO2a)*tAOM)/2) - pi/4).*Ew(w,wLa,sigmaL).*symm01a(w,w01As,C01As,gammaA)).^2) .* cos(t21*(w01As-wR)) + ...
    ((-sin(phi+pi/2).*sin(((2*pi*(AO1a-AO2a)*tAOM)/2) - pi/4).*Ew(w,wLa,sigmaL).*symm01b(w,w01Bs,C01Bs,gammaA)).^2) .* cos(t21*(w01Bs-wR)) );


symmModY = @(w,tAOM,t21,phi) 2.*exp(-abs(t21)/TauDephase).*( ...
    ((cos(phi+pi/2).*cos(((2*pi*(AO1a-AO2a)*tAOM)/2) - pi/4).*Ew(w,wLa,sigmaL).*symm00(w,w00s,C00s,gammaA)).^2 ).* sin(t21*(w00s-wR)) + ...
    ((cos(phi+pi/2).*cos(((2*pi*(AO1a-AO2a)*tAOM)/2) - pi/4).*Ew(w,wLa,sigmaL).*symm01a(w,w01As,C01As,gammaA)).^2) .* sin(t21*(w01As-wR)) + ...
    ((cos(phi+pi/2).*cos(((2*pi*(AO1a-AO2a)*tAOM)/2) - pi/4).*Ew(w,wLa,sigmaL).*symm01b(w,w01Bs,C01Bs,gammaA)).^2) .* sin(t21*(w01Bs-wR)) + ...
    ((-sin(phi+pi/2).*sin(((2*pi*(AO1a-AO2a)*tAOM)/2) - pi/4).*Ew(w,wLa,sigmaL).*symm00(w,w00s,C00s,gammaA)).^2 ).* sin(t21*(w00s-wR)) + ...
    ((-sin(phi+pi/2).* sin(((2*pi*(AO1a-AO2a)*tAOM)/2)- pi/4).*Ew(w,wLa,sigmaL).*symm01a(w,w01As,C01As,gammaA)).^2) .* sin(t21*(w01As-wR)) + ...
    ((-sin(phi+pi/2).*sin(((2*pi*(AO1a-AO2a)*tAOM)/2) - pi/4).*Ew(w,wLa,sigmaL).*symm01b(w,w01Bs,C01Bs,gammaA)).^2) .* sin(t21*(w01Bs-wR)) );


antiModX = @(w,tAOM,t21,phi) 2.*exp(-abs(t21)/TauDephase).*( ...
    ((-sin(phi).*sin(((2*pi*(AO1a-AO2a)*tAOM)/2) - pi/4).*Ew(w,wLa,sigmaL).*anti00(w,w00a,C00a,gammaA)).^2) .* cos(t21*(w00a-wR)) + ...
    ((-sin(phi).*sin(((2*pi*(AO1a-AO2a)*tAOM)/2) - pi/4).*Ew(w,wLa,sigmaL).*anti01a(w,w01Aa,C01Aa,gammaA)).^2) .* cos(t21*(w01Aa-wR)) + ...
    ((-sin(phi).*sin(((2*pi*(AO1a-AO2a)*tAOM)/2) - pi/4).*Ew(w,wLa,sigmaL).*anti01b(w,w01Ba,C01Ba,gammaA)).^2) .* cos(t21*(w01Ba-wR)) + ...
    ((cos(phi).*cos(((2*pi*(AO1a-AO2a)*tAOM)/2) - pi/4).*Ew(w,wLa,sigmaL).*anti00(w,w00a,C00a,gammaA)).^2) .* cos(t21*(w00a-wR)) + ...
    ((cos(phi).*cos(((2*pi*(AO1a-AO2a)*tAOM)/2) - pi/4).*Ew(w,wLa,sigmaL).*anti01a(w,w01Aa,C01Aa,gammaA)).^2) .* cos(t21*(w01Aa-wR)) + ...
    ((cos(phi).*cos(((2*pi*(AO1a-AO2a)*tAOM)/2) - pi/4).*Ew(w,wLa,sigmaL).*anti01b(w,w01Ba,C01Ba,gammaA)).^2) .* cos(t21*(w01Ba-wR)) ) ;


antiModY = @(w,tAOM,t21,phi) 2.*exp(-abs(t21)/TauDephase).*( ...
    ((-sin(phi).*sin(((2*pi*(AO1a-AO2a)*tAOM)/2) - pi/4).*Ew(w,wLa,sigmaL).*anti00(w,w00a,C00a,gammaA)).^2) .* sin(t21*(w00a-wR)) + ...
    ((-sin(phi).*sin(((2*pi*(AO1a-AO2a)*tAOM)/2) - pi/4).*Ew(w,wLa,sigmaL).*anti01a(w,w01Aa,C01Aa,gammaA)).^2) .* sin(t21*(w01Aa-wR)) + ...
    ((-sin(phi).*sin(((2*pi*(AO1a-AO2a)*tAOM)/2) - pi/4).*Ew(w,wLa,sigmaL).*anti01b(w,w01Ba,C01Ba,gammaA)).^2) .* sin(t21*(w01Ba-wR)) + ...
    ((cos(phi).*cos(((2*pi*(AO1a-AO2a)*tAOM)/2) - pi/4).*Ew(w,wLa,sigmaL).*anti00(w,w00a,C00a,gammaA)).^2) .* sin(t21*(w00a-wR)) + ...
    ((cos(phi).*cos(((2*pi*(AO1a-AO2a)*tAOM)/2) - pi/4).*Ew(w,wLa,sigmaL).*anti01a(w,w01Aa,C01Aa,gammaA)).^2) .* sin(t21*(w01Aa-wR)) + ...
    ((cos(phi).*cos(((2*pi*(AO1a-AO2a)*tAOM)/2) - pi/4).*Ew(w,wLa,sigmaL).*anti01b(w,w01Ba,C01Ba,gammaA)).^2) .* sin(t21*(w01Ba-wR)) ) ;

%     interference terms WITHOUT polz rotation 
symmModXNR = @(w,tAOM,t21) 2.*exp(-abs(t21)/TauDephase).*( ...
    ((Ew(w,wLa,sigmaL).*symm00(w,w00s,C00s,gammaA)).^2) .* cos(t21*(w00s-wR)) + ...
    ((Ew(w,wLa,sigmaL).*symm01a(w,w01As,C01As,gammaA)).^2) .* cos(t21*(w01As-wR)) + ...
    ((Ew(w,wLa,sigmaL).*symm01b(w,w01Bs,C01Bs,gammaA)).^2) .* cos(t21*(w01Bs-wR)) ) ;

antiModXNR = @(w,tAOM,t21) 2.*exp(-abs(t21)/TauDephase).*( ...
    ((Ew(w,wLa,sigmaL).*anti00(w,w00a,C00a,gammaA)).^2) .* cos(t21*(w00a-wR)) + ...
    ((Ew(w,wLa,sigmaL).*anti01a(w,w01Aa,C01Aa,gammaA)).^2) .* cos(t21*(w01Aa-wR)) + ...
    ((Ew(w,wLa,sigmaL).*anti01b(w,w01Ba,C01Ba,gammaA)).^2) .* cos(t21*(w01Ba-wR)) ) ;

symmModYNR = @(w,tAOM,t21) 2.*exp(-abs(t21)/TauDephase).*( ...
    ((Ew(w,wLa,sigmaL).*symm00(w,w00s,C00s,gammaA)).^2) .* sin(t21*(w00s-wR)) + ...
    ((Ew(w,wLa,sigmaL).*symm01a(w,w01As,C01As,gammaA)).^2) .* sin(t21*(w01As-wR)) + ...
    ((Ew(w,wLa,sigmaL).*symm01b(w,w01Bs,C01Bs,gammaA)).^2) .* sin(t21*(w01Bs-wR)) ) ;

antiModYNR = @(w,tAOM,t21) 2.*exp(-abs(t21)/TauDephase).*( ...
    ((Ew(w,wLa,sigmaL).*anti00(w,w00a,C00a,gammaA)).^2) .* sin(t21*(w00a-wR)) + ...
    ((Ew(w,wLa,sigmaL).*anti01a(w,w01Aa,C01Aa,gammaA)).^2) .* sin(t21*(w01Aa-wR)) + ...
    ((Ew(w,wLa,sigmaL).*anti01b(w,w01Ba,C01Ba,gammaA)).^2) .* sin(t21*(w01Ba-wR)) ) ;


%     background terms - 8/17/21 : going to include a sin(phi0) and
%     cos(phi0) for molecualr orientation in both the Y and X quads
%     respectively upon advice from Andy 
symmBkgdDemodX = @(w,tAOM,phi) 2.*abs(cos(phi)).* (...
            ((cos(phi+pi/2).*cos(((2*pi*(AO1a-AO2a)*tAOM)/2) - pi/4).*Ew(w,wLa,sigmaL).*symm00(w,w00s,C00s,gammaA)).^2 ) + ...
            ((cos(phi+pi/2).*cos(((2*pi*(AO1a-AO2a)*tAOM)/2) - pi/4).*Ew(w,wLa,sigmaL).*symm01a(w,w01As,C01As,gammaA)).^2)  + ...
            ((cos(phi+pi/2).*cos(((2*pi*(AO1a-AO2a)*tAOM)/2) - pi/4).*Ew(w,wLa,sigmaL).*symm01b(w,w01Bs,C01Bs,gammaA)).^2)  + ...
            ((-sin(phi+pi/2).*sin(((2*pi*(AO1a-AO2a)*tAOM)/2) - pi/4).*Ew(w,wLa,sigmaL).*symm00(w,w00s,C00s,gammaA)).^2 ) + ...
            ((-sin(phi+pi/2).*sin(((2*pi*(AO1a-AO2a)*tAOM)/2) - pi/4).*Ew(w,wLa,sigmaL).*symm01a(w,w01As,C01As,gammaA)).^2) + ...
            ((-sin(phi+pi/2).*sin(((2*pi*(AO1a-AO2a)*tAOM)/2) - pi/4).*Ew(w,wLa,sigmaL).*symm01b(w,w01Bs,C01Bs,gammaA)).^2) );
        
antiBkgdDemodX =  @(w,tAOM,phi) 2.*abs(cos(phi)).*( ...
             ((-sin(phi).*sin(((2*pi*(AO1a-AO2a)*tAOM)/2) - pi/4).*Ew(w,wLa,sigmaL).*anti00(w,w00a,C00a,gammaA)).^2)  + ...
            ((-sin(phi).*sin(((2*pi*(AO1a-AO2a)*tAOM)/2) - pi/4).*Ew(w,wLa,sigmaL).*anti01a(w,w01Aa,C01Aa,gammaA)).^2) + ...
            ((-sin(phi).*sin(((2*pi*(AO1a-AO2a)*tAOM)/2) - pi/4).*Ew(w,wLa,sigmaL).*anti01b(w,w01Ba,C01Ba,gammaA)).^2) + ...
            ((cos(phi).*cos(((2*pi*(AO1a-AO2a)*tAOM)/2) - pi/4).*Ew(w,wLa,sigmaL).*anti00(w,w00a,C00a,gammaA)).^2) + ...
            ((cos(phi).*cos(((2*pi*(AO1a-AO2a)*tAOM)/2) - pi/4).*Ew(w,wLa,sigmaL).*anti01a(w,w01Aa,C01Aa,gammaA)).^2) + ...
            ((cos(phi).*cos(((2*pi*(AO1a-AO2a)*tAOM)/2) - pi/4).*Ew(w,wLa,sigmaL).*anti01b(w,w01Ba,C01Ba,gammaA)).^2) ) ;
        
symmBkgdDemodY = @(w,tAOM,phi) 2.*abs(sin(phi)).* (...
            ((cos(phi+pi/2).*cos(((2*pi*(AO1a-AO2a)*tAOM)/2) - pi/4).*Ew(w,wLa,sigmaL).*symm00(w,w00s,C00s,gammaA)).^2 ) + ...
            ((cos(phi+pi/2).*cos(((2*pi*(AO1a-AO2a)*tAOM)/2) - pi/4).*Ew(w,wLa,sigmaL).*symm01a(w,w01As,C01As,gammaA)).^2)  + ...
            ((cos(phi+pi/2).*cos(((2*pi*(AO1a-AO2a)*tAOM)/2) - pi/4).*Ew(w,wLa,sigmaL).*symm01b(w,w01Bs,C01Bs,gammaA)).^2)  + ...
            ((-sin(phi+pi/2).*sin(((2*pi*(AO1a-AO2a)*tAOM)/2) - pi/4).*Ew(w,wLa,sigmaL).*symm00(w,w00s,C00s,gammaA)).^2 ) + ...
            ((-sin(phi+pi/2).*sin(((2*pi*(AO1a-AO2a)*tAOM)/2) - pi/4).*Ew(w,wLa,sigmaL).*symm01a(w,w01As,C01As,gammaA)).^2) + ...
            ((-sin(phi+pi/2).*sin(((2*pi*(AO1a-AO2a)*tAOM)/2) - pi/4).*Ew(w,wLa,sigmaL).*symm01b(w,w01Bs,C01Bs,gammaA)).^2) );
        
antiBkgdDemodY =  @(w,tAOM,phi) 2.*abs(sin(phi)).*( ...
             ((-sin(phi).*sin(((2*pi*(AO1a-AO2a)*tAOM)/2) - pi/4).*Ew(w,wLa,sigmaL).*anti00(w,w00a,C00a,gammaA)).^2)  + ...
            ((-sin(phi).*sin(((2*pi*(AO1a-AO2a)*tAOM)/2) - pi/4).*Ew(w,wLa,sigmaL).*anti01a(w,w01Aa,C01Aa,gammaA)).^2) + ...
            ((-sin(phi).*sin(((2*pi*(AO1a-AO2a)*tAOM)/2) - pi/4).*Ew(w,wLa,sigmaL).*anti01b(w,w01Ba,C01Ba,gammaA)).^2) + ...
            ((cos(phi).*cos(((2*pi*(AO1a-AO2a)*tAOM)/2) - pi/4).*Ew(w,wLa,sigmaL).*anti00(w,w00a,C00a,gammaA)).^2) + ...
            ((cos(phi).*cos(((2*pi*(AO1a-AO2a)*tAOM)/2) - pi/4).*Ew(w,wLa,sigmaL).*anti01a(w,w01Aa,C01Aa,gammaA)).^2) + ...
            ((cos(phi).*cos(((2*pi*(AO1a-AO2a)*tAOM)/2) - pi/4).*Ew(w,wLa,sigmaL).*anti01b(w,w01Ba,C01Ba,gammaA)).^2) ) ;
        
        %     background terms WITHOUT  polz rotation - these should not
        %     appear in the demodulation, but will carry it through for now
symmBkgdDemodXNR = @(w,tAOM) 2.* ...
            ((Ew(w,wLa,sigmaL).*symm00(w,w00s,C00s,gammaA)).^2 ...
            + (Ew(w,wLa,sigmaL).*symm01a(w,w01As,C01As,gammaA)).^2 ...
            + (Ew(w,wLa,sigmaL).*symm01b(w,w01Bs,C01Bs,gammaA)).^2) ;
        
antiBkgdDemodXNR =  @(w,tAOM) 2.* ...
             ((Ew(w,wLa,sigmaL).*anti00(w,w00a,C00a,gammaA)).^2 ...
            + (Ew(w,wLa,sigmaL).*anti01a(w,w01Aa,C01Aa,gammaA)).^2 ...
            + (Ew(w,wLa,sigmaL).*anti01b(w,w01Ba,C01Ba,gammaA)).^2) ;
        
symmBkgdDemodYNR = @(w,tAOM) 2.* ...
            ((Ew(w,wLa,sigmaL).*symm00(w,w00s,C00s,gammaA)).^2 ...
            + (Ew(w,wLa,sigmaL).*symm01a(w,w01As,C01As,gammaA)).^2 ...
            + (Ew(w,wLa,sigmaL).*symm01b(w,w01Bs,C01Bs,gammaA)).^2) ;
        
antiBkgdDemodYNR =  @(w,tAOM) 2.* ...
             ((Ew(w,wLa,sigmaL).*anti00(w,w00a,C00a,gammaA)).^2 ...
            + (Ew(w,wLa,sigmaL).*anti01a(w,w01Aa,C01Aa,gammaA)).^2 ...
            + (Ew(w,wLa,sigmaL).*anti01b(w,w01Ba,C01Ba,gammaA)).^2) ;


numPlots=length(t21Axis);
    DimRow= floor(sqrt(numPlots)); 
    DimCol= floor(sqrt(numPlots)); 
    
    if ((DimRow*DimCol)< numPlots)
        DimRow=DimRow+1;
        if(DimRow*DimCol)< numPlots
            DimCol=DimCol+1; 
        end
    end
    
symmInfXquadNRAlt=zeros(1,numel(t21Axis));
symmInfYquadNRAlt=zeros(1,numel(t21Axis));

%  with polz rotation   
symmInfXquad=zeros(1,numel(t21Axis));
antiInfXquad=zeros(1,numel(t21Axis));
symmBkgdXquad=zeros(1,numel(t21Axis));
antiBkgdXquad=zeros(1,numel(t21Axis));

symmInfYquad=zeros(1,numel(t21Axis));
antiInfYquad=zeros(1,numel(t21Axis));
symmBkgdYquad=zeros(1,numel(t21Axis));
antiBkgdYquad=zeros(1,numel(t21Axis));

% no polz rotation terms
symmInfXquadNR=zeros(1,numel(t21Axis));
antiInfXquadNR=zeros(1,numel(t21Axis));
symmBkgdXquadNR=zeros(1,numel(t21Axis));
antiBkgdXquadNR=zeros(1,numel(t21Axis));

symmInfYquadNR=zeros(1,numel(t21Axis));
antiInfYquadNR=zeros(1,numel(t21Axis));
symmBkgdYquadNR=zeros(1,numel(t21Axis));
antiBkgdYquadNR=zeros(1,numel(t21Axis));

% with polz rotation
symmModAmp=zeros(numel(t21Axis),numel(Taxis));
antiModAmp=zeros(numel(t21Axis),numel(Taxis));

symmBkgdAmp=zeros(numel(t21Axis),numel(Taxis));
antiBkgdAmp=zeros(numel(t21Axis),numel(Taxis));

symmTotAmp=zeros(numel(t21Axis),numel(Taxis));
antiTotAmp=zeros(numel(t21Axis),numel(Taxis));

% arrays to hodl the value under no polarization rotation
symmModAmpNR=zeros(numel(t21Axis),numel(Taxis));
antiModAmpNR=zeros(numel(t21Axis),numel(Taxis));

symmBkgdAmpNR=zeros(numel(t21Axis),numel(Taxis));
antiBkgdAmpNR=zeros(numel(t21Axis),numel(Taxis));

symmTotAmpNR=zeros(numel(t21Axis),numel(Taxis));
antiTotAmpNR=zeros(numel(t21Axis),numel(Taxis));

for m=1:numel(t21Axis)
%     interference terms with polz rotation 
        
% % with rotation
% symmInfXquadAlt(m)= (integral2(symmModXalt,100e12,1000e12,0,inf));
curT21=t21Axis(m); 

if sumComp
% symmInfYquadNRAlt(m)= sum(symmModYNR(wAxis,Taxis)); 
% symmInfXquadNRAlt(m)= sum(symmModXNR(wAxis,Taxis)); 
symmInfXquad(m)= sum(symmModX(wAxis,Taxis,curT21,phi0));
antiInfXquad(m)= sum(antiModX(wAxis,Taxis,curT21,phi0));
symmBkgdXquad(m)= sum(symmBkgdDemodX(wAxis,Taxis,phi0));
antiBkgdXquad(m)= sum(antiBkgdDemodX(wAxis,Taxis,phi0));

symmInfYquad(m)=sum(symmModY(wAxis,Taxis,curT21,phi0));
antiInfYquad(m)=sum(antiModY(wAxis,Taxis,curT21,phi0));
symmBkgdYquad(m)=sum(symmBkgdDemodY(wAxis,Taxis,phi0));
antiBkgdYquad(m)=sum(antiBkgdDemodY(wAxis,Taxis,phi0));

% without rotation
symmInfXquadNR(m)= sum(symmModXNR(wAxis,Taxis,curT21));
antiInfXquadNR(m)= sum(antiModXNR(wAxis,Taxis,curT21));
symmBkgdXquadNR(m)= sum(symmBkgdDemodXNR(wAxis,Taxis));
antiBkgdXquadNR(m)= sum(antiBkgdDemodXNR(wAxis,Taxis));

symmInfYquadNR(m)=sum(symmModYNR(wAxis,Taxis,curT21));
antiInfYquadNR(m)=sum(antiModYNR(wAxis,Taxis,curT21));
symmBkgdYquadNR(m)=sum(symmBkgdDemodYNR(wAxis,Taxis));
antiBkgdYquadNR(m)=sum(antiBkgdDemodYNR(wAxis,Taxis));
end

if integrateComp
% integrate THE ABOVE FUNCTIONS ARE 1 VALUE OF T21 OVER THE AOM CYCLE AND
% RELEVANT FREQUENCIES 

% NOTE: the use of integral2() may be problematic to acheive this sum. IT
% is much faster than the discrete sum() but has shown differences in peak
% position upon FFT compared to sum(). 
symmInfXquad(m)= (integral2(symmModX,400e12,800e12,0,inf));
% symmInfXquad(m)= (integral2(symmModX,0,inf,100e12,1000e12));
antiInfXquad(m)= (integral2(antiModX,100e12,1000e12,0,inf));
symmBkgdXquad(m)= (integral2(symmBkgdDemodX,100e12,1000e12,0,inf));
antiBkgdXquad(m)= (integral2(antiBkgdDemodX,100e12,1000e12,0,inf));

symmInfYquad(m)=(integral2(symmModY,100e12,1000e12,0,inf));
antiInfYquad(m)=(integral2(antiModY,100e12,1000e12,0,inf));
symmBkgdYquad(m)=(integral2(symmBkgdDemodY,100e12,1000e12,0,inf));
antiBkgdYquad(m)=(integral2(antiBkgdDemodY,100e12,1000e12,0,inf));

% without rotation
symmInfXquadNR(m)= (integral2(symmModXNR,100e12,1000e12,0,inf));
antiInfXquadNR(m)= (integral2(antiModXNR,100e12,1000e12,0,inf));
symmBkgdXquadNR(m)= (integral2(symmBkgdDemodXNR,100e12,1000e12,0,inf));
antiBkgdXquadNR(m)= (integral2(antiBkgdDemodXNR,100e12,1000e12,0,inf));

symmInfYquadNR(m)=(integral2(symmModYNR,100e12,1000e12,0,inf));
antiInfYquadNR(m)=(integral2(antiModYNR,100e12,1000e12,0,inf));
symmBkgdYquadNR(m)=(integral2(symmBkgdDemodYNR,100e12,1000e12,0,inf));
antiBkgdYquadNR(m)=(integral2(antiBkgdDemodYNR,100e12,1000e12,0,inf));
end

for j=1:numel(Taxis)
    % symmetric manifold background term (we'll assume along the V axis, so
    % assign the cos(A01-AO2) to it and the sin() to the H axis
    symmBkgd =  2*((cos(((2*pi*(AO1a-AO2a)*Taxis(j))/2) - pi/4).*Ew(wAxis,wLa,sigmaL).*symm00(wAxis,w00s,C00s,gammaA)).^2 ...
            + (cos(((2*pi*(AO1a-AO2a)*Taxis(j))/2) - pi/4).*Ew(wAxis,wLa,sigmaL).*symm01a(wAxis,w01As,C01As,gammaA)).^2 ...
            + (cos(((2*pi*(AO1a-AO2a)*Taxis(j))/2) - pi/4).*Ew(wAxis,wLa,sigmaL).*symm01b(wAxis,w01Bs,C01Bs,gammaA)).^2) ;
        
    antiBkgd =  2*((-sin(((2*pi*(AO1a-AO2a)*Taxis(j))/2) - pi/4).*Ew(wAxis,wLa,sigmaL).*anti00(wAxis,w00a,C00a,gammaA)).^2 ...
            + (-sin(((2*pi*(AO1a-AO2a)*Taxis(j))/2) - pi/4).*Ew(wAxis,wLa,sigmaL).*anti01a(wAxis,w01Aa,C01Aa,gammaA)).^2 ...
            + (-sin(((2*pi*(AO1a-AO2a)*Taxis(j))/2) - pi/4).*Ew(wAxis,wLa,sigmaL).*anti01b(wAxis,w01Ba,C01Ba,gammaA)).^2) ;
        
    symmBkgdNR =  2*((Ew(wAxis,wLa,sigmaL).*symm00(wAxis,w00s,C00s,gammaA)).^2 ...
            + (Ew(wAxis,wLa,sigmaL).*symm01a(wAxis,w01As,C01As,gammaA)).^2 ...
            + (Ew(wAxis,wLa,sigmaL).*symm01b(wAxis,w01Bs,C01Bs,gammaA)).^2) ;
        
    antiBkgdNR =  2*((Ew(wAxis,wLa,sigmaL).*anti00(wAxis,w00a,C00a,gammaA)).^2 ...
            + (Ew(wAxis,wLa,sigmaL).*anti01a(wAxis,w01Aa,C01Aa,gammaA)).^2 ...
            + (Ew(wAxis,wLa,sigmaL).*anti01b(wAxis,w01Ba,C01Ba,gammaA)).^2) ;

    % symmetric and antisymmetric modulated terms 
   symmMod = 2*exp(-abs(t21Axis(m))/TauDephase).*( ...
   (cos(((2*pi*(AO1a-AO2a)*Taxis(j))/2) - pi/4).*Ew(wAxis,wLa,sigmaL).*symm00(wAxis,w00s,C00s,gammaA)).^2 .* cos(2*pi*(AO1a-AO2a)*Taxis(j) - t21Axis(m)*w00s) + ...
    (cos(((2*pi*(AO1a-AO2a)*Taxis(j))/2) - pi/4).*Ew(wAxis,wLa,sigmaL).*symm01a(wAxis,w01As,C01As,gammaA)).^2 .* cos(2*pi*(AO1a-AO2a)*Taxis(j) - t21Axis(m)*w01As) + ...
    (cos(((2*pi*(AO1a-AO2a)*Taxis(j))/2) - pi/4).*Ew(wAxis,wLa,sigmaL).*symm01b(wAxis,w01Bs,C01Bs,gammaA)).^2 .* cos(2*pi*(AO1a-AO2a)*Taxis(j) - t21Axis(m)*w01Bs) ) ;


  antiMod = 2*exp(-abs(t21Axis(m))/TauDephase).*( ...
   (-sin(((2*pi*(AO1a-AO2a)*Taxis(j))/2) - pi/4).*Ew(wAxis,wLa,sigmaL).*anti00(wAxis,w00a,C00a,gammaA)).^2 * cos(2*pi*(AO1a-AO2a)*Taxis(j) - t21Axis(m)*w00a) + ...
    (-sin(((2*pi*(AO1a-AO2a)*Taxis(j))/2) - pi/4).*Ew(wAxis,wLa,sigmaL).*anti01a(wAxis,w01Aa,C01Aa,gammaA)).^2 * cos(2*pi*(AO1a-AO2a)*Taxis(j) - t21Axis(m)*w01Aa) + ...
    (-sin(((2*pi*(AO1a-AO2a)*Taxis(j))/2) - pi/4).*Ew(wAxis,wLa,sigmaL).*anti01b(wAxis,w01Ba,C01Ba,gammaA)).^2 * cos(2*pi*(AO1a-AO2a)*Taxis(j) - t21Axis(m)*w01Ba)) ;

  symmModNR = 2*exp(-abs(t21Axis(m))/TauDephase).*( ...
   (Ew(wAxis,wLa,sigmaL).*symm00(wAxis,w00s,C00s,gammaA)).^2 .* cos(2*pi*(AO1a-AO2a)*Taxis(j) - t21Axis(m)*w00s) + ...
    (Ew(wAxis,wLa,sigmaL).*symm01a(wAxis,w01As,C01As,gammaA)).^2 .* cos(2*pi*(AO1a-AO2a)*Taxis(j) - t21Axis(m)*w01As) + ...
    (Ew(wAxis,wLa,sigmaL).*symm01b(wAxis,w01Bs,C01Bs,gammaA)).^2 .* cos(2*pi*(AO1a-AO2a)*Taxis(j) - t21Axis(m)*w01Bs) ) ;


   antiModNR = 2*exp(-abs(t21Axis(m))/TauDephase).*( ...
   (Ew(wAxis,wLa,sigmaL).*anti00(wAxis,w00a,C00a,gammaA)).^2 * cos(2*pi*(AO1a-AO2a)*Taxis(j) - t21Axis(m)*w00a) + ...
    (Ew(wAxis,wLa,sigmaL).*anti01a(wAxis,w01Aa,C01Aa,gammaA)).^2 * cos(2*pi*(AO1a-AO2a)*Taxis(j) - t21Axis(m)*w01Aa) + ...
    (Ew(wAxis,wLa,sigmaL).*anti01b(wAxis,w01Ba,C01Ba,gammaA)).^2 * cos(2*pi*(AO1a-AO2a)*Taxis(j) - t21Axis(m)*w01Ba)) ;

    symmModAmp(m,j)= sum(symmMod);
    antiModAmp(m,j)= sum(antiMod);

    symmBkgdAmp(m,j)= sum(symmBkgd);
    antiBkgdAmp(m,j)= sum(antiBkgd);

    symmTotAmp(m,j)= (sum(symmMod)+sum(symmBkgd));
    antiTotAmp(m,j)=(sum(antiMod)+sum(antiBkgd));
    
    
    symmModAmpNR(m,j)= sum(symmModNR); 
    antiModAmpNR(m,j)= sum(antiModNR);

    symmBkgdAmpNR(m,j)= sum(symmBkgdNR);
    antiBkgdAmpNR(m,j)= sum(antiBkgdNR);

    symmTotAmpNR(m,j)= (sum(symmModNR) + sum(symmBkgdNR)); 
    antiTotAmpNR(m,j)= (sum(antiModNR) + sum(antiBkgdNR)); 
    
    
end

  

end
toc
if showTotalPlots
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
plot(Taxis,symmBkgdAmp')
title('Symmetric - Background terms');
xlabel('AOM clock cycle');
ylabel('Mag (arb)');
% legend('Modulated','Background');

figure(7)
plot(Taxis,antiBkgdAmp')
title('Antisymmetric - Background terms');
xlabel('AOM clock cycle');
ylabel('Mag (arb)');
% legend('Modulated','Background');
end

if showSubplots
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
    title(['Total Linear Signal - t21=' num2str(t21Axis(n)*1e15) 'fs']);
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

figure(16)
for n=1:numel(t21Axis)
    subplot(DimRow,DimCol,n)
    plot(Taxis,(antiBkgdAmpNR(n,:)+symmBkgdAmpNR(n,:)),'--')
    hold on
    plot(Taxis,symmBkgdAmpNR(n,:))
    plot(Taxis,antiBkgdAmpNR(n,:))
    title(['No Polz - Background Linear Signal - t21=' num2str(t21Axis(n)*1e15) 'fs']);
    xlabel('AOM clock cycle');
    ylabel('Mag (arb)');
    legend('Total Background','Symm Bkgd','Anti Bkgd');

end

figure(17)
for n=1:numel(t21Axis)
    subplot(DimRow,DimCol,n)
    plot(Taxis,(antiModAmpNR(n,:)+symmModAmpNR(n,:)+antiBkgdAmpNR(n,:)+symmBkgdAmpNR(n,:)),'--')
    hold on
    plot(Taxis,symmModAmpNR(n,:) + symmBkgdAmpNR(n,:))
    plot(Taxis,antiModAmpNR(n,:) + antiBkgdAmpNR(n,:))
    title(['No Polz - Total Linear Signal - t21=' num2str(t21Axis(n)*1e15) 'fs']);
    xlabel('AOM clock cycle');
    ylabel('Mag (arb)');
    legend('Total Signal','Symm Mod + Bkgd','Anti Mod + Bkgd');
end

figure(18)
for n=1:numel(t21Axis)
    subplot(DimRow,DimCol,n)
    plot(Taxis,(antiModAmpNR(n,:)+symmModAmpNR(n,:)),'--')
    hold on
    plot(Taxis,symmModAmpNR(n,:))
    plot(Taxis,antiModAmpNR(n,:))
    title(['Modulated Linear Signal - t21=' num2str(t21Axis(n)*1e15) 'fs']);
    xlabel('AOM clock cycle');
    ylabel('Mag (arb)');
    legend('Total Mod Signal','Symm Mod','Anti Mod');
end


end

if showTotalPlots
    
%     plots for the rotating case
figure(13)
plot(t21Axis*1e15,symmInfXquad,t21Axis*1e15,antiInfXquad,t21Axis*1e15,symmInfYquad,t21Axis*1e15,antiInfYquad)
title('Symm and Anti X and Y Quadratures - interference terms'); 
xlabel('t21 (fs)');
ylabel('Demodulation Mag (arb)'); 
legend('Symm Inf Xquad','Anti Inf Xquad','Symm Inf Yquad','Anti Inf Yquad'); 

figure(14)
plot(t21Axis*1e15,(symmBkgdXquad),t21Axis*1e15,(antiBkgdXquad),...
    t21Axis*1e15,(symmBkgdYquad),t21Axis*1e15,(antiBkgdYquad))
title('Symm and Anti X and Y Quadratures - background terms'); 
xlabel('t21 (fs)');
ylabel('Demodulation Mag (arb)'); 
legend('Symm Bkgd Xquad','Anti Bkgd Xquad','Symm Bkgd Yquad','Anti Bkgd Yquad'); 

figure(15)
plot(t21Axis*1e15,(symmInfXquad+symmBkgdXquad+antiInfXquad+antiBkgdXquad),...
    t21Axis*1e15,(symmInfYquad+symmBkgdYquad+antiInfYquad+antiBkgdYquad))
title('Total X and Y Quadratures - All terms - Symm and Anti'); 
xlabel('t21 (fs)');
ylabel('Demodulation Mag (arb)'); 
legend('Total Xquad','Total Yquad'); 

figure(16)
plot(t21Axis*1e15,(symmInfXquad+symmBkgdXquad),t21Axis*1e15,(antiInfXquad+antiBkgdXquad),...
    t21Axis*1e15,(symmInfYquad+symmBkgdYquad),t21Axis*1e15, (antiInfYquad+antiBkgdYquad))
title('X and Y Quadratures - All terms - Symm and Anti'); 
xlabel('t21 (fs)');
ylabel('Demodulation Mag (arb)'); 
legend('Symm Xquad','Anti Xquad','Symm Yquad','Anti Yquad');



% plots for the NON rotating case
figure(19)
plot(t21Axis*1e15,symmInfXquadNR,t21Axis*1e15,antiInfXquadNR,t21Axis*1e15,symmInfYquadNR,t21Axis*1e15,antiInfYquadNR)
title('No Polz - Symm and Anti X and Y Quadratures - interference terms only'); 
xlabel('t21 (fs)');
ylabel('Demodulation Mag (arb)'); 
legend('Symm Inf Xquad','Anti Inf Xquad','Symm Inf Yquad','Anti Inf Yquad'); 

figure(20)
plot(t21Axis*1e15,(symmBkgdXquadNR),t21Axis*1e15,(antiBkgdXquadNR),...
    t21Axis*1e15,(symmBkgdYquadNR),t21Axis*1e15,(antiBkgdYquadNR))
title('No Polz - Symm and Anti X and Y Quadratures - interference plus background terms'); 
xlabel('t21 (fs)');
ylabel('Demodulation Mag (arb)'); 
legend('Symm Bkgd Xquad','Anti Bkgd Xquad','Symm Bkgd Yquad','Anti Bkgd Yquad'); 

% going to not include the bkgd terms, since they shouldnt contribute in
% the non-rotating case 
figure(21)
plot(t21Axis*1e15,(symmInfXquadNR+antiInfXquadNR),...
    t21Axis*1e15,(symmInfYquadNR+antiInfYquadNR))
title('No Polz - Total X and Y Quadratures - All terms - Symm and Anti'); 
xlabel('t21 (fs)');
ylabel('Demodulation Mag (arb)'); 
legend('Total Xquad','Total Yquad'); 

end
%%
% configure the data for FFT to frequency domain 
time=t21Axis((t21Axis>=0));
totalSig= ((symmInfXquad+symmBkgdXquad+antiInfXquad+antiBkgdXquad)-1i*(symmInfYquad+symmBkgdYquad+antiInfYquad+antiBkgdYquad));
% totalSigNR= ((symmInfXquadNR+symmBkgdXquadNR+antiInfXquadNR+antiBkgdXquadNR)-1i*(symmInfYquadNR+symmBkgdYquadNR+antiInfYquadNR+antiBkgdYquadNR));
totalSigNR= ((symmInfXquadNR+antiInfXquadNR)-1i*(symmInfYquadNR+antiInfYquadNR));
symmSigNR= ((symmInfXquadNR)-1i*(symmInfYquadNR));
symmSigNRAlt= ((symmInfXquadNRAlt)-1i*(symmInfYquadNRAlt));
antiSigNR= ((antiInfXquadNR)-1i*(antiInfYquadNR));

% time(1)=0;
% isolate the portion fo the signal from 0 time onward
totalSig2=totalSig; 

symmSig=((symmInfXquad+symmBkgdXquad)-1i*(symmInfYquad+symmBkgdYquad));
%antiSig=((antiInfXquad+antiBkgdXquad)-1i*(antiInfYquad+antiBkgdYquad));

% totalSig=totalSig(length(time):end); 
% totalSigNR=totalSigNR(length(time):end); 
% symmSigNR=symmSigNR(length(time):end);
% antiSigNR=antiSigNR(length(time):end);
% symmSig=symmSig(length(time):end); 
% symmSigNRAlt=symmSigNRAlt(length(time):end);

NFFT = 2^(nextpow2(length(time))+3) ;
NFFT2= 2^(nextpow2(length(t21Axis))+3) ;

window2= delayedGaussian(t21Axis, 50e-15, 100e-15).*fliplr(delayedGaussian(t21Axis, 50e-15, 100e-15));
window = delayedGaussian(time, 90e-15, 100e-15);

figure(22)
plot(time,window,time,window.*(real(totalSig)./(max(real(totalSig)) )) );

figure(24)
plot(t21Axis,window2,t21Axis,window2.*(real(totalSig2)./(max(real(totalSig2)) )) );

totalSig2=totalSig2.*window2; 
% phase and window the signals
% totalSig=totalSig.*(exp(-1i*angle(totalSig(1))));
totalSigW= window.*totalSig;
totalSigW(1) = 0.5*totalSigW(1); % weird but necessary, some NMR reference somewhere.

symmSigNRAltw=symmSigNRAlt.*window;

symmSigNRw=symmSigNR.*window;
antiSigNRw=antiSigNR.*window;
symmSigNRw(1)=0.5*symmSigNRw(1); 
antiSigNRw(1)=0.5*antiSigNRw(1); 

% totalSigNR=totalSigNR.*(exp(-1i*angle(totalSigNR(1))));
totalSigWNR= window.*totalSigNR;
totalSigWNR(1) = 0.5*totalSigWNR(1); % weird but necessary, some NMR reference somewhere.

symmSigW=window.*symmSig; 

Ts=0.5*abs(time(3)-time(2)); 
Fs=1/Ts;
falt = Fs*(0:(NFFT/2))/NFFT;

f = ceil((-NFFT/2:NFFT/2-1))/(Ts*NFFT); % freq axis in Hz
f2 = ceil((-NFFT2/2:NFFT2/2-1))/(Ts*NFFT2); % freq axis in Hz
% f = ceil((-NFFT/2:NFFT/2-1))/(Ts*NFFT); % freq axis in Hz
Faxis = Fs/2*linspace(-1,1,NFFT);
% f=Faxis;  
% f = f/c0*1e7; % optional conversion into wavenumbers.
% SFuns = fft(totalSigW, NFFT);
SF2 =(fftshift((fft(totalSig2, NFFT2,2)))); % Signal in the frequency domain, FFT shift applied 

SFsymm =(fftshift((fft(symmSigW, NFFT,2)))); % Signal in the frequency domain, FFT shift applied 

SF =(fftshift((fft(totalSigW, NFFT,2)))); % Signal in the frequency domain, FFT shift applied 
% SF = (fft(totalSigW, NFFT)); % Signal in the frequency domain, FFT shift applied 

SFNR = (fftshift(fft(totalSigWNR, NFFT,2))); % Signal in the frequency domain, FFT shift applied 
SFNRsymm = (fftshift(fft(symmSigNRw, NFFT,2))); % Signal in the frequency domain, FFT shift applied 
SFNRsymmAlt = (fftshift(fft(symmSigNRAltw, NFFT,2))); % Signal in the frequency domain, FFT shift applied 
SFNRanti = (fftshift(fft(antiSigNRw, NFFT,2))); % Signal in the frequency domain, FFT shift applied 

figure(23)
plot(f+(wR-(wR-wLa)), real(SF), f+(wR-(wR-wLa)), imag(SF), f+(wR-(wR-wLa)), abs(SF));
hold on
plot(f+(wR-(wR-wLa)), real(SFNR),'--', f+(wR-(wR-wLa)) ,imag(SFNR),'--', f+(wR-(wR-wLa)), abs(SFNR),'--');
% plot(f+(wR-(wR-wLa)), real(SFsymm),'--', f+(wR-(wR-wLa)) ,imag(SFsymm),'--', f+(wR-(wR-wLa)), abs(SFsymm),'--');
xlim([510e12 620e12]);
title('Linear Spectra - Polz and No Polz');
xlabel('Frequency (Hz)');
ylabel('Mag (arb)'); 
legend('Real Polz','Imag Polz','Abs Polz','Real No Polz','Imag No Polz','Abs No Polz'); 
hold off

figure(26)
plot(f+(wR-(wR-wLa)), real(SFNRsymm), f+(wR-(wR-wLa)), imag(SFNRsymm), f+(wR-(wR-wLa)), abs(SFNRsymm));
hold on
plot(f+(wR-(wR-wLa)), real(SFNRanti), f+(wR-(wR-wLa)) ,imag(SFNRanti), f+(wR-(wR-wLa)), abs(SFNRanti));
plot(f+(wR-(wR-wLa)), real(SFNR),'--', f+(wR-(wR-wLa)) ,imag(SFNR),'--', f+(wR-(wR-wLa)), abs(SFNR),'--');
xlim([510e12 620e12]);
title('Linear Spectra - No Polz');
xlabel('Frequency (Hz)');
ylabel('Mag (arb)'); 
legend('Real Symm','Imag Symm','Abs Symm','Real Anti','Imag Anti','Abs Anti','Real Total','Imag Total','Abs Total'); 
hold off

% figure(28)
% % plot(f+(wR-(wR-wLa)), real(SFNRsymmAlt), f+(wR-(wR-wLa)), imag(SFNRsymmAlt), f+(wR-(wR-wLa)), abs(SFNRsymmAlt));
% hold on
% plot(f+(wR-(wR-wLa)), abs(SFNRsymmAlt./(max(SFNRsymmAlt))),f+(wR-(wR-wLa)), abs(SFNRsymm./(max(SFNRsymm))));
% xlim([510e12 620e12]);
% title('Linear Spectra - No Polz - Sum not integral');
% xlabel('Frequency (Hz)');
% ylabel('Mag (arb)'); 
% legend('Real Symm','Imag Symm','Abs Symm'); 
% hold off

% figure(27)
% plot(Faxis+(wR-(wR-wLa)), real(SFNRsymm), Faxis+(wR-(wR-wLa)), imag(SFNRsymm), Faxis+(wR-(wR-wLa)), abs(SFNRsymm));
% hold on
% plot(Faxis+(wR-(wR-wLa)), real(SFNRanti), Faxis+(wR-(wR-wLa)) ,imag(SFNRanti), Faxis+(wR-(wR-wLa)), abs(SFNRanti));
% plot(Faxis+(wR-(wR-wLa)), real(SFNR),'--', Faxis+(wR-(wR-wLa)) ,imag(SFNR),'--', Faxis+(wR-(wR-wLa)), abs(SFNR),'--');
% xlim([510e12 620e12]);
% title('Linear Spectra - No Polz');
% xlabel('Frequency (Hz)');
% ylabel('Mag (arb)'); 
% legend('Real Symm','Imag Symm','Abs Symm','Real Anti','Imag Anti','Abs Anti','Real Total','Imag Total','Abs Total'); 
% hold off

% figure(25)
% plot(f2+(wR), real(SF2), f2+(wR) ,imag(SF2), f2+(wR), abs(SF2));
% xlim([540e12 620e12]);
% title('Linear Spectra - Polz - all time');
% xlabel('Frequency (Hz)');
% ylabel('Mag (arb)'); 
% legend('Real Polz','Imag Polz','Abs Polz'); 
% hold off



% % % % COLLECTION OF FUNCTIONS USED BEFORE WITH 3A AND 4A OF TAMIMI 2020 - FOR
% % % % SOME REASON THEY NEVER FUNCTIONED THE SAME OF THE 3B 4B VERSIONS WHICH DO
% % % % YIELD THE PROPER DOWN SAMPLING 
% % 
% % symmModX = @(w,tAOM) 2.*exp(-abs(t21Axis(m))/TauDephase).*cos(2*pi*(AO1a-AO2a)*tAOM - wR*t21Axis(m)).*exp(-tAOM/intTime).*(2/intTime).*( ...
% %     (cos(((2*pi*(AO1a-AO2a)*tAOM)/2) - pi/4).*Ew(w,wLa,sigmaL).*symm00(w,w00s,C00s,gammaA)).^2 .* cos(2*pi*(AO1a-AO2a)*tAOM - t21Axis(m)*w00s) + ...
% %     (cos(((2*pi*(AO1a-AO2a)*tAOM)/2) - pi/4).*Ew(w,wLa,sigmaL).*symm01a(w,w01As,C01As,gammaA)).^2 .* cos(2*pi*(AO1a-AO2a)*tAOM - t21Axis(m)*w01As) + ...
% %     (cos(((2*pi*(AO1a-AO2a)*tAOM)/2) - pi/4).*Ew(w,wLa,sigmaL).*symm01b(w,w01Bs,C01Bs,gammaA)).^2 .* cos(2*pi*(AO1a-AO2a)*tAOM - t21Axis(m)*w01Bs) ) ;
% % 
% % 
% % antiModX = @(w,tAOM) 2.*exp(-abs(t21Axis(m))/TauDephase).*cos(2*pi*(AO1a-AO2a)*tAOM - wR*t21Axis(m)).*exp(-tAOM/intTime).*(2/intTime).*( ...
% %     (-sin(((2*pi*(AO1a-AO2a)*tAOM)/2) - pi/4).*Ew(w,wLa,sigmaL).*anti00(w,w00a,C00a,gammaA)).^2 .* cos(2*pi*(AO1a-AO2a)*tAOM - t21Axis(m)*w00a) + ...
% %     (-sin(((2*pi*(AO1a-AO2a)*tAOM)/2) - pi/4).*Ew(w,wLa,sigmaL).*anti01a(w,w01Aa,C01Aa,gammaA)).^2 .* cos(2*pi*(AO1a-AO2a)*tAOM - t21Axis(m)*w01Aa) + ...
% %     (-sin(((2*pi*(AO1a-AO2a)*tAOM)/2) - pi/4).*Ew(w,wLa,sigmaL).*anti01b(w,w01Ba,C01Ba,gammaA)).^2 .* cos(2*pi*(AO1a-AO2a)*tAOM - t21Axis(m)*w01Ba) ) ;
% % 
% % symmModY = @(w,tAOM) -2.*exp(-abs(t21Axis(m))/TauDephase).*sin(2*pi*(AO1a-AO2a)*tAOM - wR*t21Axis(m)).*exp(-tAOM/intTime).*(2/intTime).*( ...
% %     (cos(((2*pi*(AO1a-AO2a)*tAOM)/2) - pi/4).*Ew(w,wLa,sigmaL).*symm00(w,w00s,C00s,gammaA)).^2 .* cos(2*pi*(AO1a-AO2a)*tAOM - t21Axis(m)*w00s) + ...
% %     (cos(((2*pi*(AO1a-AO2a)*tAOM)/2) - pi/4).*Ew(w,wLa,sigmaL).*symm01a(w,w01As,C01As,gammaA)).^2 .* cos(2*pi*(AO1a-AO2a)*tAOM - t21Axis(m)*w01As) + ...
% %     (cos(((2*pi*(AO1a-AO2a)*tAOM)/2) - pi/4).*Ew(w,wLa,sigmaL).*symm01b(w,w01Bs,C01Bs,gammaA)).^2 .* cos(2*pi*(AO1a-AO2a)*tAOM - t21Axis(m)*w01Bs) ) ;
% % 
% % antiModY = @(w,tAOM) -2.*exp(-abs(t21Axis(m))/TauDephase).*sin(2*pi*(AO1a-AO2a)*tAOM - wR*t21Axis(m)).*exp(-tAOM/intTime).*(2/intTime).*( ...
% %     (-sin(((2*pi*(AO1a-AO2a)*tAOM)/2) - pi/4).*Ew(w,wLa,sigmaL).*anti00(w,w00a,C00a,gammaA)).^2 .* cos(2*pi*(AO1a-AO2a)*tAOM - t21Axis(m)*w00a) + ...
% %     (-sin(((2*pi*(AO1a-AO2a)*tAOM)/2) - pi/4).*Ew(w,wLa,sigmaL).*anti01a(w,w01Aa,C01Aa,gammaA)).^2 .* cos(2*pi*(AO1a-AO2a)*tAOM - t21Axis(m)*w01Aa) + ...
% %     (-sin(((2*pi*(AO1a-AO2a)*tAOM)/2) - pi/4).*Ew(w,wLa,sigmaL).*anti01b(w,w01Ba,C01Ba,gammaA)).^2 .* cos(2*pi*(AO1a-AO2a)*tAOM - t21Axis(m)*w01Ba) ) ;
% % 
% % %     interference terms WITHOUT polz rotation 
% % symmModXNR = @(w,tAOM) 2.*exp(-abs(t21Axis(m))/TauDephase).*cos(2*pi*(AO1a-AO2a)*tAOM - wR*t21Axis(m)).*exp(-tAOM/intTime).*(2/intTime).*( ...
% %     (Ew(w,wLa,sigmaL).*symm00(w,w00s,C00s,gammaA)).^2 .* cos(2*pi*(AO1a-AO2a)*tAOM - t21Axis(m)*w00s) + ...
% %     (Ew(w,wLa,sigmaL).*symm01a(w,w01As,C01As,gammaA)).^2 .* cos(2*pi*(AO1a-AO2a)*tAOM - t21Axis(m)*w01As) + ...
% %     (Ew(w,wLa,sigmaL).*symm01b(w,w01Bs,C01Bs,gammaA)).^2 .* cos(2*pi*(AO1a-AO2a)*tAOM - t21Axis(m)*w01Bs) ) ;
% % 
% % antiModXNR = @(w,tAOM) 2.*exp(-abs(t21Axis(m))/TauDephase).*cos(2*pi*(AO1a-AO2a)*tAOM - wR*t21Axis(m)).*exp(-tAOM/intTime).*(2/intTime).*( ...
% %     (Ew(w,wLa,sigmaL).*anti00(w,w00a,C00a,gammaA)).^2 .* cos(2*pi*(AO1a-AO2a)*tAOM - t21Axis(m)*w00a) + ...
% %     (Ew(w,wLa,sigmaL).*anti01a(w,w01Aa,C01Aa,gammaA)).^2 .* cos(2*pi*(AO1a-AO2a)*tAOM - t21Axis(m)*w01Aa) + ...
% %     (Ew(w,wLa,sigmaL).*anti01b(w,w01Ba,C01Ba,gammaA)).^2 .* cos(2*pi*(AO1a-AO2a)*tAOM - t21Axis(m)*w01Ba) ) ;
% % 
% % symmModYNR = @(w,tAOM) -2.*exp(-abs(t21Axis(m))/TauDephase).*sin(2*pi*(AO1a-AO2a)*tAOM - wR*t21Axis(m)).*exp(-tAOM/intTime).*(2/intTime).*( ...
% %     (Ew(w,wLa,sigmaL).*symm00(w,w00s,C00s,gammaA)).^2 .* cos(2*pi*(AO1a-AO2a)*tAOM - t21Axis(m)*w00s) + ...
% %     (Ew(w,wLa,sigmaL).*symm01a(w,w01As,C01As,gammaA)).^2 .* cos(2*pi*(AO1a-AO2a)*tAOM - t21Axis(m)*w01As) + ...
% %     (Ew(w,wLa,sigmaL).*symm01b(w,w01Bs,C01Bs,gammaA)).^2 .* cos(2*pi*(AO1a-AO2a)*tAOM - t21Axis(m)*w01Bs) ) ;
% % 
% % antiModYNR = @(w,tAOM) -2.*exp(-abs(t21Axis(m))/TauDephase).*sin(2*pi*(AO1a-AO2a)*tAOM - wR*t21Axis(m)).*exp(-tAOM/intTime).*(2/intTime).*( ...
% %     (Ew(w,wLa,sigmaL).*anti00(w,w00a,C00a,gammaA)).^2 .* cos(2*pi*(AO1a-AO2a)*tAOM - t21Axis(m)*w00a) + ...
% %     (Ew(w,wLa,sigmaL).*anti01a(w,w01Aa,C01Aa,gammaA)).^2 .* cos(2*pi*(AO1a-AO2a)*tAOM - t21Axis(m)*w01Aa) + ...
% %     (Ew(w,wLa,sigmaL).*anti01b(w,w01Ba,C01Ba,gammaA)).^2 .* cos(2*pi*(AO1a-AO2a)*tAOM - t21Axis(m)*w01Ba) ) ;
% % 
% % 
% % %     background terms
% % symmBkgdDemodX = @(w,tAOM) 2.*cos(2*pi*(AO1a-AO2a)*tAOM - wR*t21Axis(m)).*exp(-tAOM/intTime).*(2/intTime).* ...
% %             ((cos(((2*pi*(AO1a-AO2a)*tAOM)/2) - pi/4).*Ew(w,wLa,sigmaL).*symm00(w,w00s,C00s,gammaA)).^2 ...
% %             + (cos(((2*pi*(AO1a-AO2a)*tAOM)/2) - pi/4).*Ew(w,wLa,sigmaL).*symm01a(w,w01As,C01As,gammaA)).^2 ...
% %             + (cos(((2*pi*(AO1a-AO2a)*tAOM)/2) - pi/4).*Ew(w,wLa,sigmaL).*symm01b(w,w01Bs,C01Bs,gammaA)).^2) ;
% %         
% % antiBkgdDemodX =  @(w,tAOM) 2.*cos(2*pi*(AO1a-AO2a)*tAOM - wR*t21Axis(m)).*exp(-tAOM/intTime).*(2/intTime).* ...
% %              ((-sin(((2*pi*(AO1a-AO2a)*tAOM)/2) - pi/4).*Ew(w,wLa,sigmaL).*anti00(w,w00a,C00a,gammaA)).^2 ...
% %             + (-sin(((2*pi*(AO1a-AO2a)*tAOM)/2) - pi/4).*Ew(w,wLa,sigmaL).*anti01a(w,w01Aa,C01Aa,gammaA)).^2 ...
% %             + (-sin(((2*pi*(AO1a-AO2a)*tAOM)/2) - pi/4).*Ew(w,wLa,sigmaL).*anti01b(w,w01Ba,C01Ba,gammaA)).^2) ;
% %         
% % symmBkgdDemodY = @(w,tAOM) -2.*sin(2*pi*(AO1a-AO2a)*tAOM - wR*t21Axis(m)).*exp(-tAOM/intTime).*(2/intTime).* ...
% %             ((cos(((2*pi*(AO1a-AO2a)*tAOM)/2) - pi/4).*Ew(w,wLa,sigmaL).*symm00(w,w00s,C00s,gammaA)).^2 ...
% %             + (cos(((2*pi*(AO1a-AO2a)*tAOM)/2) - pi/4).*Ew(w,wLa,sigmaL).*symm01a(w,w01As,C01As,gammaA)).^2 ...
% %             + (cos(((2*pi*(AO1a-AO2a)*tAOM)/2) - pi/4).*Ew(w,wLa,sigmaL).*symm01b(w,w01Bs,C01Bs,gammaA)).^2) ;
% %         
% % antiBkgdDemodY =  @(w,tAOM) -2.*sin(2*pi*(AO1a-AO2a)*tAOM - wR*t21Axis(m)).*exp(-tAOM/intTime).*(2/intTime).* ...
% %              ((-sin(((2*pi*(AO1a-AO2a)*tAOM)/2) - pi/4).*Ew(w,wLa,sigmaL).*anti00(w,w00a,C00a,gammaA)).^2 ...
% %             + (-sin(((2*pi*(AO1a-AO2a)*tAOM)/2) - pi/4).*Ew(w,wLa,sigmaL).*anti01a(w,w01Aa,C01Aa,gammaA)).^2 ...
% %             + (-sin(((2*pi*(AO1a-AO2a)*tAOM)/2) - pi/4).*Ew(w,wLa,sigmaL).*anti01b(w,w01Ba,C01Ba,gammaA)).^2) ;
% %         
% %         %     background terms WITHOUT  polz rotation
% % symmBkgdDemodXNR = @(w,tAOM) 2.*cos(2*pi*(AO1a-AO2a)*tAOM - wR*t21Axis(m)).*exp(-tAOM/intTime).*(2/intTime).* ...
% %             ((Ew(w,wLa,sigmaL).*symm00(w,w00s,C00s,gammaA)).^2 ...
% %             + (Ew(w,wLa,sigmaL).*symm01a(w,w01As,C01As,gammaA)).^2 ...
% %             + (Ew(w,wLa,sigmaL).*symm01b(w,w01Bs,C01Bs,gammaA)).^2) ;
% %         
% % antiBkgdDemodXNR =  @(w,tAOM) 2.*cos(2*pi*(AO1a-AO2a)*tAOM - wR*t21Axis(m)).*exp(-tAOM/intTime).*(2/intTime).* ...
% %              ((Ew(w,wLa,sigmaL).*anti00(w,w00a,C00a,gammaA)).^2 ...
% %             + (Ew(w,wLa,sigmaL).*anti01a(w,w01Aa,C01Aa,gammaA)).^2 ...
% %             + (Ew(w,wLa,sigmaL).*anti01b(w,w01Ba,C01Ba,gammaA)).^2) ;
% %         
% % symmBkgdDemodYNR = @(w,tAOM) -2.*sin(2*pi*(AO1a-AO2a)*tAOM - wR*t21Axis(m)).*exp(-tAOM/intTime).*(2/intTime).* ...
% %             ((Ew(w,wLa,sigmaL).*symm00(w,w00s,C00s,gammaA)).^2 ...
% %             + (Ew(w,wLa,sigmaL).*symm01a(w,w01As,C01As,gammaA)).^2 ...
% %             + (Ew(w,wLa,sigmaL).*symm01b(w,w01Bs,C01Bs,gammaA)).^2) ;
% %         
% % antiBkgdDemodYNR =  @(w,tAOM) -2.*sin(2*pi*(AO1a-AO2a)*tAOM - wR*t21Axis(m)).*exp(-tAOM/intTime).*(2/intTime).* ...
% %              ((Ew(w,wLa,sigmaL).*anti00(w,w00a,C00a,gammaA)).^2 ...
% %             + (Ew(w,wLa,sigmaL).*anti01a(w,w01Aa,C01Aa,gammaA)).^2 ...
% %             + (Ew(w,wLa,sigmaL).*anti01b(w,w01Ba,C01Ba,gammaA)).^2) ;
% %         


% symmModX = @(w,tAOM) 2.*exp(-abs(t21Axis(m))/TauDephase).*cos(2*pi*(AO1a-AO2a)*tAOM - wR*t21Axis(m)).*exp(-tAOM/intTime).*(2/intTime).* ...
%     (cos(((2*pi*(AO1a-AO2a)*tAOM)/2) - pi/4).*Ew(w,(wLa),sigmaL).* ...
%    (symm00(w,w00s,C00s,gammaA).*cos(2*pi*(AO1a-AO2a)*tAOM - t21Axis(m)*(w00s)) + ...
%     symm01a(w,w01As,C01As,gammaA).*cos(2*pi*(AO1a-AO2a)*tAOM - t21Axis(m)*(w01As)) + ...
%     symm01b(w,w01Bs,C01Bs,gammaA).*cos(2*pi*(AO1a-AO2a)*tAOM - t21Axis(m)*(w01Bs)))).^2 ;
% 
% antiModX = @(w,tAOM) 2.*exp(-abs(t21Axis(m))/TauDephase).*cos(2*pi*(AO1a-AO2a)*tAOM - wR*t21Axis(m)).*exp(-tAOM/intTime).*(2/intTime).* ...
%     (-sin(((2*pi*(AO1a-AO2a)*tAOM)/2) - pi/4).*Ew(w,wLa,sigmaL).* ...
%     (anti00(w,w00a,C00a,gammaA).*cos(2*pi*(AO1a-AO2a)*tAOM - t21Axis(m)*w00a) + ...
%     anti01a(w,w01Aa,C01Aa,gammaA).*cos(2*pi*(AO1a-AO2a)*tAOM - t21Axis(m)*w01Aa) + ...
%     anti01b(w,w01Ba,C01Ba,gammaA).*cos(2*pi*(AO1a-AO2a)*tAOM - t21Axis(m)*w01Ba))).^2 ;
% 
% % same but for the Yquad (note the minus out front to handle this case)
% symmModY = @(w,tAOM) -2.*exp(-abs(t21Axis(m))/TauDephase).*sin(2*pi*(AO1a-AO2a)*tAOM - wR*t21Axis(m)).*exp(-tAOM/intTime).*(2/intTime).* ...
%     (cos(((2*pi*(AO1a-AO2a)*tAOM)/2) - pi/4).*Ew(w,wLa,sigmaL).* ...
%    (symm00(w,w00s,C00s,gammaA).*cos(2*pi*(AO1a-AO2a)*tAOM - t21Axis(m)*w00s) + ...
%     symm01a(w,w01As,C01As,gammaA).*cos(2*pi*(AO1a-AO2a)*tAOM - t21Axis(m)*w01As) + ...
%     symm01b(w,w01Bs,C01Bs,gammaA).*cos(2*pi*(AO1a-AO2a)*tAOM - t21Axis(m)*w01Bs))).^2 ;
% 
% antiModY = @(w,tAOM) -2.*exp(-abs(t21Axis(m))/TauDephase).*sin(2*pi*(AO1a-AO2a)*tAOM - wR*t21Axis(m)).*exp(-tAOM/intTime).*(2/intTime).* ...
%     (-sin(((2*pi*(AO1a-AO2a)*tAOM)/2) - pi/4).*Ew(w,wLa,sigmaL).* ...
%     (anti00(w,w00a,C00a,gammaA).*cos(2*pi*(AO1a-AO2a)*tAOM - t21Axis(m)*w00a) + ...
%     anti01a(w,w01Aa,C01Aa,gammaA).*cos(2*pi*(AO1a-AO2a)*tAOM - t21Axis(m)*w01Aa) + ...
%     anti01b(w,w01Ba,C01Ba,gammaA).*cos(2*pi*(AO1a-AO2a)*tAOM - t21Axis(m)*w01Ba))).^2 ;