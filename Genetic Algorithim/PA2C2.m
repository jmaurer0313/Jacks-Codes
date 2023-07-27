function [C2_sim,time] = PA2C2(P,A,time,yoff,addControlMode,c2Control,K)
clockMode = 0;
plotMode = 0;
switch nargin
    case 0
        [K,A,time,~] = paramSim_8state();
        [P,~] = K2P(K,time);
        yoff = 0;
        clockMode = 1;
        plotMode = 1;
    case 3
        yoff = 0;
        clockMode = 0;
        plotMode = 0;
end
%[N,~,timesteps] = size(P);
[NumStates,~,~] = size(P);
timesteps = length(time);
% Make a row vector of the discrete probability
% x = diag(A) returns a column vector of the main diagonal elements of A.
Peq = diag(P(:,:,end))';

if length(P)~=length(time)
    [P,~,~,~] = K2P(K,time);
end
%Make a  row vector of the FRET Values
% A = linspace(.1,.9,N);

%Subtract off the mean of A
% Amean = sum(A.*Peq);
Peq = reshape(Peq, 1, length(Peq));
A = reshape(A, length(A), 1);
Amean = (sum(A.*Peq'));

Ams = A - Amean;
% disp(['Ams shape: ' num2str(size(Ams))])
%--------------------------------------------------------------------------
% Calculate 2 point TCF (Sums using loops)
%--------------------------------------------------------------------------
if clockMode == 1
    tic
end
% P=eye(length(Ams));
C2_sim = zeros(1,timesteps);
for i = 1:length(Ams)
    for j = 1:length(Ams)       
        C2temp = Ams(j)*reshape(P(j,i,:),[1 timesteps])*Ams(i)*Peq(i);
%                 C2temp = Ams(j)*Ams(i)*Peq(i);       
        C2_sim =  C2_sim + C2temp;
    end
end
C2_sim = C2_sim + yoff;

if addControlMode == 1
   C2_sim = reshape(C2_sim,1,length(C2_sim)) + reshape(c2Control,1,length(c2Control));  
end

if clockMode == 1
    elapsedTime =toc;
    disp(['Takes '  num2str(elapsedTime) ' seconds to calculate C2 for an ' num2str(NumStates) ' state model']);
end
if plotMode == 1
    figure(2)
    set(gcf,'Color','w');
    
    C2_sim_plot = plot(time,C2_sim);
    C2_sim_plot.LineWidth = 2;
    C2_sim_plot.Color = 'r';
    
    title_str = ['Two point time correlation function'];
    title(title_str,'FontSize',18);
    xlabel('\tau (sec)','fontsize',16);
    ylabel('C^{(2)}(\tau)','fontsize',16);
    set(gca,'yscale','linear');
    set(gca,'xscale','log');
    set(gca,'FontSize',14);
    grid on
    axis tight;
    
    drawnow();
end

%--------------------------------------------------------------------------
% Calculate 2 point TCF (with matrix)
%--------------------------------------------------------------------------
%%
% A = linspace(0,1,N);
% Amat = A.*eye(N);
% Amat3D = ones(N,N,length(time));
% Peq2D = ones(N,length(time));
% for i = 1:length(time)
%     Amat3D(:,:,i) = Amat;
%     Peq2D(:,i) = Peq;
% end
% % C2test = sum(Amat*P*Amat*Peq);
% C2test = sum(Amat3D*P*Amat3D*Peq2D);
%
%

