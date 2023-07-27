function [C4_sim,time] = PAK2C4(P,A,K,time,tau2,zoff,addControlMode,c4Control)
clockMode = 0;
plotMode = 0;
% MODIFICATIONS: Changed for loop for 25% increase in speed    
%This is the bottleneck of the calculations

% switch nargin
%     case 0
%         disp('Simulating paramaters for an 8 state model with paramSim_8state()');
%         [K,A,time,~] = paramSim_8state();
%         P = K2P(K,time);
%         tau2 = 0;
%         zoff = 0;
%         
%         clockMode = 1;
%         plotMode = 1;
%     case 5
%         zoff = 0;
%         
%         clockMode = 1;
%         plotMode = 1;
% end


%[N,~,timesteps] = size(P);
[NumStates,~,~] = size(P);
timesteps = length(time);
% Make a row vector of the discrete probability
% x = diag(A) returns a column vector of the main diagonal elements of A.

% if length(P)~=length(time)
% Peq = diag(P(:,:,end))';
% [P, ~, K, time] = K2P(K,time);
% else
[P, ~, K, time] = K2P(K,time); 
Peq = diag(P(:,:,end))';
% end

%Make a  row vector of the FRET Values
% A = linspace(.1,.9,N);

%Subtract off the mean of A
Peq=reshape(Peq,1,length(Peq));
A = reshape(A, length(A),1); 
Amean = sum(A.*Peq');
Ams = A - Amean;

%--------------------------------------------------------------------------
% Calculate 4point TCF (Sums using loops)
%--------------------------------------------------------------------------
if clockMode == 1
    K2P_timer = tic;
end
[P_tau2eq0] = K2P(K,tau2);%Fast because tau2 is only 1 element
% P_tau2eq0 = eye(NumStates);%Fast because tau2 is only 1 element
if clockMode == 1
    elapsedTime = toc(K2P_timer);
    disp(['Time to calculate calculate Pji(t=0) for an ' num2str(NumStates) ' state model using K2P: ' num2str(elapsedTime) ]);
end
C4_sim = zeros(timesteps,timesteps);

if clockMode == 1
    tstart_forloop = tic;
end

%now set up to calcualte the C3 surface
for i = 1:NumStates
    pi_part = Ams(i)*Peq(i);
    for j = 1:NumStates
        pji_part = reshape(P(j,i,:),[1 timesteps])*Ams(j)*pi_part;        
            for k = 1:NumStates
%                 C4temp = Ams(l)*reshape(P(l,k,:),[timesteps 1])*Ams(k)*P_tau2eq0(k,j)*Ams(j)*reshape(P(j,i,:),[1 timesteps])*Ams(i)*Peq(i);
                C4temp = Ams(k)*reshape(P(k,j,:),[timesteps 1])*pji_part;
%                 C4temp = Ams(l)*reshape(P(l,k,:),[1 timesteps])*pkji_part;
                C4_sim =  C4_sim + C4temp;
            end
    end
end

C4_sim_temp = C4_sim;
C4_sim = C4_sim_temp + zoff;

if addControlMode == 1
   C4_sim = C4_sim + c4Control;  
end
% C4_sim = C4_sim_temp ;

if clockMode == 1
    elapsedTime = toc(tstart_forloop);
    disp(['Time to calculate calculate C4 for an ' num2str(NumStates) ' state model: ' num2str(elapsedTime) ]);
end

% if plotMode == 1
    %%
    %extra stuff to make the plot for Andy more quickly
%     c4Exp = @ (x, A, tau) 2*((A*exp(-x/tau)).*(A*exp(-x/tau)).');
% 
%     c4Control = c4Exp(time,c4CtrlAmp,C2C4_ctrlTime)*C4_tau2eq0_exp(1,1);
%     figure();
%     set(gcf,'Color','w');
%     set(gcf,'Name','4-point Time Correlation Function');
%     
%     mesh(time,time,C4_sim);
%     shading interp;
%     title_str = ['4-point TCF: \tau_2 = ' num2str(tau2)];
%     title(title_str,'FontSize',12);
%     xlabel('\tau_1 (sec)','fontsize',12);
%     zlabel('C^{(4)}(\tau_1,\tau_2,\tau_3)','fontsize',12);
%     set(gca,'yscale','log');
%     set(gca,'xscale','log');
%     set(gca,'FontSize',14);
%     colormap('jet');
%     colorbar;
%     grid on
%     axis tight;
%     hold on
%     view([65 25]);
%     
% %define a plot 3 axis for each point in exp data to have an X and Y
% %     plot3timeGrid=repmat(time,1,length(time));
%     for i=1:length(time)
%         for j=1:length(time)
%         plot3(time(i),time(j),C4_tau2eq0_exp(i,j),'o','Color','red','MarkerSize',2.5);
%         end
%     end
%     
%     figure()
%     set(gcf,'Name','4-point Time Correlation Function');
%     
%     surf(time,time,C4_sim);
%     colormap('jet');
%     hold on
%     contour3(time,time,C4_sim,10, '-k', 'LineWidth',1.5);
%     shading interp;
%     title_str = ['4-point TCF: \tau_2 = ' num2str(tau2)];
%     title(title_str,'FontSize',12);
%     xlabel('\tau_1 (sec)','fontsize',12);
%     zlabel('C^{(4)}(\tau_1,\tau_2,\tau_3)','fontsize',12);
%     set(gca,'yscale','log');
%     set(gca,'xscale','log');
%     set(gca,'FontSize',14);
%     colorbar;
%     grid on
%     axis tight;
%     hold on
%     contour3(time,time,C4_tau2eq0_exp, 10, 'white', 'LineWidth',1.5)
%     view([0,90]);
% end

end



%THIS IS THE OLD C4 LOOP STRUCUTRE WHICH INCLUDED THE POSSIBILITY OF TAU2
% for i = 1:NumStates
%     pi_part = Ams(i)*Peq(i);
%     for j = 1:NumStates
%         pji_part = reshape(P(j,i,:),[1 timesteps])*Ams(j)*pi_part;
%         for k = 1:NumStates
%             pkji_part = Ams(k)*P_tau2eq0(k,j)*pji_part;
%             for l = 1:NumStates
% %                 C4temp = Ams(l)*reshape(P(l,k,:),[timesteps 1])*Ams(k)*P_tau2eq0(k,j)*Ams(j)*reshape(P(j,i,:),[1 timesteps])*Ams(i)*Peq(i);
%                 C4temp = Ams(l)*reshape(P(l,k,:),[timesteps 1])*pkji_part;
% %                 C4temp = Ams(l)*reshape(P(l,k,:),[1 timesteps])*pkji_part;
%                 C4_sim =  C4_sim + C4temp;
%             end
%         end
%     end
% end