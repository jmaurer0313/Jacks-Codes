function [chisquared,chisquared_array, chisquared_unweighted, chisquared_unweighted_array,chisquared_Weighted_array_true] = chiSqCalc_v5_ActChi2(K, A, C2time, C4time,addControlMode,c2Control,c4Control,sigmas,yoff,zoff,histCutoff)

global normalizeMode diagnoseMode
global targetHistogram weightingFactor_Amphist weightHistFunc
global C2_exp_y weightingFactor_C2 weightC2func
global C4_tau1range C4_tau2eq0_exp weightingFactor_C4_t0 wC4func
global fitHistMode fitC2Mode fitC4Mode
% global yoff zoff
global Amp_bins
global showProgressOnFit_mode
global model_name

sample_description= model_name; 
%Initialize an array to hold the chi squared values
chisquared_array = zeros(3,1);
fitHistMode=1;
fitC2Mode=1;
fitC4Mode=1;

% sigma_A = [sigma_A1; sigma_A2; sigma_A3; sigma_A4; sigma_A5; sigma_A6; sigma_A7; sigma_A8];
%------------------------------------------------------------------
% (1) Optimization Target #1: Get a single value for the entire
% hitogram: rms_array(1) (OPTIMIZATION FUNCTION:% multigoaltcf_analytical_3state)
%------------------------------------------------------------------
P = K2P(K, C2time);
if fitHistMode == 1
    
    %set index cutoffs for the ampBins of 0.01 to 0.X where X reflects the
    %data histogram (cutoff before the edge)
    [I,v]=find(Amp_bins>0.01);
    startIdx=v(1);
    [I,v]=find(Amp_bins<histCutoff); 
    finalIdx=v(end); 

    %[Peq, denom_hist_sim, hist_sim, Peq_names] = histMaker_Nstate(P, A, sigma_A, FRET_bins);
%     [~, ~, hist_sim, ~] = histMaker_Nstate(P, A, sigma_A, FRET_bins);
    [Peq, denom_hist_sim, hist_sim, Peq_names] = histMaker_Nstate_v3normCorrect(P, A, sigmas, Amp_bins);
    if size(hist_sim) == size(targetHistogram)
        chisquared_array(1) = mean((((hist_sim(startIdx:finalIdx)-targetHistogram(startIdx:finalIdx)).^2)./(targetHistogram(startIdx:finalIdx))).*weightHistFunc(startIdx:finalIdx));
    elseif size(hist_sim') == size(targetHistogram)
        chisquared_array(1) = mean((((hist_sim(startIdx:finalIdx)'-targetHistogram(startIdx:finalIdx)).^2)./targetHistogram(startIdx:finalIdx)));
    end
end

%------------------------------------------------------------------------------------------
% (2) Optimization Target #2: Get a single value for the entire 2-point TCF: rms_array(1)
%-----------------------------------------------------------------------------------------
if fitC2Mode == 1
    
    %     [P, ~, ~, ~] = K2P(K,C2_time);     %[P, V, p, time] = K2P(K,time)
    %     [time, C2_sim, ~] = P2C(P_c2, K, time, A);
    [C2_sim,~] = PA2C2(P,A,C2time,yoff,addControlMode,c2Control);
    %     C2_sim = C2_sim + yoff;
    if normalizeMode == 1
        C2_sim = C2_sim./C2_sim(1);
    end
%     C2_exp_y=reshape(C2_exp_y,1,length(C2_exp_y)); 
    %     chisquared_array(2) = sum(((reshape(C2_sim,length(C2_sim),1) - C2_exp_y).*weightC2func).^2)*weightingFactor_C2;
    %     chisquared_array(2) = sum(((C2_sim - C2_exp_y).^2).*weightC2func)*weightingFactor_C2;
    chisquared_array(2) =mean((((C2_sim - reshape(C2_exp_y,1,length(C2_exp_y))).^2)./(C2_exp_y)).*reshape(weightC2func,1,length(weightC2func)));
    
end

%--------------------------------------------------------------------------
% Optimization Target #3: Get a single value for the entire 4-point TCF: rms_array(3)
%--------------------------------------------------------------------------
if fitC4Mode == 1
    tau2 = 0;
    P_C4 = K2P(K, C4time);
    [C4_tau2eq0_sim,~] = PAK2C4(P_C4,A,K,C4time,tau2,zoff,addControlMode,c4Control);
    
    if normalizeMode == 1
        C4_tau2eq0_sim = C4_tau2eq0_sim./C4_tau2eq0_sim(1,1);
    end
    chisquared_array(3) = mean(mean((((C4_tau2eq0_sim - C4_tau2eq0_exp).^2)./(C4_tau2eq0_exp)).*wC4func));
    
end
%% Get an unweighted chisquared value
chisquared_unweighted_array = chisquared_array;
chisquared_unweighted = sum(chisquared_unweighted_array);


%% Weight each surfaces contribution
chisquared_Weighted_array_true = zeros(size(chisquared_array));
chisquared_Weighted_array_true(1) = chisquared_array(1)*weightingFactor_Amphist;
chisquared_Weighted_array_true(2) = chisquared_array(2)*weightingFactor_C2;
chisquared_Weighted_array_true(3) = chisquared_array(3)*weightingFactor_C4_t0;

%--------------------------------------------------------------------------
% Sum all of the indivual root-mean-squares to get an overall rms
%--------------------------------------------------------------------------
chisquared = sum(chisquared_Weighted_array_true);

%% Display relative contributions
if diagnoseMode == 1
    %     disp(['     The chisquared value is ' num2str(chisquared,'%.9f') ':' num2str(chisquared_Weighted_array)]);
    %Plot the relative contribution of each Chi-squared paramater
    %     histContribution = 0;
    %     C2Contriubtion = 0;
    %     C4Contribution = 0;
    if fitHistMode == 1
        histContribution = chisquared_Weighted_array_true(1)*100/chisquared;
        %         disp(['     The contribution of the Histogram to the chisquared is ' num2str(histContribution,'%.2f') '%']);
        fprintf('Hist = %.1f   , ',histContribution);
    end
    if fitC2Mode == 1
        C2Contriubtion = chisquared_Weighted_array_true(2)*100/chisquared;
        %         disp(['     The contribution of the 2ptTCF to the chisquared is ' num2str(C2Contriubtion,'%.2f') '%']);
        fprintf('C2 = %.1f   , ',C2Contriubtion);
    end
    if fitC4Mode == 1
        C4Contribution = chisquared_Weighted_array_true(3)*100/chisquared;
        %         disp(['     The contribution of the 4-point TCF to the chisquared is ' num2str(C4Contribution,'%.2f') '%']);
        fprintf('C4 = %.1f   ',C4Contribution);
    end
    fprintf('\r\n');
end


if showProgressOnFit_mode == 1
    if fitHistMode == 1
        
        % Update figure with guesses to the fit
        figure(1);
        
        clf;
        hold on;
        % Replot data Histogram
        data_hist_Plot = plot(Amp_bins,targetHistogram);
        data_hist_Plot.LineWidth = 2;
        data_hist_Plot.Color = 'blue';
        %                             x.DisplayName = 'Data';
        
        xlabel('Amp','FontSize',14);
        ylabel('Frequency','FontSize',14);
        
%         [sample_description, ~] = sample_descriptionGetter();
        title_str = [sample_description];
        title(title_str,'fontsize',14);
        
        hold on;
        
        %Plot the fit to the data
        plot(Amp_bins,hist_sim,'r-','LineWidth',3);
        hold off
        
    end
    
    if fitC2Mode == 1
        % Update figure with guesses to the fit
        figure(2);
        clf
        plot(C2time,C2_exp_y,'b.');
        hold on;
        plot(C2time,C2_sim,'r-','LineWidth',3);
        hold off;
        xlabel('\tau (sec)','fontsize',16);
        ylabel('C^{(2)}(\tau)','fontsize',16);
        set(gca,'yscale','linear');
        set(gca,'xscale','log');
        set(gca,'FontSize',14);
        
%         [sample_description, ~] = sample_descriptionGetter();
        title_str = [sample_description];
        title(title_str,'fontsize',14);
        
        
    end
    
    if fitC4Mode == 1
        % Update figure with guesses to the fit
        figure(3);
        
        %Replot the data
        clf;
        set(gcf,'Name','Four-Point TCF: C^{(4)}');
        set(gcf,'Color','w');
        
        C4dataPlot = mesh(C4_tau1range,C4_tau1range,C4_tau2eq0_exp);
        C4dataPlot.DisplayName = 'C^{(4)} Data';
        
        xlabel('\tau_1 (sec)','FontSize',14);
        ylabel('\tau_3 (sec)','FontSize',14);
        ylabel('C^{(4)}','FontSize',14);
        %         title('Experimental vs Simulated C4','FontSize',14);
        
%         [sample_description] = sample_descriptionGetter();
        title_str = [sample_description];
        title(title_str,'fontsize',14);
        xlabel('\tau_1 (sec)','fontsize',14);
        ylabel('\tau_3 (sec)','fontsize',14);
        
        zlabel('C^{(4)}(\tau_1, \tau_2, \tau_3)','fontsize',18);
        
        %----------| Clean up the plot (non-specific --> specific) |-----------
        grid on;
        set(gca,'FontSize',12);
        axis tight;
        axis square;
        colormap jet;
        colorbar;
        set(gca,'xscale','log');
        set(gca,'yscale','log');
        %         set(gca,'zscale',' log');
        %         zlim([1e-4,inf]);
        view(55.3,33.2);
        drawnow();
        hold on;
        %Plot the fit to the data
        C4_plot = surf(C4time,C4time,real(C4_tau2eq0_sim));
        %         title_str = ['Four point time correlation function'];
        %         title(title_str,'FontSize',18);
        
        set(gca,'yscale','log');
        set(gca,'xscale','log');
        set(gca,'FontSize',14);
        
        grid on
        axis tight;
        
    end
    drawnow();
    
%     disp('Press enter to continue');
%     pause();
end

