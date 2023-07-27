% Quick final figure generator for the Histogram, c2 and c3 surfaces

tau2 = 0;
addControlMode=1;
c2Exp = @ (x, A, tau) A*exp(-x/tau);
c4Exp = @ (x, A, tau) 2*((A*exp(-x/tau)).*(A*exp(-x/tau)).');
c2Control = c2Exp(C2_exp_x,c2CtrlAmp,C2C4_ctrlTime)*C2_exp_y(1);
c4Control = c4Exp(C4_tau1range,c4CtrlAmp,C2C4_ctrlTime)*C4_tau2eq0_exp(1,1);  
labelFontSize = 12;
titleFontSize = 14;
schemeFontSize = 12;

[C4_sim,C4time] = PAK2C4(P,A,K,C4_tau1range,tau2,zoff,addControlMode,c4Control);
[C2_sim,C2time] = PA2C2(P,A,C2_exp_x,yoff,addControlMode,c2Control);
FRET_bins = Amp_bins;
[Peq, denom_hist_sim, hist_sim, Peq_names] = histMaker_Nstate_v3normCorrect(P, A, sigma_A, FRET_bins);

%%
% **************Block for the histogram*******************
% ********************************************************
figure(1)
set(gcf,'Color','w');
    set(gcf,'Name','Amp Histogram');
    lineColor = [[0 1 0];[0 1 1]; [1 0 1];[0 0 0];[1 0.5 0];[0.5 0 0.8]];%'c-.';'m-.';'k-.']; % Define a list of colors to loop over
    lineStyle = char(':',':', ':',':',':',':','-.','-.','-.'); % Define a list of colors to loop over

    % Replot data Histogram
    data_hist_Plot = plot(Amp_bins,targetHistogram, 'DisplayName', 'Data');
    data_hist_Plot.LineWidth = 2.5;
    data_hist_Plot.Color = 'blue';
    %                             x.DisplayName = 'Data';
    xlabel('Visibility','FontSize',labelFontSize);
    ylabel('Frequency','FontSize',labelFontSize);

%     [sample_description, ~] = sample_descriptionGetter();
    title_str = ['Experimental vs Simulated Histograms' ...
        10 ];
    title(title_str,'fontsize',titleFontSize);
    hold on;
    % Clear plots from previous run
    if exist('histPlot','var') == 1
        delete(histPlot)
        %                             disp('histPlot deleted')
    end
%     FRET_bins = Amp_bins;
    [Peq, denom_hist_sim, hist_sim, Peq_names] = histMaker_Nstate_v3normCorrect(P, A, sigma_A, FRET_bins);

    % Plot histogram of each state
    for i = 1:length(Peq)
%                             histPlot = plot(Amp_bins, Peq(i) * exp(-((Amp_bins - A(i))/sigma_A(i)).^2)./denom_hist_sim,...
%             histPlot = plot(Amp_bins,(Peq(i)/(sqrt(2*pi)*sigma_A(i)))*exp(-((Amp_bins - A(i))/(sqrt(2)*sigma_A(i))).^2)./denom_hist_sim,...                                
%             'Color',lineColor(i,:),'LineStyle',lineStyle(i,:),'LineWidth',2,'DisplayName',[Peq_names(i,:) ' = ' num2str(Peq(i))]);
            
            histPlot = plot(Amp_bins,(Peq(i)/(sqrt(2*pi)*sigma_A(i)))*exp(-((Amp_bins - A(i))/(sqrt(2)*sigma_A(i))).^2)./denom_hist_sim,...                                
            'Color',lineColor(i,:),'LineWidth',2,'DisplayName',[Peq_names(i,:) ' = ' num2str(Peq(i))]);
        %                                 lgd = legend('show');
        % lgd.Location = 'northwest';
        hold on
    end
    set(gca, 'FontSize', labelFontSize);
    if exist('histPlotTot','var') == 1
        delete(histPlotTot)
    end
    hold on
    % Plot overall histogram
    histPlotTot = plot(Amp_bins, hist_sim,...
        '-','LineWidth',2.5,'DisplayName','Total Fit', 'Color','r');
    lgd_tot = legend('show');
    lgd_tot.FontSize = 12;
    set(lgd_tot.BoxFace, 'ColorType','truecoloralpha', 'ColorData',uint8(255*[1;1;1;.5]));  % [.5,.5,.5] is light gray; 0.8 means 20% transparent

%%    
% ************Block for the C2**********
% ********************************************************
figure(2)
plot(C2_exp_x,C2_exp_y,'b.','MarkerSize',10,'DisplayName','C^{(2)}(\tau) Data');
set(gcf,'Color','w');
hold on;
plot(C2_exp_x,C2_sim, 'r','LineWidth',2);


% [sample_description, ~] = sample_descriptionGetter();

title_str = ['Experimental vs Simulated C2' ...
    10 ];
title(title_str,'fontsize',titleFontSize);
ylim([min(C2_exp_y),max(C2_exp_y)]);
xlabel('Time (sec)','FontSize',labelFontSize);
ylabel('C^{(2)}(\tau)','FontSize',labelFontSize);
set(gca,'xscale','log');
drawnow();
set(gca, 'FontSize', labelFontSize);
grid on
set(gca,'xscale','log');
ylim([0.1e-4 0.0012]);
xlim([250e-6 2.5]); 
drawnow();



%%
%************ Block to make the plots for the C3***********
% ********************************************************

%     c4Control = c4Exp(C4time,c4CtrlAmp,C2C4_ctrlTime)*C4_tau2eq0_exp(1,1);
    figure(3);
    set(gcf,'Color','w');
    set(gcf,'Name','4-point Time Correlation Function');
    
    mesh(C4time,C4time,C4_sim);
    shading interp;
    title_str = ['4-point TCF: \tau_2 = ' num2str(tau2)];
    title(title_str,'FontSize',12);
    xlabel('\tau_1 (sec)','fontsize',12);
    zlabel('C^{(4)}(\tau_1,\tau_2,\tau_3)','fontsize',12);
    set(gca,'yscale','log');
    set(gca,'xscale','log');
    set(gca,'FontSize',14);
    colormap('jet');
    colorbar;
    grid off
    axis tight;
    hold on
    view([65 25]);
    zlim([-1e-6 2.9e-05])
    caxis([-1e-6 2.9e-05])
%define a plot 3 axis for each point in exp data to have an X and Y
%     plot3timeGrid=repmat(time,1,length(time));
    for i=1:length(C4time)
        for j=1:length(C4time)
        plot3(C4time(i),C4time(j),C4_tau2eq0_exp(i,j),'o','Color','red','MarkerSize',2.5);
        end
    end
    
    figure(4)
    set(gcf,'Name','4-point Time Correlation Function');
    
    surf(C4time,C4time,C4_sim);
    colormap('jet');
    hold on
    contour3(C4time,C4time,C4_sim,10, '-k', 'LineWidth',1.5);
    shading interp;
    title_str = ['4-point TCF: \tau_2 = ' num2str(tau2)];
    title(title_str,'FontSize',12);
    xlabel('\tau_1 (sec)','fontsize',12);
    zlabel('C^{(4)}(\tau_1,\tau_2,\tau_3)','fontsize',12);
    set(gca,'yscale','log');
    set(gca,'xscale','log');
    set(gca,'FontSize',14);
    colorbar;
    grid on
    axis tight;
    hold on
    contour3(C4time,C4time,C4_tau2eq0_exp, 10, 'white', 'LineWidth',1.5)
    view([0,90]);
    zlim([-1e-6 2.9e-05])
    caxis([-1e-6 2.9e-05])
    pbaspect([1 1 1]); 
    
%   ***************************************************************  
%     SAME SET OF PLOTS BUT SCALED COLORBAR BY THE MAXIMA OF THE DATA
    
      figure(5);
    set(gcf,'Color','w');
    set(gcf,'Name','4-point Time Correlation Function');
    
    mesh(C4time,C4time,C4_sim);
    shading interp;
    title_str = ['4-point TCF: \tau_2 = ' num2str(tau2)];
    title(title_str,'FontSize',12);
    xlabel('\tau_1 (sec)','fontsize',12);
    zlabel('C^{(4)}(\tau_1,\tau_2,\tau_3)','fontsize',12);
    set(gca,'yscale','log');
    set(gca,'xscale','log');
    set(gca,'FontSize',14);
    colormap('jet');
    colorbar;
    grid off
    axis tight;
    hold on
    view([65 25]);
    zlim([-1e-6 max([max(max(C4_tau2eq0_exp)),max(max(C4_sim))]) ])
    caxis([-1e-6 max([max(max(C4_tau2eq0_exp)),max(max(C4_sim))]) ])
%define a plot 3 axis for each point in exp data to have an X and Y
%     plot3timeGrid=repmat(time,1,length(time));
    for i=1:length(C4time)
        for j=1:length(C4time)
        plot3(C4time(i),C4time(j),C4_tau2eq0_exp(i,j),'o','Color','red','MarkerSize',2.5);
        end
    end
    
    figure(6)
    set(gcf,'Name','4-point Time Correlation Function');
    
    surf(C4time,C4time,C4_sim);
    colormap('jet');
    hold on
    contour3(C4time,C4time,C4_sim,10, '-k', 'LineWidth',1.5);
    shading interp;
    title_str = ['4-point TCF: \tau_2 = ' num2str(tau2)];
    title(title_str,'FontSize',12);
    xlabel('\tau_1 (sec)','fontsize',12);
    zlabel('C^{(4)}(\tau_1,\tau_2,\tau_3)','fontsize',12);
    set(gca,'yscale','log');
    set(gca,'xscale','log');
    set(gca,'FontSize',14);
    colorbar;
    grid on
    axis tight;
    hold on
    contour3(C4time,C4time,C4_tau2eq0_exp, 10, 'white', 'LineWidth',1.5)
    view([0,90]);
    zlim([-1e-6 max([max(max(C4_tau2eq0_exp)),max(max(C4_sim))]) ])
    caxis([-1e-6 max([max(max(C4_tau2eq0_exp)),max(max(C4_sim))]) ])
    pbaspect([1 1 1]); 