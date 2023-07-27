function [newChi2List] = postTalapasPlotter_v2(numModels2Plot, myTableSorted, startFold, programName,genAlgMode,PatSea_version_Num,outputNumber)

model2plot_list = table2array(myTableSorted(1:numModels2Plot,1));

c2Exp = @ (x, A, tau) A*exp(-x/tau);
c4Exp = @ (x, A, tau) 2*((A*exp(-x/tau)).*(A*exp(-x/tau)).');
newChi2List=[];
histCutoff=0.4; 
for ii = 1:length(model2plot_list) 
    
    if genAlgMode == 1
    cd([startFold filesep() programName  model2plot_list{ii} filesep() 'lowestChiSquare']); 
    load([startFold filesep() programName  model2plot_list{ii} filesep() 'lowestChiSquare' filesep() 'fitInputData.mat']);
    load([startFold filesep() programName  model2plot_list{ii} filesep() 'lowestChiSquare' filesep() 'BestFitResults.mat']);
%     cd([startFold filesep() model2plot_list{ii} filesep() 'lowestChiSquare']); 
%     load([startFold filesep() model2plot_list{ii} filesep() 'lowestChiSquare' filesep() 'fitInputData.mat']);
%     load([startFold filesep() model2plot_list{ii} filesep() 'lowestChiSquare' filesep() 'BestFitResults.mat']);
    sigma_A=sigmas;
    else
    % change for patternsearch
    cd([startFold filesep()  model2plot_list{ii} filesep()  'PatSea_Nstate_polz_outputs_v' num2str(PatSea_version_Num) filesep() 'output' num2str(outputNumber)]); 
    load([startFold filesep()  model2plot_list{ii} filesep()  'PatSea_Nstate_polz_outputs_v' num2str(PatSea_version_Num) filesep() 'output' num2str(outputNumber) filesep() 'BestFitResults_PatSea_Nstate_polz.mat']);
    load([startFold filesep()  model2plot_list{ii} filesep()  'PatSea_Nstate_polz_outputs_v' num2str(PatSea_version_Num) filesep() 'output' num2str(outputNumber) filesep() 'fitInputData.mat']);
    
    end
    addControlMode = 1;
    
    wd = pwd;
    extension = 'pdf';

%     labelFontSize = 18;
%     titleFontSize = 20;
%     schemeFontSize = 16;
    
labelFontSize = 12;
titleFontSize = 14;
schemeFontSize = 12;
if genAlgMode
c2Control = c2Exp(C2_exp_x,c2ctrlAmp,c2ctrlTau)*C2_exp_y(1);
c4Control = c4Exp(C4_tau1range,c4ctrlAmp,c2ctrlTau)*C4_tau2eq0_exp(1,1);
else
c2Control = c2Exp(C2_exp_x,c2CtrlAmp,C2C4_ctrlTime)*C2_exp_y(1);
c4Control = c4Exp(C4_tau1range,c4CtrlAmp,C2C4_ctrlTime)*C4_tau2eq0_exp(1,1);    
end

saveFigMode = 1;
figure(1)
clf
talapasHistPlotter(Amp_bins,targetHistogram, P, A, sigma_A, labelFontSize,titleFontSize,saveFigMode,wd)%,sample_description)
figure(2)
clf
talapasC2plotter(C2_exp_x, C2_exp_y, P, A, yoff, addControlMode,c2Control,labelFontSize,titleFontSize,saveFigMode,wd)
figure(3)
clf
tau2 = 0;
talapasC4plotter(C4_tau1range,C4_tau3range,C4_tau2eq0_exp,P,A,K,tau2,zoff,addControlMode,c4Control,labelFontSize,titleFontSize,saveFigMode,wd)
figure(4)
clf
talapasSchemePlotter(K, model2plot_list{ii}, labelFontSize,titleFontSize,schemeFontSize,saveFigMode,wd)


saveFigMode = 0;

figure(111)
clf
subplot(2,3,1)

talapasHistPlotter(Amp_bins,targetHistogram, P, A, sigma_A, labelFontSize,titleFontSize,saveFigMode,wd)%,sample_description)

figure(111)
subplot(2,3,2)

talapasC2plotter(C2_exp_x, C2_exp_y, P, A, yoff, addControlMode,c2Control,labelFontSize,titleFontSize,saveFigMode,wd)

figure(111)
subplot(2,3,3)

talapasC4plotter(C4_tau1range,C4_tau3range,C4_tau2eq0_exp,P,A,K,tau2,zoff,addControlMode,c4Control,labelFontSize,titleFontSize,saveFigMode,wd)

figure(111)
ps = subplot(2,3,[4,6]);
talapasSchemePlotter(K, model2plot_list{ii}, labelFontSize,titleFontSize,schemeFontSize,saveFigMode,wd)
set(gca, 'position', [0.1300    0.1100    0.7750    0.3412] );

figure(111)
%         set(gcf, 'Position',  [100, 100, 2000, 1700])
set(gcf,'Position',[9          48        1230         749])

pngOutputLoc = [startFold filesep()  'pdfSummaries_' programName];
if exist(pngOutputLoc,'dir') ~= 7
    mkdir(pngOutputLoc)
end
vectorFigSaver(['pdfSummary_' model2plot_list{ii} ],pngOutputLoc, 'pdf')

%calculate the chi2 from thee 3 surfaces above using the "actual chi2"
%version
if ~genAlgMode
C2time = C2_exp_x;
C4time = C4_tau1range; 
weightHistFunc=ones(1,length(targetHistogram)); 
[chisquared,chisquared_array, chisquared_unweighted, chisquared_unweighted_array,chisquared_Weighted_array_true] = chiSqCalc_v5_ActChi2(K, A, C2time, C4time,addControlMode,c2Control,c4Control,sigma_A,yoff,zoff,histCutoff);
newChi2List=[newChi2List ; chisquared_Weighted_array_true'];
end

end
%%

end

function talapasHistPlotter(Amp_bins,targetHistogram, P, A, sigma_A, labelFontSize, titleFontSize,saveFigMode,wd)%, sample_description)
% figure(1)    
% clf
%%
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
    FRET_bins = Amp_bins;
    [Peq, denom_hist_sim, hist_sim, Peq_names] = histMaker_Nstate_v3normCorrect(P, A, sigma_A, FRET_bins);

    % Plot histogram of each state
    for i = 1:length(Peq)
%                             histPlot = plot(Amp_bins, Peq(i) * exp(-((Amp_bins - A(i))/sigma_A(i)).^2)./denom_hist_sim,...
            histPlot = plot(Amp_bins,(Peq(i)/(sqrt(2*pi)*sigma_A(i)))*exp(-((Amp_bins - A(i))/(sqrt(2)*sigma_A(i))).^2)./denom_hist_sim,...                                
            'Color',lineColor(i,:),'LineStyle',lineStyle(i,:),'LineWidth',2,'DisplayName',[Peq_names(i,:) ' = ' num2str(Peq(i))]);
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


    if saveFigMode == 1
    vectorFigSaver('hist_bestFit',wd,'pdf')
    end
    
end


function talapasC2plotter(C2_exp_x, C2_exp_y, P, A, yoff, addControlMode,c2Control,labelFontSize,titleFontSize,saveFigMode,wd)
C2time = C2_exp_x;
[C2_sim,C2time] = PA2C2(P,A,C2time,yoff,addControlMode,c2Control);


%Plot the correlation function
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
drawnow();

if saveFigMode == 1
vectorFigSaver('C2_bestFit',wd,'pdf')
end

end



function talapasC4plotter(C4_tau1range,C4_tau3range,C4_tau2eq0_exp,P,A,K,tau2,zoff,addControlMode,c4Control,labelFontSize,titleFontSize,saveFigMode,wd)
tau2 = 0;   
set(gcf,'Name','Comparing all surfaces');
set(gcf,'Color','w');

C4dataPlot = mesh(C4_tau1range,C4_tau3range,C4_tau2eq0_exp);
C4dataPlot.DisplayName = 'C^{(4)} Data';

xlabel('\tau_1 (sec)','FontSize',labelFontSize);
ylabel('\tau_3 (sec)','FontSize',labelFontSize);
ylabel('C^{(4)}','FontSize',14);
%         title('Experimental vs Simulated C4','FontSize',14);
tau2usec = 0;
% [sample_description, ~] = sample_descriptionGetter();
%             [sample_description, save_prefix] = sample_descriptionGetter();
title_str = ['Experimental vs Simulated '...
    10 'C^{(4)}(\tau_1, \tau_2 = ' num2str(tau2usec) ', \tau_3)' ...
    10 ];
title(title_str,'fontsize',titleFontSize);
set(gca, 'FontSize', labelFontSize);
xlabel('\tau_1 (sec)','fontsize',labelFontSize);
ylabel('\tau_3 (sec)','fontsize',labelFontSize);

zlabel('C^{(4)}(\tau_1, \tau_2, \tau_3)','fontsize',labelFontSize);
grid on;
%                     set(gca,'FontSize',12);
axis tight;
axis square;
colormap jet;
colorbar;
set(gca,'xscale','log');
set(gca,'yscale','log');
%         set(gca,'zscale','log');
%         zlim([1e-4,inf]);
view(55.3,33.2);
drawnow();

%                     surf(C4_tau1range,C4_tau1range,C4_tau2eq0_exp)
%                         time fix
%                         [P, V, K, C4time] = K2P(K,C4time);
tau2 = 0;
    [C4_sim,C4time] = PAK2C4(P,A,K,C4_tau1range,tau2,zoff,addControlMode,c4Control);

C4_sim = C4_sim; %Add the y-offset^2 to C4_sim

hold on;
if isreal(C4_sim) == 0
    disp(['C4 surf is complex, plotting only the real part'])
    C4_plot = surf(C4_tau1range,C4_tau1range,real(C4_sim));
else
    C4_plot = surf(C4_tau1range,C4_tau1range,C4_sim);
end
set(gca,'xscale','log');
set(gca,'yscale','log');
drawnow();

if saveFigMode == 1
vectorFigSaver('C4_bestFit',wd,'pdf')
end

end

 


function talapasSchemePlotter(K, model_name, labelFontSize,titleFontSize,schemeFontSize,saveFigMode,wd)
[G, labels] = model_schematic_plotter(K,model_name);

p = plot(G,'EdgeLabel', G.Edges.Power,'interpreter','latex');
p.NodeFontSize = 25;
p.NodeLabelColor = 'r';
p.EdgeFontSize = schemeFontSize;
p.MarkerSize = 5;
p.ArrowSize = 15;
p.LineWidth = 2; %(1/(tijs_wDB));
% titles = title(['Model: ' strrep(model_list{ii},'_',' ') ], 'FontSize',titleFontSize);
titles = title(['Model: ' strrep(model_name,'_',' ') ], 'FontSize',titleFontSize);
set(gcf,'Color','w');
set(gcf, 'PaperSize', [6 6]);
set(gca, 'FontSize', labelFontSize);
%                     box off
set(gca,'xcolor','w','ycolor','w','xtick',[],'ytick',[])
gcf
ylimits = ylim;
xlimits = xlim;
ymax = ylimits(2);
xmin = xlimits(1);
vert_spacing = ymax/length(labels); 
for p = 1:length(labels)
%     text(-1.8,2.0-(0.2*p), [strrep(labels{p},'t','t_{') '}'], 'FontSize',labelFontSize,'Color','k','interpreter','tex')
%     text(-2,1.2-(0.2*p), [strrep(strrep(labels{p},'t','t_{'),'$','') '}'], 'FontSize',labelFontSize,'Color','k','interpreter','tex')
    text(xmin*(1+0.2), ymax-1.3*vert_spacing*p, [strrep(labels{p},'$','')], 'FontSize',labelFontSize,'Color','k','interpreter','tex')
end
  
if saveFigMode == 1
vectorFigSaver('modelNetwork_bestFit',wd,'pdf')
end

end


