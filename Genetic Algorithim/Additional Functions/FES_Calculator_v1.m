% wrtie a plotter/calcualtor for FES of optmized SM data sets from PatSea
%%
%first call and run a histogram of the data, then invert to the free energy
%surface
saveFigMode=0;
wd='';
labelFontSize=10;
titleFontSize=14;
talapasHistPlotter(Amp_bins,targetHistogram, P, A, sigma_A, labelFontSize,titleFontSize,saveFigMode,wd)


%%
[Peq, denom_hist_sim, hist_sim, Peq_names] = histMaker_Nstate_v3normCorrect(P, A, sigma_A, Amp_bins);
Kb=1.38e-23;
T=294.8; 
KbT=Kb*T; 
barWidth=0.015;
FESbins=linspace(-0.75,0.75,300); 

%generate exhaustive string list for N=4 rate constants, then compare it
%against the param_strings from a particular model to see if that entry
%should be null or not
allRateStrings = {'k12' 'k21' 'k13' 'k31' 'k14' 'k41' 'k23' 'k32' ...
                'k24' 'k42' 'k34' 'k43'} ;
%             For 5 States 
% allRateStrings = {'k12' 'k21' 'k13' 'k31' 'k14' 'k41' 'k23' 'k32' ...
%                 'k24' 'k42' 'k34' 'k43'} ;
rateLabels=param_strings(1:length(param_strings)-11);
rates=1./tijs; 
kmax=max(max(K)); 
transBarriers=repmat({'-'},1,length(allRateStrings));
%barriers given by the sum of state 1 E + transistion energy from 1-> M
absBarriers=repmat({'-'},1,length(allRateStrings)); 
%the order of the rate matrix when singly index is as follows:
%k11, k12,k13,k14,k21,k22,k23,k24,.... So going to make a hard coded index
%list to match to the tables outout (for N=4)
invRatesIdxList=[2,5,3,9,4,13,7,10,8,14,12,15];
fullInvRates=zeros(1,length(invRatesIdxList)); 

FEScurves=zeros(4,length(FESbins));
for i=1:length(A)
% FESCurves(i,:)=-log((Peq(i)/(sqrt(2*pi)*sigma_A(i)))*exp(-((FESbins - A(i))/sigma_A(i)).^2)./denom_hist_sim);
FESCurves(i,:)=-log((Peq(i))*exp(-((FESbins - A(i))/sigma_A(i)).^2)./denom_hist_sim);
end
FESCurves=FESCurves-min(min(FESCurves)); 

for m=1:length(invRatesIdxList)
    %make the units in milliseconds    
    fullInvRates(m)=floor(((1/(K(invRatesIdxList(m))))*1000)*1000)/1000;     
end



for j=1:length(allRateStrings)
  [bool,idx]=(find(contains(rateLabels,allRateStrings(j))>0));
  if bool>0
      %leaving out the prefactor of KbT to have the units be in terms of
      %KbT
      transBarriers(j)={floor(-log(rates(idx)/kmax)*1000)/1000};
  end
end

%calculate the params for a histogram table
% coeffs=Peq./(sqrt(2*pi)*sigma_A'); 
%due to intial error in plot, this facotr of sqrt is needed
sigmas=sigma_A./sqrt(2); 
coeffs=Peq./(sqrt(2*pi)*sigmas'); 

outputLabels=[]; 
%leaving out the prefactor of KbT to have the units be in terms of
 %KbT
% freeEnergies=-log(Peq./(sqrt(2*pi)*sigmas'));
freeEnergies=-log(Peq);
freeErescaled=(floor((freeEnergies-min(freeEnergies))*1000)./1000)';
barrierHeights=[]; 
for i=1:length(sigmas)
%     myCoef=num2str(coeffs(i));
%     myA=num2str(A(i));
%     mySig=num2str(sigmas(i));
%     myPeq= num2str(Peq(i));
    outputLabels=[outputLabels,coeffs(i),A(i),sigmas(i),Peq(i)];
end
format short
outputLabels=floor(outputLabels*1000)./1000;

%now calculate the abs barriers, acccount for possible differences in the
%heights going forward/backward between two states 
for j=1:length(allRateStrings)
  [bool,idx]=(find(contains(rateLabels,allRateStrings(j))>0));
  if bool>0
      %leaving out the prefactor of KbT to have the units be in terms of
      %KbT
      absBarriers(j)={floor(-log(rates(idx)/kmax)*1000)/1000 + freeErescaled(str2num(allRateStrings{j}(2)))};
  end
end

%draw a figure with the absolute free energies of each minima. the
%transitiuon barrers from "absBarriers" and then put directional arrows to
%show the transitions
figure(3)
plot(A,freeErescaled,'square','LineWidth',2);
hold on
%define the detialed balance condotion for each model
includeDetBal=1;
startState=2;
endState=1;
detEpureBarrier=floor(-log(K(endState,startState)/kmax)*1000)/1000 ;
detBalE={floor(-log(K(endState,startState)/kmax)*1000)/1000 + freeErescaled(startState)};
for j=1:length(allRateStrings)
  [bool,idx]=(find(contains(rateLabels,allRateStrings(j))>0));
  if bool>0
      %first calculate the half-way distance along x of the two states
      %involved
      state1=str2num(allRateStrings{j}(2));
      state2=str2num(allRateStrings{j}(3));
%       [bool2,minIdx]=find((abs(FESCurves(state1,:)-FESCurves(state2,:)))==min(abs(FESCurves(state1,:)-FESCurves(state2,:))))
      x1=A(str2num(allRateStrings{j}(2)));
      x2=A(str2num(allRateStrings{j}(3)));
      if x1>x2
          halfDistMark=((x1-x2)/2)+x2; 
          span=(x2-x1);
      else
          halfDistMark=((x2-x1)/2)+x1; 
          span=(x2-x1);
      end
      
%       plot(FESbins(minIdx), (absBarriers{j}),'*','LineWidth',1);
        plot(halfDistMark, (absBarriers{j}),'*','LineWidth',1);
%       x=[A(str2num(allRateStrings{j}(2))),A(str2num(allRateStrings{j}(3)))];
%       y=[freeErescaled(str2num(allRateStrings{j}(2))),(absBarriers{j})];
      x=[A(str2num(allRateStrings{j}(2))),A(str2num(allRateStrings{j}(2)))+span/2];
      y=[freeErescaled(str2num(allRateStrings{j}(2))),(absBarriers{j})];
      text(x(1)-0.025,y(1),['G_' allRateStrings{j}(2) '^0']);
      text(halfDistMark-0.025,y(2),['G_' allRateStrings{j}(2) '^0' '_' allRateStrings{j}(3)]);
      plot([x(1)-barWidth,x(1)+barWidth],[y(1),y(1)],[halfDistMark-barWidth,halfDistMark+barWidth],[y(2),y(2)],'LineWidth',2); 
%       annotation('textarrow',x./(max(x)),y./(max(y)),'String',allRateStrings{j},headle);
%       quiver(A(str2num(allRateStrings{j}(2))),freeErescaled(str2num(allRateStrings{j}(2))),...
%              span/2,transBarriers{j},0 );
%       absBarriers(j)={floor(-log(rates(idx)/kmax)*1000)/1000 + freeErescaled(str2num(allRateStrings{j}(2)))};
  end
end
%tack on the detailed balance condition
% [bool2,minIdx]=find((abs(FESCurves(startState,:)-FESCurves(endState,:)))==min(abs(FESCurves(startState,:)-FESCurves(endState,:))))
% x=FESbins(minIdx);
x1=A(startState);
x2=A(endState);
      if x1>x2
          halfDistMark=((x1-x2)/2)+x2; 
          span=(x2-x1);
      else
          halfDistMark=((x2-x1)/2)+x1; 
          span=(x2-x1);
      end
   if includeDetBal
y=[detBalE{1}];
plot(halfDistMark, detBalE{1},'*','LineWidth',1);
plot([halfDistMark-barWidth,halfDistMark+barWidth],[y(1),y(1)],'LineWidth',2);
text(halfDistMark-0.025,y,['G_' num2str(startState) '^0' '_' num2str(num2str(endState))]);
% set(gca,'xscale','log'); 
    end
hold on
plot(FESbins,FESCurves(1,:),FESbins,FESCurves(2,:),FESbins,FESCurves(3,:),FESbins,FESCurves(4,:),'LineWidth',1.5)
plot([-1,1],[0,0],'--','Color','black')
ylim([-1, 25])
xlim([-0.1 0.5])
xlabel('Visibility','FontSize',18);
ylabel('Free Energy (K_bT)','FontSize',18);



%%
function talapasHistPlotter(Amp_bins,targetHistogram, P, A, sigma_A, labelFontSize, titleFontSize,saveFigMode,wd)%, sample_description)
figure(1)    
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
    
    
    figure(2)    
% clf
%%
    set(gcf,'Color','w');
    set(gcf,'Name','Visb FES');
    lineColor = [[0 1 0];[0 1 1]; [1 0 1];[0 0 0];[1 0.5 0];[0.5 0 0.8]];%'c-.';'m-.';'k-.']; % Define a list of colors to loop over
    lineStyle = char(':',':', ':',':',':',':','-.','-.','-.'); % Define a list of colors to loop over

    % Replot data Histogram
%     data_hist_Plot = plot(Amp_bins,targetHistogram, 'DisplayName', 'Data');
%     data_hist_Plot.LineWidth = 2.5;
%     data_hist_Plot.Color = 'blue';
    %                             x.DisplayName = 'Data';
    xlabel('Visibility','FontSize',labelFontSize);
    ylabel('Free Energy G (KbT)','FontSize',labelFontSize);

%     [sample_description, ~] = sample_descriptionGetter();
    title_str = ['Optimized Free Energy Surface' ...
        10 ];
    title(title_str,'fontsize',titleFontSize);
    hold on;
    
    % Clear plots from previous run
    if exist('histPlot','var') == 1
        delete(histPlot)
        %                             disp('histPlot deleted')
    end
    FRET_bins = Amp_bins;
    FESbins=linspace(-0.75,0.75,300); 
%     FEScurves=[];
    [Peq, denom_hist_sim, hist_sim, Peq_names] = histMaker_Nstate_v3normCorrect(P, A, sigma_A, FRET_bins);

    % Plot histogram of each state
    for i = 1:length(Peq)
            figure(2)
            if i==3
            plot(FESbins,-log((Peq(i)/(sqrt(2*pi)*sigma_A(i)))*exp(-((FESbins + A(i))/sigma_A(i)).^2)./denom_hist_sim));
            else
            plot(FESbins,-log((Peq(i)/(sqrt(2*pi)*sigma_A(i)))*exp(-((FESbins - A(i))/sigma_A(i)).^2)./denom_hist_sim));
%                             histPlot = plot(Amp_bins, Peq(i) * exp(-((Amp_bins - A(i))/sigma_A(i)).^2)./denom_hist_sim,...
            end
            
            figure(1)
            histPlot = plot(Amp_bins,(Peq(i)/(sqrt(2*pi)*sigma_A(i)))*exp(-((Amp_bins - A(i))/sigma_A(i)).^2)./denom_hist_sim,...                                
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

    figure(2)
    xlim([-0.5 0.5]); 
    ylim([0 40]);

    if saveFigMode == 1
    vectorFigSaver('hist_bestFit',wd,'pdf')
    end
    
end