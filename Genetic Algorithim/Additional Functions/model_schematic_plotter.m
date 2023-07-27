function [G, Power] = model_schematic_plotter(K,model_name)
% function [G, Power] = model_schematic_plotter(K,Nstates,tijs, rates_str_DB, rates_idx_DB, param_strings, chisquared)

plotMode = 0;

G = digraph(K>0);

if plotMode == 1
figure(11)
plot(G)
end

edges = table2array(G.Edges);

% build tijs from K
tijs_wDB = [];
tijs_str = [];
for i = 1:length(edges)
%     tij = 1/K(edges(i,1), edges(i,2));
    tij = 1/K(edges(i,2), edges(i,1));
    tijs_wDB = [tijs_wDB; tij];
%     tijs_str = [tijs_str; 't_{' num2str(edges(i,1)) num2str(edges(i,2)) '}'];
%     tijs_str = [tijs_str; '$t_{' num2str(edges(i,1)) num2str(edges(i,2)) '}$']; 
    tijs_str = [tijs_str; 't' num2str(edges(i,1)) num2str(edges(i,2)) ];

end


tijs_label = [];
for k = 1:length(tijs_wDB)
tijs_label = [tijs_label; tijs_str(k,:) ' = ' num2str(tijs_wDB(k),'%.3e')];
end

Power = cellstr(tijs_label);
G = digraph(edges(:,1),edges(:,2), table(Power));


if plotMode == 1
% figure(11)
% clf
p = plot(G,'EdgeLabel', G.Edges.Power,'interpreter','latex');
p.NodeFontSize = 25;
p.NodeLabelColor = 'r';
p.EdgeFontSize = 15;
p.MarkerSize = 5;
p.ArrowSize = 15;
p.LineWidth = 2; %(1/(tijs_wDB));
title(['Model: ' strrep(model_name,'_',' ')], 'FontSize',20)
set(gcf,'Color','w');
set(gcf, 'PaperSize', [6 6]);
set(gca, 'FontSize', 15);
%                     box off
set(gca,'xcolor','w','ycolor','w','xtick',[],'ytick',[])

ylimits = ylim;
xlimits = xlim;
ymax = ylimits(2);
xmin = xlimits(1);
vert_spacing = ymax/length(tijs_label); 
tijs_label = string(tijs_label);
for p = 1:length(tijs_label)
%     text(-1.8,2.0-(0.2*p), [strrep(labels{p},'t','t_{') '}'], 'FontSize',labelFontSize,'Color','k','interpreter','tex')
%     text(-2,1.2-(0.2*p), [strrep(strrep(labels{p},'t','t_{'),'$','') '}'], 'FontSize',labelFontSize,'Color','k','interpreter','tex')
%     text(xmin*(1+0.2), ymax-1.3*vert_spacing*p, [strrep(tijs_label{p},'$','')], 'FontSize',labelFontSize,'Color','k','interpreter','tex')
    text(xmin*(1+0.2), ymax-1.3*vert_spacing*p, [strrep(tijs_label(p,:),'$','')], 'FontSize',15,'Color','k','interpreter','tex')
end

% for p = 1:size(tijs_label,1)
% %     text(-2.9,0.5-(0.05*p), tijs_label(p,:), 'FontSize',labelFontSize,'Color','k')%,'interpreter','latex')
%     text(-2.2,1-(0.08*p), strrep(tijs_label(p,:),'$',''), 'FontSize',labelFontSize,'Color','k','interpreter','tex')
% end
set(gcf,'Position',[500    500    1000    800])
end

end





%% Below this line is old junk

% p = plot(G,'EdgeLabel', G.Edges.Power);
% p.NodeFontSize = 20;
% p.NodeLabelColor='r';
% p.EdgeFontSize = 14;
% p.MarkerSize = 5;
% p.ArrowSize = 15;
% p.LineWidth = 2; %(1/(tijs_wDB));
% set(gcf,'Color','w');



% set(gcf, 'PaperSize', [6 6]);

% T = table(tijs_str,tijs_wDB);%,Weight,'RowNames',LastName)
% T = table(tijs_wDB,'RowNames',string(tijs_str));%,Weight,'RowNames',LastName)

% text(0,0,tijs_str


% fig = uifigure;
% subplot(1,2,2)
% uitable(fig,'Data',T)
% 
% uitable(fig,'Data',T,'ColumnName',T.Properties.VariableNames,...           
% 'RowName',T.Properties.RowNames,'Units', 'Normalized');%, 'Position',hAx.Position);% set table on top




% text(0.4,1.1,['\chi^2 = ' num2str(chisquared)],'FontSize',12);




% if length(rates_idx_DB) == 1
%     rates_idx_DB = [rates_idx_DB];
%     rates_str_DB = [rates_str_DB];
% end

% tijs_wDB = zeros(size(edges));
% tijs_str = [];

% for i = 1:length(rates_idx_DB)
% %     [row,col] = ind2sub([Nstates Nstates], rates_idx_DB(i));
% 
%     db_label = char(rates_str_DB(i));
%     db_label = [str2num(db_label(2)) str2num(db_label(3))];
% 
%     % row for detailed balance
%     db_loc = find(edges == db_label,1);
% 
%     db_val = 1/K(rates_idx_DB);
%     
%     tijs_str = strrep(param_strings(1:length(tijs)), 'k','t_{');
%     tijs_str = char(reshape(tijs_str, length(tijs_str),1));
%     tijs_str_arr = [];
%     for j = 1:length(tijs_str)
%         tijs_str_arr = [tijs_str_arr; tijs_str(j,:), '}'];
%     end
%     tijs_str = tijs_str_arr;
% %     tijs_str = strrep(tijs_str, 'k','t');
%     db_str = strrep(char(reshape(rates_str_DB(i),1,1)), 'k','t_{');
%     db_str = [db_str '}'];
%     tijs_str = [tijs_str(1:db_loc-1,:); db_str ;tijs_str(db_loc:end,:)];
%     
%     tijs_wDB = [tijs(1:db_loc-1), db_val, tijs(db_loc:end)];
% end
