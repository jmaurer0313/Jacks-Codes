%AUTHOR: Jack Maurer 
% Goal: generic scatter plot code fr power dependence
printOpt=1; 
% LinPreAtt=dir('*LinPreAtt.mat'); 
% load('LinPreAtt'); 

fig1 = figure(1);
clf(1)
fig1pos = get(fig1, 'Position');
set(fig1, 'PaperUnits','centimeters')
set(fig1, 'PaperSize', [18, 18/fig1pos(3)*fig1pos(4)])
hold on
scatter(LinPreAtt(:,1),LinPreAtt(:,2),'Filled');
scatter(LinPreAtt(:,1),LinPreAtt(:,3),'Filled');
title('Linear Signal Pre-attenuation');
xlabel('Transmission');
ylabel('Magnitude (arb)');
xlim([0.0 1.1]); 
ylim([0 0.3]); 
grid on
legend('X Quad', 'Y Quad'); 
if printOpt
print -dpdf -r300 -bestfit LinPreAtt.pdf
end 

fig2 = figure(2);
clf(2)
fig2pos = get(fig2, 'Position');
set(fig2, 'PaperUnits','centimeters')
set(fig2, 'PaperSize', [18, 18/fig2pos(3)*fig2pos(4)])
hold on
scatter(LinPostAtt(:,1),LinPostAtt(:,2),'Filled');
scatter(LinPostAtt(:,1),LinPostAtt(:,3),'Filled');
title('Linear Signal Post-attenuation');
xlabel('Transmission');
ylabel('Magnitude (arb)');
xlim([0.0 1.1]); 
ylim([0 0.3]); 
grid on
legend('X Quad', 'Y Quad'); 
if printOpt
print -dpdf -r300 -bestfit LinPostAtt.pdf
end 


fig3 = figure(3);
clf(3)
fig3pos = get(fig3, 'Position');
set(fig3, 'PaperUnits','centimeters')
set(fig3, 'PaperSize', [18, 18/fig3pos(3)*fig3pos(4)])
hold on
scatter(NonlinPreAtt(:,1),NonlinPreAtt(:,2),'Filled');
scatter(NonlinPreAtt(:,1),NonlinPreAtt(:,3),'Filled');
title('Nonlinear Signal Pre-attenuation');
xlabel('Transmission');
ylabel('Magnitude (arb)');
xlim([0.0 1.1]); 
ylim([0 0.015]); 
grid on
legend('X Quad', 'Y Quad'); 
if printOpt
print -dpdf -r300 -bestfit NonlinPreAtt.pdf
end 


fig4 = figure(4);
clf(4)
fig4pos = get(fig4, 'Position');
set(fig4, 'PaperUnits','centimeters')
set(fig4, 'PaperSize', [18, 18/fig4pos(3)*fig4pos(4)])
hold on
scatter(NonlinPostAtt(:,1),NonlinPostAtt(:,2),'Filled');
scatter(NonlinPostAtt(:,1),NonlinPostAtt(:,3),'Filled');
title('Nonlinear Signal Post-attenuation');
xlabel('Transmission');
ylabel('Magnitude (arb)');
xlim([0.0 1.1]); 
ylim([0 0.015]); 
grid on
legend('X Quad', 'Y Quad'); 
if printOpt
print -dpdf -r300 -bestfit NonlinePostAtt.pdf
end 
