function [] = formatedboxplot(x)
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here
M={'M1', 'M2', 'M3', 'M4','M5', 'M6','M7','M8','M9','M10', 'M11','M12'};

boxplot(x,'Widths',0.8,'Notch','off','Labels',M,'Whisker',1)
lines = findobj(gcf, 'type', 'line', 'Tag', 'Median');
set(lines, 'Color', 'b'); % changes the median line colour
% Change the boxplot color from blue to green
a = get(get(gca,'children'),'children');   % Get the handles of all the objects
t = get(a,'tag');   % List the names of all the objects 
set(a(1:12),'MarkerSize',6)
set(a,'LineWidth',1);
set(a,'Color','b') %set ouliers to black colour
set(a([25 27 29 31 33 35]),'Color','m') %set boxes to magenta colour
set(a([13 15 17 19 21 23]),'Color','m') %set median to magenta colour
set(a([37 39 41 43 45 47 61 63 65 67 69 71]),'Color','m') %set lower whiskers to magenta colour
set(a([49 51 53 55 57 59 73 75 77 79 81 83]),'Color','m') %set upper whiskers to magenta colour

hold on
x1=ones(size(x)).*(1+(rand(size(x))-0.5)/2);
x2=ones(size(x)).*(1+(rand(size(x))-0.5)/4);
x3=ones(size(x)).*(1+(rand(size(x))-0.5)/6);
x4=ones(size(x)).*(1+(rand(size(x))-0.5)/8);
x5=ones(size(x)).*(1+(rand(size(x))-0.5)/10);
x6=ones(size(x)).*(1+(rand(size(x))-0.5)/12);
x7=ones(size(x)).*(1+(rand(size(x))-0.5)/14);
x8=ones(size(x)).*(1+(rand(size(x))-0.5)/16);
x9=ones(size(x)).*(1+(rand(size(x))-0.5)/18);
x10=ones(size(x)).*(1+(rand(size(x))-0.5)/20);
x11=ones(size(x)).*(1+(rand(size(x))-0.5)/22);
x12=ones(size(x)).*(1+(rand(size(x))-0.5)/24);

f1=scatter(x1(:,1),x(:,1),'k','filled');f1.MarkerFaceAlpha = 0.5;hold on 
f2=scatter(x2(:,2).*2,x(:,2),'k','filled');f2.MarkerFaceAlpha = f1.MarkerFaceAlpha;
f3=scatter(x3(:,3).*3,x(:,3),'k','filled');f3.MarkerFaceAlpha = f1.MarkerFaceAlpha;
f4=scatter(x4(:,4).*4,x(:,4),'k','filled');f4.MarkerFaceAlpha = f1.MarkerFaceAlpha;
f5=scatter(x5(:,5).*5,x(:,5),'k','filled');f5.MarkerFaceAlpha = f1.MarkerFaceAlpha;
f6=scatter(x6(:,6).*6,x(:,6),'k','filled');f6.MarkerFaceAlpha = f1.MarkerFaceAlpha;
f7=scatter(x7(:,7).*7,x(:,7),'k','filled');f7.MarkerFaceAlpha = f1.MarkerFaceAlpha;
f8=scatter(x8(:,8).*8,x(:,8),'k','filled');f8.MarkerFaceAlpha = f1.MarkerFaceAlpha;
f9=scatter(x9(:,9).*9,x(:,9),'k','filled');f9.MarkerFaceAlpha = f1.MarkerFaceAlpha;
f10=scatter(x10(:,10).*10,x(:,10),'k','filled');f10.MarkerFaceAlpha = f1.MarkerFaceAlpha;
f11=scatter(x11(:,11).*11,x(:,11),'k','filled');f11.MarkerFaceAlpha = f1.MarkerFaceAlpha;
f12=scatter(x12(:,12).*12,x(:,12),'k','filled');f12.MarkerFaceAlpha = f1.MarkerFaceAlpha;

ylabel('RMS (^o)','FontWeight','bold')
xticks([1.5,3.5,5.5,7.5,9.5,11.5])
xticklabels(M)
xtickangle(45)
grid minor
%%create a second axes.
% ax1 = gca; % the first axes
% inter=([2.5,4.5,6.5,8.5,10.5,12.5]); %intervals to show a line between two boxes
% ax2 = axes('Position',ax1.Position,...
%   'XAxisLocation','bottom',...
%   'YAxisLocation','left',...
%   'Color','none',... 
%   'Ylim',ax1.YLim,...
%   'XLim',ax1.XLim,...
%   'TickLength',[0 0],...
%   'XTick', inter,  ...
%   'YTickLabel', [],  ...
%   'XTickLabel', []  );
% linkaxes([ax1 ax2],'xy')
% grid on
% set(gcf,'CurrentAxes',ax1); % return to ax1
end

