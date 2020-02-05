%% cell speed distribution

n = 1;
m= 2*n;
eps = 10;


%sim1 = 'speedD5CiLonly.csv';%speedNObiasm12eps1.txt';%'.txt'
%sim1 = 'Every1minspeedD5CoACiLgrowing.csv'
%sim1 = 'MayorEvery1minspeedD1eps10CiLonly.csv';
sim1 = 'InvestigateEvery7minspeedD1eps10CiLonly.csv';
M1 = csvread(sim1);

M1 = M1/7;

%sim1 = 'speedD5CoACiL.csv';%'speedNObiasm2eps1.txt';%'speedNObiasm12eps1LessRandom.txt'
%sim1 = 'Every3minspeedD5CoACiLgrowing.csv'
%sim1 = 'MayorEvery5minspeedD1eps10CiLonly.csv';
sim1 = 'InvestigateEvery7minspeedD1eps10CoACiL.csv';
M2 = csvread(sim1);

M2 = M2/7;

%sim1 = 'ALLcellsspeedD10CoACiLgrowing.csv';%'speedNObiasm2eps10.txt';%'speedNObiasm12eps1LessRandom.txt'
%sim1 = 'Every5minspeedD5CoACiLgrowing.csv'
%sim1 = 'MayorEvery7minspeedD1eps10CiLonly.csv'
sim1 = 'InvestigateEvery7minspeedD1eps10CiLonlyGROWINGDOMAIN.csv';
M3 = csvread(sim1);

M3 = M3/7;

%M3 = M3old(M3old<5);


%sim1 = 'speedD5CoACiL.csv';%'speedNObiasm2eps1.txt';%'speedNObiasm12eps1LessRandom.txt'
%sim1 = 'speedD10CoACiLfixed.csv'
sim1 = 'InvestigateEvery7minspeedD1eps10CoACiLGROWINGDOMAIN.csv';
M4 = csvread(sim1);

M4 = M4/7;
% 
% %sim1 = 'speedD5CoACiL.csv';%'speedNObiasm2eps1.txt';%'speedNObiasm12eps1LessRandom.txt'
% %sim1 = 'speedD10CoACiLfixed.csv'
% sim1 = 'Every9minspeedD5CoACiLgrowing.csv'
% M5 = csvread(sim1);
% 
% M5 = M5/9;


%figure
%hist(M1,150)

% group = [    ones(size(M1));
%          3 * ones(size(M2));
%          6 * ones(size(M3));
%          9 * ones(size(M4));
%          11 * ones(size(M5))];
% figure
% boxplot([M1; M2; M3;M4;M5],group)
% set(gca,'XTickLabel',{'(a)','(b)','(c)','(d)','(e)'})
%   

%% only for three

group = [    ones(size(M1));
         3 * ones(size(M2));
         6 * ones(size(M3));
         9 * ones(size(M4))];
figure
boxplot([M1; M2; M3; M4],group)
set(gca,'XTickLabel',{'(a)','(b)','(c)','(d)'})
  



set(findobj(gca,'type','line'),'linew',3)

ylabel(['Cell speed, ',char(181),'m/min'])


%
 %title(['m = ' num2str(m) ',\epsilon = ' num2str(eps)]) % only repulsion
 %xlabel(['Cell speed, ',char(181),'m/min'])
 %ylabel('Frequency')
 
  %xlabel(['Cell speed, ',char(181),'m/min'])
 %ylabel('Frequency')
 
 
 %xlim([0,6])
 ylim([0,6])
 set(gca,'FontSize',36)
 ax=gca;
 box on
set(gca,'linewidth',4) 

% figure 
% len = size(M1);
% x = [1:len(1)];
% scatter(x,M1,'b')
% % figure
% % plot(M2,'r')
% 
% 
% coefficients = polyfit(x, M1', 1)
% xFit = linspace(min(x), max(x), 1000);
% yFit = polyval(coefficients , xFit);
% hold on;
% plot(xFit, yFit, 'r-', 'LineWidth', 2);
% grid on;


figure 
len = size(M1);
%sim1 = 'ALLcellspositionsD10CoACiLgrowing.csv';%'positionsD10CoACiLfixed.csv';%
sim1 = 'InvestigateEvery7minpositionsD1eps10CiLonly.csv';
x = csvread(sim1);
%x =x (M3old<5);
scatter(x,M1,'b','filled')
% figure
% plot(M2,'r')


coefficients = polyfit(x, M1, 1)
xFit = linspace(min(x), max(x), 1000);
yFit = polyval(coefficients , xFit);
hold on;
plot(xFit, yFit, 'r-', 'LineWidth', 2);
grid on;

set(findobj(gca,'type','line'),'linew',3)

ylabel(['Cell speed, ',char(181),'m/min'],'FontSize',14)


%
 %title(['m = ' num2str(m) ',\epsilon = ' num2str(eps)]) % only repulsion
 xlabel(['Distance from the neural tube, ',char(181),'m'],'FontSize',14)

 
 
 xlim([0,500])
 ylim([0,2])
 set(gca,'FontSize',36)
 ax=gca;
 box on
set(gca,'linewidth',4) 
