% plot velocity, change in x and change in y

sim1 = 'InvestigateEvery7minVelocityXD10eps10CiLonly.csv';
X = csvread(sim1);
X = X/7;


sim1 = 'InvestigateEvery7minVelocityYD10eps10CiLonly.csv';
Y = csvread(sim1);
Y = Y/7;


sim1 = 'InvestigateEvery7minVelocityXD10eps10CiLonlyGROWINGDOMAIN.csv';
XG = csvread(sim1);
XG = XG/7;


sim1 = 'InvestigateEvery7minVelocityYD10eps10CiLonlyGROWINGDOMAIN.csv';
YG = csvread(sim1);
YG = YG/7;

figure
scatter(X,Y,'filled')
hold on
plot([0,0],[-6,6],'k','linewidth',2)
plot([-6,6],[0,0],'k','linewidth',2)

% hold on
% 
% scatter(XG,YG)

set(gca,'linewidth',2)
xlabel('\Delta X','FontSize',14)
 
ylabel('\Delta Y','FontSize',14)
set(gca,'FontSize',36)
 ax = gca;
%legend([b2,b1],'CiL only','CoA+CiL');
 box on
 
 set(gca,'linewidth',4)   

