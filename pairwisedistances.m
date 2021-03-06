%% plot distances

x=[1:5];

sigma = 40;

%sim1 = 'GOOD DATA February 08/pairwiseDistance1D0p5eps.csv'; % pairwiseDistance0p05BIAS1D0p5eps
sim1 = 'convergence/SimNr20Distanceseps1D0p01.csv'; % pairwiseDistance0p05BIAS1D0p5eps
M1 = csvread(sim1);

averageCoACiL(1) = mean(M1)
stdsCoACiL(1) =std(M1)


%sim1 = 'GOOD DATA February 08/pairwiseDistance1D1eps.csv';
sim1 = 'convergence/SimNr20Distanceseps1D10.csv'

M1 = csvread(sim1);

averageCoACiL(2) = mean(M1);
stdsCoACiL(2) =std(M1);


%sim1 = 'GOOD DATA February 08/pairwiseDistance1D5eps.csv';
sim1 = 'convergence/SimNr20Distanceseps200D0p01.csv';
M1 = csvread(sim1);

averageCoACiL(3) = mean(M1);
stdsCoACiL(3) =std(M1);

%sim1 = 'GOOD DATA February 08/pairwiseDistance1D10eps.csv';
sim1 = 'convergence/SimNr20Distanceseps200D10.csv';
M1 = csvread(sim1);

averageCoACiL(4) = mean(M1);
stdsCoACiL(4) =std(M1);

%sim1 = 'GOOD DATA February 08/pairwiseDistance10D1eps.csv';
sim1 = 'convergence/SimNr20Distanceseps200D0p01-3cells.csv'
M1 = csvread(sim1);

averageCoACiL(5) = mean(M1);
stdsCoACiL(5) =std(M1);

%sim1 = 'GOOD DATA February 08/pairwiseDistance0p01D1eps.csv'; % pairwiseDistance0p05BIAS1D0p5eps
% sim1 = 'Distanceseps200D0p01-3cells.csv'
% M1 = csvread(sim1);
% 
% averageCoACiL(6) = mean(M1)
% stdsCoACiL(6) =std(M1)


figure

potentialMin = 2^(1/4)*sigma;

xline = [0,6];
yline = [potentialMin,potentialMin];


b1= bar(x,averageCoACiL,0.6);   
hold on
plot(xline,yline,'linewidth',4)

er1 = errorbar(x,averageCoACiL,stdsCoACiL,'k.','linewidth',2)   
%er.Color = [0 0 0];      

 xticks([1,2,3,4,5]);%,6])
 xticklabels({'(a)','(b)','(c)','(d)','(e)'});%,'(f)'})
 xlim([0,6])

set(gca,'linewidth',2)
%xlabel(['Bias and interaction type'],'FontSize',14)
 
ylabel(['Pairwise distance, ',char(181),'m'],'FontSize',14)
set(gca,'FontSize',36)
 ax = gca;
%legend([b2,b1],'CiL only','CoA+CiL');
 box on
 
 set(gca,'linewidth',4)   