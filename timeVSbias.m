%% plot time to invade versus bias, for Xenopus model

xCoACiL = [2,5,8,11,14]; 
xCiLonly =[1,4,7,10,13];

total=1080;

%% first CoACiL

sim1 = 'TOSAVETimeCoACiL0p01betaNarrowDomain.csv';


M1 = csvread(sim1);
M1 = M1(:,1)*18/total;

averageCoACiL(1) = mean(M1)
stdsCoACiL(1) =std(M1)

%% first CiLonly

sim1 = 'TOSAVETimeCiLonly0p01betaNarrowDomain.csv';


M1 = csvread(sim1);
M1 = M1(:,1)*18/total;

averageCiLonly(1) = mean(M1)
stdsCiLonly(1) =std(M1)



%% second CoA CiL
sim1 = 'TOSAVETimeCoACiL0p02betaNarrowDomain.csv';


M1 = csvread(sim1);
M1 = M1(:,1)*18/total;

averageCoACiL(2) = mean(M1)
stdsCoACiL(2) =std(M1)

%% second CiL only
sim1 = 'TOSAVETimeCiLonly0p02betaNarrowDomain.csv';


M1 = csvread(sim1);
M1 = M1(:,1)*18/total;

averageCiLonly(2) = mean(M1)
stdsCiLonly(2) =std(M1)


%% third CoA CiL

sim1 = 'TOSAVETimeCoACiL0p03betaNarrowDomain.csv';


M1 = csvread(sim1);
M1 = M1(:,1)*18/total;

averageCoACiL(3) = mean(M1)
stdsCoACiL(3) =std(M1)


%% third CiL only

sim1 = 'TOSAVETimeCiLonly0p03betaNarrowDomain.csv';


M1 = csvread(sim1);
M1 = M1(:,1)*18/total;

averageCiLonly(3) = mean(M1)
stdsCiLonly(3) =std(M1)

%% fourth CoA CiL

sim1 = 'TOSAVETimeCoACiL0p04betaNarrowDomain.csv';


M1 = csvread(sim1);
M1 = M1(:,1)*18/total;

averageCoACiL(4) = mean(M1)
stdsCoACiL(4) =std(M1)


%% fourth CiL only

sim1 = 'TOSAVETimeCiLonly0p04betaNarrowDomain.csv';


M1 = csvread(sim1);
M1 = M1(:,1)*18/total;

averageCiLonly(4) = mean(M1)
stdsCiLonly(4) =std(M1)


%% fifth CoA CiL

sim1 = 'TOSAVETimeCoACiL0p05betaNarrowDomain.csv';


M1 = csvread(sim1);
M1 = M1(:,1)*18/total;

averageCoACiL(5) = mean(M1)
stdsCoACiL(5) =std(M1)


%% fifth CiL only

sim1 = 'TOSAVETimeCiLonly0p05betaNarrowDomain.csv';


M1 = csvread(sim1);
M1 = M1(:,1)*18/total;

averageCiLonly(5) = mean(M1)
stdsCiLonly(5) =std(M1)


figure
b1= bar(xCoACiL,averageCoACiL,0.3)   
hold on
b2 = bar(xCiLonly,averageCiLonly,0.3) 


er1 = errorbar(xCoACiL,averageCoACiL,stdsCoACiL,'k.','linewidth',2)   
er2 = errorbar(xCiLonly, averageCiLonly, stdsCiLonly,'k.','linewidth',2)
%er.Color = [0 0 0];      

xticks([1.5,4.5,7.5,10.5,13.5])
xticklabels({'0.01','0.02','0.03','0.04','0.05'})
xlim([0.3,15.3])

set(gca,'linewidth',2)
xlabel(['Bias, \beta'],'FontSize',14)
ylabel(['Time to invasion, hrs'],'FontSize',14)
set(gca,'FontSize',36)
 ax = gca;
legend([b2,b1],'CiL only','CoA+CiL');
 box on
 grid on
 ylim([0,40])
 set(gca,'linewidth',4)   

