% %% check whether it converges by comparing the ratio
% N = 50;
% 
% % figure 
% % 20 simulations
% % TimeConv = 19; % check convergence at this time, 10 max
% % %% dt = 0.010, D = 0
% % 
% % filename = sprintf('convergence/CoACiLeps200D0masscentr0p01dt.csv');
% % Mass_Centre = load(filename);
% % Mass_Centreh = Mass_Centre(1:20);
% % 
% % %% dt = 0.0050, D = 0
% % 
% % filename = sprintf('convergence/CoACiLeps200D0masscentr0p005dt.csv');
% % Mass_Centre = load(filename);
% % Mass_Centreh2 = Mass_Centre(1:20);
% % 
% % %% dt = 0.0025, D = 0
% % 
% % filename = sprintf('convergence/CoACiLeps200D0masscentr0p0025dt.csv');
% % Mass_Centre = load(filename);
% % Mass_Centreh4 = Mass_Centre(1:20);
% % 
% % 
% % Ratio = (Mass_Centreh(TimeConv) - Mass_Centreh2(TimeConv))/(Mass_Centreh2(TimeConv) - Mass_Centreh4(TimeConv))
% % 
% % 
% 
% 
% 
% 
% 
% %  xlabel('Time, hrs','FontSize',36)
% % % ylabel('positions, ',14)
% %  pldt0p1 = plot(time,Mass_Centre,'-.g','LineWidth',3)
% %  %ylim([1,6]);
% %   box on
% %  set(gca,'FontSize',36)
% %  ylabel(['Distance from the neural tube, ',char(181),'m'],'FontSize',34)
% %  ax = gca;
% %  set(gca,'linewidth',4) 
% 
% 
% %%
%  % Models with stochasticity, averages over twenty
%  %%
% 
%  %% dt = 0.02, D = 5
%  
% filename = sprintf('convergence/CoACiLeps200D5masscentr100sim0p02dt.csv');
% Mass_Centre = load(filename);
% 
% 
% MatrixforAll = zeros(20,N);
% 
% for i =1:N
%     MatrixforAll(:,i) = Mass_Centre((i-1)*20+1:i*(20));
% end
% 
% Meandt0p02 = mean(MatrixforAll');
% % 
%  
% figure
% time = [1:20];
% 
% pldt0p02 = plot(time,mean(MatrixforAll'),'-m','LineWidth',3)
% hold on
% 
% errorbar(time,mean(MatrixforAll'),std(MatrixforAll'),'-.m','LineWidth',2)
%  %ylim([1,6]);
% 
%  set(gca,'FontSize',36)
%   xlabel('Time, hrs','FontSize',36)
%  ylabel(['Distance from the neural tube, ',char(181),'m'],'FontSize',34)
%  ax = gca;
%  set(gca,'linewidth',4) 
%  %legend(pldt0p01,'\Delta t = 0.010')
%  
%  
%  %% dt = 0.01, D = 5
%  
% filename = sprintf('convergence/CoACiLeps200D5masscentr100sim0p01dt.csv');
% Mass_Centre = load(filename);
% 
% 
% MatrixforAll = zeros(20,N);
% 
% for i =1:N
%     MatrixforAll(:,i) = Mass_Centre((i-1)*20+1:i*(20));
% end
% 
% Meandt0p01 = mean(MatrixforAll');
% % 
%  
% hold on
% time = [1:20];
% 
% pldt0p01 = plot(time,mean(MatrixforAll'),'-k','LineWidth',3)
% hold on
% 
% errorbar(time,mean(MatrixforAll'),std(MatrixforAll'),'-.k','LineWidth',2)
%  %ylim([1,6]);
% 
%  set(gca,'FontSize',36)
%   xlabel('Time, hrs','FontSize',36)
%  ylabel(['Distance from the neural tube, ',char(181),'m'],'FontSize',34)
%  ax = gca;
%  set(gca,'linewidth',4) 
%  %legend(pldt0p01,'\Delta t = 0.010')
%  
%   %% dt = 0.005, D = 5
%  
% filename = sprintf('convergence/CoACiLeps200D5masscentr100sim0p005dt.csv');
% Mass_Centre = load(filename);
% 
% 
% MatrixforAll = zeros(20,N);
% 
% for i =1:N
%     MatrixforAll(:,i) = Mass_Centre((i-1)*20+1:i*(20));
% end
% 
% % 
% % 
% Meandt0p005 = mean(MatrixforAll');
% 
% hold on
% 
%     pldt0p005 = plot(time,mean(MatrixforAll'),'-.b','LineWidth',3)
% errorbar(time,mean(MatrixforAll'),std(MatrixforAll'),'-b','LineWidth',2)
%  %ylim([1,6]);
% 
%  
%   %% dt = 0.0025, D = 5
%  
% filename = sprintf('convergence/CoACiLeps200D5masscentr100sim0p0025dt.csv');
% Mass_Centre = load(filename);
% 
% 
% MatrixforAll = zeros(20,N);
% 
% for i =1:N
%     MatrixforAll(:,i) = Mass_Centre((i-1)*20+1:i*(20));
% end
% 
% % 
% % 
% Meandt0p0025 = mean(MatrixforAll');
% 
% hold on
% 
%     pldt0p0025 = plot(time,mean(MatrixforAll'),':g','LineWidth',3)
% errorbar(time,mean(MatrixforAll'),std(MatrixforAll'),'-g','LineWidth',2)
%  %ylim([1,6]);
% 
%   
% %   %% dt = 0.00125, D = 5
% %  
% % filename = sprintf('convergence/CoACiLeps200D5masscentr100sim0p00125dt.csv');
% % Mass_Centre = load(filename);
% % 
% % 
% % MatrixforAll = zeros(20,50);
% % 
% % for i =1:20
% %     MatrixforAll(:,i) = Mass_Centre((i-1)*20+1:i*(20));
% % end
% % 
% % % 
% % % 
% % Meandt0p00125 = mean(MatrixforAll');
% % 
% % hold on
% % 
% %     pldt0p00125 = plot(time,mean(MatrixforAll'),':c','LineWidth',3)
% % errorbar(time,mean(MatrixforAll'),std(MatrixforAll'),'-c','LineWidth',2)
%  
% 
% %  %% dt = 0.001, D = 5
% %  
% % filename = sprintf('convergence/CoACiLeps200D5masscentr20sim0p001dt.csv');
% % Mass_Centre = load(filename);
% % 
% % 
% % MatrixforAll = zeros(20,20);
% % 
% % for i =1:20
% %     MatrixforAll(:,i) = Mass_Centre((i-1)*20+1:i*(20));
% % end
% % 
% % % 
% % % 
% % hold on
% % 
% %     pldt0p001 = plot(time,mean(MatrixforAll'),'--r','LineWidth',3)
% % errorbar(time,mean(MatrixforAll'),std(MatrixforAll'),'-r','LineWidth',2)
% %  %ylim([1,6]);
% 
%  set(gca,'FontSize',36)
%   xlabel('Time, hrs','FontSize',36)
%  ylabel(['Distance from the neural tube, ',char(181),'m'],'FontSize',34)
%  ax = gca;
%  set(gca,'linewidth',4) 
% %legend([pldt0p02,pldt0p01,pldt0p005,pldt0p0025,pldt0p00125],'\Delta t = 0.02000','\Delta t = 0.01000','\Delta t = 0.00500','\Delta t = 0.00250','\Delta t = 0.00125') 
% legend([pldt0p02,pldt0p01,pldt0p005,pldt0p0025],'\Delta t = 0.02000','\Delta t = 0.01000','\Delta t = 0.00500','\Delta t = 0.00250') 
% 
%  %ylim([200,700])
% %  
% %  
% 
% TimeConv = 20;
% error1 = zeros(1,TimeConv);
% error2 = zeros(1,TimeConv);
% for tim = 1:TimeConv
%     error1(tim) = abs(Meandt0p01(tim) - Meandt0p005(tim))/Meandt0p01(tim);
%     error2(tim) = abs(Meandt0p005(tim) - Meandt0p0025(tim))/Meandt0p005(tim);
%     error3(tim) = abs(Meandt0p02(tim)-Meandt0p01(tim))/Meandt0p02(tim);
% end
% 
% figure 
% 
% plot(time,error1,'-k','Linewidth',2)
% hold on
% plot(time,error2,':r','Linewidth',2)
% plot(time,error3,'-.g','Linewidth',2)
%  set(gca,'FontSize',36)
%   xlabel('Time, hrs','FontSize',36)
%  ylabel(['Error'],'FontSize',34)
%  ax = gca;
%  set(gca,'linewidth',4) 
 
 

%% 1000 simulations

%% check whether it converges by comparing the ratio
N = 100;


 %% dt = 0.04, D = 5
 
filename = sprintf('convergence/CentreMasseps200D51000simdt0p02N5.csv');
Mass_Centre = load(filename);


MatrixforAll = zeros(20,N);

for i =1:N
    MatrixforAll(:,i) = Mass_Centre((i-1)*20+1:i*(20));
end

Meandt0p04 = mean(MatrixforAll');
% 
 
figure
time = [1:20];

pldt0p04 = plot(time,mean(MatrixforAll'),'-m','LineWidth',3)
hold on

errorbar(time,mean(MatrixforAll'),std(MatrixforAll'),'-.m','LineWidth',2)
 %ylim([1,6]);

 set(gca,'FontSize',36)
  xlabel('Time, hrs','FontSize',36)
 ylabel(['Distance from the neural tube, ',char(181),'m'],'FontSize',34)
 ax = gca;
 set(gca,'linewidth',4) 
 %legend(pldt0p01,'\Delta t = 0.010')
 
 
 %% dt = 0.02, D = 5
 
filename = sprintf('convergence/CentreMasseps200D51000simdt0p02N5.csv');
Mass_Centre = load(filename);


MatrixforAll = zeros(20,N);

for i =1:N
    MatrixforAll(:,i) = Mass_Centre((i-1)*20+1:i*(20));
end

Meandt0p02 = mean(MatrixforAll');
% 
 
hold on
time = [1:20];

pldt0p02 = plot(time,mean(MatrixforAll'),'-g','LineWidth',3)
hold on

errorbar(time,mean(MatrixforAll'),std(MatrixforAll'),'-.g','LineWidth',2)
 %ylim([1,6]);

 set(gca,'FontSize',36)
  xlabel('Time, hrs','FontSize',36)
 ylabel(['Distance from the neural tube, ',char(181),'m'],'FontSize',34)
 ax = gca;
 set(gca,'linewidth',4) 
 %legend(pldt0p01,'\Delta t = 0.010')
 
  
 %% dt = 0.01, D = 5
 
filename = sprintf('convergence/CentreMasseps200D51000simdt0p01N5.csv');
Mass_Centre = load(filename);


MatrixforAll = zeros(20,N);

for i =1:N
    MatrixforAll(:,i) = Mass_Centre((i-1)*20+1:i*(20));
end

Meandt0p01 = mean(MatrixforAll');
% 
 
hold on
time = [1:20];

pldt0p01 = plot(time,mean(MatrixforAll'),'-k','LineWidth',3)
hold on

errorbar(time,mean(MatrixforAll'),std(MatrixforAll'),'-.k','LineWidth',2)
 %ylim([1,6]);

 set(gca,'FontSize',36)
  xlabel('Time, hrs','FontSize',36)
 ylabel(['Distance from the neural tube, ',char(181),'m'],'FontSize',34)
 ax = gca;
 set(gca,'linewidth',4) 
 %legend(pldt0p01,'\Delta t = 0.010')
 

  %% dt = 0.005, D = 5
 
filename = sprintf('convergence/CentreMasseps200D5100simdt0p005N5.csv');
Mass_Centre = load(filename);


MatrixforAll = zeros(20,N);

for i =1:N
    MatrixforAll(:,i) = Mass_Centre((i-1)*20+1:i*(20));
end

% 
% 
Meandt0p005 = mean(MatrixforAll');

hold on

    pldt0p005 = plot(time,mean(MatrixforAll'),'-.r','LineWidth',3)
errorbar(time,mean(MatrixforAll'),std(MatrixforAll'),'-r','LineWidth',2)
 %ylim([1,6]);


  %% dt = 0.0025, D = 5
 
filename = sprintf('convergence/CentreMasseps200D51000simdt0p0025N5.csv');
Mass_Centre = load(filename);


MatrixforAll = zeros(20,N);

for i =1:N
    MatrixforAll(:,i) = Mass_Centre((i-1)*20+1:i*(20));
end

% 
% 
Meandt0p0025 = mean(MatrixforAll');

hold on

    pldt0p0025 = plot(time,mean(MatrixforAll'),':b','LineWidth',3)
errorbar(time,mean(MatrixforAll'),std(MatrixforAll'),'-b','LineWidth',2)
 %ylim([1,6]);

  %% dt = 0.00125, D = 5
 
filename = sprintf('convergence/CentreMasseps200D51000simdt0p00125N5.csv');
Mass_Centre = load(filename);


MatrixforAll = zeros(20,N);

for i =1:N
    MatrixforAll(:,i) = Mass_Centre((i-1)*20+1:i*(20));
end

% 
% 
Meandt0p00125 = mean(MatrixforAll');

hold on

    pldt0p00125 = plot(time,mean(MatrixforAll'),':c','LineWidth',3)
errorbar(time,mean(MatrixforAll'),std(MatrixforAll'),'-c','LineWidth',2)
 

%  %% dt = 0.001, D = 5
%  
% filename = sprintf('convergence/CoACiLeps200D5masscentr20sim0p001dt.csv');
% Mass_Centre = load(filename);
% 
% 
% MatrixforAll = zeros(20,20);
% 
% for i =1:20
%     MatrixforAll(:,i) = Mass_Centre((i-1)*20+1:i*(20));
% end
% 
% % 
% % 
% hold on
% 
%     pldt0p001 = plot(time,mean(MatrixforAll'),'--r','LineWidth',3)
% errorbar(time,mean(MatrixforAll'),std(MatrixforAll'),'-r','LineWidth',2)
%  %ylim([1,6]);

 set(gca,'FontSize',36)
  xlabel('Time, hrs','FontSize',36)
 ylabel(['Distance from the neural tube, ',char(181),'m'],'FontSize',34)
 ax = gca;
 set(gca,'linewidth',4) 
%legend([pldt0p02,pldt0p01,pldt0p005,pldt0p0025,pldt0p00125],'\Delta t = 0.02000','\Delta t = 0.01000','\Delta t = 0.00500','\Delta t = 0.00250','\Delta t = 0.00125') 
legend([pldt0p04,pldt0p02,pldt0p01,pldt0p005,pldt0p0025,pldt0p00125],'\Delta t = 0.04000','\Delta t = 0.02000','\Delta t = 0.01000','\Delta t = 0.00500','\Delta t = 0.00250','\Delta t = 0.00125') 

 %ylim([200,700])
%  
%  

TimeConv = 20;
error1 = zeros(1,TimeConv);
error2 = zeros(1,TimeConv);
error3 = zeros(1,TimeConv);
error4 = zeros(1,TimeConv);
error5= zeros(1,TimeConv);
for tim = 1:TimeConv
    error1(tim) = abs(Meandt0p04(tim) - Meandt0p00125(tim))/Meandt0p04(tim);
    error2(tim) = abs(Meandt0p02(tim) - Meandt0p00125(tim))/Meandt0p02(tim);
    error3(tim) = abs(Meandt0p01(tim)-Meandt0p00125(tim))/Meandt0p01(tim);
    error4(tim) = abs(Meandt0p005(tim) - Meandt0p00125(tim))/Meandt0p005(tim);
    error5(tim) = abs(Meandt0p0025(tim) - Meandt0p00125(tim))/Meandt0p0025(tim);
end

% for tim = 1:TimeConv
%     error1(tim) = (Meandt0p04(tim) - Meandt0p02(tim))/(Meandt0p02(tim) - Meandt0p01(tim));
%     error2(tim) = (Meandt0p02(tim) - Meandt0p01(tim))/(Meandt0p01(tim) - Meandt0p005(tim));
%     error3(tim) = (Meandt0p01(tim) - Meandt0p005(tim))/(Meandt0p005(tim) - Meandt0p0025(tim));
%     error4(tim) =  (Meandt0p005(tim) - Meandt0p0025(tim))/(Meandt0p0025(tim) - Meandt0p00125(tim));
% 
% end



figure 

er1 = plot(time,error1,'-m','Linewidth',2)
hold on
er2 = plot(time,error2,':g','Linewidth',2)
er3 = plot(time,error3,'-.k','Linewidth',2)
er4 = plot(time,error4,'--r','Linewidth',2)
er5 = plot(time,error5,'-b','Linewidth',2)
 set(gca,'FontSize',36)
  xlabel('Time, hrs','FontSize',36)
 ylabel(['Error'],'FontSize',34)
 ax = gca;
 set(gca,'linewidth',4) 
 legend([er1,er2,er3,er4,er5],'\Delta t = 0.04000','\Delta t = 0.02000','\Delta t = 0.01000','\Delta t = 0.00500','\Delta t = 0.00250') 

