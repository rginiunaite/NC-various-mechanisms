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
 
 

% %% 1000 simulations
% 
% %% check whether it converges by comparing the ratio
% N = 100;
% 
% 
%  %% dt = 0.04, D = 5
%  
% filename = sprintf('convergence/CentreMasseps200D51000simdt0p02N5.csv');
% Mass_Centre = load(filename);
% 
% 
% MatrixforAll = zeros(20,N);
% 
% for i =1:N
%     MatrixforAll(:,i) = Mass_Centre((i-1)*20+1:i*(20));
% end
% 
% Meandt0p04 = mean(MatrixforAll');
% % 
%  
% figure
% time = [1:20];
% 
% pldt0p04 = plot(time,mean(MatrixforAll'),'-m','LineWidth',3)
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
%  %% dt = 0.02, D = 5
%  
% filename = sprintf('convergence/CentreMasseps200D51000simdt0p02N5.csv');
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
% hold on
% time = [1:20];
% 
% pldt0p02 = plot(time,mean(MatrixforAll'),'-g','LineWidth',3)
% hold on
% 
% errorbar(time,mean(MatrixforAll'),std(MatrixforAll'),'-.g','LineWidth',2)
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
% filename = sprintf('convergence/CentreMasseps200D51000simdt0p01N5.csv');
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
% 
%   %% dt = 0.005, D = 5
%  
% filename = sprintf('convergence/CentreMasseps200D5100simdt0p005N5.csv');
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
%     pldt0p005 = plot(time,mean(MatrixforAll'),'-.r','LineWidth',3)
% errorbar(time,mean(MatrixforAll'),std(MatrixforAll'),'-r','LineWidth',2)
%  %ylim([1,6]);
% 
% 
%   %% dt = 0.0025, D = 5
%  
% filename = sprintf('convergence/CentreMasseps200D51000simdt0p0025N5.csv');
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
%     pldt0p0025 = plot(time,mean(MatrixforAll'),':b','LineWidth',3)
% errorbar(time,mean(MatrixforAll'),std(MatrixforAll'),'-b','LineWidth',2)
%  %ylim([1,6]);
% 
%   %% dt = 0.00125, D = 5
%  
% filename = sprintf('convergence/CentreMasseps200D51000simdt0p00125N5.csv');
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
% Meandt0p00125 = mean(MatrixforAll');
% 
% hold on
% 
%     pldt0p00125 = plot(time,mean(MatrixforAll'),':c','LineWidth',3)
% errorbar(time,mean(MatrixforAll'),std(MatrixforAll'),'-c','LineWidth',2)
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
% legend([pldt0p04,pldt0p02,pldt0p01,pldt0p005,pldt0p0025,pldt0p00125],'\Delta t = 0.04000','\Delta t = 0.02000','\Delta t = 0.01000','\Delta t = 0.00500','\Delta t = 0.00250','\Delta t = 0.00125') 
% 
%  %ylim([200,700])
% %  
% %  
% 
% TimeConv = 20;
% error1 = zeros(1,TimeConv);
% error2 = zeros(1,TimeConv);
% error3 = zeros(1,TimeConv);
% error4 = zeros(1,TimeConv);
% error5= zeros(1,TimeConv);
% for tim = 1:TimeConv
%     error1(tim) = abs(Meandt0p04(tim) - Meandt0p00125(tim))/Meandt0p04(tim);
%     error2(tim) = abs(Meandt0p02(tim) - Meandt0p00125(tim))/Meandt0p02(tim);
%     error3(tim) = abs(Meandt0p01(tim)-Meandt0p00125(tim))/Meandt0p01(tim);
%     error4(tim) = abs(Meandt0p005(tim) - Meandt0p00125(tim))/Meandt0p005(tim);
%     error5(tim) = abs(Meandt0p0025(tim) - Meandt0p00125(tim))/Meandt0p0025(tim);
% end
% 
% % for tim = 1:TimeConv
% %     error1(tim) = (Meandt0p04(tim) - Meandt0p02(tim))/(Meandt0p02(tim) - Meandt0p01(tim));
% %     error2(tim) = (Meandt0p02(tim) - Meandt0p01(tim))/(Meandt0p01(tim) - Meandt0p005(tim));
% %     error3(tim) = (Meandt0p01(tim) - Meandt0p005(tim))/(Meandt0p005(tim) - Meandt0p0025(tim));
% %     error4(tim) =  (Meandt0p005(tim) - Meandt0p0025(tim))/(Meandt0p0025(tim) - Meandt0p00125(tim));
% % 
% % end
% 
% 
% 
% figure 
% 
% er1 = plot(time,error1,'-m','Linewidth',2)
% hold on
% er2 = plot(time,error2,':g','Linewidth',2)
% er3 = plot(time,error3,'-.k','Linewidth',2)
% er4 = plot(time,error4,'--r','Linewidth',2)
% er5 = plot(time,error5,'-b','Linewidth',2)
%  set(gca,'FontSize',36)
%   xlabel('Time, hrs','FontSize',36)
%  ylabel(['Error'],'FontSize',34)
%  ax = gca;
%  set(gca,'linewidth',4) 
%  legend([er1,er2,er3,er4,er5],'\Delta t = 0.04000','\Delta t = 0.02000','\Delta t = 0.01000','\Delta t = 0.00500','\Delta t = 0.00250') 


 
%  %% 2000 simulations, ordered timesteps
% 
% %% check whether it converges by comparing the ratio
% N = 2000;
% end_time = 11;
% 
%  
%  %% dt = 0.02, D = 5
%  
% %filename = sprintf('convergence/CoACiLeps200D5masscentr500simdt0p02ORDER2.csv');
% %filename = sprintf('convergence/CentreMasseps200D51000simdt0p02.csv');
% %filename = sprintf('convergence/CentreMassD5eps200sim2000dt0p02.csv');
% %filename = sprintf('convergence/CentreMassD5eps200sim2000dt0p02cells25.csv');
% %filename = sprintf('convergence/CentreMassD5eps200sim2000dt0p02cells10.csv');
% %filename = sprintf('convergence/CentreMassD5eps200sim500dt0p02cells25.csv');
% filename = sprintf('convergence/CentreMassD1eps200sim2000dt0p02.csv');
% 
% Mass_Centre = load(filename);
% 
% 
% MatrixforAll0p02 = zeros(end_time,N);
% 
% for i =1:N
%     MatrixforAll0p02(:,i) = Mass_Centre((i-1)*end_time+1:i*(end_time));
% end
% 
% Meandt0p02 = mean(MatrixforAll0p02');
% Stddt0p02 = std(MatrixforAll0p02');
% % 
%  
% hold on
% time = [0:end_time-1];
% 
% pldt0p02 = plot(time,Meandt0p02,'-g','LineWidth',3)
% hold on
% 
% errorbar(time,Meandt0p02,Stddt0p02,'-.g','LineWidth',2)
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
% %filename = sprintf('convergence/CoACiLeps200D5masscentr500simdt0p01ORDER2.csv');
% %filename = sprintf('convergence/CentreMasseps200D51000simdt0p01.csv');
% %filename = sprintf('convergence/CentreMassD5eps200sim2000dt0p01.csv');
% %filename = sprintf('convergence/CentreMassD5eps200sim2000dt0p01cells25.csv');
% %filename = sprintf('convergence/CentreMassD5eps200sim2000dt0p01cells10.csv');
% %filename = sprintf('convergence/CentreMassD5eps200sim500dt0p01cells25.csv');
% filename = sprintf('convergence/CentreMassD1eps200sim2000dt0p01.csv');
% 
% 
% 
% Mass_Centre = load(filename);
% 
% 
% MatrixforAll = zeros(end_time,N);
% 
% for i =1:N
%     MatrixforAll(:,i) =Mass_Centre((i-1)*end_time+1:i*(end_time));
% end
% 
% Meandt0p01 = mean(MatrixforAll');
% Stddt0p01 = std(MatrixforAll');
% % 
%  
% hold on
% 
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
% 
%   %% dt = 0.005, D = 5
%  
% %filename = sprintf('convergence/CoACiLeps200D5masscentr500simdt0p005ORDERV2.csv');
% %filename = sprintf('convergence/CentreMasseps200D51000simdt0p005.csv');
% %filename = sprintf('convergence/CentreMassD5eps200sim2000dt0p005.csv');
% %filename = sprintf('convergence/CentreMassD5eps200sim2000dt0p005cells25.csv');
% %filename = sprintf('convergence/CentreMassD5eps200sim2000dt0p005cells10.csv');
% %filename = sprintf('convergence/CentreMassD5eps200sim500dt0p005cells25.csv');
% filename = sprintf('convergence/CentreMassD1eps200sim2000dt0p005.csv');
% 
% 
% Mass_Centre = load(filename);
% 
% 
% MatrixforAll = zeros(end_time,N);
% 
% for i =1:N
%     MatrixforAll(:,i) = Mass_Centre((i-1)*end_time+1:i*(end_time));
% end
% 
% % 
% % 
% Meandt0p005 = mean(MatrixforAll');
% Stddt0p005 = std(MatrixforAll');
% 
% hold on
% 
%     pldt0p005 = plot(time,mean(MatrixforAll'),'-.r','LineWidth',3)
% errorbar(time,mean(MatrixforAll'),std(MatrixforAll'),'-r','LineWidth',2)
%  %ylim([1,6]);
% 
% % 
%   %% dt = 0.0025, D = 5
%  
% %filename = sprintf('convergence/CoACiLeps200D5masscentr500simdt0p0025ORDERV2.csv');
% %filename = sprintf('convergence/CentreMasseps200D51000simdt0p0025.csv');
% %filename = sprintf('convergence/CentreMassD5eps200sim2000dt0p0025.csv');
% %filename = sprintf('convergence/CentreMassD5eps200sim2000dt0p0025cells25.csv');
% %filename = sprintf('convergence/CentreMassD5eps200sim2000dt0p0025cells10.csv');
% %filename = sprintf('convergence/CentreMassD5eps200sim500dt0p0025cells25.csv');
% filename = sprintf('convergence/CentreMassD1eps200sim2000dt0p0025.csv');
% 
% Mass_Centre = load(filename);
% 
% 
% MatrixforAll = zeros(end_time,N);
% 
% for i =1:N
%     MatrixforAll(:,i) = Mass_Centre((i-1)*end_time+1:i*(end_time));
% end
% 
% % 
% % 
% Meandt0p0025 = mean(MatrixforAll');
% Stddt0p0025 = std(MatrixforAll');
% 
% hold on
% 
%     pldt0p0025 = plot(time,mean(MatrixforAll'),':b','LineWidth',3)
% errorbar(time,mean(MatrixforAll'),std(MatrixforAll'),'-b','LineWidth',2)
%  %ylim([1,6]);
% 
%   %% dt = 0.00125, D = 5
%  
% %filename = sprintf('convergence/CoACiLeps200D5masscentr500simdt0p00125ORDER.csv');
% %filename = sprintf('convergence/CentreMasseps200D51000simdt0p00125.csv');
% %filename = sprintf('convergence/CentreMassD5eps200sim2000dt0p00125.csv');
% %filename = sprintf('convergence/CentreMassD5eps200sim2000dt0p00125cells25.csv');
% %filename = sprintf('convergence/CentreMassD5eps200sim2000dt0p00125cells10.csv');
% %filename = sprintf('convergence/CentreMassD5eps200sim500dt0p00125cells25.csv');
% filename = sprintf('convergence/CentreMassD1eps200sim2000dt0p00125.csv');
% 
% Mass_Centre = load(filename);
% 
% 
% MatrixforAll0p00125 = zeros(end_time,N);
% 
% for i =1:N
%     MatrixforAll0p00125(:,i) = Mass_Centre((i-1)*end_time+1:i*(end_time));
% end
% 
% % 
% % 
% Meandt0p00125 = mean(MatrixforAll0p00125');
% Stddt0p00125 = std(MatrixforAll0p00125');
% 
% hold on
% 
%     pldt0p00125 = plot(time,Meandt0p00125,':c','LineWidth',3)
% errorbar(time,Meandt0p00125, Stddt0p00125,'-c','LineWidth',2)
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
%  xticks([0,2,4,6,8, 10]);%,6])
%  xticklabels({'0','1','2','3','4','5'}); % this is because I saved the data every 30min for five hours
%  ax = gca;
%  set(gca,'linewidth',4) 
% %legend([pldt0p02,pldt0p01,pldt0p005,pldt0p0025,pldt0p00125],'\Delta t = 0.02000','\Delta t = 0.01000','\Delta t = 0.00500','\Delta t = 0.00250','\Delta t = 0.00125') 
% legend([pldt0p02,pldt0p01,pldt0p005,pldt0p0025,pldt0p00125],'\Delta t = 0.02000','\Delta t = 0.01000','\Delta t = 0.00500','\Delta t = 0.00250','\Delta t = 0.00125') 
% 
%  %ylim([200,700])
% %  
% %  
% 
% TimeConv = end_time;
% stderror2 = zeros(1,TimeConv);
% stderror3 = zeros(1,TimeConv);
% stderror4 = zeros(1,TimeConv);
% stderror5 = zeros(1,TimeConv);
% error2 = zeros(1,TimeConv);
% error3 = zeros(1,TimeConv);
% error4 = zeros(1,TimeConv);
% error5= zeros(1,TimeConv);
% for tim = 1:TimeConv
%     error2(tim) = abs(Meandt0p02(tim) - Meandt0p00125(tim))/Meandt0p02(tim);
%     stderror2(tim) = sqrt((Stddt0p02(tim)^2 + Stddt0p00125(tim)^2 )/Meandt0p02(tim)^2 -  + (Meandt0p02(tim)-Meandt0p00125(tim))^2/Meandt0p02(tim)^4 * Stddt0p02(tim)^2   );    
%     error3(tim) = abs(Meandt0p01(tim)-Meandt0p00125(tim))/Meandt0p01(tim);
%     stderror3(tim) = sqrt((Stddt0p01(tim)^2 + Stddt0p00125(tim)^2 )/Meandt0p01(tim)^2 -  + (Meandt0p01(tim)-Meandt0p00125(tim))^2/Meandt0p01(tim)^4 * Stddt0p01(tim)^2   );    
% 
%     error4(tim) = abs(Meandt0p005(tim) - Meandt0p00125(tim))/Meandt0p005(tim);
%     stderror4(tim) = sqrt((Stddt0p005(tim)^2 + Stddt0p00125(tim)^2 )/Meandt0p005(tim)^2 -  + (Meandt0p005(tim)-Meandt0p00125(tim))^2/Meandt0p005(tim)^4 * Stddt0p005(tim)^2   );    
% 
%     error5(tim) = abs(Meandt0p0025(tim) - Meandt0p00125(tim))/Meandt0p0025(tim);
%     stderror5(tim) = sqrt((Stddt0p0025(tim)^2 + Stddt0p00125(tim)^2 )/Meandt0p0025(tim)^2 -  + (Meandt0p0025(tim)-Meandt0p00125(tim))^2/Meandt0p0025(tim)^4 * Stddt0p0025(tim)^2   );    
% 
% end
% 
% % for tim = 1:TimeConv
% %     error1(tim) = (Meandt0p04(tim) - Meandt0p02(tim))/(Meandt0p02(tim) - Meandt0p01(tim));
% %     error2(tim) = (Meandt0p02(tim) - Meandt0p01(tim))/(Meandt0p01(tim) - Meandt0p005(tim));
% %     error3(tim) = (Meandt0p01(tim) - Meandt0p005(tim))/(Meandt0p005(tim) - Meandt0p0025(tim));
% %     error4(tim) =  (Meandt0p005(tim) - Meandt0p0025(tim))/(Meandt0p0025(tim) - Meandt0p00125(tim));
% % 
% % end
% 
% 
% 
% figure 
% 
% %er1 = plot(time,error1,'-m','Linewidth',2)
% hold on
% er2 = plot(time,error2,':g','Linewidth',2)
% %errorbar(time,error2, stderror2/sqrt(N),':g','LineWidth',2)
% 
% 
% er3 = plot(time,error3,'-.k','Linewidth',2)
% %errorbar(time,error3,stderror3/sqrt(N),'-.k','Linewidth',2)
% er4 = plot(time,error4,'--r','Linewidth',2)
% %errorbar(time,error4,stderror4/sqrt(N),'--r','Linewidth',2)
% 
% er5 = plot(time,error5,'-b','Linewidth',2)
% %errorbar(time,error5,stderror5/sqrt(N),'-b','Linewidth',2)
% 
%  set(gca,'FontSize',36)
%   xlabel('Time, hrs','FontSize',36)
%  ylabel(['Error'],'FontSize',34)
%   xticks([0,2,4,6,8, 10]);%,6])
%   xticklabels({'0','1','2','3','4','5'}); % this is because I saved the data every 30min for five hours
%  yticks([0,0.002,0.004,0.006,0.008, 0.010]);%,6])
%  yticklabels({'0','0.002','0.004','0.006','0.008','0.010'}); % no standard
%  %error
% %   yticks([-0.01,0.00,0.01, 0.02]);%,6])
% %  yticklabels({'-0.01','0.00','0.01','0.02'}); % with standard error
% %  ylim([-0.015,0.025])
%  ax = gca;
%  set(gca,'linewidth',4) 
%  legend([er2,er3,er4,er5],'\Delta t = 0.02000','\Delta t = 0.01000','\Delta t = 0.00500','\Delta t = 0.00250') 
% 
%  

%  %% 10000 simulations, not ordered timesteps
% 
% %% check whether it converges by comparing the ratio
% N = 10000;
% end_time = 11;
% 
%  figure
%  % dt = 0.02, D = 5
%  
% %filename = sprintf('convergence/CoACiLeps200D5masscentr500simdt0p02ORDER2.csv');
% %filename = sprintf('convergence/CentreMasseps200D51000simdt0p02.csv');
% %filename = sprintf('convergence/CentreMassD5eps200sim2000dt0p02.csv');
% %filename = sprintf('convergence/CentreMassD5eps200sim2000dt0p02cells25.csv');
% %filename = sprintf('convergence/CentreMassD5eps200sim2000dt0p02cells10.csv');
% %filename = sprintf('convergence/CentreMassD5eps200sim500dt0p02cells25.csv');
% %filename = sprintf('convergence/CentreMassD5eps200sim10000dt0p02.csv');
% %filename = sprintf('convergence/CentremassD5eps200dt0p02sim30000.csv');
% %filename = sprintf('RepOnlyCentreMassD5eps200dt0p02sim20000.csv');
% %filename = sprintf('convergence/CentreMassD1eps150dt0p02sim20000.csv');
% %filename = sprintf('ORDERCentreMassD5eps200dt0p02sim10.csv');
% %filename = sprintf('CentreMasseps200D5dt0p02sim20000.csv');
% %filename = sprintf('convergence/CentreMassD5eps200dt0p02BIAS1000sim.csv');
% filename = sprintf('convergence/CentreMassD5eps100dt0p02sim20000.csv');
% 
% 
% 
% Mass_Centre = load(filename);
% 
% for i =1:length(Mass_Centre)
%     if Mass_Centre(i) == 200
%         if rem(i-1,11) ~= 0
%            vector(i) =  i;
%         end
%     end
% end
% 
% 
% MatrixforAll0p02 = zeros(end_time,N);
% 
% for i =1:N
%     MatrixforAll0p02(:,i) = Mass_Centre((i-1)*end_time+1:i*(end_time));
% end
% 
% Meandt0p02 = mean(MatrixforAll0p02');
% Stddt0p02 = std(MatrixforAll0p02');
% % 
%  
% hold on
% time = [0:end_time-1];
% 
% pldt0p02 = plot(time,Meandt0p02,'-g','LineWidth',3)
% hold on
% 
% errorbar(time,Meandt0p02,Stddt0p02,'-.g','LineWidth',2)
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
% %filename = sprintf('convergence/CoACiLeps200D5masscentr500simdt0p01ORDER2.csv');
% %filename = sprintf('convergence/CentreMasseps200D51000simdt0p01.csv');
% %filename = sprintf('convergence/CentreMassD5eps200sim2000dt0p01.csv');
% %filename = sprintf('convergence/CentreMassD5eps200sim2000dt0p01cells25.csv');
% %filename = sprintf('convergence/CentreMassD5eps200sim2000dt0p01cells10.csv');
% %filename = sprintf('convergence/CentreMassD5eps200sim500dt0p01cells25.csv');
% %filename = sprintf('convergence/CentreMassD5eps200sim10000dt0p01.csv');
% %filename = sprintf('convergence/CentremassD5eps200dt0p01sim30000.csv');
% %filename = sprintf('RepOnlyCentreMassD5eps200dt0p01sim20000.csv');
% %filename = sprintf('convergence/CentreMassD1eps150dt0p01sim20000.csv');
% %filename = sprintf('ORDERCentreMassD5eps200dt0p01sim10.csv');
% %filename = sprintf('CentreMasseps200D5dt0p01sim20000.csv');
% %filename = sprintf('convergence/CentreMassD5eps200dt0p01BIAS1000sim.csv');
% filename = sprintf('convergence/CentreMassD5eps100dt0p01sim20000.csv');
% 
% 
% 
% Mass_Centre = load(filename);
% for i =1:length(Mass_Centre)
%     if Mass_Centre(i) == 200
%         if rem(i-1,11) ~= 0
%            vector0p01(i) =  i;
%         end
%     end
% end
% 
% 
% MatrixforAll = zeros(end_time,N);
% 
% for i =1:N
%     MatrixforAll(:,i) =Mass_Centre((i-1)*end_time+1:i*(end_time));
% end
% 
% Meandt0p01 = mean(MatrixforAll');
% Stddt0p01 = std(MatrixforAll');
% % 
%  
% hold on
% 
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
% 
%   %% dt = 0.005, D = 5
%  
% %filename = sprintf('convergence/CoACiLeps200D5masscentr500simdt0p005ORDERV2.csv');
% %filename = sprintf('convergence/CentreMasseps200D51000simdt0p005.csv');
% %filename = sprintf('convergence/CentreMassD5eps200sim2000dt0p005.csv');
% %filename = sprintf('convergence/CentreMassD5eps200sim2000dt0p005cells25.csv');
% %filename = sprintf('convergence/CentreMassD5eps200sim2000dt0p005cells10.csv');
% %filename = sprintf('convergence/CentreMassD5eps200sim500dt0p005cells25.csv');
% %filename = sprintf('convergence/CentreMassD5eps200sim10000dt0p005.csv');
% %filename = sprintf('convergence/CentremassD5eps200dt0p005sim30000.csv');
% %filename = sprintf('RepOnlyCentreMassD5eps200dt0p005sim20000.csv');
% %filename = sprintf('convergence/CentreMassD1eps150dt0p005sim20000.csv');
% %filename = sprintf('ORDERCentreMassD5eps200dt0p005sim10.csv');
% %filename = sprintf('CentreMasseps200D5dt0p005sim20000.csv');
% %filename = sprintf('convergence/CentreMassD5eps200dt0p005BIAS1000sim.csv');
% filename = sprintf('convergence/CentreMassD5eps100dt0p005sim20000.csv');
% 
% 
% Mass_Centre = load(filename);
% 
% 
% MatrixforAll = zeros(end_time,N);
% 
% for i =1:N
%     MatrixforAll(:,i) = Mass_Centre((i-1)*end_time+1:i*(end_time));
% end
% 
% % 
% % 
% Meandt0p005 = mean(MatrixforAll');
% Stddt0p005 = std(MatrixforAll');
% 
% hold on
% 
%     pldt0p005 = plot(time,mean(MatrixforAll'),'-.r','LineWidth',3)
% errorbar(time,mean(MatrixforAll'),std(MatrixforAll'),'-r','LineWidth',2)
%  %ylim([1,6]);
% 
% % 
%   %% dt = 0.0025, D = 5
% 
% %filename = sprintf('convergence/CoACiLeps200D5masscentr500simdt0p0025ORDERV2.csv');
% %filename = sprintf('convergence/CentreMasseps200D51000simdt0p0025.csv');
% %filename = sprintf('convergence/CentreMassD5eps200sim2000dt0p0025.csv');
% %filename = sprintf('convergence/CentreMassD5eps200sim2000dt0p0025cells25.csv');
% %filename = sprintf('convergence/CentreMassD5eps200sim2000dt0p0025cells10.csv');
% %filename = sprintf('convergence/CentreMassD5eps200sim500dt0p0025cells25.csv');
% %filename = sprintf('convergence/CentreMassD5eps200sim10000dt0p0025.csv');
% %filename = sprintf('convergence/CentremassD5eps200dt0p0025sim30000.csv');
% %filename = sprintf('RepOnlyCentreMassD5eps200dt0p0025sim20000.csv');
% %filename = sprintf('convergence/CentreMassD1eps150dt0p0025sim20000.csv');
% %filename = sprintf('ORDERCentreMassD5eps200dt0p0025sim10.csv');
% %filename = sprintf('CentreMasseps200D5dt0p0025sim20000.csv');
% %filename = sprintf('convergence/CentreMassD5eps200dt0p0025BIAS1000sim.csv');
% filename = sprintf('convergence/CentreMassD5eps100dt0p0025sim20000.csv');
% 
% 
% Mass_Centre = load(filename);
% 
% 
% MatrixforAll = zeros(end_time,N);
% 
% for i =1:N
%     MatrixforAll(:,i) = Mass_Centre((i-1)*end_time+1:i*(end_time));
% end
% 
% % 
% % 
% Meandt0p0025 = mean(MatrixforAll');
% Stddt0p0025 = std(MatrixforAll');
% 
% hold on
% 
%     pldt0p0025 = plot(time,mean(MatrixforAll'),':b','LineWidth',3)
% errorbar(time,mean(MatrixforAll'),std(MatrixforAll'),'-b','LineWidth',2)
%  %ylim([1,6]);
% 
% 
%   %% dt = 0.00125, D = 5
% 
% %filename = sprintf('convergence/CoACiLeps200D5masscentr500simdt0p00125ORDER.csv');
% %filename = sprintf('convergence/CentreMasseps200D51000simdt0p00125.csv');
% %filename = sprintf('convergence/CentreMassD5eps200sim2000dt0p00125.csv');
% %filename = sprintf('convergence/CentreMassD5eps200sim2000dt0p00125cells25.csv');
% %filename = sprintf('convergence/CentreMassD5eps200sim2000dt0p00125cells10.csv');
% %filename = sprintf('convergence/CentreMassD5eps200sim500dt0p00125cells25.csv');
% %filename = sprintf('convergence/CentreMassD1eps200sim2000dt0p00125.csv');
% %filename = sprintf('convergence/CentreMassD5eps200sim10000dt0p00125.csv');
% %filename = sprintf('convergence/CentremassD5eps200dt0p00125sim30000.csv');
% %filename = sprintf('convergence/trial0p00125sim24466.csv');
% %filename = sprintf('RepOnlyCentreMassD5eps200dt0p00125sim20000.csv');
% %filename = sprintf('convergence/CentreMassD1eps150dt0p00125sim20000.csv');
% %filename = sprintf('ORDERCentreMassD5eps200dt0p00125sim10.csv');
% %filename = sprintf('CentreMasseps200D5dt0p00125sim20000.csv');
% %filename = sprintf('convergence/CentreMassD5eps200dt0p00125BIAS1000sim.csv');
% filename = sprintf('convergence/CentreMassD5eps100dt0p00125sim20000.csv');
% 
% 
% 
% Mass_Centre = load(filename);
% 
% MatrixforAll0p00125 = zeros(end_time,N);
% 
% for i =1:N
% 
%     MatrixforAll0p00125(:,i) = Mass_Centre((i-1)*end_time+1:i*(end_time));
% end
% 
% % 
% % 
% Meandt0p00125 = mean(MatrixforAll0p00125');
% Stddt0p00125 = std(MatrixforAll0p00125');
% 
% 
% hold on
% 
%     pldt0p00125 = plot(time,Meandt0p00125,':c','LineWidth',3)
% errorbar(time,Meandt0p00125, Stddt0p00125,'-c','LineWidth',2)
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
%  grid on
%   xlabel('Time, hrs','FontSize',36)
%  ylabel(['Distance from the neural tube, ',char(181),'m'],'FontSize',34)
%  xticks([0,2,4,6,8, 10]);%,6])
%  xticklabels({'0','1','2','3','4','5'}); % this is because I saved the data every 30min for five hours
%  ax = gca;
%  set(gca,'linewidth',4) 
% %legend([pldt0p02,pldt0p01,pldt0p005,pldt0p0025,pldt0p00125],'\Delta t = 0.02000','\Delta t = 0.01000','\Delta t = 0.00500','\Delta t = 0.00250','\Delta t = 0.00125') 
% legend([pldt0p02,pldt0p01,pldt0p005,pldt0p0025,pldt0p00125],'\Delta t = 0.02000','\Delta t = 0.01000','\Delta t = 0.00500','\Delta t = 0.00250','\Delta t = 0.00125') 
% box on
%  %ylim([200,700])
% %  
% %  
% 
% TimeConv = end_time;
% stderror2 = zeros(1,TimeConv);
% stderror3 = zeros(1,TimeConv);
% stderror4 = zeros(1,TimeConv);
% stderror5 = zeros(1,TimeConv);
% error2 = zeros(1,TimeConv);
% error3 = zeros(1,TimeConv);
% error4 = zeros(1,TimeConv);
% error5= zeros(1,TimeConv);
% diff2 = zeros(1,TimeConv);
% diff3 = zeros(1,TimeConv);
% diff4 = zeros(1,TimeConv);
% diff5= zeros(1,TimeConv);
% 
% % relative error
% for tim = 1:TimeConv
%     error2(tim) = abs(Meandt0p02(tim) - Meandt0p00125(tim))/Meandt0p00125(tim);
%     diff2(tim) = (Meandt0p02(tim) - Meandt0p00125(tim))/Meandt0p00125(tim);
%     stderror2(tim) = sqrt((Stddt0p02(tim)^2 + Stddt0p00125(tim)^2 )/Meandt0p00125(tim)^2 - 2*(Meandt0p02(tim)-Meandt0p00125(tim))/Meandt0p00125(tim)^3 * Stddt0p00125(tim)^2  + (Meandt0p02(tim)-Meandt0p00125(tim))^2/Meandt0p00125(tim)^4 * Stddt0p00125(tim)^2    + diff2(tim)^2 - error2(tim)^2);    
%     
%     error3(tim) = abs(Meandt0p01(tim)-Meandt0p00125(tim))/Meandt0p00125(tim);
%     diff3(tim) = (Meandt0p01(tim)-Meandt0p00125(tim))/Meandt0p00125(tim);
%     stderror3(tim) = sqrt((Stddt0p01(tim)^2 + Stddt0p00125(tim)^2 )/Meandt0p00125(tim)^2 - 2*(Meandt0p01(tim)-Meandt0p00125(tim))/Meandt0p00125(tim)^3 * Stddt0p00125(tim)^2 + (Meandt0p01(tim)-Meandt0p00125(tim))^2/Meandt0p00125(tim)^4 * Stddt0p00125(tim)^2   + diff3(tim)^2 - error3(tim)^2); 
% 
%     error4(tim) = abs(Meandt0p005(tim) - Meandt0p00125(tim))/Meandt0p00125(tim);
%     diff4(tim) = (Meandt0p005(tim) - Meandt0p00125(tim))/Meandt0p00125(tim);
%     stderror4(tim) = sqrt((Stddt0p005(tim)^2 + Stddt0p00125(tim)^2 )/Meandt0p00125(tim)^2 - 2*(Meandt0p005(tim)-Meandt0p00125(tim))/Meandt0p00125(tim)^3 * Stddt0p00125(tim)^2 + (Meandt0p005(tim)-Meandt0p00125(tim))^2/Meandt0p00125(tim)^4 * Stddt0p00125(tim)^2  + diff4(tim)^2 - error4(tim)^2); 
% 
%     error5(tim) = abs(Meandt0p0025(tim) - Meandt0p00125(tim))/Meandt0p00125(tim);
%     diff5(tim) = (Meandt0p0025(tim) - Meandt0p00125(tim))/Meandt0p00125(tim);
%     stderror5(tim) = sqrt((Stddt0p0025(tim)^2 + Stddt0p00125(tim)^2 )/Meandt0p00125(tim)^2 -2*(Meandt0p0025(tim)-Meandt0p00125(tim))/Meandt0p00125(tim)^3 * Stddt0p00125(tim)^2  + (Meandt0p0025(tim)-Meandt0p00125(tim))^2/Meandt0p00125(tim)^4 * Stddt0p00125(tim)^2   + diff5(tim)^2 - error5(tim)^2);   
% 
% end
% 
% % % absolute error
% % for tim = 1:TimeConv
% %     error2(tim) = abs(Meandt0p02(tim) - Meandt0p00125(tim));
% %     diff2(tim) = (Meandt0p02(tim) - Meandt0p00125(tim));
% %     stderror2(tim) = sqrt((Stddt0p02(tim)^2 + Stddt0p00125(tim)^2 ));
% %     
% %     error3(tim) = abs(Meandt0p01(tim)-Meandt0p00125(tim));
% %     diff3(tim) = (Meandt0p01(tim)-Meandt0p00125(tim));
% %     stderror3(tim) = sqrt(Stddt0p01(tim)^2 + Stddt0p00125(tim)^2 );
% % 
% %     error4(tim) = abs(Meandt0p005(tim) - Meandt0p00125(tim));
% %     diff4(tim) = (Meandt0p005(tim) - Meandt0p00125(tim));
% %     stderror4(tim) = sqrt(Stddt0p005(tim)^2 + Stddt0p00125(tim)^2 );
% % 
% %     error5(tim) = abs(Meandt0p0025(tim) - Meandt0p00125(tim));
% %     diff5(tim) = (Meandt0p0025(tim) - Meandt0p00125(tim));
% %     stderror5(tim) = sqrt((Stddt0p0025(tim)^2 + Stddt0p00125(tim)^2 ));
% % 
% % end
% 
% % for tim = 1:TimeConv
% %     error1(tim) = (Meandt0p04(tim) - Meandt0p02(tim))/(Meandt0p02(tim) - Meandt0p01(tim));
% %     error2(tim) = (Meandt0p02(tim) - Meandt0p01(tim))/(Meandt0p01(tim) - Meandt0p005(tim));
% %     error3(tim) = (Meandt0p01(tim) - Meandt0p005(tim))/(Meandt0p005(tim) - Meandt0p0025(tim));
% %     error4(tim) =  (Meandt0p005(tim) - Meandt0p0025(tim))/(Meandt0p0025(tim) - Meandt0p00125(tim));
% % 
% % end
% 
% figure 
% 
% %er1 = plot(time,error1,'-m','Linewidth',2)
% hold on
% er2 = plot(time,error2,':g','Linewidth',2)
% errorbar(time,error2, stderror2/sqrt(N),':g','LineWidth',2)
% 
% 
% er3 = plot(time,error3,'-.k','Linewidth',2)
% errorbar(time,error3,stderror3/sqrt(N),'-.k','Linewidth',2)
% er4 = plot(time,error4,'--r','Linewidth',2)
% errorbar(time,error4,stderror4/sqrt(N),'--r','Linewidth',2)
% 
% er5 = plot(time,error5,'-b','Linewidth',2)
% errorbar(time,error5,stderror5/sqrt(N),'-b','Linewidth',2)
% 
%  set(gca,'FontSize',36)
%   xlabel('Time, hrs','FontSize',36)
%  ylabel(['Error'],'FontSize',34)
%   xticks([0,2,4,6,8, 10]);%,6])
%   xticklabels({'0','1','2','3','4','5'}); % this is because I saved the data every 30min for five hours
%  %yticks([0,0.002,0.004,0.006,0.008, 0.010]);%,6])
%  %yticklabels({'0','0.002','0.004','0.006','0.008','0.010'}); % no standard
%  %error
% %  yticks([0.000,0.002,0.004, 0.006,0.008,0.010]);%,6])
% % yticklabels({'0.000','0.002','0.004', '0.006','0.008','0.010'}); % with standard error
% %  ylim([-0.015,0.025])
%  ax = gca;
%  grid on
%  set(gca,'linewidth',4) 
%  legend([er2,er3,er4,er5],'\Delta t = 0.02000','\Delta t = 0.01000','\Delta t = 0.00500','\Delta t = 0.00250') 
% box on
 

%% Ordered timesteps, 30min, data every 5min

%% check whether it converges by comparing the ratio
N = 1000;
end_time = 7;

 figure
 % dt = 0.02, D = 5
 

%filename = sprintf('convergence/CentreMassD5eps200dt0p02CORRECT1000sim.csv');
filename = sprintf('convergence/5CellsCentreMAssD5eps200dt0p02time30min1000sim.csv');



Mass_Centre = load(filename);

for i =1:length(Mass_Centre)
    if Mass_Centre(i) == 200
        if rem(i-1,11) ~= 0
           vector(i) =  i;
        end
    end
end


MatrixforAll0p02 = zeros(end_time,N);

for i =1:N
    MatrixforAll0p02(:,i) = Mass_Centre((i-1)*end_time+1:i*(end_time));
end

Meandt0p02 = mean(MatrixforAll0p02');
Stddt0p02 = std(MatrixforAll0p02');
% 
 
hold on
time = [0:end_time-1];

pldt0p02 = plot(time,Meandt0p02,'-g','LineWidth',3);
hold on

errorbar(time,Meandt0p02,Stddt0p02,'-.g','LineWidth',2);
 %ylim([1,6]);

 set(gca,'FontSize',36)
  xlabel('Time, hrs','FontSize',36)
 ylabel(['Distance from the neural tube, ',char(181),'m'],'FontSize',34)
 ax = gca;
 set(gca,'linewidth',4) 
 %legend(pldt0p01,'\Delta t = 0.010')
 
  
 %% dt = 0.01, D = 5
 

%filename = sprintf('convergence/CentreMassD5eps200dt0p01CORRECT1000sim.csv');
filename = sprintf('convergence/5CellsCentreMAssD5eps200dt0p01time30min1000sim.csv');



Mass_Centre = load(filename);
for i =1:length(Mass_Centre)
    if Mass_Centre(i) == 200
        if rem(i-1,11) ~= 0
           vector0p01(i) =  i;
        end
    end
end


MatrixforAll0p01 = zeros(end_time,N);

for i =1:N
    MatrixforAll0p01(:,i) =Mass_Centre((i-1)*end_time+1:i*(end_time));
end

Meandt0p01 = mean(MatrixforAll0p01');
Stddt0p01 = std(MatrixforAll0p01');
% 
 
hold on


pldt0p01 = plot(time,mean(MatrixforAll0p01'),'-k','LineWidth',3);
hold on

errorbar(time,mean(MatrixforAll0p01'),std(MatrixforAll0p01'),'-.k','LineWidth',2);
 %ylim([1,6]);

 set(gca,'FontSize',36)
  xlabel('Time, hrs','FontSize',36)
 ylabel(['Distance from the neural tube, ',char(181),'m'],'FontSize',34)
 ax = gca;
 set(gca,'linewidth',4) 
 %legend(pldt0p01,'\Delta t = 0.010')
 

  %% dt = 0.005, D = 5
  
%filename = sprintf('convergence/CentreMassD5eps200dt0p005CORRECT1000sim.csv');
filename = sprintf('convergence/5CellsCentreMAssD5eps200dt0p005time30min1000sim.csv');


Mass_Centre = load(filename);


MatrixforAll0p005 = zeros(end_time,N);

for i =1:N
    MatrixforAll0p005(:,i) = Mass_Centre((i-1)*end_time+1:i*(end_time));
end

% 
% 
Meandt0p005 = mean(MatrixforAll0p005');
Stddt0p005 = std(MatrixforAll0p005');

hold on

    pldt0p005 = plot(time,mean(MatrixforAll0p005'),'-.r','LineWidth',3);
errorbar(time,mean(MatrixforAll0p005'),std(MatrixforAll0p005'),'-r','LineWidth',2);
 %ylim([1,6]);

% 
  %% dt = 0.0025, D = 5


%filename = sprintf('convergence/CentreMassD5eps200dt0p0025CORRECT1000sim.csv');
filename = sprintf('convergence/5CellsCentreMAssD5eps200dt0p0025time30min1000sim.csv');


Mass_Centre = load(filename);


MatrixforAll0p0025 = zeros(end_time,N);

for i =1:N
    MatrixforAll0p0025(:,i) = Mass_Centre((i-1)*end_time+1:i*(end_time));
end

% 
% 
Meandt0p0025 = mean(MatrixforAll0p0025');
Stddt0p0025 = std(MatrixforAll0p0025');

hold on

    pldt0p0025 = plot(time,mean(MatrixforAll0p0025'),':b','LineWidth',3);
errorbar(time,mean(MatrixforAll0p0025'),std(MatrixforAll0p0025'),'-b','LineWidth',2);
 %ylim([1,6]);


  %% dt = 0.00125, D = 5


%filename = sprintf('convergence/CentreMassD5eps200dt0p00125CORRECT1000sim.csv');
filename = sprintf('convergence/5CellsCentreMAssD5eps200dt0p00125time30min1000sim.csv');



Mass_Centre = load(filename);

MatrixforAll0p00125 = zeros(end_time,N);

for i =1:N

    MatrixforAll0p00125(:,i) = Mass_Centre((i-1)*end_time+1:i*(end_time));
end

% 
% 
Meandt0p00125 = mean(MatrixforAll0p00125');
Stddt0p00125 = std(MatrixforAll0p00125');


hold on

    pldt0p00125 = plot(time,Meandt0p00125,':c','LineWidth',3);
errorbar(time,Meandt0p00125, Stddt0p00125,'-c','LineWidth',2);
 

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
 grid on
  xlabel('Time, hrs','FontSize',36)
 ylabel(['Distance from the neural tube, ',char(181),'m'],'FontSize',34)
 xticks([0,1,2,3,4, 5,6]);%,6])
  xticklabels({'0','5','10','15','20','25','30'}); % this is because I saved the data every 30min for five hours
 ax = gca;
 set(gca,'linewidth',4) 
%legend([pldt0p02,pldt0p01,pldt0p005,pldt0p0025,pldt0p00125],'\Delta t = 0.02000','\Delta t = 0.01000','\Delta t = 0.00500','\Delta t = 0.00250','\Delta t = 0.00125') 
legend([pldt0p02,pldt0p01,pldt0p005,pldt0p0025,pldt0p00125],'\Delta t = 0.02000','\Delta t = 0.01000','\Delta t = 0.00500','\Delta t = 0.00250','\Delta t = 0.00125') 
box on
 %ylim([200,700])
%  
%  

TimeConv = end_time;
stderror2 = zeros(1,TimeConv);
stderror3 = zeros(1,TimeConv);
stderror4 = zeros(1,TimeConv);
stderror5 = zeros(1,TimeConv);
error2 = zeros(1,TimeConv);
error3 = zeros(1,TimeConv);
error4 = zeros(1,TimeConv);
error5= zeros(1,TimeConv);
diff2 = zeros(1,TimeConv);
diff3 = zeros(1,TimeConv);
diff4 = zeros(1,TimeConv);
diff5= zeros(1,TimeConv);


% I will calculate pairwise differences and use their standard deviation

difference0p02 = zeros(TimeConv,N);
meandiff0p02 = zeros(1,TimeConv);
stddiff0p02 =zeros(1,TimeConv);

difference0p01 = zeros(TimeConv,N);
meandiff0p01 = zeros(1,TimeConv);
stddiff0p01 =zeros(1,TimeConv);

difference0p005 = zeros(TimeConv,N);
meandiff0p005 = zeros(1,TimeConv);
stddiff0p005 =zeros(1,TimeConv);

difference0p0025 = zeros(TimeConv,N);
meandiff0p0025 = zeros(1,TimeConv);
stddiff0p0025 =zeros(1,TimeConv);

difference0p00125 = zeros(TimeConv,N);
meandiff0p00125 = zeros(1,TimeConv);
stddiff0p00125 =zeros(1,TimeConv);

for i=1:N
    for tim =1:TimeConv
        difference0p02(tim,i) = abs(MatrixforAll0p02(tim,i)  - MatrixforAll0p00125(tim))/MatrixforAll0p00125(tim);
        difference0p01(tim,i) = abs(MatrixforAll0p01(tim,i)  - MatrixforAll0p00125(tim))/MatrixforAll0p00125(tim);
        difference0p005(tim,i) = abs(MatrixforAll0p005(tim,i)  - MatrixforAll0p00125(tim))/MatrixforAll0p00125(tim);
        difference0p0025(tim,i) = abs(MatrixforAll0p0025(tim,i)  - MatrixforAll0p00125(tim))/MatrixforAll0p00125(tim);
        difference0p00125(tim,i) = abs(MatrixforAll0p00125(tim,i)  - MatrixforAll0p00125(tim))/MatrixforAll0p00125(tim);
    
    end
end

for tim=1:TimeConv
    meandiff0p02(tim) = mean(difference0p02(tim,:));
    stddiff0p02(tim) = std(difference0p02(tim,:));
    
    meandiff0p01(tim) = mean(difference0p01(tim,:));
    stddiff0p01(tim) = std(difference0p01(tim,:));
    
    meandiff0p005(tim) = mean(difference0p005(tim,:));
    stddiff0p005(tim) = std(difference0p005(tim,:));
    
    meandiff0p0025(tim) = mean(difference0p0025(tim,:));
    stddiff0p0025(tim) = std(difference0p0025(tim,:));
    
    
end


% % relative error
% for tim = 1:TimeConv
%     error2(tim) = abs(Meandt0p02(tim) - Meandt0p00125(tim))/Meandt0p00125(tim);
%     diff2(tim) = (Meandt0p02(tim) - Meandt0p00125(tim))/Meandt0p00125(tim);
%     stderror2(tim) = sqrt((Stddt0p02(tim)^2 + Stddt0p00125(tim)^2 )/Meandt0p00125(tim)^2 - 2*(Meandt0p02(tim)-Meandt0p00125(tim))/Meandt0p00125(tim)^3 * Stddt0p00125(tim)^2  + (Meandt0p02(tim)-Meandt0p00125(tim))^2/Meandt0p00125(tim)^4 * Stddt0p00125(tim)^2    + diff2(tim)^2 - error2(tim)^2);    
%     
%     error3(tim) = abs(Meandt0p01(tim)-Meandt0p00125(tim))/Meandt0p00125(tim);
%     diff3(tim) = (Meandt0p01(tim)-Meandt0p00125(tim))/Meandt0p00125(tim);
%     stderror3(tim) = sqrt((Stddt0p01(tim)^2 + Stddt0p00125(tim)^2 )/Meandt0p00125(tim)^2 - 2*(Meandt0p01(tim)-Meandt0p00125(tim))/Meandt0p00125(tim)^3 * Stddt0p00125(tim)^2 + (Meandt0p01(tim)-Meandt0p00125(tim))^2/Meandt0p00125(tim)^4 * Stddt0p00125(tim)^2   + diff3(tim)^2 - error3(tim)^2); 
% 
%     error4(tim) = abs(Meandt0p005(tim) - Meandt0p00125(tim))/Meandt0p00125(tim);
%     diff4(tim) = (Meandt0p005(tim) - Meandt0p00125(tim))/Meandt0p00125(tim);
%     stderror4(tim) = sqrt((Stddt0p005(tim)^2 + Stddt0p00125(tim)^2 )/Meandt0p00125(tim)^2 - 2*(Meandt0p005(tim)-Meandt0p00125(tim))/Meandt0p00125(tim)^3 * Stddt0p00125(tim)^2 + (Meandt0p005(tim)-Meandt0p00125(tim))^2/Meandt0p00125(tim)^4 * Stddt0p00125(tim)^2  + diff4(tim)^2 - error4(tim)^2); 
% 
%     error5(tim) = abs(Meandt0p0025(tim) - Meandt0p00125(tim))/Meandt0p00125(tim);
%     diff5(tim) = (Meandt0p0025(tim) - Meandt0p00125(tim))/Meandt0p00125(tim);
%     stderror5(tim) = sqrt((Stddt0p0025(tim)^2 + Stddt0p00125(tim)^2 )/Meandt0p00125(tim)^2 -2*(Meandt0p0025(tim)-Meandt0p00125(tim))/Meandt0p00125(tim)^3 * Stddt0p00125(tim)^2  + (Meandt0p0025(tim)-Meandt0p00125(tim))^2/Meandt0p00125(tim)^4 * Stddt0p00125(tim)^2   + diff5(tim)^2 - error5(tim)^2);   
% 
% end


figure 

%er1 = plot(time,error1,'-m','Linewidth',2)
hold on
er2 = plot(time,meandiff0p02,':g','Linewidth',2);
errorbar(time,meandiff0p02, stddiff0p02/sqrt(N),':g','LineWidth',2);

er3 = plot(time,meandiff0p01,'-.k','Linewidth',2);
errorbar(time,meandiff0p01, stddiff0p01/sqrt(N),'-.k','Linewidth',2);

er4 = plot(time,meandiff0p005,'--r','Linewidth',2);
errorbar(time,meandiff0p005, stddiff0p005/sqrt(N),'--r','Linewidth',2);

er5 = plot(time,meandiff0p0025,'-b','Linewidth',2);
errorbar(time,meandiff0p0025, stddiff0p0025/sqrt(N),'-b','Linewidth',2);


 set(gca,'FontSize',36)
  xlabel('Time, min','FontSize',36)
 ylabel(['Error'],'FontSize',34)
  xticks([0,1,2,3,4, 5,6]);%,6])
  xticklabels({'0','5','10','15','20','25','30'}); % this is because I saved the data every 30min for five hours
 %yticks([0,0.002,0.004,0.006,0.008, 0.010]);%,6])
 %yticklabels({'0','0.002','0.004','0.006','0.008','0.010'}); % no standard
 %error
%  yticks([0.000,0.002,0.004, 0.006,0.008,0.010]);%,6])
% yticklabels({'0.000','0.002','0.004', '0.006','0.008','0.010'}); % with standard error
%  ylim([-0.015,0.025])
 ax = gca;
 grid on
 set(gca,'linewidth',4) 
 legend([er2,er3,er4,er5],'\Delta t = 0.02000','\Delta t = 0.01000','\Delta t = 0.00500','\Delta t = 0.00250') 
box on
