% centre of mass

time = [1:20];
 figure 
% %% dt = 0.1, D = 0
% 
% filename = sprintf('convergence/CoACiLeps200D0masscentr0p1dt.csv');
% Mass_Centre = load(filename);
% Mass_Centre = Mass_Centre(1:20);
% 
%  xlabel('Time, hrs','FontSize',36)
% % ylabel('positions, ',14)
%  
% 
%  
%  pldt0p1 = plot(time,Mass_Centre,'-.g','LineWidth',3)
%  %ylim([1,6]);
%   box on
%  set(gca,'FontSize',36)
%  ylabel(['Distance from the neural tube, ',char(181),'m'],'FontSize',34)
%  ax = gca;
%  set(gca,'linewidth',4) 
% 
% 
% %% dt = 0.01, D = 0
% 
% filename = sprintf('convergence/CoACiLeps200D0masscentr0p01dt.csv');
% Mass_Centre = load(filename);
% Mass_Centre = Mass_Centre(1:20);
% 
%  xlabel('Time, hrs','FontSize',36)
% % ylabel('positions, ',14)
%  
% hold on
%  
%  pldt0p01 = plot(time,Mass_Centre,'-b','LineWidth',3)
%  %ylim([1,6]);
%   box on
%  set(gca,'FontSize',36)
%  ylabel(['Distance from the neural tube, ',char(181),'m'],'FontSize',34)
%  ax = gca;
%  set(gca,'linewidth',4) 
%  
% 
% %% dt = 0.001, D = 0
% 
% filename = sprintf('convergence/CoACiLeps200D0masscentr0p001dt.csv');
% Mass_Centre = load(filename);
% Mass_Centre = Mass_Centre(1:20);
% 
%  xlabel('Time, hrs','FontSize',36)
% % ylabel('positions, ',14)
%  
% hold on
%  
%  pldt0p001 = plot(time,Mass_Centre,'r--','LineWidth',3)
%  %ylim([1,6]);
%   box on
%  set(gca,'FontSize',36)
%  ylabel(['Distance from the neural tube, ',char(181),'m'],'FontSize',34)
%  ax = gca;
%  set(gca,'linewidth',4) 
%  
% %legend([pldt0p1,pldt0p01,pldt0p001],'\Delta t = 0.100','\Delta t = 0.010','\Delta t = 0.001')
% 
%  %%
%  % Models with stochasticity
%  %%
%  
% %% dt = 0.1, D = 5
%  
% filename = sprintf('convergence/CoACiLeps200D5masscentr0p1dt.csv');
% Mass_Centre = load(filename);
% Mass_Centre = Mass_Centre(1:20);
% 
%  xlabel('Time, hrs','FontSize',36)
% % ylabel('positions, ',14)
%  
% 
% figure
% 
% pldt0p1 = plot(time,Mass_Centre,'-.g','LineWidth',3)
%  %ylim([1,6]);
% 
%  set(gca,'FontSize',36)
%  ylabel(['Distance from the neural tube, ',char(181),'m'],'FontSize',34)
%  ax = gca;
%  set(gca,'linewidth',4) 
%  
%  
%  %% dt = 0.01, D = 5
%  
% filename = sprintf('convergence/CoACiLeps200D5masscentr0p01dt.csv');
% Mass_Centre = load(filename);
% Mass_Centre = Mass_Centre(1:20);
% 
%  xlabel('Time, hrs','FontSize',36)
% % ylabel('positions, ',14)
%  
% hold on
%  
% pldt0p01 = plot(time,Mass_Centre,'-b','LineWidth',3)
%  %ylim([1,6]);
%   box on
%  set(gca,'FontSize',36)
%  ylabel(['Distance from the neural tube, ',char(181),'m'],'FontSize',34)
%  ax = gca;
%  set(gca,'linewidth',4) 
% 
%  %% dt = 0.01, D = 5
%  
% filename = sprintf('convergence/CoACiLeps200D5masscentr0p001dt.csv');
% Mass_Centre = load(filename);
% Mass_Centre = Mass_Centre(1:20);
% 
%  xlabel('Time, hrs','FontSize',36)
% % ylabel('positions, ',14)
%  
% hold on
%  
% pldt0p001 = plot(time,Mass_Centre,'r--','LineWidth',3)
%  %ylim([1,6]);
%   box on
%  set(gca,'FontSize',36)
%  ylabel(['Distance from the neural tube, ',char(181),'m'],'FontSize',34)
%  ax = gca;
%  set(gca,'linewidth',4) 
% 
%  
%  
%  legend([pldt0p1,pldt0p01,pldt0p001],'\Delta t = 0.100','\Delta t = 0.010','\Delta t = 0.001')
% %  


 %%
 % Models with stochasticity, averages over twenty
 %%
 
% %% dt = 0.1, D = 5
%  
% filename = sprintf('convergence/CoACiLeps200D5masscentr20sim0p1dt.csv');
% Mass_Centre = load(filename);
% 
% 
% MatrixforAll = zeros(20,20);
% 
% for i =1:20
%     MatrixforAll(:,i) = Mass_Centre((i-1)*20+1:i*(20));
% end
% 
% 
%  xlabel('Time, hrs','FontSize',36)
%  
% 
% for i = 1:20
%     hold on
%     pldt0p1 = plot(time,MatrixforAll(:,i),'-k','LineWidth',3)
% end
%  ylim([200,700]);
% 
%  set(gca,'FontSize',36)
%  ylabel(['Distance from the neural tube, ',char(181),'m'],'FontSize',34)
%  ax = gca;
%  set(gca,'linewidth',4) 
%  legend(pldt0p1,'\Delta t = 0.100')
%  
%  
%  
%  %% dt = 0.01, D = 5
%  
% filename = sprintf('convergence/CoACiLeps200D5masscentr20sim0p01dt.csv');
% Mass_Centre = load(filename);
% 
% 
% MatrixforAll = zeros(20,20);
% 
% for i =1:20
%     MatrixforAll(:,i) = Mass_Centre((i-1)*20+1:i*(20));
% end
% 
% 
% 
%  
% figure
% for i = 1:20
%     hold on
%     pldt0p01 = plot(time,MatrixforAll(:,i),'-k','LineWidth',3)
% end
%  %ylim([1,6]);
% 
%  set(gca,'FontSize',36)
%   xlabel('Time, hrs','FontSize',36)
%  ylabel(['Distance from the neural tube, ',char(181),'m'],'FontSize',34)
%  ax = gca;
%  set(gca,'linewidth',4) 
%  legend(pldt0p01,'\Delta t = 0.010')
%  ylim([200,700])
%  
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
% 
% 
%  
% figure 
% 
% for i = 1:20
%     hold on
%     pldt0p001 = plot(time,MatrixforAll(:,i),'-k','LineWidth',3)
% end
%  %ylim([1,6]);
% 
%  set(gca,'FontSize',36)
%   xlabel('Time, hrs','FontSize',36)
%  ylabel(['Distance from the neural tube, ',char(181),'m'],'FontSize',34)
%  ax = gca;
%  set(gca,'linewidth',4) 
%  legend(pldt0p001,'\Delta t = 0.001')
%  ylim([200,700])
 
 
 
 
  %%
 % D=1 \epsilon = 200, sho how different it can be.
 %%
 
 %% nseed = 10, fast dynamics
 
filename = sprintf('convergence/CoACiLeps200D1nseed10fast.csv');
Mass_Centre = load(filename);


length = size(Mass_Centre);

time1 = [1: length(1)];

 xlabel('Time, hrs','FontSize',36)
 

    pfast = plot(time1,Mass_Centre,'-b','LineWidth',3)
%ylim([200,700]);

 set(gca,'FontSize',36)
 ylabel(['Distance from the neural tube, ',char(181),'m'],'FontSize',34)
 ax = gca;
 set(gca,'linewidth',4) 
% legend(pldt0p1,'\Delta t = 0.100')
 
   %% nseed = 15, slow dynamics
 
filename = sprintf('convergence/CoACiLeps200D1nseed16slow.csv');
Mass_Centre = load(filename);

length = size(Mass_Centre);

time2 = [1: length(1)];


    hold on
    pslow = plot(time2,Mass_Centre,'--r','LineWidth',3)

 %ylim([1,6]);

 set(gca,'FontSize',36)
  xlabel('Time, hrs','FontSize',36)
 ylabel(['Distance from the neural tube, ',char(181),'m'],'FontSize',34)
 ax = gca;
 set(gca,'linewidth',4) 
% legend(pldt0p01,'\Delta t = 0.010')
 %ylim([200,700])
 
 
 
 %% nseed = 9, get stuck, steady state
 
filename = sprintf('convergence/CoACiLeps200D1nseed9slow.csv');
Mass_Centre = load(filename);

length = size(Mass_Centre);

time2 = [1: length(1)];


    hold on
    pstuck= plot(time2,Mass_Centre,':g','LineWidth',3)

 %ylim([1,6]);

 set(gca,'FontSize',36)
  xlabel('Time, hrs','FontSize',36)
 ylabel(['Distance from the neural tube, ',char(181),'m'],'FontSize',34)
 ax = gca;
 set(gca,'linewidth',4) 

 
 legend([pfast,pslow,pstuck],'sim 1','sim 2','sim 3')
 %ylim([200,700])
 
 xlim([0,100]);
