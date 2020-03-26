%% check whether it converges by comparing the ratio


%figure 

% TimeConv = 10; % check convergence at this time, 10 max
% %% dt = 0.010, D = 0
% 
% filename = sprintf('convergence/CoACiLeps200D0masscentr0p01dt.csv');
% Mass_Centre = load(filename);
% Mass_Centreh = Mass_Centre(1:20);
% 
% %% dt = 0.0050, D = 0
% 
% filename = sprintf('convergence/CoACiLeps200D0masscentr0p005dt.csv');
% Mass_Centre = load(filename);
% Mass_Centreh2 = Mass_Centre(1:20);
% 
% %% dt = 0.0025, D = 0
% 
% filename = sprintf('convergence/CoACiLeps200D0masscentr0p0025dt.csv');
% Mass_Centre = load(filename);
% Mass_Centreh4 = Mass_Centre(1:20);
% 
% 
% Ratio = (Mass_Centreh(TimeConv) - Mass_Centreh2(TimeConv))/(Mass_Centreh2(TimeConv) - Mass_Centreh4(TimeConv))
% 
% 




%  xlabel('Time, hrs','FontSize',36)
% % ylabel('positions, ',14)
%  pldt0p1 = plot(time,Mass_Centre,'-.g','LineWidth',3)
%  %ylim([1,6]);
%   box on
%  set(gca,'FontSize',36)
%  ylabel(['Distance from the neural tube, ',char(181),'m'],'FontSize',34)
%  ax = gca;
%  set(gca,'linewidth',4) 


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
 %% dt = 0.01, D = 5
 
filename = sprintf('convergence/CoACiLeps200D5masscentr20sim0p01dt.csv');
Mass_Centre = load(filename);


MatrixforAll = zeros(20,20);

for i =1:20
    MatrixforAll(:,i) = Mass_Centre((i-1)*20+1:i*(20));
end


% 
 
figure
time = [1:20];

pldt0p01 = plot(time,mean(MatrixforAll'),'-k','LineWidth',3)
errorbar(time,mean(MatrixforAll'),std(MatrixforAll'),'-k','LineWidth',2)
 %ylim([1,6]);

 set(gca,'FontSize',36)
  xlabel('Time, hrs','FontSize',36)
 ylabel(['Distance from the neural tube, ',char(181),'m'],'FontSize',34)
 ax = gca;
 set(gca,'linewidth',4) 
 %legend(pldt0p01,'\Delta t = 0.010')
 ylim([200,700])
 
 %% dt = 0.001, D = 5
 
filename = sprintf('convergence/CoACiLeps200D5masscentr20sim0p001dt.csv');
Mass_Centre = load(filename);


MatrixforAll = zeros(20,20);

for i =1:20
    MatrixforAll(:,i) = Mass_Centre((i-1)*20+1:i*(20));
end

% 
% 
hold on

    pldt0p001 = plot(time,mean(MatrixforAll'),'--r','LineWidth',3)
errorbar(time,mean(MatrixforAll'),std(MatrixforAll'),'-r','LineWidth',2)
 %ylim([1,6]);

 set(gca,'FontSize',36)
  xlabel('Time, hrs','FontSize',36)
 ylabel(['Distance from the neural tube, ',char(181),'m'],'FontSize',34)
 ax = gca;
 set(gca,'linewidth',4) 
legend([pldt0p01,pldt0p001],'\Delta t = 0.010','\Delta t = 0.001')

 ylim([200,700])
%  
%  
%  