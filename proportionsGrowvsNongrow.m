%% plot average of twenty simulations, for the model with fixed and growing domain
% I have for variance 1 DATAnon-growing-width120-normal1var.csv, and variance 1p5 DATAnon-growing-width120-normal1p5var.csv
% also any comparisons for different proportions

N = 20; % number of simulations

% % non growing
%sim1 = 'DATA for 24hrs experiments Basic model/DATAnon-growing-width120-normal1p5var.csv'; %% for 24hrs
%sim1 = 'DATA for 18hrs experiments for Basic model/DATA18hrsnon-growing-width120-normal1p0var.csv'; %% for 18hrs
sim1 = 'DATAD20eps1m2CoACiLGrowingDomain.csv' %%% if growing domain add to the end of all files: GrowingDomain



M1 = csvread(sim1);
%M12 = csvread(sim12);
%M1 = (M1(:,1) + M12(:,1))*0.5;
M1 = M1(:,1)/20;

for i = 2:N
    
  % filename =
  % sprintf('DATA for 24hrs experiments Basic model/sepdata-non-growing-width120-normal1p5var%i.csv',i-1); for
  % 24hrs
  % filename = sprintf('DATA for 18hrs experiments for Basic model/sepdata18hrs-non-growing-width120-normal1p0var%i.csv',i-1); % for 18hrs
    filename = sprintf('sepdataD20eps1m2CoACiLGrowingDomain%i.csv',i-1);
    sepdata = load(filename);
        alldata(:,i) = sepdata; 

end

std1 = std(alldata');

%% data set number 2

sim2 = 'DATAD20eps5m2CoACiLGrowingDomain.csv'

M2 = csvread(sim2);
M2 = M2(:,1)/20;

for i = 2:N
    
    filename = sprintf('sepdataD20eps5m2CoACiLGrowingDomain%i.csv',i-1);
    sepdata = load(filename);
        alldata(:,i) = sepdata; 

end

std2 = std(alldata');

% 
% 
%% dataset number 3

sim3 = 'DATAD20eps10m2CoACiLGrowingDomain.csv'

M3 = csvread(sim3);
M3 = M3(:,1)/20;

for i = 2:N
    
    filename = sprintf('sepdataD20eps10m2CoACiLGrowingDomain%i.csv',i-1);
    sepdata = load(filename);
        alldata(:,i) = sepdata; 

end

std3 = std(alldata');
% % 
% %% dataset number 4
% 
% sim4 = 'DATAD5eps1m2CoACiLGrowingDomain.csv'
% 
% M4 = csvread(sim4);
% M4 = M4(:,1)/20;
% 
% for i = 2:N
%     
%    filename = sprintf('sepdataD5eps1m2CoACiLGrowingDomain%i.csv',i-1);
%     sepdata = load(filename);
%         alldata(:,i) = sepdata; 
% 
% end
% 
% std4 = std(alldata');
% 
% 



% 
% %% growing
% %sim2 = 'DATA for 24hrs experiments Basic model/DATAgrowing-width120-normal1p5var.csv'; % 24hrs
% sim2 = 'DATA for 18hrs experiments for Basic model/DATA18hrsgrowing-width120-normal1p0var.csv'; % 18hrs
% 
% M2 = csvread(sim2);
% M2 = M2(:,1)/20;
% 
% % 
%  for i = 2:N
% %    
%    %filename = sprintf('DATA for 24hrs experiments Basic model/sepdata-growing-width120-normal1p5var%i.csv',i-1);
%    %% 24hrs
%    filename = sprintf('DATA for 18hrs experiments for Basic model/sepdata18hrs-growing-width120-normal1p0var%i.csv',i-1); % 18hrs
% 
%     sepdata2 = load(filename);
%         alldata2(:,i) = sepdata2; 
% % 
%  end
% % 
%  std2 = std(alldata2');
% % 
% 

%% ploting for growth no growth
% figure
% 
% x = [0:55:1098];
% 
% 
%  h1 = plot(x,M1,'linewidth',4)
%  hold on
%    errorbar (x,M1, std1,'b.','linewidth',2)
%     
%  h2 = plot(x,M2,'--','linewidth',4)
%     errorbar (x,M2, std2,'k.','linewidth',2)    
% 
%  xlabel(['Distance from the neural tube, ',char(181),'m'],'FontSize',14)
% set(gca,'linewidth',2)
% ylabel(['Number of cells (per 55 ',char(181),'m)'],'FontSize',14)
% set(gca,'FontSize',36)
%  ax = gca;
% 
% 
%  box on
% legend([h1,h2],'No growth','Domain grows');
% 
%  set(gca,'linewidth',4)   


%% plotting many others

figure

x = [0:55:1098];


 h1 = plot(x,M1,'linewidth',4)
 hold on
 % errorbar (x,M1, std1,'k.','linewidth',2)
    
 h2 = plot(x,M2,'--','linewidth',4)
   % errorbar (x,M2, std2,'k.','linewidth',2)   
%     
 h3 = plot(x,M3,'-.','linewidth',4)
   % errorbar (x,M3, std3,'k.','linewidth',2)  
%     
%   h4 = plot(x,M4,':','linewidth',4)
% %  % errorbar (x,M4, std4,'k.','linewidth',2)  

 xlabel(['Distance from the neural tube, ',char(181),'m'],'FontSize',14)
set(gca,'linewidth',2)
ylabel(['Number of cells (per 55 ',char(181),'m)'],'FontSize',14)
set(gca,'FontSize',36)
 ax = gca;

%ylim([0,30])
 box on
%legend([h1,h2,h3,h4],'D=0.5','D=1.0','D=2.0','D=5.0');
%legend([h1,h2,h3,h4],'\epsilon=0.5','\epsilon=1.0','\epsilon=2.0','\epsilon=5.0');
%legend([h2,h3],'m=2','m=12')
legend([h1,h2,h3],'\epsilon=1.0','\epsilon=5.0','\epsilon=10.0')
 set(gca,'linewidth',4)   
