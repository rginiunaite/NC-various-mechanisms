% %% plot mean and variance of twenty simulations, varying variance of normal distribution, whether the 
% %%filopodia go only inside the domain or everywhere
% 
% 
sim1 = 'fixed-width120-uni2.csv';
M1 = csvread(sim1);
M1 = M1(:,1);
x = [0:55:1099];



% d=dir('fixed-width120-uni*.csv');   % files starting like this
% n=length(d);        % how many there were
% data=zeros(length(M1),n);     % preallocate a cell array to hold results
% temp = zeros (length(M1),2); % need this because data is in two columns
% for i=1:n
%     temp = csvread(d(i).name);
%   data(:,i)= temp(:,1);  % read each file
% end
% 
% average = mean(data'); %
% 
% std_average = std(average(1:7));
% 
% standarDev = std(data');
% 
% 
% 
% figure 
% % 
% %h1 = plot(x,average,'-')
% 
% err = standarDev;
% h1 = errorbar(x,average,err)
% 
% 
% h1.LineWidth =4;
% h2.LineWidth =4;
%  h3.LineWidth =4;
% % h4.LineWidth =6;
% % h5.LineWidth =6;
% % h6.LineWidth =6;
% 
% xlabel('Distance from the neural tube, \mu m','FontSize',14)
% set(gca,'linewidth',2)
% ylabel('Number of cells','FontSize',14)
% set(gca,'FontSize',36)
%  ax = gca;
%  
%  box on
% 
%  set(gca,'linewidth',4)
%  
%  
%  %% add another curve
%  
% 
% d=dir('growing-width120-uni*.csv');   % files starting like this
% n=length(d);        % how many there were
% data=zeros(length(M1),n);     % preallocate a cell array to hold results
% temp = zeros (length(M1),2); % need this because data is in two columns
% for i=1:n
%     temp = csvread(d(i).name);
%   data(:,i)= temp(:,1);  % read each file
% end
% 
% average = mean(data'); %
% 
% std_average1 = std(average);
% 
% standarDev = std(data');
% 
% x = [0:55:1099];
% 
% hold on 
% % 
% %h1 = plot(x,average,'-')
% 
% err = standarDev;
% h1 = errorbar(x,average,err,'--')
% 
% 
% h1.LineWidth =4;
% h2.LineWidth =4;
%  h3.LineWidth =4;
% % h4.LineWidth =6;
% % h5.LineWidth =6;
% % h6.LineWidth =6;
% 
% xlabel('Distance from the neural tube, \mu m','FontSize',14)
% set(gca,'linewidth',2)
% ylabel('Number of cells','FontSize',14)
% set(gca,'FontSize',36)
%  ax = gca;
%  
%  box on
% 
%  set(gca,'linewidth',4)
%  
%  
%   %% add another curve
% 
% d=dir('growing-width120-uni-inside*.csv');   % files starting like this
% n=length(d);        % how many there were
% data=zeros(length(M1),n);     % preallocate a cell array to hold results
% temp = zeros (length(M1),2); % need this because data is in two columns
% for i=1:n
%     temp = csvread(d(i).name);
%   data(:,i)= temp(:,1);  % read each file
% end
% 
% average = mean(data'); %
% 
% std_average2 = std(average);
% 
% standarDev = std(data');
% 
% 
% hold on 
% % 
% %h1 = plot(x,average,'-')
% 
% err = standarDev;
% h1 = errorbar(x,average,err,'-.')
% 
% 
% h1.LineWidth =4;
% h2.LineWidth =4;
%  h3.LineWidth =4;
% % h4.LineWidth =6;
% % h5.LineWidth =6;
% % h6.LineWidth =6;
% 
% xlabel('Distance from the neural tube, \mu m','FontSize',14)
% set(gca,'linewidth',2)
% ylabel('Number of cells','FontSize',14)
% set(gca,'FontSize',36)
%  ax = gca;
%  
%  box on
% 
%  set(gca,'linewidth',4)
%  
%  
%  

% %% add another curve



% d=dir('growing-width120-normal1point5var*.csv');   % files starting like this
% n=length(d);        % how many there were
% data=zeros(length(M1),n);     % preallocate a cell array to hold results
% temp = zeros (length(M1),2); % need this because data is in two columns
% for i=1:n
%     temp = csvread(d(i).name);
%   data(:,i)= temp(:,1);  % read each file
% end
% 
% average = mean(data'); %
% 
% std_average2 = std(average);
% 
% standarDev = std(data');
% 
% 
% figure% 
% %h1 = plot(x,average,'-')
% 
% err = standarDev;
% h1 = errorbar(x,average,err,'-')
% 
% 
% 
% h1.LineWidth =4;
% h2.LineWidth =4;
%  h3.LineWidth =4;
% % h4.LineWidth =6;
% % h5.LineWidth =6;
% % h6.LineWidth =6;
% 
% xlabel('Distance from the neural tube, \mu m','FontSize',14)
% set(gca,'linewidth',2)
% ylabel('Number of cells','FontSize',14)
% set(gca,'FontSize',36)
%  ax = gca;
%  
%  box on
% 
%  set(gca,'linewidth',4)
 
%  
%  
% 
% 
   %% add another curve

d=dir('growing-width120-normal2var*.csv');   % files starting like this
n=length(d);        % how many there were
data=zeros(length(M1),n);     % preallocate a cell array to hold results
temp = zeros (length(M1),2); % need this because data is in two columns
for i=1:n
    temp = csvread(d(i).name);
  data(:,i)= temp(:,1);  % read each file
end

average = mean(data'); %

std_average2 = std(average);

standarDev = std(data');


figure
% 
%h1 = plot(x,average,'-')

err = standarDev;
h1 = errorbar(x,average,err,':')


h1.LineWidth =4;
h2.LineWidth =4;
 h3.LineWidth =4;
% h4.LineWidth =6;
% h5.LineWidth =6;
% h6.LineWidth =6;

xlabel('Distance from the neural tube, \mu m','FontSize',14)
set(gca,'linewidth',2)
ylabel('Number of cells','FontSize',14)
set(gca,'FontSize',36)
 ax = gca;
 
 box on

 set(gca,'linewidth',4)
 
 


%% add another curve


d=dir('growing-width180-normal2var*.csv');   % files starting like this
%d=dir('growing-width120-normal2point5var*.csv');   % files starting like this
n=length(d);        % how many there were
data=zeros(length(M1),n);     % preallocate a cell array to hold results
temp = zeros (length(M1),2); % need this because data is in two columns
for i=1:n
    temp = csvread(d(i).name);
  data(:,i)= temp(:,1);  % read each file
end

average = mean(data'); %

std_average2 = std(average);

standarDev = std(data');


hold on 
% 
%h1 = plot(x,average,'-')

err = standarDev;
h1 = errorbar(x,average,err,'-.')


h1.LineWidth =4;
h2.LineWidth =4;
 h3.LineWidth =4;
% h4.LineWidth =6;
% h5.LineWidth =6;
% h6.LineWidth =6;
xlabel('Distance from the neural tube, \mu m','FontSize',14)
set(gca,'linewidth',2)
ylabel('Number of cells','FontSize',14)
set(gca,'FontSize',36)
 ax = gca;
 
 box on

 set(gca,'linewidth',4)
 
%  
%  
 legend ('\sigma^2 = 1.5', '\sigma^2 = 2.0','\sigma^2 = 2.5')


ylim([0 20]);

% legend for different width domains
 legend ('domain width = 120 \mu m', 'domain width = 180 \mu m')



%% legend for the first comparison

%legend ('Domain - fixed, movement - random','Domain - growing, movement - random', 'Domain - growing, movement - filopodia inside','Domain - growing, filopodia - biased')

