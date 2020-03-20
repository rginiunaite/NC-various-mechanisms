%% convergence

%% I do it for D=0, no stochasticity, eps 200.

%% dt =0.1

filename = sprintf('convergence/CoACiLeps200D0positions0p1dt.csv');
data = load(filename);
xcoord = data(:,1);

sizevec = size(xcoord);

cell=zeros(sizevec(1)/5,5);

for i=1:5
cell(:,i) = xcoord(i:5:end);
end

time =[1:20];

figure
for i=1:5
    hold on
    h1 = plot(time,cell(1:20,i),'-.g','LineWidth',3);
end

%% dt = 0.01

filename = sprintf('convergence/CoACiLeps200D0positions0p01dt.csv');
data = load(filename);
xcoord = data(:,1);

sizevec = size(xcoord);

cell=zeros(sizevec(1)/5,5);

for i=1:5
cell(:,i) = xcoord(i:5:end);
end

time =[1:20];


for i=1:5
    hold on
    h2 = plot(time,cell(1:20,i),'-b','LineWidth',3);
end


 xlabel('Time, hrs','FontSize',36)
% ylabel('positions, ',14)
 
 
 %ylim([1,6]);
  box on
 set(gca,'FontSize',36)
 ylabel(['Distance from the neural tube, ',char(181),'m'],'FontSize',34)
 ax = gca;
 set(gca,'linewidth',4) 
 
 
%% dt = 0.001 
 
filename = sprintf('convergence/CoACiLeps200D0positions0p001dt.csv');
data = load(filename);
xcoord = data(:,1);

sizevec2 = size(xcoord);

cell2=zeros(sizevec2(1)/5,5);

for i=1:5
cell2(:,i) = xcoord(i:5:end);
end

time =[1:20];

for i=1:5
    hold on
    h3= plot(time,cell2(1:20,i),'r--','LineWidth',3);
end

legend([h1,h2,h3],'\Delta t = 0.100','\Delta t = 0.010','\Delta t = 0.001')
ylim([0,600]);
 