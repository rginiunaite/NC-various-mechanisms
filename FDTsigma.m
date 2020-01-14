% FDT versus sigma, from Basic.cpp

%cells1 = 'DATA for 24hrs experiments Basic model/sigma05width120nongrow.txt' % 24hrs, sigma = 0.5
cells1 = 'DATA for 18hrs experiments for Basic model/FDT18hrssigma05width120.txt' % 18hrs
cellpos1 = csvread(cells1);

xcoord1 = cellpos1(:,1);

average(1) = mean(xcoord1);
stddev(1) = std(xcoord1);


%cells1 = 'DATA for 24hrs experiments Basic model/sigma10width120nongrow.txt' % 24hrs, sigma = 1.0
cells1 = 'DATA for 18hrs experiments for Basic model/FDT18hrssigma10width120.txt' % 18hrs

cellpos1 = csvread(cells1);

xcoord1 = cellpos1(:,1);

average(2) = mean(xcoord1);
stddev(2) = std(xcoord1);

%cells1 = 'DATA for 24hrs experiments Basic model/sigma15width120nongrow.txt' % 24hrs, sigma = 1.5 
cells1 = 'DATA for 18hrs experiments for Basic model/FDT18hrssigma15width120.txt' % 18hrs
cellpos1 = csvread(cells1);

xcoord1 = cellpos1(:,1);

average(3) = mean(xcoord1);
stddev(3) = std(xcoord1);

%cells1 = 'DATA for 24hrs experiments Basic model/sigma20width120nongrow.txt' % 24hrs, sigma = 2.0
cells1 = 'DATA for 18hrs experiments for Basic model/FDT18hrssigma20width120.txt' % 18hrs

cellpos1 = csvread(cells1);

xcoord1 = cellpos1(:,1);

average(4) = mean(xcoord1);
stddev(4) = std(xcoord1);


sigma = [0.5, 1.0, 1.5, 2.0];

FDT = average

figure
hbD = bar(sigma, FDT,'y'); % physical
hold off
figure
p(1) = plot(sigma,FDT,'-.','linewidth',4)
hold on
% For each set of bars, find the centers of the bars, and write error bars
pause(0.1); %pause allows the figure to be created

for ib = 1:numel(hbD)
    %XData property is the tick labels/group centers; XOffset is the offset
    %of each distinct group
    xData = hbD(ib).XData+hbD(ib).XOffset;
    errorbar(xData,FDT(ib,:),stddev(ib,:),'k.','linewidth',2)
end



set(gca,'FontSize',30)
ax = gca;

xlabel(['\sigma^2'],'FontSize',14)
ylabel(['Furthest distance travelled, ',char(181),'m'],'FontSize',14)

 set(gca,'FontSize',30)
ax = gca;


 box on

 set(gca,'linewidth',4)

xticks([0.0,0.5,1.0,1.5,2.0])%,2.0])
xticklabels({'0.0','0.5', '1.0','1.5','2.0'});%, '2.0'})
xlim([0,2.5])

% %% width 240
% 
% 
% 
% %cells1 = 'DATA for 24hrs experiments Basic model/sigma05width240nongrow.txt' % 
% cells1 = 'DATA for 18hrs experiments for Basic model/FDT18hrssigma05width240.txt' % 18hrs
% 
% cellpos1 = csvread(cells1);
% 
% xcoord1 = cellpos1(:,1);
% 
% average(1) = mean(xcoord1);
% stddev(1) = std(xcoord1);
% 
% 
% %cells1 = 'DATA for 24hrs experiments Basic model/sigma10width240nongrow.txt' % 
% cells1 = 'DATA for 18hrs experiments for Basic model/FDT18hrssigma10width240.txt' % 18hrs
% 
% cellpos1 = csvread(cells1);
% 
% xcoord1 = cellpos1(:,1);
% 
% average(2) = mean(xcoord1);
% stddev(2) = std(xcoord1);
% 
% %cells1 = 'DATA for 24hrs experiments Basic model/sigma15width240nongrow.txt' % 
% cells1 = 'DATA for 18hrs experiments for Basic model/FDT18hrssigma15width240.txt' % 18hrs
% cellpos1 = csvread(cells1);
% 
% xcoord1 = cellpos1(:,1);
% 
% average(3) = mean(xcoord1);
% stddev(3) = std(xcoord1);
% 
% %cells1 = 'DATA for 24hrs experiments Basic model/sigma20width240nongrow.txt' % 
% cells1 = 'DATA for 18hrs experiments for Basic model/FDT18hrssigma20width240.txt' % 18hrs
% 
% cellpos1 = csvread(cells1);
% 
% xcoord1 = cellpos1(:,1);
% 
% average(4) = mean(xcoord1);
% stddev(4) = std(xcoord1);
% 
% 
% sigma = [0.5, 1.0, 1.5, 2.0];
% 
% FDT = average
% 
% 
% hold on
% 
% p(2) = plot(sigma,FDT,'linewidth',4)
% %hbD = bar(sigma, FDT,'y'); % physical
% hold on
% % For each set of bars, find the centers of the bars, and write error bars
% pause(0.1); %pause allows the figure to be created
% 
% for ib = 1:numel(hbD)
%     %XData property is the tick labels/group centers; XOffset is the offset
%     %of each distinct group
%     xData = hbD(ib).XData+hbD(ib).XOffset;
%     errorbar(xData,FDT(ib,:),stddev(ib,:),'k.','linewidth',2)
% end
% 
% %% width 60
% 
% 
% 
% %cells1 = 'DATA for 24hrs experiments Basic model/sigma05width60nongrow.txt' % 
% cells1 = 'DATA for 18hrs experiments for Basic model/FDT18hrssigma05width60.txt' % 18hrs
% cellpos1 = csvread(cells1);
% 
% xcoord1 = cellpos1(:,1);
% 
% average(1) = mean(xcoord1);
% stddev(1) = std(xcoord1);
% 
% 
% %cells1 = 'DATA for 24hrs experiments Basic model/sigma10width60nongrow.txt' % control, change D3 P6
% cells1 = 'DATA for 18hrs experiments for Basic model/FDT18hrssigma10width60.txt' % 18hrs
% cellpos1 = csvread(cells1);
% 
% xcoord1 = cellpos1(:,1);
% 
% average(2) = mean(xcoord1);
% stddev(2) = std(xcoord1);
% 
% %cells1 = 'DATA for 24hrs experiments Basic model/sigma15width60nongrow.txt' % 
% cells1 = 'DATA for 18hrs experiments for Basic model/FDT18hrssigma15width60.txt' % 18hrs
% cellpos1 = csvread(cells1);
% 
% xcoord1 = cellpos1(:,1);
% 
% average(3) = mean(xcoord1);
% stddev(3) = std(xcoord1);
% 
% %cells1 = 'DATA for 24hrs experiments Basic model/sigma20width60nongrow.txt' % 
% cells1 = 'DATA for 18hrs experiments for Basic model/FDT18hrssigma20width60.txt' % 18hrs
% cellpos1 = csvread(cells1);
% 
% xcoord1 = cellpos1(:,1);
% 
% average(4) = mean(xcoord1);
% stddev(4) = std(xcoord1);
% 
% 
% sigma = [0.5, 1.0, 1.5, 2.0];
% 
% FDT = average
% 
% 
% hold on
% 
% p(3) = plot(sigma,FDT,'--','linewidth',4)
% %hbD = bar(sigma, FDT,'y'); % physical
% hold on
% % For each set of bars, find the centers of the bars, and write error bars
% pause(0.1); %pause allows the figure to be created
% 
% for ib = 1:numel(hbD)
%     %XData property is the tick labels/group centers; XOffset is the offset
%     %of each distinct group
%     xData = hbD(ib).XData+hbD(ib).XOffset;
%     errorbar(xData,FDT(ib,:),stddev(ib,:),'k.','linewidth',2)
% end
% %legend(h(2:3),'This two','This three');
% 
%  legend(p,['L_y = 120',char(181),'m'],['L_y = 240',char(181),'m'],['L_y = 60',char(181),'m']);
% %legend('Width=120','','Width=240','','Width=60')


% %% compare with growing domain
% 
% %cells1 = 'DATA for 24hrs experiments Basic model/sigma05width120domaingrow.txt' % 24hrs, sigma = 0.5
% cells1 = 'DATA for 18hrs experiments for Basic model/FDT18hrssigma05width120grows.txt' % 18hrs
% cellpos1 = csvread(cells1);
% 
% xcoord1 = cellpos1(:,1);
% 
% average(1) = mean(xcoord1);
% stddev(1) = std(xcoord1);
% 
% 
% %cells1 = 'DATA for 24hrs experiments Basic model/sigma10width120domaingrow.txt' % 24hrs, sigma = 1.0
% cells1 = 'DATA for 18hrs experiments for Basic model/FDT18hrssigma10width120grows.txt' % 18hrs
% cellpos1 = csvread(cells1);
% 
% xcoord1 = cellpos1(:,1);
% 
% average(2) = mean(xcoord1);
% stddev(2) = std(xcoord1);
% 
% %cells1 = 'DATA for 24hrs experiments Basic model/sigma15width120domaingrow.txt' % 24hrs, sigma =1.5
% cells1 = 'DATA for 18hrs experiments for Basic model/FDT18hrssigma15width120grows.txt' % 18hrs
% 
% cellpos1 = csvread(cells1);
% 
% xcoord1 = cellpos1(:,1);
% 
% average(3) = mean(xcoord1);
% stddev(3) = std(xcoord1);
% 
% % cells1 = 'DATA for 24hrs experiments Basic model/sigma20width120domaingrow.txt' % 24hrs, sigma =2.0
% cells1 = 'DATA for 18hrs experiments for Basic model/FDT18hrssigma20width120grows.txt' % 18hrs
% cellpos1 = csvread(cells1);
% 
% xcoord1 = cellpos1(:,1);
% 
% average(4) = mean(xcoord1);
% stddev(4) = std(xcoord1);
% 
% 
% sigma = [0.5, 1.0, 1.5, 2.0];
% 
% FDT = average
% 
% 
% hold on
% 
% p(2) = plot(sigma,FDT,'linewidth',4)
% %hbD = bar(sigma, FDT,'y'); % physical
% hold on
% % For each set of bars, find the centers of the bars, and write error bars
% pause(0.1); %pause allows the figure to be created
% 
% for ib = 1:numel(hbD)
%     %XData property is the tick labels/group centers; XOffset is the offset
%     %of each distinct group
%     xData = hbD(ib).XData+hbD(ib).XOffset;
%     errorbar(xData,FDT(ib,:),stddev(ib,:),'k.','linewidth',2)
% end
% 
% legend(p,'No growth','Domain grows');


%% sampling everywhere, and a cell simply does not move if it is out of bounds


cells1 = 'DATA for 18hrs experiments for Basic model/FDT18hrssigma05width120REJECT.txt'
cellpos1 = csvread(cells1);

xcoord1 = cellpos1(:,1);

average(1) = mean(xcoord1);
stddev(1) = std(xcoord1);


cells1 = 'DATA for 18hrs experiments for Basic model/FDT18hrssigma10width120REJECT.txt'
cellpos1 = csvread(cells1);

xcoord1 = cellpos1(:,1);

average(2) = mean(xcoord1);
stddev(2) = std(xcoord1);

cells1 = 'DATA for 18hrs experiments for Basic model/FDT18hrssigma15width120REJECT.txt'

cellpos1 = csvread(cells1);

xcoord1 = cellpos1(:,1);

average(3) = mean(xcoord1);
stddev(3) = std(xcoord1);

cells1 = 'DATA for 18hrs experiments for Basic model/FDT18hrssigma20width120REJECT.txt'
cellpos1 = csvread(cells1);

xcoord1 = cellpos1(:,1);

average(4) = mean(xcoord1);
stddev(4) = std(xcoord1);


sigma = [0.5, 1.0, 1.5, 2.0];

FDT = average


hold on

p(2) = plot(sigma,FDT,'r','linewidth',4)
%hbD = bar(sigma, FDT,'y'); % physical
hold on
% For each set of bars, find the centers of the bars, and write error bars
pause(0.1); %pause allows the figure to be created

for ib = 1:numel(hbD)
    %XData property is the tick labels/group centers; XOffset is the offset
    %of each distinct group
    xData = hbD(ib).XData+hbD(ib).XOffset;
    a =2
    errorbar(xData,FDT(ib,:),stddev(ib,:),'r.','linewidth',2)
end

legend(p,'Sample inside','Reject steps');
