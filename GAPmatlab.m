% gap statistic in matlab

% load fisheriris;
% rng('default');  % For reproducibility
% eva = evalclusters(meas,'kmeans','gap','KList',[1:6])
% figure;
% plot(eva);


%cellpos1 = 'positionsBreak.csv';
cellpos1 = 'positionsVerycont.csv';

cellpos1 = csvread(cellpos1);

cellposx = cellpos1(:,2);
cellposy = cellpos1(:,3);

cellposfinal(:,1) = cellposx;
cellposfinal(:,2) = cellposy;

eva = evalclusters(cellposfinal,'kmeans','gap','KList',[1:6])

figure;
plot(eva);

eva.StdLogW;

standarderrpr = eva.SE


figure
PetalLength = cellposx;
PetalWidth = cellposy;
ClusterGroup = eva.OptimalY;
gscatter(PetalLength,PetalWidth,ClusterGroup,'rbgkc','xod^*');