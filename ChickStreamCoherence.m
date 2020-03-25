%% Chick model, we split the domain into sections of length 100\mu m 
%and check how many of them have less than 5 cells
clear all
N = 20;
threshold = 5;

% sim1 = 'ChickCiLOnlyeps1D1.csv';
% 
% M1 = csvread(sim1);
% M1 = M1(:,1)/20;

epsvalues = [1,50,100,150,200,250];
Dvalues = [1,12,9,6,3,1];

TotalMatrix = zeros(6,6);

alldata = zeros(11,N); % 11 because there are 11 sections, the domain is split 1100/10
logicalmatrix = zeros(11,N); % 11 because there are 11 sections, the domain is split 1100/10
logicalcol = zeros(11,1);

ascendMatrix = zeros(6,6); % this matrix will check whether the logic matrix decreasing, 
%if it is not decreasing it means that there are some gaps, we define this as a broken stream
for j =1:1
    D = Dvalues(j);
    for k = 1:6
        eps = epsvalues(k);
        for i = 1:N

           filename = sprintf('sepdataChicCiLOnlyeps%dD%dnvalue%d.csv',eps, D, i-1);
           sepdata = load(filename);
           alldata(:,i) = sepdata;
           logicalcol = alldata(:,i)>threshold;
            if  issorted(logicalcol,'descend') == false
                ascendMatrix(j,k) = ascendMatrix(j,k) +1;
               
            end
           logicalmatrix(:,i) = logicalcol;
           
            
        end
        TotalMatrix(j,k) = nnz(logicalmatrix);
    end
end

% x=[100:-5:0] %an example
% issorted(x,'strictdescend'issorted(str,2,'descend'))  %to check if they are sorted
% abs(diff(x))==5



TotalMatrix = TotalMatrix/220; % we use 220 because 20 simulations and 11 sections in each, we look how many sections do not have any


figure
  imagesc(TotalMatrix);        % draw image and scale colormap to values range
  colorbar;          % show color scale
 % caxis([10, 50]);
 xticks([1,2,3,4,5, 6]);%,6])
 xticklabels({'1','50','100','150','200','250'});

 xlabel('\epsilon','FontSize',14)
 ylabel('D','FontSize',14)
  yticks([1,2,3,4,5, 6]);%,6])
yticklabels({'15','12','9','6','3','1'});
 %ylim([1,6]);
  box on
 set(gca,'FontSize',36)
 ax = gca;
 set(gca,'linewidth',4) 
 
figure
  imagesc(ascendMatrix);        % draw image and scale colormap to values range
  colorbar;          % show color scale
 % caxis([10, 50]);
 xticks([1,2,3,4,5, 6]);%,6])
 xticklabels({'1','50','100','150','200','250'});

 xlabel('\epsilon','FontSize',14)
 ylabel('D','FontSize',14)
  yticks([1,2,3,4,5, 6]);%,6])
yticklabels({'15','12','9','6','3','1'});
 %ylim([1,6]);
  box on
 set(gca,'FontSize',36)
 ax = gca;
 set(gca,'linewidth',4) 