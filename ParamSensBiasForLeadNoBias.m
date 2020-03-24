%% Bias for leaders and CoA CiL with no bias


%% Vary D


filename = sprintf('Parameter Sensitivity/LeadOnlyvaryDandEpsMATRIX.csv');
matrix = load(filename);
figure
  imagesc(matrix);        % draw image and scale colormap to values range
  colorbar;          % show color scale
 % caxis([10, 50]);
 xticks([1,2,3,4,5, 6]);%,6])
 xticklabels({'1','50','100','150','200','250'});

 
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
 
%  figure 
%  
%  plot (alldata(3,:))
 
 
 
 
 
 %%  Vary beta
 
 
 
 
 
%  
%  
%  
%  xlabel('\epsilon','FontSize',14)
%  ylabel('\beta','FontSize',14)
%   yticks([1,2,3,4,5, 6]);%,6])
% yticklabels({'0.05','0.04','0.03','0.02','0.01','1'});
%  %ylim([1,6]);
%   box on
%  set(gca,'FontSize',36)
%  ax = gca;
%  set(gca,'linewidth',4) 