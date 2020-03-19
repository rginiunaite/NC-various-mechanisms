% parameter sensitivity analysis

N = 5;
alldata = zeros(N,N);

% eps = 50
for i = 1:N
    count =i;
    if i==1
        k=1;
    else
        k = 3*(i-1);
    end
    i=k;
    %filename = sprintf('Parameter Sensitivity/CoACiLeps50D%itime.csv', i);
    filename = sprintf('Parameter Sensitivity/CiLOnlyeps1D%itime.csv', i);
    sepdata = load(filename);
    meanvalue = std(sepdata);
    alldata(count,1) = meanvalue; 

end

% eps = 100
for i = 1:N
    count =i;
    if i==1
        k=1;
    else
        k = 3*(i-1);
    end
    i=k;
   %  filename = sprintf('Parameter Sensitivity/CoACiLeps100D%itime.csv', i);
   filename = sprintf('Parameter Sensitivity/CiLOnlyeps50D%itime.csv', i);
    sepdata = load(filename);
    meanvalue = std(sepdata);
    alldata(count,2) = meanvalue; 

end

% eps = 150
for i = 1:N
    count =i;
    if i==1
        k=1;
    else
        k = 3*(i-1);
    end
    i=k;
    %filename = sprintf('Parameter Sensitivity/CoACiLeps150D%itime.csv', i);
    filename = sprintf('Parameter Sensitivity/CiLOnlyeps100D%itime.csv', i);
    sepdata = load(filename);
    meanvalue = std(sepdata);
    alldata(count,3) = meanvalue; 

end


% eps = 200
for i = 1:N
    count =i;
    if i==1
        k=1;
    else
        k = 3*(i-1);
    end
    i=k;
   %  filename = sprintf('Parameter Sensitivity/CoACiLeps200D%itime.csv', i);
   filename = sprintf('Parameter Sensitivity/CiLOnlyeps150D%itime.csv', i);
    sepdata = load(filename);
    meanvalue = std(sepdata);
    alldata(count,4) = meanvalue; 

end

% eps = 250
for i = 1:N
    count =i;
    if i==1
        k=1;
    else
        k = 3*(i-1);
    end
    i=k;
    %filename = sprintf('Parameter Sensitivity/CoACiLeps250D%itime.csv', i);
    filename = sprintf('Parameter Sensitivity/CiLOnlyeps200D%itime.csv', i);
    sepdata = load(filename);
    meanvalue = std(sepdata);
    alldata(count,5) = meanvalue; 

end

% % eps = 250
% for i = 1:N
%     count =i;
%     if i==1
%         k=1;
%     else
%         k = 3*(i-1);
%     end
%     i=k;
%     filename = sprintf('Parameter Sensitivity/CoACiLeps250D%itime.csv', i);
%     %filename = sprintf('Parameter Sensitivity/CiLOnlyeps250D%itime.csv', i);
%     sepdata = load(filename);
%     meanvalue = mean(sepdata);
%     alldata(count,6) = meanvalue; 
% 
% end

figure

alldata = alldata/60;
    %colormap('hot');   % set colormap
    imagesc(alldata);        % draw image and scale colormap to values range
    colorbar;          % show color scale
   % caxis([10, 50]);
 xticks([1,2,3,4,5, 6]);%,6])
 xticklabels({'1','50','100','150','200','250'});

 xlabel('\epsilon','FontSize',14)
 ylabel('D','FontSize',14)
  yticks([1,2,3,4,5, 6]);%,6])
 yticklabels({'1','3','6','9','12','15'});
 %ylim([1,6]);
  box on
 set(gca,'FontSize',36)
 ax = gca;
 set(gca,'linewidth',4) 
 
%  figure 
%  
%  plot (alldata(3,:))
 
    