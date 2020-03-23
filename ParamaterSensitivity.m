% parameter sensitivity analysis

N = 5;
alldata = zeros(N,N);

%% Vary D and eps


% % eps = 1
% for i = 1:N
%     count =i;
%     if i==6
%         k= 1;
%     else
%         k = 15- 3*(i-1);
%     end
%     i=k;
%     filename = sprintf('Parameter Sensitivity/CoACiLeps1D%itime.csv', i);
%     %filename = sprintf('Parameter Sensitivity/CiLOnlyeps1D%itime.csv', i);
%     sepdata = load(filename);
%     stdvalue = mean(sepdata);
%     alldata(count,1) = stdvalue; 
% 
% end
% 
% % eps = 50
% for i = 1:N
%     count =i;
%     if i==6
%         k= 1;
%     else
%         k = 15- 3*(i-1);
%     end
%     i=k;
%     filename = sprintf('Parameter Sensitivity/CoACiLeps50D%itime.csv', i);
%    %filename = sprintf('Parameter Sensitivity/CiLOnlyeps50D%itime.csv', i);
%     sepdata = load(filename);
%     stdvalue = mean(sepdata);
%     alldata(count,2) = stdvalue; 
% 
% end
% 
% % eps = 100
% for i = 1:N
%     count =i;
%     if i==6
%         k= 1;
%     else
%         k = 15- 3*(i-1);
%     end
%     i=k;
%     filename = sprintf('Parameter Sensitivity/CoACiLeps100D%itime.csv', i);
%     %filename = sprintf('Parameter Sensitivity/CiLOnlyeps100D%itime.csv', i);
%     sepdata = load(filename);
%     stdvalue = mean(sepdata);
%     alldata(count,3) = stdvalue; 
% 
% end
% 
% 
% % eps = 150
% for i = 1:N
%     count =i;
%     if i==6
%         k= 1;
%     else
%         k = 15- 3*(i-1);
%     end
%     i=k;
%     filename = sprintf('Parameter Sensitivity/CoACiLeps150D%itime.csv', i);
%   %filename = sprintf('Parameter Sensitivity/CiLOnlyeps150D%itime.csv', i);
%     sepdata = load(filename);
%     stdvalue = mean(sepdata);
%     alldata(count,4) = stdvalue; 
% 
% end
% 
% % eps = 200
% for i = 1:N
%     count =i;
%     if i==6
%         k= 1;
%     else
%         k = 15- 3*(i-1);
%     end
%     i=k;
%     filename = sprintf('Parameter Sensitivity/CoACiLeps200D%itime.csv', i);
%     %filename = sprintf('Parameter Sensitivity/CiLOnlyeps200D%itime.csv', i);
%     sepdata = load(filename);
%     stdvalue = mean(sepdata);
%     alldata(count,5) = stdvalue; 
% 
% end
% 
% % eps = 250
% for i = 1:N
%     count =i;
%     if i==6
%         k= 1;
%     else
%         k = 15- 3*(i-1);
%     end
%     i=k;
%     filename = sprintf('Parameter Sensitivity/CoACiLeps250D%itime.csv', i);
%     %filename = sprintf('Parameter Sensitivity/CiLOnlyeps250D%itime.csv', i);
%     sepdata = load(filename);
%     stdvalue = mean(sepdata);
%     alldata(count,6) = stdvalue; 
% 
% end
% 
% figure
% 
% alldata = alldata/60;
%     %colormap('hot');   % set colormap
%     imagesc(alldata);        % draw image and scale colormap to values range
%     colorbar;          % show color scale
%    % caxis([10, 50]);
%  xticks([1,2,3,4,5, 6]);%,6])
%  xticklabels({'1','50','100','150','200','250'});
% 
%  xlabel('\epsilon','FontSize',14)
%  ylabel('D','FontSize',14)
%   yticks([1,2,3,4,5, 6]);%,6])
% yticklabels({'15','12','9','6','3','1'});
%  %ylim([1,6]);
%   box on
%  set(gca,'FontSize',36)
%  ax = gca;
%  set(gca,'linewidth',4) 
%  
% %  figure 
% %  
% %  plot (alldata(3,:))
%  
%     






%% Vary beta and


% eps = 1
for i = 1:N
    count =i;
    if i==6
        k= 1;
    else
        k = 6-i;
    end
    i=k;
    filename = sprintf('Parameter Sensitivity/CoACiLeps1beta0p0%i.csv', i);
    %filename = sprintf('Parameter Sensitivity/CiLOnlyeps1D%itime.csv', i);
    sepdata = load(filename);
    stdvalue = mean(sepdata);
    alldata(count,1) = stdvalue; 

end

% eps = 50
for i = 1:N
    count =i;
    if i==6
        k= 1;
    else
        k = 6-i;
    end
    i=k;
      filename = sprintf('Parameter Sensitivity/CoACiLeps50beta0p0%i.csv', i);
   %filename = sprintf('Parameter Sensitivity/CiLOnlyeps50D%itime.csv', i);
    sepdata = load(filename);
    stdvalue = mean(sepdata);
    alldata(count,2) = stdvalue; 

end

% eps = 100
for i = 1:N
    count =i;
    if i==6
        k= 1;
    else
        k = 6-i;
    end
    i=k;
       filename = sprintf('Parameter Sensitivity/CoACiLeps100beta0p0%i.csv', i);
    %filename = sprintf('Parameter Sensitivity/CiLOnlyeps100D%itime.csv', i);
    sepdata = load(filename);
    stdvalue = mean(sepdata);
    alldata(count,3) = stdvalue; 

end


% eps = 150
for i = 1:N
    count =i;
    if i==6
        k= 1;
    else
        k = 6-i;
    end
    i=k;
      filename = sprintf('Parameter Sensitivity/CoACiLeps150beta0p0%i.csv', i);
  %filename = sprintf('Parameter Sensitivity/CiLOnlyeps150D%itime.csv', i);
    sepdata = load(filename);
    stdvalue = mean(sepdata);
    alldata(count,4) = stdvalue; 

end

% eps = 200
for i = 1:N
    count =i;
    if i==6
        k= 1;
    else
          k = 6-i;
    end
    i=k;
      filename = sprintf('Parameter Sensitivity/CoACiLeps200beta0p0%i.csv', i);
    %filename = sprintf('Parameter Sensitivity/CiLOnlyeps200D%itime.csv', i);
    sepdata = load(filename);
    stdvalue = mean(sepdata);
    alldata(count,5) = stdvalue; 

end

% eps = 250
for i = 1:N
    count =i;
    if i==6
        k= 1;
    else
        k = 6-i;
    end
    i=k;
      filename = sprintf('Parameter Sensitivity/CoACiLeps250beta0p0%i.csv', i);
    %filename = sprintf('Parameter Sensitivity/CiLOnlyeps250D%itime.csv', i);
    sepdata = load(filename);
    stdvalue = mean(sepdata);
    alldata(count,6) = stdvalue; 

end

figure

alldata = alldata/60;
    %colormap('hot');   % set colormap
    imagesc(alldata);        % draw image and scale colormap to values range
    colorbar;          % show color scale
   % caxis([10, 50]);
 xticks([1,2,3,4,5, 6]);%,6])
 xticklabels({'1','50','100','150','200','250'});

 xlabel('\epsilon','FontSize',14)
 ylabel('\beta','FontSize',14)
  yticks([1,2,3,4,5, 6]);%,6])
yticklabels({'0.05','0.04','0.03','0.02','0.01','1'});
 %ylim([1,6]);
  box on
 set(gca,'FontSize',36)
 ax = gca;
 set(gca,'linewidth',4) 
 
%  figure 
%  
%  plot (alldata(3,:))
 
    
