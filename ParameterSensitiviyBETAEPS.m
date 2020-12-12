% parameter sensitivity analysis
% Vary beta and epsilon
% D = 7

model = 2; % 0 - AttrRepALLBIASEd, 1 - RepOnlyALLBIASED, 2 - AttrRepBiasedLeaders, not in this on: 3 AttrRep


N = 5;
alldata = zeros(N,N);
alldatastd = zeros(N,N);
nanmatrix = zeros(N,N);

Dvalues  = [13,10,7,4,1];
betavalues = [13,10,7,4,1];
epsvalues = [200,150,100,50,1];


% %% Vary beta and eps
% 
% 
% eps = 1
for count = 1:N
    i = betavalues(count);
    % old
    %filename = sprintf('Parameter Sensitivity/CoACiLeps1beta0p0%i.csv', i);
    %filename = sprintf('Parameter Sensitivity/CiLOnlyeps1D5beta0p0%i.csv', i);
    % new
    %filename = sprintf('NEW DATA XENOPUS/Attrrepeps1D6beta0p%i.csv', i);
    %filename = sprintf('NEW DATA XENOPUS/RepOnlyeps1D6beta0p%i.csv', i);
    
    
     
    if model == 0
    filename = sprintf('XENOPUS DATA FINAL/AttrRepALLBIASEDeps1D7beta0p%i.csv', i);
    end
   if model == 1
     filename = sprintf('XENOPUS DATA FINAL/RepOnlyALLBIASEDeps1D7beta0p%i.csv', i);
   end 
        
    if model == 2
    filename = sprintf('XENOPUS DATA FINAL/AttrRepBIASEDLEADERSeps1D7beta0p%i.csv', i);
    end
%    if model == 3
%      filename = sprintf('XENOPUS DATA FINAL/AttrRepeps150D7beta0p%i.csv', i);
%    end

    sepdata = load(filename);  
     indices = find(abs(sepdata)>2999); % only values smaller than 3000, which is 50h
    sepdata(indices) = NaN;
    
    % count the number of nan
    nanmatrix(count,1) = sum(isnan(sepdata));
    
    meanvalue = nanmean(sepdata);  
    alldata(count,1) = meanvalue; 
    stdvalue = nanstd(sepdata);
    alldatastd(count,1) = stdvalue; 
end

% 
% eps = 50
for count = 1:N
    i = betavalues(count);
    %  filename = sprintf('Parameter Sensitivity/CoACiLeps50beta0p0%i.csv', i);
%   filename = sprintf('Parameter Sensitivity/CiLOnlyeps50D5beta0p0%i.csv', i);

    % new
    %filename = sprintf('NEW DATA XENOPUS/Attrrepeps50D6beta0p%i.csv', i);
    %filename = sprintf('NEW DATA XENOPUS/RepOnlyeps50D6beta0p%i.csv', i);
   % filename = sprintf('NEW DATA XENOPUS/AttrRepOnlyBiasedLeaderseps50D6beta0p%i.csv', i);

     
    if model == 0
    filename = sprintf('XENOPUS DATA FINAL/AttrRepALLBIASEDeps50D7beta0p%i.csv', i);
    end
   if model == 1
     filename = sprintf('XENOPUS DATA FINAL/RepOnlyALLBIASEDeps50D7beta0p%i.csv', i);
   end
    if model == 2
    filename = sprintf('XENOPUS DATA FINAL/AttrRepBIASEDLEADERSeps50D7beta0p%i.csv', i);
    end
    sepdata = load(filename);  
     indices = find(abs(sepdata)>2999); % only values smaller than 3000, which is 50h
    sepdata(indices) = NaN;
    
    % count the number of nan
    nanmatrix(count,2) = sum(isnan(sepdata));
    
    meanvalue = nanmean(sepdata);  
    alldata(count,2) = meanvalue; 
    stdvalue = nanstd(sepdata);
    alldatastd(count,2) = stdvalue; 

end

% eps = 100
for count = 1:N
    i = betavalues(count);
      % filename = sprintf('Parameter Sensitivity/CoACiLeps100beta0p0%i.csv', i);
%    filename = sprintf('Parameter Sensitivity/CiLOnlyeps100D5beta0p0%i.csv', i);

    % new
   % filename = sprintf('NEW DATA XENOPUS/Attrrepeps100D6beta0p%i.csv', i);
   % filename = sprintf('NEW DATA XENOPUS/RepOnlyeps100D6beta0p%i.csv', i);
   %   filename = sprintf('NEW DATA XENOPUS/AttrRepOnlyBiasedLeaderseps100D6beta0p%i.csv', i);
   
     
    if model == 0
    filename = sprintf('XENOPUS DATA FINAL/AttrRepALLBIASEDeps100D7beta0p%i.csv', i);
    end
   if model == 1
     filename = sprintf('XENOPUS DATA FINAL/RepOnlyALLBIASEDeps100D7beta0p%i.csv', i);
   end
    if model == 2
    filename = sprintf('XENOPUS DATA FINAL/AttrRepBIASEDLEADERSeps100D7beta0p%i.csv', i);
    end
    sepdata = load(filename);  
     indices = find(abs(sepdata)>2999); % only values smaller than 3000, which is 50h
    sepdata(indices) = NaN;
    
    % count the number of nan
    nanmatrix(count,3) = sum(isnan(sepdata));
    
    meanvalue = nanmean(sepdata);  
    alldata(count,3) = meanvalue; 
    stdvalue = nanstd(sepdata);
    alldatastd(count,3) = stdvalue; 

end


% eps = 150
for count = 1:N
    i = betavalues(count);
     % filename = sprintf('Parameter Sensitivity/CoACiLeps150beta0p0%i.csv', i);
%  filename = sprintf('Parameter Sensitivity/CiLOnlyeps150D5beta0p0%i.csv', i);
  
      % new
   % filename = sprintf('NEW DATA XENOPUS/Attrrepeps150D6beta0p%i.csv', i);
    %  filename = sprintf('NEW DATA XENOPUS/RepOnlyeps150D6beta0p%i.csv', i);
    %filename = sprintf('NEW DATA XENOPUS/AttrRepOnlyBiasedLeaderseps150D6beta0p%i.csv', i);
     
    if model == 0
    filename = sprintf('XENOPUS DATA FINAL/AttrRepALLBIASEDeps150D7beta0p%i.csv', i);
    end
   if model == 1
     filename = sprintf('XENOPUS DATA FINAL/RepOnlyALLBIASEDeps150D7beta0p%i.csv', i);
   end
        if model == 2
    filename = sprintf('XENOPUS DATA FINAL/AttrRepBIASEDLEADERSeps150D7beta0p%i.csv', i);
    end  
      sepdata = load(filename);  
     indices = find(abs(sepdata)>2999); % only values smaller than 3000, which is 50h
    sepdata(indices) = NaN;
    
    % count the number of nan
    nanmatrix(count,4) = sum(isnan(sepdata));
    
    meanvalue = nanmean(sepdata);  
    alldata(count,4) = meanvalue; 
    stdvalue = nanstd(sepdata);
    alldatastd(count,4) = stdvalue; 

end

% eps = 200
for count = 1:N
    i = betavalues(count);
    % filename = sprintf('Parameter Sensitivity/CoACiLeps200beta0p0%i.csv', i);
   %  filename = sprintf('Parameter Sensitivity/CiLOnlyeps200D5beta0p0%i.csv', i);
     
         % new
    %filename = sprintf('NEW DATA XENOPUS/Attrrepeps200D6beta0p%i.csv', i);
     %    filename = sprintf('NEW DATA XENOPUS/RepOnlyeps200D6beta0p%i.csv', i);
     % filename = sprintf('NEW DATA XENOPUS/AttrRepOnlyBiasedLeaderseps200D6beta0p%i.csv', i);

     
    if model == 0
    filename = sprintf('XENOPUS DATA FINAL/AttrRepALLBIASEDeps200D7beta0p%i.csv', i);
    end
   if model == 1
     filename = sprintf('XENOPUS DATA FINAL/RepOnlyALLBIASEDeps200D7beta0p%i.csv', i);
   end
         if model == 2
    filename = sprintf('XENOPUS DATA FINAL/AttrRepBIASEDLEADERSeps200D7beta0p%i.csv', i);
    end 
    sepdata = load(filename);  
     indices = find(abs(sepdata)>2999); % only values smaller than 3000, which is 50h
    sepdata(indices) = NaN;
    
    % count the number of nan
    nanmatrix(count,5) = sum(isnan(sepdata));
    
    meanvalue = nanmean(sepdata);  
    alldata(count,5) = meanvalue; 
    stdvalue = nanstd(sepdata);
    alldatastd(count,5) = stdvalue; 
end

% % eps = 250
% for count = 1:N
%     i = betavalues(count);
%     % filename = sprintf('Parameter Sensitivity/CoACiLeps250beta0p0%i.csv', i);
%     % filename = sprintf('Parameter Sensitivity/CiLOnlyeps250D5beta0p0%i.csv', i);
%      
%     % new
%     %filename = sprintf('NEW DATA XENOPUS/Attrrepeps250D6beta0p%i.csv', i);
%      % filename = sprintf('NEW DATA XENOPUS/RepOnlyeps250D6beta0p%i.csv', i);
%       filename = sprintf('NEW DATA XENOPUS/AttrRepOnlyBiasedLeaderseps250D6beta0p%i.csv', i);
% 
%     sepdata = load(filename);  
%      indices = find(abs(sepdata)>2999); % only values smaller than 3000, which is 50h
%     sepdata(indices) = NaN;
%     
%     % count the number of nan
%     nanmatrix(count,6) = sum(isnan(sepdata));
%     
%     meanvalue = nanmean(sepdata);  
%     alldata(count,6) = meanvalue; 
%     stdvalue = nanstd(sepdata);
%     alldatastd(count,6) = stdvalue; 
% end

figure

alldata = alldata/60;
alldatastd = alldatastd/60;

%    
% indices = find(abs(alldata)>50);
% alldata(indices) = NaN;
% alldatastd(indices) = NaN;

    %colormap('hot');   % set colormap
   b =  imagesc(alldata);        % draw image and scale colormap to values range
   set(b,'AlphaData',~isnan(alldata)) % set to white NaN values
    colorbar;          % show color scale
    oldcmap = colormap;
    colormap(flipud(oldcmap))
    
   % caxis([10, 50]);
 xticks([1,2,3,4,5, 6]);%,6])
 xticklabels({'1','50','100','150','200'});

 xlabel('\epsilon','FontSize',14)
 ylabel('\beta','FontSize',14)
  yticks([1,2,3,4,5, 6]);%,6])
%yticklabels({'0.05','0.04','0.03','0.02','0.01','1'});
yticklabels({'1.3','1.0','0.7','0.4','0.1'});

box on
 set(gca,'FontSize',36)
 ax = gca;
 set(gca,'linewidth',4) 
 
 % plots std
  figure
% 
    %colormap('hot');   % set colormap
   s= imagesc(alldatastd);        % draw image and scale colormap to values range
    set(s,'AlphaData',~isnan(alldatastd)) % set to white NaN values

 C =   colorbar;          % show color scale
    L=cellfun(@(x)sprintf('%.1f',x),num2cell(get(C,'xtick')),'Un',0);
set(C,'xticklabel',L)
    
    
       oldcmap = colormap;
    colormap(flipud(oldcmap))
    
    
    
 %caxis([10, 70]); % for mean
 %caxis([0, 17]); % for std
 xticks([1,2,3,4,5, 6]);%,6])
 xticklabels({'1','50','100','150','200'});


 xlabel('\epsilon','FontSize',14)
 ylabel('\beta','FontSize',14)
  yticks([1,2,3,4,5, 6]);%,6])
%yticklabels({'0.05','0.04','0.03','0.02','0.01','1'});
yticklabels({'1.3','1.0','0.7','0.4','0.1'});

 %ylim([1,6]);
  box on
 set(gca,'FontSize',36)
 ax = gca;
 set(gca,'linewidth',4) 
