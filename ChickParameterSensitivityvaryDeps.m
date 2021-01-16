% Chick parameter sensitivity analysis
% Vary D and epsilon
%beta =0.7

model = 0; % 0 - AttrRepALLBiasedeps, 1 - RepOnlyALLBiased, 2 - AttrRepBiasedLeaders,  3 AttrRep, 4 Rep Only
final = 1200; % 1440 if 24 hours, 1080-if 18h


N = 5;
alldata = zeros(N,N);
alldatastd = zeros(N,N);
nanmatrix = zeros(N,N);

Dvalues  = [13,10,7,4,1];
betavalues = [13,10,7,4,1];
epsvalues = [0,19,38,56,75];


% %% Vary beta and eps
% 
% 
% eps = 0.4
for count = 1:N
    i = Dvalues(count);
    % old
    %filename = sprintf('Parameter Sensitivity/CoACiLeps1beta0p0%i.csv', i);
    %filename = sprintf('Parameter Sensitivity/CiLOnlyeps1D5beta0p0%i.csv', i);
    % new
    %filename = sprintf('NEW DATA CHICK/Attrrepeps1D6beta0p%i.csv', i);
    %filename = sprintf('NEW DATA CHICK/RepOnlyeps1D6beta0p%i.csv', i);
    
    
     
    if model == 0
    filename = sprintf('CHICK DATA FINAL3/AttrRepALLBiasedeps0D%ibeta0p4.csv', i);
    end
   if model == 1
     filename = sprintf('CHICK DATA FINAL3/RepOnlyALLBiasedeps0D%ibeta0p4.csv', i);
   end 
        
    if model == 2
    filename = sprintf('CHICK DATA FINAL3/AttrRepBiasedLeaderseps0D%ibeta0p4.csv', i);
    end
   if model == 3
     filename = sprintf('CHICK DATA FINAL3/AttrRepeps0D%ibeta0p0.csv', i);
   end
   if model == 4
     filename = sprintf('CHICK DATA FINAL3/RepOnlyeps0D%ibeta0p0.csv', i);
   end
    sepdata = load(filename);  
     indices = find(abs(sepdata)>final); % only values smaller than final, which is 18h
    sepdata(indices) = NaN;
    
    % count the number of nan
    nanmatrix(count,1) = sum(isnan(sepdata));
    
    meanvalue = nanmean(sepdata);  
    alldata(count,1) = meanvalue; 
    stdvalue = nanstd(sepdata);
    alldatastd(count,1) = stdvalue; 
end

% 
% eps = 19
for count = 1:N
    i = Dvalues(count);
    %  filename = sprintf('Parameter Sensitivity/CoACiLeps50beta0p0%i.csv', i);
%   filename = sprintf('Parameter Sensitivity/CiLOnlyeps50D5beta0p0%i.csv', i);

    % new
    %filename = sprintf('NEW DATA CHICK/Attrrepeps50D6beta0p%i.csv', i);
    %filename = sprintf('NEW DATA CHICK/RepOnlyeps50D6beta0p%i.csv', i);
   % filename = sprintf('NEW DATA CHICK/AttrRepOnlyBiasedLeaderseps50D6beta0p%i.csv', i);

     
    if model == 0
    filename = sprintf('CHICK DATA FINAL3/AttrRepALLBiasedeps19D%ibeta0p4.csv', i);
    end
   if model == 1
     filename = sprintf('CHICK DATA FINAL3/RepOnlyALLBiasedeps19D%ibeta0p4.csv', i);
   end 
        
    if model == 2
    filename = sprintf('CHICK DATA FINAL3/AttrRepBiasedLeaderseps19D%ibeta0p4.csv', i);
    end
   if model == 3
     filename = sprintf('CHICK DATA FINAL3/AttrRepeps19D%ibeta0p0.csv', i);
   end
      if model == 4
     filename = sprintf('CHICK DATA FINAL3/RepOnlyeps19D%ibeta0p0.csv', i);
   end
   
    sepdata = load(filename);  
     indices = find(abs(sepdata)>final); % only values smaller than 18h
    sepdata(indices) = NaN;
    
    % count the number of nan
    nanmatrix(count,2) = sum(isnan(sepdata));
    
    meanvalue = nanmean(sepdata);  
    alldata(count,2) = meanvalue; 
    stdvalue = nanstd(sepdata);
    alldatastd(count,2) = stdvalue; 

end

% eps = 38
for count = 1:N
    i = Dvalues(count);
      % filename = sprintf('Parameter Sensitivity/CoACiLeps100beta0p0%i.csv', i);
%    filename = sprintf('Parameter Sensitivity/CiLOnlyeps100D5beta0p0%i.csv', i);

    % new
   % filename = sprintf('NEW DATA CHICK/Attrrepeps100D6beta0p%i.csv', i);
   % filename = sprintf('NEW DATA CHICK/RepOnlyeps100D6beta0p%i.csv', i);
   %   filename = sprintf('NEW DATA CHICK/AttrRepOnlyBiasedLeaderseps100D6beta0p%i.csv', i);
   
    if model == 0
    filename = sprintf('CHICK DATA FINAL3/AttrRepALLBiasedeps38D%ibeta0p4.csv', i);
    end
   if model == 1
     filename = sprintf('CHICK DATA FINAL3/RepOnlyALLBiasedeps38D%ibeta0p4.csv', i);
   end 
        
    if model == 2
    filename = sprintf('CHICK DATA FINAL3/AttrRepBiasedLeaderseps38D%ibeta0p4.csv', i);
    end
   if model == 3
     filename = sprintf('CHICK DATA FINAL3/AttrRepeps38D%ibeta0p0.csv', i);
   end
      if model == 4
     filename = sprintf('CHICK DATA FINAL3/RepOnlyeps38D%ibeta0p0.csv', i);
   end
    sepdata = load(filename);  
     indices = find(abs(sepdata)>final); % only values smaller than 18h
    sepdata(indices) = NaN;
    
    % count the number of nan
    nanmatrix(count,3) = sum(isnan(sepdata));
    
    meanvalue = nanmean(sepdata);  
    alldata(count,3) = meanvalue; 
    stdvalue = nanstd(sepdata);
    alldatastd(count,3) = stdvalue; 

end


% eps = 56
for count = 1:N
    i = Dvalues(count);
     % filename = sprintf('Parameter Sensitivity/CoACiLeps150beta0p0%i.csv', i);
%  filename = sprintf('Parameter Sensitivity/CiLOnlyeps150D5beta0p0%i.csv', i);
  
      % new
   % filename = sprintf('NEW DATA CHICK/Attrrepeps150D6beta0p%i.csv', i);
    %  filename = sprintf('NEW DATA CHICK/RepOnlyeps150D6beta0p%i.csv', i);
    %filename = sprintf('NEW DATA CHICK/AttrRepOnlyBiasedLeaderseps150D6beta0p%i.csv', i);
     
        if model == 0
    filename = sprintf('CHICK DATA FINAL3/AttrRepALLBiasedeps56D%ibeta0p4.csv', i);
    end
   if model == 1
     filename = sprintf('CHICK DATA FINAL3/RepOnlyALLBiasedeps56D%ibeta0p4.csv', i);
   end 
        
    if model == 2
    filename = sprintf('CHICK DATA FINAL3/AttrRepBiasedLeaderseps56D%ibeta0p4.csv', i);
    end
   if model == 3
     filename = sprintf('CHICK DATA FINAL3/AttrRepeps56D%ibeta0p0.csv', i);
   end
     if model == 4
     filename = sprintf('CHICK DATA FINAL3/RepOnlyeps56D%ibeta0p0.csv', i);
   end  
      sepdata = load(filename);  
     indices = find(abs(sepdata)>final); % only values smaller than 18h
    sepdata(indices) = NaN;
    
    % count the number of nan
    nanmatrix(count,4) = sum(isnan(sepdata));
    
    meanvalue = nanmean(sepdata);  
    alldata(count,4) = meanvalue; 
    stdvalue = nanstd(sepdata);
    alldatastd(count,4) = stdvalue; 

end

% eps = 75
for count = 1:N
    i = Dvalues(count);
    % filename = sprintf('Parameter Sensitivity/CoACiLeps200beta0p0%i.csv', i);
   %  filename = sprintf('Parameter Sensitivity/CiLOnlyeps200D5beta0p0%i.csv', i);
     
         % new
    %filename = sprintf('NEW DATA CHICK/Attrrepeps200D6beta0p%i.csv', i);
     %    filename = sprintf('NEW DATA CHICK/RepOnlyeps200D6beta0p%i.csv', i);
     % filename = sprintf('NEW DATA CHICK/AttrRepOnlyBiasedLeaderseps200D6beta0p%i.csv', i);

     if model == 0
    filename = sprintf('CHICK DATA FINAL3/AttrRepALLBiasedeps75D%ibeta0p4.csv', i);
    end
   if model == 1
     filename = sprintf('CHICK DATA FINAL3/RepOnlyALLBiasedeps75D%ibeta0p4.csv', i);
   end 
        
    if model == 2
    filename = sprintf('CHICK DATA FINAL3/AttrRepBiasedLeaderseps75D%ibeta0p4.csv', i);
    end
   if model == 3
     filename = sprintf('CHICK DATA FINAL3/AttrRepeps75D%ibeta0p0.csv', i);
   end
     if model == 4
     filename = sprintf('CHICK DATA FINAL3/RepOnlyeps75D%ibeta0p0.csv', i);
   end
       
    sepdata = load(filename);  
     indices = find(abs(sepdata)>final); % only values smaller than 18h
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
%     %filename = sprintf('NEW DATA CHICK/Attrrepeps250D6beta0p%i.csv', i);
%      % filename = sprintf('NEW DATA CHICK/RepOnlyeps250D6beta0p%i.csv', i);
%       filename = sprintf('NEW DATA CHICK/AttrRepOnlyBiasedLeaderseps250D6beta0p%i.csv', i);
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

    C =   colorbar;          % show color scale
%     L=cellfun(@(x)sprintf('%.1f',x),num2cell(get(C,'xtick')),'Un',0);
% set(C,'xticklabel',L)
   
    oldcmap = colormap;
    colormap(flipud(oldcmap))
    
   % caxis([10, 50]);
 xticks([1,2,3,4,5, 6]);%,6])
 xticklabels({'0.4','19.0','38.0','56.0','75.0'});

 xlabel('\epsilon','FontSize',14)
 ylabel('D','FontSize',14)
  yticks([1,2,3,4,5, 6]);%,6])
%yticklabels({'0.05','0.04','0.03','0.02','0.01','1'});
yticklabels({'13','10','7','4','1'});

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
 xticklabels({'0.4','19.0','38.0','56.0','75.0'});


 xlabel('\epsilon','FontSize',14)
 ylabel('D','FontSize',14)
  yticks([1,2,3,4,5, 6]);%,6])
%yticklabels({'0.05','0.04','0.03','0.02','0.01','1'});
yticklabels({'13','10','7','4','1'});

 %ylim([1,6]);
  box on
 set(gca,'FontSize',36)
 ax = gca;
 set(gca,'linewidth',4) 