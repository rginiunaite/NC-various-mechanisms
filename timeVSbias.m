%% plot time to invade versus bias, for Xenopus model

% xCoACiL = [2,5,8,11,14]; 
% xCiLonly =[1,4,7,10,13];
animal = 1; % 0-Xenopus, 1-chick
xCoACiL = [3,7,11]; 
xCiLonly =[2,6,10];
xVolExc = [1,5,9]

total=1080;

%% first CoACiL

%sim1 = 'TOSAVETimeCoACiL0p01betaNarrowDomain.csv';
%sim1 = 'TOSAVEXenopusonlybias0p01D1.csv';
%sim1 = 'BIAS data/AttractionRepulsionD5beta0p01.csv';

if animal ==0
    sim1 = 'XENOPUS DATA FINAL2/AttrRepALLBIASEDto850eps200D7beta0p4.csv';
else
    sim1 = 'CHICK DATA FINAL/AttrRepALLBiasedeps75D7beta0p4.csv';
end


M1 = csvread(sim1);
M1 = M1(:,1)/60;

averageCoACiL(1) = mean(M1)
stdsCoACiL(1) =std(M1)

%% first CiLonly

%sim1 = 'TOSAVETimeCiLonly0p01betaNarrowDomain.csv';
% sim1 = 'BIAS data/RepulsionD5beta0p01.csv';

if animal ==0
    sim1 = 'XENOPUS DATA FINAL2/RepOnlyALLBIASEDto850eps200D7beta0p4.csv';
else
    sim1 = 'CHICK DATA FINAL/RepOnlyALLBiasedeps75D7beta0p4.csv';
end


M1 = csvread(sim1);
M1 = M1(:,1)/60;

averageCiLonly(1) = mean(M1)
stdsCiLonly(1) =std(M1)

%% first VolExc

%sim1 = 'TOSAVETimeCiLonly0p01betaNarrowDomain.csv';
%sim1 = 'BIAS data/VolumeExclusionD5beta0p01.csv';

if animal ==0
    sim1 = 'XENOPUS DATA FINAL2/RepOnlyALLBIASEDto850eps1D7beta0p4.csv';
else
    sim1 = 'CHICK DATA FINAL/RepOnlyALLBiasedeps0D7beta0p4.csv';
end

M1 = csvread(sim1);
M1 = M1(:,1)/60;


averageVolExc(1) = mean(M1)
stdsVolExc(1) =std(M1)


%% second CoA CiL
%sim1 = 'TOSAVETimeCoACiL0p02betaNarrowDomain.csv';
%sim1 = 'TOSAVEXenopusonlybias0p02D1.csv';
%sim1 = 'BIAS data/AttractionRepulsionD5beta0p03.csv';

if animal ==0
    sim1 = 'XENOPUS DATA FINAL2/AttrRepALLBIASEDto850eps200D7beta0p7.csv';
else
    sim1 = 'CHICK DATA FINAL/AttrRepALLBiasedeps75D7beta0p7.csv';
end

M1 = csvread(sim1);
M1 = M1(:,1)/60;

averageCoACiL(2) = mean(M1)
stdsCoACiL(2) =std(M1)

%% second CiL only
%sim1 = 'TOSAVETimeCiLonly0p02betaNarrowDomain.csv';
%sim1 = 'BIAS data/RepulsionD5beta0p03.csv';

if animal == 0
    sim1 = 'XENOPUS DATA FINAL2/RepOnlyALLBIASEDto850eps200D7beta0p7.csv';
else
    sim1 = 'CHICK DATA FINAL/RepOnlyALLBiasedeps75D7beta0p7.csv';
end
M1 = csvread(sim1);
M1 = M1(:,1)/60;

averageCiLonly(2) = mean(M1)
stdsCiLonly(2) =std(M1)


%% second VolExc
%sim1 = 'TOSAVETimeCiLonly0p02betaNarrowDomain.csv';
%sim1 = 'BIAS data/VolumeExclusionD5beta0p03.csv';

if animal == 0
    sim1 = 'XENOPUS DATA FINAL2/RepOnlyALLBIASEDto850eps1D7beta0p7.csv';
else
    sim1 = 'CHICK DATA FINAL/RepOnlyALLBiasedeps0D7beta0p7.csv';
end
M1 = csvread(sim1);
M1 = M1(:,1)/60;

averageVolExc(2) = mean(M1)
stdsVolExc(2) =std(M1)

%% third CoA CiL

%sim1 = 'TOSAVETimeCoACiL0p03betaNarrowDomain.csv';
%sim1 = 'TOSAVEXenopusonlybias0p03D1.csv';
%sim1 = 'BIAS data/AttractionRepulsionD5beta0p05.csv';

if animal == 0
    sim1 = 'XENOPUS DATA FINAL2/AttrRepALLBIASEDto850eps200D7beta0p10.csv';
else
    sim1 = 'CHICK DATA FINAL/AttrRepALLBiasedeps75D7beta0p10.csv';
end

M1 = csvread(sim1);
M1 = M1(:,1)/60;

averageCoACiL(3) = mean(M1)
stdsCoACiL(3) =std(M1)


%% third CiL only

%sim1 = 'TOSAVETimeCiLonly0p03betaNarrowDomain.csv';
%sim1 = 'BIAS data/RepulsionD5beta0p05.csv';

if animal == 0
    sim1 = 'XENOPUS DATA FINAL2/RepOnlyALLBIASEDto850eps200D7beta0p10.csv';
else
    sim1 = 'CHICK DATA FINAL/RepOnlyALLBiasedeps75D7beta0p10.csv';
end

M1 = csvread(sim1);
M1 = M1(:,1)/60;

averageCiLonly(3) = mean(M1)
stdsCiLonly(3) =std(M1)

%% third VolExc

%sim1 = 'TOSAVETimeCiLonly0p03betaNarrowDomain.csv';
%sim1 = 'BIAS data/VolumeExclusionD5beta0p05.csv';

if animal == 0
    sim1 = 'XENOPUS DATA FINAL2/RepOnlyALLBIASEDto850eps1D7beta0p10.csv';
else
    sim1 = 'CHICK DATA FINAL/RepOnlyALLBiasedeps0D7beta0p10.csv';
end

M1 = csvread(sim1);
M1 = M1(:,1)/60;

averageVolExc(3) = mean(M1)
stdsVolExc(3) =std(M1)

% %% fourth CoA CiL
% 
% %sim1 = 'TOSAVETimeCoACiL0p04betaNarrowDomain.csv';
% %sim1 = 'TOSAVEXenopusonlybias0p04D1.csv';
% 
% M1 = csvread(sim1);
% M1 = M1(:,1)*18/total;
% 
% averageCoACiL(4) = mean(M1)
% stdsCoACiL(4) =std(M1)
% 
% 
% %% fourth CiL only
% 
% %sim1 = 'TOSAVETimeCiLonly0p04betaNarrowDomain.csv';
% 
% 
% M1 = csvread(sim1);
% M1 = M1(:,1)*18/total;
% 
% averageCiLonly(4) = mean(M1)
% stdsCiLonly(4) =std(M1)
% 
% 
% %% fifth CoA CiL
% 
% %sim1 = 'TOSAVETimeCoACiL0p05betaNarrowDomain.csv';
% %sim1 = 'TOSAVEXenopusonlybias0p05D1.csv';
% 
% M1 = csvread(sim1);
% M1 = M1(:,1)*18/total;
% 
% averageCoACiL(5) = mean(M1)
% stdsCoACiL(5) =std(M1)
% 
% 
% %% fifth CiL only
% 
% %sim1 = 'TOSAVETimeCiLonly0p05betaNarrowDomain.csv';
% 
% 
% M1 = csvread(sim1);
% M1 = M1(:,1)*18/total;
% 
% averageCiLonly(5) = mean(M1)
% stdsCiLonly(5) =std(M1)


figure
b1= bar(xCoACiL,averageCoACiL,0.15) 
 
hold on
b2 = bar(xCiLonly,averageCiLonly,0.15) 
b3 = bar(xVolExc,averageVolExc,0.15)
 
er1 = errorbar(xCoACiL,averageCoACiL,stdsCoACiL,'k.','linewidth',2)   
er2 = errorbar(xCiLonly, averageCiLonly, stdsCiLonly,'k.','linewidth',2)
er3 = errorbar(xVolExc, averageVolExc, stdsVolExc,'k.','linewidth',2)    

%xticks ([1 5 9]); 


%xticks([1.5,4.5,7.5,10.5,13.5])
%xticklabels({'0.01','0.02','0.03','0.04','0.05'})

xlim([0.3,12.3])
XTick = [2,6,10];
set(gca,'xtick',XTick)
xticklabels({'0.4','0.7','1.0'})
set(gca,'linewidth',2)
xlabel(['Bias, \beta'],'FontSize',14)
ylabel(['Time to invasion, hrs'],'FontSize',14)
set(gca,'FontSize',36)
 ax = gca;
%legend([b2,b1],'CiL only','CoA+CiL');
legend([b3,b2,b1],'Repulsion-only, \epsilon = 1','Repulsion-only, \epsilon = 200','Attraction and repulsion');% Xenopus
%legend([b3,b2,b1],'Repulsion-only, \epsilon = 0.4','Repulsion-only, \epsilon = 75.0','Attraction and repulsion');

 box on
 grid on
 %ylim([0,40])
 set(gca,'linewidth',4)   

