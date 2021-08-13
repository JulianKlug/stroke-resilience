% may2020: creates figures after FDA ananylsis, spontaneous changes
% Mitsouko van Assche

clear all
mainpath = '/Users/jk1/unige_onedrive/OneDrive - unige.ch/BCT/atlas_BNA/BNA_240_flipped_N32_retroicor_SBB4_prop_bin';
fpath = fullfile('/Users/jk1/stroke_research/resilience_stroke/longitudinal_analysis/figures');
STpath = fullfile(mainpath,'ST');
HCpath = fullfile(mainpath,'HC');
propthr = 0.3:0.1:1; %since below graph are disconnected

% ibm color palette 
tp1_color = [4/255,155/255,154/255];
tp2_color = [1/255,45/255,152/255];
tp3_color = [168/255,109/255,254/255];

faceAlpha = 0.05;

%% GLOBAL EFFICIENCY
cd(STpath)
datafile_ST = 'CharPath240_prop_bin.mat';%%%%
load(datafile_ST)
data_TP1 = GlobEfficiency_bin.ST01;
data_TP2 = GlobEfficiency_bin.ST02;
data_TP3 = GlobEfficiency_bin.ST03;
cd(HCpath)
datafile_HC = 'CharPath240_bin_HC.mat';%%%%
load(datafile_HC)
data_HC = GlobEfficiency;

cd(fpath)
for p=1:numel(propthr)
    propname = ['top' num2str(propthr(p)*100)];
    data1 = data_TP1.(propname);
    mean_data1(p) = mean(data1);
    std_data1(p) = std(data1);
    
    data2 = data_TP2.(propname);
    mean_data2(p) = mean(data2);
    std_data2(p) = std(data2);
    
    data3 = data_TP3.(propname);
    mean_data3(p) = mean(data3);
    std_data3(p) = std(data3);
    
    datahc = data_HC.(propname);
    mean_datahc(p) = mean(datahc);
    std_datahc(p) = std(datahc);
end

xaxis = 1:1:8;
data1_sup = mean_data1 + std_data1;
data1_low = mean_data1 - std_data1;
data2_sup = mean_data2 + std_data2;
data2_low = mean_data2 - std_data2;
data3_sup = mean_data3 + std_data3;
data3_low = mean_data3 - std_data3;
datahc_low = mean_datahc + std_datahc;
datahc_sup = mean_datahc - std_datahc;

figure
p1 = plot(mean_data1,'Color',tp1_color,'LineWidth',2)
hold on
fill([xaxis fliplr(xaxis)], [data1_sup, fliplr(data1_low)],tp1_color,'FaceAlpha',faceAlpha,'EdgeColor','none')
p = plot(mean_datahc,'Color','k','LineWidth',2)
fill([xaxis fliplr(xaxis)], [datahc_sup, fliplr(datahc_low)],'k','FaceAlpha',faceAlpha,'EdgeColor','none')
legend([p1 p],{'Stroke TP1', 'HC'},'Location','northwest');
% title('FDA analysis of global efficiency','FontWeight','bold')
xticklabels(propthr)
xlabel('Edge density','FontWeight','bold')
ylabel('Global efficiency','FontWeight','bold')
set(gca,'FontSize',13)
ylim([0.64 0.84])
saveas(gcf,'FDA_GlobEff_transv1')
print(gcf, '-dtiff', 'FDA_GlobEff_transv1.tiff');

figure
p2 = plot(mean_data2,'Color',tp2_color,'LineWidth',2)
hold on
fill([xaxis fliplr(xaxis)], [data2_sup, fliplr(data2_low)],tp2_color,'FaceAlpha',faceAlpha,'EdgeColor','none')
p = plot(mean_datahc,'Color','k','LineWidth',2)
fill([xaxis fliplr(xaxis)], [datahc_sup, fliplr(datahc_low)],'k','FaceAlpha',faceAlpha,'EdgeColor','none')
legend([p2 p],{'Stroke TP2', 'HC'},'Location','northwest');
% title('FDA analysis of global efficiency')
xticklabels(propthr)
xlabel('Edge density','FontWeight','bold')
ylabel('Global efficiency','FontWeight','bold')
set(gca,'FontSize',13)
ylim([0.64 0.84])
saveas(gcf,'FDA_GlobEff_transv2')
print(gcf, '-dtiff', 'FDA_GlobEff_transv2.tiff');

figure
p3 = plot(mean_data3,'Color',tp3_color,'LineWidth',2)
hold on
fill([xaxis fliplr(xaxis)], [data3_sup, fliplr(data3_low)],tp3_color,'FaceAlpha',faceAlpha,'EdgeColor','none')
p = plot(mean_datahc,'Color','k','LineWidth',2)
fill([xaxis fliplr(xaxis)], [datahc_sup, fliplr(datahc_low)],'k','FaceAlpha',faceAlpha,'EdgeColor','none')
legend([p3 p],{'Stroke TP3', 'HC'},'Location','northwest');
% title('FDA analysis of global efficiency')
xticklabels(propthr)
xlabel('Edge density','FontWeight','bold')
ylabel('Global efficiency','FontWeight','bold')
set(gca,'FontSize',13)
ylim([0.64 0.84])
saveas(gcf,'FDA_GlobEff_transv3')
print(gcf, '-dtiff', 'FDA_GlobEff_transv3.tiff');

%% PATH LENGTH
data_TP1 = Charpathlength_bin.ST01;
data_TP2 = Charpathlength_bin.ST02;
data_TP3 = Charpathlength_bin.ST03;
data_HC = Charpathlength;
for p=1:numel(propthr)
    propname = ['top' num2str(propthr(p)*100)];
    data1 = data_TP1.(propname);
    mean_data1(p) = mean(data1);
    std_data1(p) = std(data1);
    
    data2 = data_TP2.(propname);
    mean_data2(p) = mean(data2);
    std_data2(p) = std(data2);
    
    data3 = data_TP3.(propname);
    mean_data3(p) = mean(data3);
    std_data3(p) = std(data3);
    
    datahc = data_HC.(propname);
    mean_datahc(p) = mean(datahc);
    std_datahc(p) = std(datahc);
end

xaxis = 1:1:8;
data1_sup = mean_data1 + std_data1;
data1_low = mean_data1 - std_data1;
data2_sup = mean_data2 + std_data2;
data2_low = mean_data2 - std_data2;
data3_sup = mean_data3 + std_data3;
data3_low = mean_data3 - std_data3;
datahc_low = mean_datahc + std_datahc;
datahc_sup = mean_datahc - std_datahc;


figure
p1 = plot(mean_data1,'Color',tp1_color,'LineWidth',2)
hold on
fill([xaxis fliplr(xaxis)], [data1_sup, fliplr(data1_low)],tp1_color,'FaceAlpha',faceAlpha,'EdgeColor','none')
p = plot(mean_datahc,'Color','k','LineWidth',2)
fill([xaxis fliplr(xaxis)], [datahc_sup, fliplr(datahc_low)],'k','FaceAlpha',faceAlpha,'EdgeColor','none')
legend([p1 p],{'Stroke TP1', 'HC'},'Location','northwest');
% title('FDA analysis of path length')
xticklabels(propthr)
xlabel('Edge density','FontWeight','bold')
ylabel('Path Length','FontWeight','bold')
set(gca,'FontSize',13)
ylim([1.35 1.8])
saveas(gcf,'FDA_Pathlength_transv1')
print(gcf, '-dtiff', 'FDA_Pathlength_transv1.tiff');

figure
p2 = plot(mean_data2,'Color',tp2_color,'LineWidth',2)
hold on
fill([xaxis fliplr(xaxis)], [data2_sup, fliplr(data2_low)],tp2_color,'FaceAlpha',faceAlpha,'EdgeColor','none')
p = plot(mean_datahc,'Color','k','LineWidth',2)
fill([xaxis fliplr(xaxis)], [datahc_sup, fliplr(datahc_low)],'k','FaceAlpha',faceAlpha,'EdgeColor','none')
legend([p2 p],{'Stroke TP2', 'HC'},'Location','northwest');
% title('FDA analysis of path length')
xticklabels(propthr)
xlabel('Edge density','FontWeight','bold')
ylabel('Path Length','FontWeight','bold')
set(gca,'FontSize',13)
ylim([1.35 1.8])
saveas(gcf,'FDA_Pathlength_transv2')
print(gcf, '-dtiff', 'FDA_Pathlength_transv2.tiff');

figure
p3 = plot(mean_data3,'Color',tp3_color,'LineWidth',2)
hold on
fill([xaxis fliplr(xaxis)], [data3_sup, fliplr(data3_low)],tp3_color,'FaceAlpha',faceAlpha,'EdgeColor','none')
p = plot(mean_datahc,'Color','k','LineWidth',2)
fill([xaxis fliplr(xaxis)], [datahc_sup, fliplr(datahc_low)],'k','FaceAlpha',faceAlpha,'EdgeColor','none')
legend([p3 p],{'Stroke TP3', 'HC'},'Location','northwest');
% title('FDA analysis of path length')
xticklabels(propthr)
xlabel('Edge density','FontWeight','bold')
ylabel('Path Length','FontWeight','bold')
set(gca,'FontSize',13)
ylim([1.35 1.8])
saveas(gcf,'FDA_Pathlength_transv3')
print(gcf, '-dtiff', 'FDA_Pathlength_transv3.tiff');

%% CC
cd(STpath)
datafile_ST = 'CC240_prop_bin.mat';%%%%
load(datafile_ST)
data_TP1 = CC_all_mean_bin.ST01;
data_TP2 = CC_all_mean_bin.ST02;
data_TP3 = CC_all_mean_bin.ST03;
cd(HCpath)
datafile_HC = 'CC240_bin_HC.mat';%%%%
load(datafile_HC)
data_HC = CC_all_mean;
cd(fpath)
for p=1:numel(propthr)
    propname = ['top' num2str(propthr(p)*100)];
    data1 = data_TP1.(propname);
    mean_data1(p) = mean(data1);
    std_data1(p) = std(data1);
    
    data2 = data_TP2.(propname);
    mean_data2(p) = mean(data2);
    std_data2(p) = std(data2);
    
    data3 = data_TP3.(propname);
    mean_data3(p) = mean(data3);
    std_data3(p) = std(data3);
    
    datahc = data_HC.(propname);
    mean_datahc(p) = mean(datahc);
    std_datahc(p) = std(datahc);
end

xaxis = 1:1:8;
data1_sup = mean_data1 + std_data1;
data1_low = mean_data1 - std_data1;
data2_sup = mean_data2 + std_data2;
data2_low = mean_data2 - std_data2;
data3_sup = mean_data3 + std_data3;
data3_low = mean_data3 - std_data3;
datahc_low = mean_datahc + std_datahc;
datahc_sup = mean_datahc - std_datahc;

figure
p1 = plot(mean_data1,'Color',tp1_color,'LineWidth',2)
hold on
fill([xaxis fliplr(xaxis)], [data1_sup, fliplr(data1_low)],tp1_color,'FaceAlpha',faceAlpha,'EdgeColor','none')
p = plot(mean_datahc,'Color','k','LineWidth',2)
fill([xaxis fliplr(xaxis)], [datahc_sup, fliplr(datahc_low)],'k','FaceAlpha',faceAlpha,'EdgeColor','none')
legend([p1 p],{'Stroke TP1', 'HC'},'Location','northwest');
% title('FDA analysis of clustering coefficient')
xticklabels(propthr)
xlabel('Edge density','FontWeight','bold')
ylabel('Clustering coefficient','FontWeight','bold')
set(gca,'FontSize',13)
saveas(gcf,'FDA_CC_transv1')
print(gcf, '-dtiff', 'FDA_CC_transv1.tiff');

figure
p2 = plot(mean_data2,'Color',tp2_color,'LineWidth',2)
hold on
fill([xaxis fliplr(xaxis)], [data2_sup, fliplr(data2_low)],tp2_color,'FaceAlpha',faceAlpha,'EdgeColor','none')
p = plot(mean_datahc,'Color','k','LineWidth',2)
fill([xaxis fliplr(xaxis)], [datahc_sup, fliplr(datahc_low)],'k','FaceAlpha',faceAlpha,'EdgeColor','none')
legend([p2 p],{'Stroke TP2', 'HC'},'Location','northwest');
% title('FDA analysis of clustering coefficient')
xticklabels(propthr)
xlabel('Edge density','FontWeight','bold')
ylabel('Clustering coefficient','FontWeight','bold')
set(gca,'FontSize',13)
saveas(gcf,'FDA_CC_transv2')
print(gcf, '-dtiff', 'FDA_CC_transv2.tiff');

figure
p3 = plot(mean_data3,'Color',tp3_color,'LineWidth',2)
hold on
fill([xaxis fliplr(xaxis)], [data3_sup, fliplr(data3_low)],tp3_color,'FaceAlpha',faceAlpha,'EdgeColor','none')
p = plot(mean_datahc,'Color','k','LineWidth',2)
fill([xaxis fliplr(xaxis)], [datahc_sup, fliplr(datahc_low)],'k','FaceAlpha',faceAlpha,'EdgeColor','none')
legend([p3 p],{'Stroke TP3', 'HC'},'Location','northwest');
% title('FDA analysis of clustering coefficient')
xticklabels(propthr)
xlabel('Edge density','FontWeight','bold')
ylabel('Clustering coefficient','FontWeight','bold')
set(gca,'FontSize',13)
saveas(gcf,'FDA_CC_transv3')
print(gcf, '-dtiff', 'FDA_CC_transv3.tiff');

%% BC
cd(STpath)
datafile_ST = 'Betweenness240_bin.mat';%%%%
load(datafile_ST)
data_TP1 = betweenness.ST01;
data_TP2 = betweenness.ST02;
data_TP3 = betweenness.ST03;
cd(HCpath)
datafile_HC = 'Betweenness240_bin_HC.mat';%%%%
load(datafile_HC)
data_HC = betweenness;
cd(fpath)
for p=1:numel(propthr)
    propname = ['top' num2str(propthr(p)*100)];
    data1 = mean(data_TP1.(propname),2);
    mean_data1(p) = mean(data1);
    std_data1(p) = std(data1);
    
    data2 = mean(data_TP2.(propname),2);
    mean_data2(p) = mean(data2);
    std_data2(p) = std(data2);
    
    data3 = mean(data_TP3.(propname),2);
    mean_data3(p) = mean(data3);
    std_data3(p) = std(data3);
    
    datahc = mean(data_HC.(propname),2);
    mean_datahc(p) = mean(datahc);
    std_datahc(p) = std(datahc);
end

xaxis = 1:1:8;
data1_sup = mean_data1 + std_data1;
data1_low = mean_data1 - std_data1;
data2_sup = mean_data2 + std_data2;
data2_low = mean_data2 - std_data2;
data3_sup = mean_data3 + std_data3;
data3_low = mean_data3 - std_data3;
datahc_low = mean_datahc + std_datahc;
datahc_sup = mean_datahc - std_datahc;


figure
p1 = plot(mean_data1,'Color',tp1_color,'LineWidth',2)
hold on
fill([xaxis fliplr(xaxis)], [data1_sup, fliplr(data1_low)],tp1_color,'FaceAlpha',faceAlpha,'EdgeColor','none')
p = plot(mean_datahc,'Color','k','LineWidth',2)
fill([xaxis fliplr(xaxis)], [datahc_sup, fliplr(datahc_low)],'k','FaceAlpha',faceAlpha,'EdgeColor','none')
legend([p1 p],{'Stroke TP1', 'HC'},'Location','northeast');
% title('FDA analysis of betweenness centrality')
xticklabels(propthr)
xlabel('Edge density','FontWeight','bold')
ylabel('Betweenness centrality','FontWeight','bold')
set(gca,'FontSize',13)
saveas(gcf,'FDA_BC_transv1')
print(gcf, '-dtiff', 'FDA_BC_transv1.tiff');

figure
p2 = plot(mean_data2,'Color',tp2_color,'LineWidth',2)
hold on
fill([xaxis fliplr(xaxis)], [data2_sup, fliplr(data2_low)],tp2_color,'FaceAlpha',faceAlpha,'EdgeColor','none')
p = plot(mean_datahc,'Color','k','LineWidth',2)
fill([xaxis fliplr(xaxis)], [datahc_sup, fliplr(datahc_low)],'k','FaceAlpha',faceAlpha,'EdgeColor','none')
legend([p2 p],{'Stroke TP2', 'HC'},'Location','northeast');
% title('FDA analysis of betweenness centrality')
xticklabels(propthr)
xlabel('Edge density','FontWeight','bold')
ylabel('Betweenness centrality','FontWeight','bold')
set(gca,'FontSize',13)
saveas(gcf,'FDA_BC_transv2')
print(gcf, '-dtiff', 'FDA_BC_transv2.tiff');

figure
p3 = plot(mean_data3,'Color',tp3_color,'LineWidth',2)
hold on
fill([xaxis fliplr(xaxis)], [data3_sup, fliplr(data3_low)],tp3_color,'FaceAlpha',faceAlpha,'EdgeColor','none')
p = plot(mean_datahc,'Color','k','LineWidth',2)
fill([xaxis fliplr(xaxis)], [datahc_sup, fliplr(datahc_low)],'k','FaceAlpha',faceAlpha,'EdgeColor','none')
legend([p3 p],{'Stroke TP3', 'HC'},'Location','northeast');
% title('FDA analysis of betweenness centrality')
xticklabels(propthr)
xlabel('Edge density','FontWeight','bold')
ylabel('Betweenness centrality','FontWeight','bold')
set(gca,'FontSize',13)
saveas(gcf,'FDA_BC_transv3')
print(gcf, '-dtiff', 'FDA_BC_transv3.tiff');

%% SIGMA
% cf
% 'D:\STROKE\MRIdata\BCT\atlas_BNA\BNA_240_flipped_N32_retroicor_SBB4_prop_bin\HC'

%% DEGREE
cd(STpath)
datafile_ST = 'Degrees240_prop_bin.mat';%%%%
load(datafile_ST)
data_TP1 = mymeandegree.ST01;
data_TP2 = mymeandegree.ST02;
data_TP3 = mymeandegree.ST03;
cd(HCpath)
datafile_HC = 'Degrees240_bin_HC.mat';%%%%
load(datafile_HC)
data_HC = mymeandegree_bin;

cd(fpath)
for p=1:numel(propthr)
    propname = ['top' num2str(propthr(p)*100)];
    data1 = data_TP1.(propname);
    mean_data1(p) = mean(data1);
    std_data1(p) = std(data1);
    
    data2 = data_TP2.(propname);
    mean_data2(p) = mean(data2);
    std_data2(p) = std(data2);
    
    data3 = data_TP3.(propname);
    mean_data3(p) = mean(data3);
    std_data3(p) = std(data3);
    
    datahc = data_HC.(propname);
    mean_datahc(p) = mean(datahc);
    std_datahc(p) = std(datahc);
end

xaxis = 1:1:8;
data1_sup = mean_data1 + std_data1;
data1_low = mean_data1 - std_data1;
data2_sup = mean_data2 + std_data2;
data2_low = mean_data2 - std_data2;
data3_sup = mean_data3 + std_data3;
data3_low = mean_data3 - std_data3;
datahc_low = mean_datahc + std_datahc;
datahc_sup = mean_datahc - std_datahc;

figure
p1 = plot(mean_data1,'Color',tp1_color,'LineWidth',2)
hold on
fill([xaxis fliplr(xaxis)], [data1_sup, fliplr(data1_low)],tp1_color,'FaceAlpha',faceAlpha,'EdgeColor','none')
p = plot(mean_datahc,'Color','k','LineWidth',2)
fill([xaxis fliplr(xaxis)], [datahc_sup, fliplr(datahc_low)],'k','FaceAlpha',faceAlpha,'EdgeColor','none')
legend([p1 p],{'Stroke TP1', 'HC'},'Location','northwest');
% title('FDA analysis of degrees')
xticklabels(propthr)
xlabel('Edge density','FontWeight','bold')
ylabel('Degree','FontWeight','bold')
set(gca,'FontSize',13)
ylim([70 160])
saveas(gcf,'FDA_degrees_transv1')
print(gcf, '-dtiff', 'FDA_degrees_transv1.tiff');

figure
p2 = plot(mean_data2,'Color',tp2_color,'LineWidth',2)
hold on
fill([xaxis fliplr(xaxis)], [data2_sup, fliplr(data2_low)],tp2_color,'FaceAlpha',faceAlpha,'EdgeColor','none')
p = plot(mean_datahc,'Color','k','LineWidth',2)
fill([xaxis fliplr(xaxis)], [datahc_sup, fliplr(datahc_low)],'k','FaceAlpha',faceAlpha,'EdgeColor','none')
legend([p2 p],{'Stroke TP2', 'HC'},'Location','northwest');
% title('FDA analysis of degrees')
xticklabels(propthr)
xlabel('Edge density','FontWeight','bold')
ylabel('Degree','FontWeight','bold')
set(gca,'FontSize',13)
ylim([70 160])
saveas(gcf,'FDA_degrees_transv2')
print(gcf, '-dtiff', 'FDA_degrees_transv2.tiff');

figure
p3 = plot(mean_data3,'Color',tp3_color,'LineWidth',2)
hold on
fill([xaxis fliplr(xaxis)], [data3_sup, fliplr(data3_low)],tp3_color,'FaceAlpha',faceAlpha,'EdgeColor','none')
p = plot(mean_datahc,'Color','k','LineWidth',2)
fill([xaxis fliplr(xaxis)], [datahc_sup, fliplr(datahc_low)],'k','FaceAlpha',faceAlpha,'EdgeColor','none')
legend([p3 p],{'Stroke TP3', 'HC'},'Location','northwest');
% title('FDA analysis of degrees')
xticklabels(propthr)
xlabel('Edge density','FontWeight','bold')
ylabel('Degree','FontWeight','bold')
set(gca,'FontSize',13)
ylim([70 160])
saveas(gcf,'FDA_degrees_transv3')
print(gcf, '-dtiff', 'FDA_degrees_transv3.tiff');