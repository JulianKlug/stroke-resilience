clear all
mainpath = '/Users/jk1/unige_onedrive/OneDrive - unige.ch/BCT/attacks240';
respath_st01 = fullfile(mainpath,'attack_ST01');
respath_st02 = fullfile(mainpath,'attack_ST02');
respath_st03 = fullfile(mainpath,'attack_ST03');
res_st01 = 'clin_attack_ST01_06-10-2020 16-06v4.mat';
res_st02 = 'clin_attack_ST02_06-11-2020 16-18v4.mat';
res_st03 = 'clin_attack_ST03_06-11-2020 16-26v4.mat';
TPcompa = 'TP2vsTP3';  %'TP2vsTP3'; 'TP1vsTP2'; 

finalpath = fullfile('/Users/jk1/stroke_research/resilience_stroke/attack_analysis/figures/clinical_attack');
propthr = 0.3:0.1:1;

cd (respath_st01)
data_st01 = load(res_st01);
cd (respath_st02)
data_st02 = load(res_st02);
cd (respath_st03)
data_st03 = load(res_st03);
cd(finalpath)

propnamelist={};
for p = 1:numel(propthr)
    propname = ['top' num2str(propthr(p)*100)];
    propnamelist{p} =propname;
end

%% for GE
GE_ST01_pre = data_st01.GlobEff_bin;
GE_ST02_pre = data_st02.GlobEff_bin;
GE_ST03_pre = data_st03.GlobEff_bin;

GE_ST01_post_corti = data_st01.GlobEff_bin_new.C;
GE_ST02_post_corti = data_st02.GlobEff_bin_new.C;
GE_ST03_post_corti = data_st03.GlobEff_bin_new.C;

GE_ST01_post_subcorti = data_st01.GlobEff_bin_new.SUB;
GE_ST02_post_subcorti = data_st02.GlobEff_bin_new.SUB;
GE_ST03_post_subcorti = data_st03.GlobEff_bin_new.SUB;

GE_ST01_post_csub = data_st01.GlobEff_bin_new.CSUB;
GE_ST02_post_csub = data_st02.GlobEff_bin_new.CSUB;
GE_ST03_post_csub = data_st03.GlobEff_bin_new.CSUB;

%% remove subjects with missing data
switch TPcompa  
    case 'TP1vsTP2'
        for p = 1:numel(propthr)
            propname = ['top' num2str(propthr(p)*100)];
            GE_ST_pre.(propname)(4) = [];%%%
        end
        GE_ST01_post_corti(4,:)=[];
        GE_ST01_post_subcorti(4,:)=[];
        GE_ST01_post_csub(4,:)=[];
        
        GE_ST02_post_corti(1,:)=[];
        GE_ST02_post_subcorti(1,:)=[];
        GE_ST02_post_csub(1,:)=[];
        
        % standardize matrix types
        for p = 1:numel(propthr)
            propname = ['top' num2str(propthr(p)*100)];
            GE_pre(:,p) = GE_ST_pre.(propname)';
        end
        
        % standardize variable names
        GE_gp1_post_corti = GE_ST01_post_corti;
        GE_gp2_post_corti = GE_ST02_post_corti;

        GE_gp1_post_subcorti = GE_ST01_post_subcorti;
        GE_gp2_post_subcorti = GE_ST02_post_subcorti;

        GE_gp1_post_csub = GE_ST01_post_csub;
        GE_gp2_post_csub = GE_ST02_post_csub;
        
        clear GE_ST01_pre GE_ST02_pre 
        clear GE_ST01_post_corti GE_ST02_post_corti GE_ST01_post_subcorti GE_ST02_post_subcorti GE_ST01_post_csub GE_ST02_post_csub
        
    case 'TP2vsTP3'
        for p = 1:numel(propthr)
            propname = ['top' num2str(propthr(p)*100)];
            GE_ST02_pre.(propname)(12) = [];%%%
        end
        for p = 1:numel(propthr)
            propname = ['top' num2str(propthr(p)*100)];
            GE_ST03_pre.(propname)(5) = [];%%%
        end
        GE_ST02_post_corti(12,:)=[];
        GE_ST02_post_subcorti(12,:)=[];
        GE_ST02_post_csub(12,:)=[];
        
        GE_ST03_post_corti(5,:)=[];
        GE_ST03_post_subcorti(5,:)=[];
        GE_ST03_post_csub(5,:)=[];
        
         % standardize matrix types
        for p = 1:numel(propthr)
            propname = ['top' num2str(propthr(p)*100)];
            GE_gp1_pre(:,p) = GE_ST02_pre.(propname)';
            GE_gp2_pre(:,p) = GE_ST03_pre.(propname)';
        end
        
        % standardize variable names
        GE_gp1_post_corti = GE_ST02_post_corti;
        GE_gp2_post_corti = GE_ST03_post_corti;
        
        GE_gp1_post_subcorti = GE_ST02_post_subcorti;
        GE_gp2_post_subcorti = GE_ST03_post_subcorti;
        
        GE_gp1_post_csub = GE_ST02_post_csub;
        GE_gp2_post_csub = GE_ST03_post_csub;
        
        clear  GE_ST02_pre GE_ST03_pre
        clear  GE_ST02_post_corti GE_ST03_post_corti  GE_ST02_post_subcorti GE_ST03_post_subcorti GE_ST01_post_csub GE_ST02_post_csub GE_ST03_post_csub
end

%% calculate means each group
mean_GE_gp1_pre = mean(GE_gp1_pre,1); 
mean_GE_gp2_pre = mean(GE_gp2_pre,1);

mean_GE_gp1_post_corti = mean(GE_gp1_post_corti,1);
mean_GE_gp2_post_corti = mean(GE_gp2_post_corti,1);
mean_GE_gp1_post_subcorti = mean(GE_gp1_post_subcorti,1);
mean_GE_gp2_post_subcorti = mean(GE_gp2_post_subcorti,1);
mean_GE_gp1_post_csub = mean(GE_gp1_post_csub,1);
mean_GE_gp2_post_csub = mean(GE_gp2_post_csub,1);

% calculate std each group
std_GE_gp1_pre = std(GE_gp1_pre,1); 
std_GE_gp2_pre = std(GE_gp2_pre,1);

std_GE_gp1_post_corti = std(GE_gp1_post_corti,1);
std_GE_gp2_post_corti = std(GE_gp2_post_corti,1);
std_GE_gp1_post_subcorti = std(GE_gp1_post_subcorti,1);
std_GE_gp2_post_subcorti = std(GE_gp2_post_subcorti,1);
std_GE_gp1_post_csub = std(GE_gp1_post_csub,1);
std_GE_gp2_post_csub = std(GE_gp2_post_csub,1);

% plot mean corti TP2vsTP3 CORTI
xaxis = 1:1:numel(propthr);

% ibm color palette 
tp1_color = [4/255,155/255,154/255];
tp2_color = [1/255,45/255,152/255];
tp3_color = [168/255,109/255,254/255];

faceAlpha = 0.2;

figure
p1 = plot(mean_GE_gp1_pre,'Color',tp2_color,'LineWidth',2)
hold on
fill([xaxis fliplr(xaxis)], [(mean_GE_gp1_pre + std_GE_gp1_pre), fliplr(mean_GE_gp1_pre - std_GE_gp1_pre)],tp2_color,'FaceAlpha', faceAlpha,'EdgeColor','none')
p2 = plot(mean_GE_gp2_pre,'Color',tp3_color,'LineWidth',2)
fill([xaxis fliplr(xaxis)], [(mean_GE_gp2_pre + std_GE_gp2_pre), fliplr(mean_GE_gp2_pre - std_GE_gp2_pre)],tp3_color,'FaceAlpha', faceAlpha,'EdgeColor','none')
p3 = plot(mean_GE_gp1_post_corti,'Linestyle','--','Color',tp2_color, 'LineWidth',2)
fill([xaxis fliplr(xaxis)], [(mean_GE_gp1_post_corti + std_GE_gp1_post_corti), fliplr(mean_GE_gp1_post_corti - std_GE_gp1_post_corti)],tp2_color,'FaceAlpha', faceAlpha,'EdgeColor','none')
p4 = plot(mean_GE_gp2_post_corti,'Linestyle','--','Color',tp3_color, 'LineWidth',2)
fill([xaxis fliplr(xaxis)], [(mean_GE_gp2_post_corti + std_GE_gp2_post_corti), fliplr(mean_GE_gp2_post_corti - std_GE_gp2_post_corti)],tp3_color,'FaceAlpha', faceAlpha,'EdgeColor','none')

xticklabels(propthr)
xlabel('Edge density','FontWeight','bold')
ylabel('Global efficiency','FontWeight','bold')
ylim([0.35 0.85])
title('Cortical MCA territory','FontWeight','bold')
switch TPcompa  
    case 'TP1vsTP2'
        legend ([p1 p3 p2 p4],{'Stroke TP1 pre', 'Stroke TP1 post', 'Stroke TP2 post','Stroke TP2 post'}, 'Location','southeast')
    case 'TP2vsTP3'
        legend ([p1 p3 p2 p4],{'Stroke TP2 pre', 'Stroke TP2 post', 'Stroke TP3 pre','Stroke TP3 post'}, 'Location','southeast')
end
set(gca,'FontSize',13)
plotname = [TPcompa '_cortical_GE_paper'];
saveas(gcf,plotname)
print(gcf, '-dtiff', [plotname '.tiff']);

% plot mean subcorti TP2vsTP3 SUBCORTI
figure
p1 = plot(mean_GE_gp1_pre,'Color',tp2_color,'LineWidth',2)
hold on
fill([xaxis fliplr(xaxis)], [(mean_GE_gp1_pre + std_GE_gp1_pre), fliplr(mean_GE_gp1_pre - std_GE_gp1_pre)],tp2_color,'FaceAlpha', faceAlpha,'EdgeColor','none')
p2 = plot(mean_GE_gp2_pre,'Color',tp3_color,'LineWidth',2)
fill([xaxis fliplr(xaxis)], [(mean_GE_gp2_pre + std_GE_gp2_pre), fliplr(mean_GE_gp2_pre - std_GE_gp2_pre)],tp3_color,'FaceAlpha', faceAlpha,'EdgeColor','none')
p3 = plot(mean_GE_gp1_post_subcorti,'Linestyle','--','Color',tp2_color, 'LineWidth',2)
fill([xaxis fliplr(xaxis)], [(mean_GE_gp1_post_subcorti + std_GE_gp1_post_subcorti), fliplr(mean_GE_gp1_post_subcorti - std_GE_gp1_post_subcorti)],tp2_color,'FaceAlpha', faceAlpha,'EdgeColor','none')
p4 = plot(mean_GE_gp2_post_subcorti,'Linestyle','--','Color',tp3_color, 'LineWidth',2)
fill([xaxis fliplr(xaxis)], [(mean_GE_gp2_post_subcorti + std_GE_gp2_post_subcorti), fliplr(mean_GE_gp2_post_subcorti - std_GE_gp2_post_subcorti)],tp3_color,'FaceAlpha', faceAlpha,'EdgeColor','none')

xticklabels(propthr)
xlabel('Edge density','FontWeight','bold')
ylabel('Global efficiency','FontWeight','bold')
ylim([0.35 0.85])
title('Subcortical MCA territory','FontWeight','bold')
switch TPcompa  
    case 'TP1vsTP2'
        legend ([p1 p3 p2 p4],{'Stroke TP1 pre', 'Stroke TP1 post', 'Stroke TP2 post','Stroke TP2 post'}, 'Location','southeast')
    case 'TP2vsTP3'
        legend ([p1 p3 p2 p4],{'Stroke TP2 pre', 'Stroke TP2 post', 'Stroke TP3 pre','Stroke TP3 post'}, 'Location','southeast')
end
set(gca,'FontSize',13)
plotname = [TPcompa '_subcortical_GE_paper'];
saveas(gcf,plotname)
print(gcf, '-dtiff', [plotname '.tiff']);

% plot mean subcorti TP2vsTP3 CORTICO-SUBCORTI
figure
p1 = plot(mean_GE_gp1_pre,'Color',tp2_color,'LineWidth',2)
hold on
fill([xaxis fliplr(xaxis)], [(mean_GE_gp1_pre + std_GE_gp1_pre), fliplr(mean_GE_gp1_pre - std_GE_gp1_pre)],tp2_color,'FaceAlpha', faceAlpha,'EdgeColor','none')
p2 = plot(mean_GE_gp2_pre,'Color',tp3_color,'LineWidth',2)
fill([xaxis fliplr(xaxis)], [(mean_GE_gp2_pre + std_GE_gp2_pre), fliplr(mean_GE_gp2_pre - std_GE_gp2_pre)],tp3_color,'FaceAlpha', faceAlpha,'EdgeColor','none')
p3 = plot(mean_GE_gp1_post_csub,'Linestyle','--','Color',tp2_color, 'LineWidth',2)
fill([xaxis fliplr(xaxis)], [(mean_GE_gp1_post_csub + std_GE_gp1_post_csub), fliplr(mean_GE_gp1_post_csub - std_GE_gp1_post_csub)],tp2_color,'FaceAlpha', faceAlpha,'EdgeColor','none')
p4 = plot(mean_GE_gp2_post_csub,'Linestyle','--','Color',tp3_color, 'LineWidth',2)
fill([xaxis fliplr(xaxis)], [(mean_GE_gp2_post_csub + std_GE_gp2_post_csub), fliplr(mean_GE_gp2_post_csub - std_GE_gp2_post_csub)],tp3_color,'FaceAlpha', faceAlpha,'EdgeColor','none')

xticklabels(propthr)
xlabel('Edge density','FontWeight','bold')
ylabel('Global efficiency','FontWeight','bold')
ylim([0.35 0.85])
title('Cortico-subortical MCA territory','FontWeight','bold')
switch TPcompa  
    case 'TP1vsTP2'
        legend ([p1 p3 p2 p4],{'Stroke TP1 pre', 'Stroke TP1 post', 'Stroke TP2 post','Stroke TP2 post'}, 'Location','southeast')
    case 'TP2vsTP3'
        legend ([p1 p3 p2 p4],{'Stroke TP2 pre', 'Stroke TP2 post', 'Stroke TP3 pre','Stroke TP3 post'}, 'Location','southeast')
end
set(gca,'FontSize',13)
plotname = [TPcompa '_csub_GE_paper'];
saveas(gcf,plotname)
print(gcf, '-dtiff', [plotname '.tiff']);


