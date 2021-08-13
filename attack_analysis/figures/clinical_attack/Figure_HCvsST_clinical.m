%perform FDA applied to graph metrics for stroke connect project
%v1: nov 2019: Mitsouko van Assche; application to BNA atlas 28 rois
%v1.1: nov 2019: Mitsouko van Assche; because there are some missing data
%in the patients (P1 at TP1, P5 at TP2, P17 at TP3), the exact matched
%control group is used to calculate reference values before permutation
%v3: may20: adapted for density range 30-100%
clear all
mainpath = '/Users/jk1/unige_onedrive/OneDrive - unige.ch/BCT/attacks240';
respath_st01 = fullfile(mainpath,'attack_ST01');
respath_st02 = fullfile(mainpath,'attack_ST02');
respath_st03 = fullfile(mainpath,'attack_ST03');
respath_hc = fullfile(mainpath,'attack_HC');
finalpath = fullfile('/Users/jk1/stroke_research/resilience_stroke/attack_analysis/figures/clinical_attack');
res_st01 = 'clin_attack_ST01_06-10-2020 16-06v4.mat';
res_st02 = 'clin_attack_ST02_06-11-2020 16-18v4.mat';
res_st03 = 'clin_attack_ST03_06-11-2020 16-26v4.mat';
res_hc = 'clin_attack_HC_06-11-2020 16-34v4.mat';

compa1 = 'TP1'; 
compa2 = 'TP2'; 
compa3 = 'TP3'; 
CV = 'LOSGO'; % LOSO, LOSGO
TP = {'ST01','ST02','ST03'};%%%%
nPerms=10000; %%
baseFileName = 'clinical attack';
propthr = 0.3:0.1:1;

cd (respath_st01)
data_st01 = load(res_st01);
cd (respath_st02)
data_st02 = load(res_st02);
cd (respath_st03)
data_st03 = load(res_st03);
cd (respath_hc)
data_hc = load(res_hc);
cd(finalpath)

propnamelist={};
for p = 1:numel(propthr)
    propname = ['top' num2str(propthr(p)*100)];
    propnamelist{p} =propname;
end

%% for GE
GE_ST01_pre_corti = data_st01.GlobEff_bin;
GE_ST02_pre_corti = data_st02.GlobEff_bin;
GE_ST03_pre_corti = data_st03.GlobEff_bin;
GE_ST01_post_corti = data_st01.GlobEff_bin_new.C;
GE_ST02_post_corti = data_st02.GlobEff_bin_new.C;
GE_ST03_post_corti = data_st03.GlobEff_bin_new.C;

GE_ST01_pre_subcorti = data_st01.GlobEff_bin;
GE_ST02_pre_subcorti = data_st02.GlobEff_bin;
GE_ST03_pre_subcorti = data_st03.GlobEff_bin;
GE_ST01_post_subcorti = data_st01.GlobEff_bin_new.SUB;
GE_ST02_post_subcorti = data_st02.GlobEff_bin_new.SUB;
GE_ST03_post_subcorti = data_st03.GlobEff_bin_new.SUB;

GE_ST01_pre_csub = data_st01.GlobEff_bin;
GE_ST02_pre_csub = data_st02.GlobEff_bin;
GE_ST03_pre_csub = data_st03.GlobEff_bin;
GE_ST01_post_csub = data_st01.GlobEff_bin_new.CSUB;
GE_ST02_post_csub = data_st02.GlobEff_bin_new.CSUB;
GE_ST03_post_csub = data_st03.GlobEff_bin_new.CSUB;

GE_HC01_pre_corti = data_hc.GlobEff_bin;
GE_HC02_pre_corti = data_hc.GlobEff_bin;
GE_HC03_pre_corti = data_hc.GlobEff_bin;
GE_HC01_post_corti = data_hc.GlobEff_bin_new.C;
GE_HC02_post_corti = data_hc.GlobEff_bin_new.C;
GE_HC03_post_corti = data_hc.GlobEff_bin_new.C;

GE_HC01_pre_subcorti = data_hc.GlobEff_bin;
GE_HC02_pre_subcorti = data_hc.GlobEff_bin;
GE_HC03_pre_subcorti = data_hc.GlobEff_bin;
GE_HC01_post_subcorti = data_hc.GlobEff_bin_new.SUB;
GE_HC02_post_subcorti = data_hc.GlobEff_bin_new.SUB;
GE_HC03_post_subcorti = data_hc.GlobEff_bin_new.SUB;

GE_HC01_pre_csub = data_hc.GlobEff_bin;
GE_HC02_pre_csub = data_hc.GlobEff_bin;
GE_HC03_pre_csub = data_hc.GlobEff_bin;
GE_HC01_post_csub = data_hc.GlobEff_bin_new.CSUB;
GE_HC02_post_csub = data_hc.GlobEff_bin_new.CSUB;
GE_HC03_post_csub = data_hc.GlobEff_bin_new.CSUB;

% remove controls not paired with any patient
for p = 1:numel(propthr)
    propname = ['top' num2str(propthr(p)*100)];
    GE_HC01_pre_corti.(propname)(9)=[]; %cortical
    GE_HC02_pre_corti.(propname)(5)=[]; 
    GE_HC03_pre_corti.(propname)(8)=[];  
    GE_HC01_pre_subcorti.(propname)(9)=[];%subcortical
    GE_HC02_pre_subcorti.(propname)(5)=[]; 
    GE_HC03_pre_subcorti.(propname)(8)=[];  
    GE_HC01_pre_csub.(propname)(9)=[];%corticosubcortical
    GE_HC02_pre_csub.(propname)(5)=[]; 
    GE_HC03_pre_csub.(propname)(8)=[];  
    
end
GE_HC01_post_corti(9,:)=[];
GE_HC02_post_corti(5,:)=[];
GE_HC03_post_corti(8,:)=[];
GE_HC01_post_subcorti(9,:)=[];
GE_HC02_post_subcorti(5,:)=[]; 
GE_HC03_post_subcorti(8,:)=[];
GE_HC01_post_csub(9,:)=[];
GE_HC02_post_csub(5,:)=[];
GE_HC03_post_csub(8,:)=[];

% standardize matrix types
for p = 1:numel(propthr)
    propname = ['top' num2str(propthr(p)*100)];
    GE_TP1_pre_corti(:,p) = GE_ST01_pre_corti.(propname)';
    GE_TP2_pre_corti(:,p) = GE_ST02_pre_corti.(propname)';
    GE_TP3_pre_corti(:,p) = GE_ST03_pre_corti.(propname)';
    GE_HC1_pre_corti(:,p) = GE_HC01_pre_corti.(propname)';
    GE_HC2_pre_corti(:,p) = GE_HC02_pre_corti.(propname)';
    GE_HC3_pre_corti(:,p) = GE_HC03_pre_corti.(propname)';
end
for p = 1:numel(propthr)
    propname = ['top' num2str(propthr(p)*100)];
    GE_TP1_pre_subcorti(:,p) = GE_ST01_pre_subcorti.(propname)';
    GE_TP2_pre_subcorti(:,p) = GE_ST02_pre_subcorti.(propname)';
    GE_TP3_pre_subcorti(:,p) = GE_ST03_pre_subcorti.(propname)';
    GE_HC1_pre_subcorti(:,p) = GE_HC01_pre_subcorti.(propname)';
    GE_HC2_pre_subcorti(:,p) = GE_HC02_pre_subcorti.(propname)';
    GE_HC3_pre_subcorti(:,p) = GE_HC03_pre_subcorti.(propname)';
end
for p = 1:numel(propthr)
    propname = ['top' num2str(propthr(p)*100)];
    GE_TP1_pre_csub(:,p) = GE_ST01_pre_csub.(propname)';
    GE_TP2_pre_csub(:,p) = GE_ST02_pre_csub.(propname)';
    GE_TP3_pre_csub(:,p) = GE_ST03_pre_csub.(propname)';
    GE_HC1_pre_csub(:,p) = GE_HC01_pre_csub.(propname)';
    GE_HC2_pre_csub(:,p) = GE_HC02_pre_csub.(propname)';
    GE_HC3_pre_csub(:,p) = GE_HC03_pre_csub.(propname)';
end
% standardize variable names
GE_TP1_post_corti = GE_ST01_post_corti;
GE_TP2_post_corti = GE_ST02_post_corti;
GE_TP3_post_corti = GE_ST03_post_corti;
GE_HC1_post_corti = GE_HC01_post_corti;
GE_HC2_post_corti = GE_HC02_post_corti;
GE_HC3_post_corti = GE_HC03_post_corti;

GE_TP1_post_subcorti = GE_ST01_post_subcorti;
GE_TP2_post_subcorti = GE_ST02_post_subcorti;
GE_TP3_post_subcorti = GE_ST03_post_subcorti;
GE_HC1_post_subcorti = GE_HC01_post_subcorti;
GE_HC2_post_subcorti = GE_HC02_post_subcorti;
GE_HC3_post_subcorti = GE_HC03_post_subcorti;

GE_TP1_post_csub = GE_ST01_post_csub;
GE_TP2_post_csub = GE_ST02_post_csub;
GE_TP3_post_csub = GE_ST03_post_csub;
GE_HC1_post_csub = GE_HC01_post_csub;
GE_HC2_post_csub = GE_HC02_post_csub;
GE_HC3_post_csub = GE_HC03_post_csub;

clear GE_ST01_pre_corti GE_ST02_pre_corti GE_ST03_pre_corti GE_ST01_pre_subcorti GE_ST02_pre_subcorti GE_ST03_pre_subcorti GE_ST01_pre_csub GE_ST02_pre_csub GE_ST03_pre_csub
clear GE_ST01_post_corti GE_ST02_post_corti GE_ST03_post_corti GE_ST01_post_subcorti GE_ST02_post_subcorti GE_ST03_post_subcorti GE_ST01_post_csub GE_ST02_post_csub GE_ST03_post_csub
clear GE_HC01_pre_corti GE_HC02_pre_corti GE_HC03_pre_corti GE_HC01_pre_subcorti GE_HC02_pre_subcorti GE_HC03_pre_subcorti GE_HC01_pre_csub GE_HC02_pre_csub GE_HC03_pre_csub
clear GE_HC01_post_corti GE_HC02_post_corti GE_HC03_post_corti GE_HC01_post_subcorti GE_HC02_post_subcorti GE_HC03_post_subcorti GE_HC01_post_csub GE_HC02_post_csub GE_HC03_post_csub

%% calculate means each group
mean_GE_TP1_pre_corti = mean(GE_TP1_pre_corti,1); 
mean_GE_TP2_pre_corti = mean(GE_TP2_pre_corti,1);
mean_GE_TP3_pre_corti = mean(GE_TP3_pre_corti,1);
mean_GE_TP1_pre_subcorti = mean(GE_TP1_pre_subcorti,1);
mean_GE_TP2_pre_subcorti = mean(GE_TP2_pre_subcorti,1);
mean_GE_TP3_pre_subcorti = mean(GE_TP3_pre_subcorti,1);
mean_GE_TP1_pre_csub = mean(GE_TP1_pre_csub,1);
mean_GE_TP2_pre_csub = mean(GE_TP2_pre_csub,1);
mean_GE_TP3_pre_csub = mean(GE_TP3_pre_csub,1);

mean_GE_HC1_pre_corti = mean(GE_HC1_pre_corti,1);
mean_GE_HC2_pre_corti = mean(GE_HC2_pre_corti,1);
mean_GE_HC3_pre_corti = mean(GE_HC3_pre_corti,1);
mean_GE_HC1_pre_subcorti = mean(GE_HC1_pre_subcorti,1);
mean_GE_HC2_pre_subcorti = mean(GE_HC2_pre_subcorti,1);
mean_GE_HC3_pre_subcorti = mean(GE_HC3_pre_subcorti,1);
mean_GE_HC1_pre_csub = mean(GE_HC1_pre_csub,1);
mean_GE_HC2_pre_csub = mean(GE_HC2_pre_csub,1);
mean_GE_HC3_pre_csub = mean(GE_HC3_pre_csub,1);

mean_GE_TP1_post_corti = mean(GE_TP1_post_corti,1);
mean_GE_TP2_post_corti = mean(GE_TP2_post_corti,1);
mean_GE_TP3_post_corti = mean(GE_TP3_post_corti,1);
mean_GE_TP1_post_subcorti = mean(GE_TP1_post_subcorti,1);
mean_GE_TP2_post_subcorti = mean(GE_TP2_post_subcorti,1);
mean_GE_TP3_post_subcorti = mean(GE_TP3_post_subcorti,1);
mean_GE_TP1_post_csub = mean(GE_TP1_post_csub,1);
mean_GE_TP2_post_csub = mean(GE_TP2_post_csub,1);
mean_GE_TP3_post_csub = mean(GE_TP3_post_csub,1);

mean_GE_HC1_post_corti = mean(GE_HC1_post_corti,1);
mean_GE_HC2_post_corti = mean(GE_HC2_post_corti,1);
mean_GE_HC3_post_corti = mean(GE_HC3_post_corti,1);
mean_GE_HC1_post_subcorti = mean(GE_HC1_post_subcorti,1);
mean_GE_HC2_post_subcorti = mean(GE_HC2_post_subcorti,1);
mean_GE_HC3_post_subcorti = mean(GE_HC3_post_subcorti,1);
mean_GE_HC1_post_csub = mean(GE_HC1_post_csub,1);
mean_GE_HC2_post_csub = mean(GE_HC2_post_csub,1);
mean_GE_HC3_post_csub = mean(GE_HC3_post_csub,1);

% calculate std each group
std_GE_TP1_pre_corti = std(GE_TP1_pre_corti,1); 
std_GE_TP2_pre_corti = std(GE_TP2_pre_corti,1);
std_GE_TP3_pre_corti = std(GE_TP3_pre_corti,1);
std_GE_TP1_pre_subcorti = std(GE_TP1_pre_subcorti,1);
std_GE_TP2_pre_subcorti = std(GE_TP2_pre_subcorti,1);
std_GE_TP3_pre_subcorti = std(GE_TP3_pre_subcorti,1);
std_GE_TP1_pre_csub = std(GE_TP1_pre_csub,1);
std_GE_TP2_pre_csub = std(GE_TP2_pre_csub,1);
std_GE_TP3_pre_csub = std(GE_TP3_pre_csub,1);

std_GE_HC1_pre_corti = std(GE_HC1_pre_corti,1);
std_GE_HC2_pre_corti = std(GE_HC2_pre_corti,1);
std_GE_HC3_pre_corti = std(GE_HC3_pre_corti,1);
std_GE_HC1_pre_subcorti = std(GE_HC1_pre_subcorti,1);
std_GE_HC2_pre_subcorti = std(GE_HC2_pre_subcorti,1);
std_GE_HC3_pre_subcorti = std(GE_HC3_pre_subcorti,1);
std_GE_HC1_pre_csub = std(GE_HC1_pre_csub,1);
std_GE_HC2_pre_csub = std(GE_HC2_pre_csub,1);
std_GE_HC3_pre_csub = std(GE_HC3_pre_csub,1);

std_GE_TP1_post_corti = std(GE_TP1_post_corti,1);
std_GE_TP2_post_corti = std(GE_TP2_post_corti,1);
std_GE_TP3_post_corti = std(GE_TP3_post_corti,1);
std_GE_TP1_post_subcorti = std(GE_TP1_post_subcorti,1);
std_GE_TP2_post_subcorti = std(GE_TP2_post_subcorti,1);
std_GE_TP3_post_subcorti = std(GE_TP3_post_subcorti,1);
std_GE_TP1_post_csub = std(GE_TP1_post_csub,1);
std_GE_TP2_post_csub = std(GE_TP2_post_csub,1);
std_GE_TP3_post_csub = std(GE_TP3_post_csub,1);

std_GE_HC1_post_corti = std(GE_HC1_post_corti,1);
std_GE_HC2_post_corti = std(GE_HC2_post_corti,1);
std_GE_HC3_post_corti = std(GE_HC3_post_corti,1);
std_GE_HC1_post_subcorti = std(GE_HC1_post_subcorti,1);
std_GE_HC2_post_subcorti = std(GE_HC2_post_subcorti,1);
std_GE_HC3_post_subcorti = std(GE_HC3_post_subcorti,1);
std_GE_HC1_post_csub = std(GE_HC1_post_csub,1);
std_GE_HC2_post_csub = std(GE_HC2_post_csub,1);
std_GE_HC3_post_csub = std(GE_HC3_post_csub,1);

% plot mean corti TP3
% bluecol = [167/255, 109/255, 254/255];
% bluecol2 = [117/255, 76/255 178/255];
% bluecol3 = [67/255, 44/255, 102/255];

bluecol = [168/255,109/255,254/255];
bluecol2 = [168/255,109/255,254/255];
bluecol3 = [168/255,109/255,254/255];

faceAlpha = 0.2;

xaxis = 1:1:numel(propthr);

figure
p1 = plot(mean_GE_TP3_pre_corti,'Color',bluecol,'LineWidth',2) % TP3 pre
hold on
fill([xaxis fliplr(xaxis)], [(mean_GE_TP3_pre_corti + std_GE_TP3_pre_corti), fliplr(mean_GE_TP3_pre_corti - std_GE_TP3_pre_corti)],bluecol,'FaceAlpha', faceAlpha,'EdgeColor','none') % TP3 pre
p2 = plot(mean_GE_TP3_post_corti,'Linestyle','--','Color',bluecol,'LineWidth',2) % TP3 post 
fill([xaxis fliplr(xaxis)], [(mean_GE_TP3_post_corti + std_GE_TP3_post_corti), fliplr(mean_GE_TP3_post_corti - std_GE_TP3_post_corti)],bluecol,'FaceAlpha', faceAlpha,'EdgeColor','none') % TP3 post 
p3 = plot(mean_GE_HC3_pre_corti,'Color','k','LineWidth',2)
fill([xaxis fliplr(xaxis)], [(mean_GE_HC3_pre_corti + std_GE_HC3_pre_corti), fliplr(mean_GE_HC3_pre_corti - std_GE_HC3_pre_corti)],'k','FaceAlpha', faceAlpha,'EdgeColor','none') % TP3 post 
p4 = plot(mean_GE_HC3_post_corti,'Linestyle','--','Color','k','LineWidth',2)
fill([xaxis fliplr(xaxis)], [(mean_GE_HC3_post_corti + std_GE_HC3_post_corti), fliplr(mean_GE_HC3_post_corti - std_GE_HC3_post_corti)],'k','FaceAlpha', faceAlpha,'EdgeColor','none') % TP3 post 
legend ([p1 p3 p2 p4],{'Stroke TP3 pre', 'HC pre','Stroke TP3 post','HC post'},'Location','southeast')
xticklabels(propthr)
xlabel('Edge density','FontWeight','bold')
ylabel('Global efficiency','FontWeight','bold')
ylim([0.35 0.85])
title('Cortical MCA territory','FontWeight','bold')
plotname = ['HCvsTP3_cortical_GE'];
ylim([0.35 0.85])
set(gca,'FontSize',13)
saveas(gcf,plotname)
print(gcf, '-dtiff', [plotname '.tiff']);

% plot mean subcorti TP3
figure
p1 = plot(mean_GE_TP3_pre_subcorti,'Color',bluecol2,'LineWidth',2) % TP3 pre
hold on
fill([xaxis fliplr(xaxis)], [(mean_GE_TP3_pre_subcorti + std_GE_TP3_pre_subcorti), fliplr(mean_GE_TP3_pre_subcorti - std_GE_TP3_pre_subcorti)],bluecol2,'FaceAlpha', faceAlpha,'EdgeColor','none') % TP3 pre
p2 = plot(mean_GE_TP3_post_subcorti,'Linestyle','--','Color',bluecol2,'LineWidth',2) % TP3 post 
fill([xaxis fliplr(xaxis)], [(mean_GE_TP3_post_subcorti + std_GE_TP3_post_subcorti), fliplr(mean_GE_TP3_post_subcorti - std_GE_TP3_post_subcorti)],bluecol2,'FaceAlpha', faceAlpha,'EdgeColor','none') % TP3 post 
p3 = plot(mean_GE_HC3_pre_subcorti,'Color','k','LineWidth',2)
fill([xaxis fliplr(xaxis)], [(mean_GE_HC3_pre_subcorti + std_GE_HC3_pre_subcorti), fliplr(mean_GE_HC3_pre_subcorti - std_GE_HC3_pre_subcorti)],'k','FaceAlpha', faceAlpha,'EdgeColor','none') % TP3 post 
p4 = plot(mean_GE_HC3_post_subcorti,'Linestyle','--','Color','k','LineWidth',2)
fill([xaxis fliplr(xaxis)], [(mean_GE_HC3_post_subcorti + std_GE_HC3_post_subcorti), fliplr(mean_GE_HC3_post_subcorti - std_GE_HC3_post_subcorti)],'k','FaceAlpha', faceAlpha,'EdgeColor','none') % TP3 post 
legend ([p1 p3 p2 p4],{'Stroke TP3 pre', 'HC pre','Stroke TP3 post','HC post'},'Location','southeast')
xticklabels(propthr)
xlabel('Edge density','FontWeight','bold')
ylabel('Global efficiency','FontWeight','bold')
ylim([0.35 0.85])
title('Subcortical MCA territory','FontWeight','bold')
plotname = ['HCvsTP3_subcortical_GE'];
ylim([0.35 0.85])
set(gca,'FontSize',13)
saveas(gcf,plotname)
print(gcf, '-dtiff', [plotname '.tiff']);

% plot mean cortico-subcorti
figure
p1 = plot(mean_GE_TP3_pre_csub,'Color',bluecol3,'LineWidth',2) % TP3 pre
hold on
fill([xaxis fliplr(xaxis)], [(mean_GE_TP3_pre_csub + std_GE_TP3_pre_csub), fliplr(mean_GE_TP3_pre_csub - std_GE_TP3_pre_csub)],bluecol3,'FaceAlpha', faceAlpha,'EdgeColor','none') % TP3 pre
p2 = plot(mean_GE_TP3_post_csub,'Linestyle','--','Color',bluecol3,'LineWidth',2) % TP3 post 
fill([xaxis fliplr(xaxis)], [(mean_GE_TP3_post_csub + std_GE_TP3_post_csub), fliplr(mean_GE_TP3_post_csub - std_GE_TP3_post_csub)],bluecol3,'FaceAlpha', faceAlpha,'EdgeColor','none') % TP3 post 
p3 = plot(mean_GE_HC3_pre_csub,'Color','k','LineWidth',2)
fill([xaxis fliplr(xaxis)], [(mean_GE_HC3_pre_csub + std_GE_HC3_pre_csub), fliplr(mean_GE_HC3_pre_csub - std_GE_HC3_pre_csub)],'k','FaceAlpha', faceAlpha,'EdgeColor','none') % TP3 post 
p4 = plot(mean_GE_HC3_post_csub,'Linestyle','--','Color','k','LineWidth',2)
fill([xaxis fliplr(xaxis)], [(mean_GE_HC3_post_csub + std_GE_HC3_post_csub), fliplr(mean_GE_HC3_post_csub - std_GE_HC3_post_csub)],'k','FaceAlpha', faceAlpha,'EdgeColor','none') % TP3 post 
legend ([p1 p3 p2 p4],{'Stroke TP3 pre', 'HC pre','Stroke TP3 post','HC post'},'Location','southeast')
xticklabels(propthr)
xlabel('Edge density','FontWeight','bold')
ylabel('Global efficiency','FontWeight','bold')
ylim([0.35 0.85])
title('Cortico-subcortical MCA territory','FontWeight','bold')
plotname = ['HCvsTP3_csub_GE'];
ylim([0.35 0.85])
set(gca,'FontSize',13)
saveas(gcf,plotname)
print(gcf, '-dtiff', [plotname '.tiff']);


