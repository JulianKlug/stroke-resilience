
% illustrate random failure like in Lo et al PNAS 2015
% May2020; Mitsouko van Assche

clear all

% set paths
path1 = '/Users/jk1/unige_onedrive/OneDrive - unige.ch/BCT/attacks240/attack_ST01';
data1 = 'rnd_attack_degrees_ST01_05-12-2020 17-46v4.mat';
% data1 = 'trg_attack_meandegrees_ST01_05-06-2020 16-15v4.mat';

path2 = '/Users/jk1/unige_onedrive/OneDrive - unige.ch/BCT/attacks240/attack_ST02';
data2 = 'rnd_attack_degrees_ST02_05-12-2020 20-16v4.mat';
% data2 = 'trg_attack_meandegrees_ST02_05-06-2020 18-35v4.mat';

path3 = '/Users/jk1/unige_onedrive/OneDrive - unige.ch/BCT/attacks240/attack_ST03';
data3 = 'rnd_attack_degrees_ST03_05-13-2020 18-19v4.mat';
% data3 = 'trg_attack_meandegrees_ST03_05-07-2020 00-36v4.mat';

pathhc = '/Users/jk1/unige_onedrive/OneDrive - unige.ch/BCT/attacks240/attack_HC';
datahc = 'rnd_attack_degrees_HC_05-12-2020 16-06v4';
% datahc = 'trg_attack_meandegrees_HC_05-06-2020 14-14v4';

fpath = '/Users/jk1/stroke_research/resilience_stroke/attack_analysis/figures/random_attack';
nroi = 240;
npatients = 15;
nhc = 16;
%
cd(path1)
load(data1)
GE1 = GlobEff_bin_new;
cd(path2)
load(data2)
GE2 = GlobEff_bin_new;
cd(path3)
load(data3)
GE3 = GlobEff_bin_new;
cd(pathhc)
load(datahc)
GEhc = GlobEff_bin_new;
cd(fpath)
clearvars -except GE1 GE2 GE3 GEhc data1 data2 data3 datahc nroi npatients nhc
%
% remove thresholds 0, 0.1, 0.2 as there are unconnected graphs
propthr = 0.3:0.1:1;
for attacknb = 1:nroi
    itername = ['iter' num2str(attacknb)];
    thisGE = GE1.(itername);
    GEff1.(itername) = thisGE(:,4:11);
    thisGE = GE2.(itername);
    GEff2.(itername) = thisGE(:,4:11);
    thisGE = GE3.(itername);
    GEff3.(itername) = thisGE(:,4:11);
    thisGE = GEhc.(itername);
    GEffhc.(itername) = thisGE(:,4:11);
end
% basic parameters
% for attacknb = 1:nroi
%     itername = ['iter' num2str(attacknb)];
%     meanGE1(attacknb,:) = mean(GEff1.(itername),1);
%     stdGE1(attacknb,:) = std(GEff1.(itername),1);
%     meanGE2(attacknb,:) = mean(GEff2.(itername),1);
%     stdGE2(attacknb,:) = std(GEff2.(itername),1);
%     meanGE3(attacknb,:) = mean(GEff3.(itername),1);
%     stdGE3(attacknb,:) = std(GEff3.(itername),1);
%     meanGEhc(attacknb,:) = mean(GEffhc.(itername),1);
%     stdGEhc(attacknb,:) = std(GEffhc.(itername),1);
% end
%% normalized version !
% for attacknb = 1:nroi 
%     itername = ['iter' num2str(attacknb)];
%     normGEff1.(itername) = GEff1.(itername) ./ GEff1.iter1;
%     normGEff2.(itername) = GEff2.(itername) ./ GEff2.iter1;
%     normGEff3.(itername) = GEff3.(itername) ./ GEff3.iter1;
%     normGEffhc.(itername) = GEffhc.(itername) ./ GEffhc.iter1;
% end
% % un-normalized version !
for attacknb = 1:nroi 
    itername = ['iter' num2str(attacknb)];
    normGEff1.(itername) = GEff1.(itername);
    normGEff2.(itername) = GEff2.(itername);
    normGEff3.(itername) = GEff3.(itername);
    normGEffhc.(itername) = GEffhc.(itername);
end
%rearrange matrices
for attacknb = 1:nroi
    for s = 1:npatients
        for p = 1:numel(propthr)
            itername = ['iter' num2str(attacknb)];
            normGEff1bis(s,attacknb,p) = normGEff1.(itername)(s,p);
            normGEff2bis(s,attacknb,p) = normGEff2.(itername)(s,p);
            normGEff3bis(s,attacknb,p) = normGEff3.(itername)(s,p);
        end
    end
end
for attacknb = 1:nroi
    for s = 1:nhc
        for p = 1:numel(propthr)
            itername = ['iter' num2str(attacknb)];
            normGEffhcbis(s,attacknb,p) = normGEffhc.(itername)(s,p);
        end
    end
end
%AUC per subjects
for s = 1:npatients
    for p = 1:numel(propthr)
        temp = cumtrapz(normGEff1bis(s,:,p));
        areaGE1(s,p) = temp(end)/100;
        temp = cumtrapz(normGEff2bis(s,:,p));
        areaGE2(s,p) = temp(end)/100;
        temp = cumtrapz(normGEff3bis(s,:,p));
        areaGE3(s,p) = temp(end)/100;
    end
end
for s = 1:nhc
    for p = 1:numel(propthr)
        temp = cumtrapz(normGEffhcbis(s,:,p));
        areaGEhc(s,p) = temp(end)/100;
    end
end
%std AUC and mean AUC
for p = 1:numel(propthr)
    std_areaGE1(p) = std(areaGE1(:,p),1);
    std_areaGE2(p) = std(areaGE2(:,p),1);
    std_areaGE3(p) = std(areaGE3(:,p),1);
    std_areaGEhc(p) = std(areaGEhc(:,p),1);
    mean_areaGE1(p) = mean(areaGE1(:,p),1);
    mean_areaGE2(p) = mean(areaGE2(:,p),1);
    mean_areaGE3(p) = mean(areaGE3(:,p),1);
    mean_areaGEhc(p) = mean(areaGEhc(:,p),1);
end

% prepare plots
SUP1 = mean_areaGE1 + std_areaGE1;
LOW1 = mean_areaGE1 - std_areaGE1;
SUP2 = mean_areaGE2 + std_areaGE2;
LOW2 = mean_areaGE2 - std_areaGE2; 
SUP3 = mean_areaGE3 + std_areaGE3;
LOW3 = mean_areaGE3 - std_areaGE3;
SUP = mean_areaGEhc + std_areaGEhc;
LOW = mean_areaGEhc - std_areaGEhc;
xaxis = 1:1:numel(propthr);

% original colors
% tp1_color = [0.227, 0.51, 0.29];
% tp2_color = [0.929 0.463 0.024];
% tp3_color = [0.929 0.024 0.643];

% % rocket color palette
% tp1_color = [0.38092887, 0.12061482, 0.32506528];
% tp2_color = [0.7965014, 0.10506637, 0.31063031];
% tp3_color = [0.95922872, 0.53307513, 0.3748895];

% ibm color palette 
tp1_color = [4/255,155/255,154/255];
tp2_color = [1/255,45/255,152/255];
tp3_color = [168/255,109/255,254/255];

faceAlpha = 0.05;

% plots
figure
plot1 = plot(mean_areaGE1,'Color',tp1_color,'LineWidth',2)
hold on
fill([xaxis fliplr(xaxis)], [SUP1, fliplr(LOW1)],tp1_color,'FaceAlpha',faceAlpha,'EdgeColor','none')
plothc = plot(mean_areaGEhc,'Color','k','LineWidth',2)
fill([xaxis fliplr(xaxis)], [SUP, fliplr(LOW)],'k','FaceAlpha',faceAlpha,'EdgeColor','none')
legend([plot1 plothc],{'Stroke TP1','HC'},'Location','northwest');
ylim([0.47 0.67])
ylabel('Global efficiency','FontWeight','bold')
xlabel('Edge density','FontWeight','bold')
xticklabels(propthr)
% title('Targeted attack')
title('Random attack','FontWeight','bold')
% plotname = ['FDA_GlobEff_degrees_AUC_TP1'];
plotname = ['FDA_GlobEff_rnd_degrees_AUC_TP1'];
set(gca,'FontSize',13)
saveas(gcf,plotname)
print(gcf, '-dtiff', [plotname '.tiff']);

figure
plot2 = plot(mean_areaGE2,'Color',tp2_color,'LineWidth',2)
hold on
fill([xaxis fliplr(xaxis)], [SUP2, fliplr(LOW2)],tp2_color,'FaceAlpha',faceAlpha,'EdgeColor','none')
plothc = plot(mean_areaGEhc,'Color','k','LineWidth',2)
fill([xaxis fliplr(xaxis)], [SUP, fliplr(LOW)],'k','FaceAlpha',faceAlpha,'EdgeColor','none')
legend([plot2 plothc],{'Stroke TP2','HC'},'Location','northwest');
ylim([0.47 0.67])
ylabel('Global efficiency','FontWeight','bold')
xlabel('Edge density','FontWeight','bold')
xticklabels(propthr)
% title('Targeted attack')
title('Random attack','FontWeight','bold')
% plotname = ['FDA_GlobEff_degrees_AUC_TP2'];
plotname = ['FDA_GlobEff_rnd_degrees_AUC_TP2'];
set(gca,'FontSize',13)
saveas(gcf,plotname)
print(gcf, '-dtiff', [plotname '.tiff']);

figure
plot3 = plot(mean_areaGE3,'Color',tp3_color,'LineWidth',2)
hold on
fill([xaxis fliplr(xaxis)], [SUP3, fliplr(LOW3)],tp3_color,'FaceAlpha',faceAlpha,'EdgeColor','none')
plothc = plot(mean_areaGEhc,'Color','k','LineWidth',2)
fill([xaxis fliplr(xaxis)], [SUP, fliplr(LOW)],'k','FaceAlpha',faceAlpha,'EdgeColor','none')
legend([plot3 plothc],{'Stroke TP3','HC'},'Location','northwest');
ylim([0.47 0.67])
ylabel('Global efficiency','FontWeight','bold')
xlabel('Edge density','FontWeight','bold')
xticklabels(propthr)
% title('Targeted attack')
title('Random attack','FontWeight','bold')
% plotname = ['FDA_GlobEff_degrees_AUC_TP3'];
plotname = ['FDA_GlobEff_rnd_degrees_AUC_TP3'];
set(gca,'FontSize',13)
saveas(gcf,plotname)
print(gcf, '-dtiff', [plotname '.tiff']);

figure
plot1 = plot(mean_areaGE1,'Color',tp1_color,'LineWidth',2)
hold on
fill([xaxis fliplr(xaxis)], [SUP1, fliplr(LOW1)],tp1_color,'FaceAlpha',faceAlpha,'EdgeColor','none')
plot2 = plot(mean_areaGE2,'Color',tp2_color,'LineWidth',2)
fill([xaxis fliplr(xaxis)], [SUP2, fliplr(LOW2)],tp2_color,'FaceAlpha',faceAlpha,'EdgeColor','none')
legend([plot1 plot2],{'TP1','TP2'},'Location','northwest');
% ylim([50 65])
ylabel('Global efficiency')
xlabel('Connection density')
xticklabels(propthr)
% title('Targeted attack')
title('Random attack')
% plotname = ['FDA_GlobEff_degrees_AUC_TP1TP2'];
plotname = ['FDA_GlobEff_rnd_degrees_AUC_TP1TP2'];
saveas(gcf,plotname)
print(gcf, '-dtiff', [plotname '.tiff']);

figure
plot2 = plot(mean_areaGE1,'Color',tp2_color,'LineWidth',2)
hold on
fill([xaxis fliplr(xaxis)], [SUP2, fliplr(LOW2)],tp2_color,'FaceAlpha',faceAlpha,'EdgeColor','none')
plot3 = plot(mean_areaGE2,'Color',tp3_color,'LineWidth',2)
fill([xaxis fliplr(xaxis)], [SUP3, fliplr(LOW3)],tp3_color,'FaceAlpha',faceAlpha,'EdgeColor','none')
legend([plot2 plot3],{'TP2','TP3'},'Location','northwest');
% ylim([50 65])
ylabel('Global efficiency')
xlabel('Connection density')
xticklabels(propthr)
% title('Targeted attack')
title('Random attack')
% plotname = ['FDA_GlobEff_degrees_AUC_TP2TP3'];
plotname = ['FDA_GlobEff_rnd_degrees_AUC_TP2TP3'];
saveas(gcf,plotname)
print(gcf, '-dtiff', [plotname '.tiff']);


    

