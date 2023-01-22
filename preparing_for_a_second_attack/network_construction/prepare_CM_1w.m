%prepare_CM_1

%1: based on prepa_CM_2.m to retrieve CMs used for previous analyses,
%   threshold matrices with function W = threshold_absolute(W, thr) using
%   thr=0, threshold matrices with function W = threshold_proportional(W,
%   p) using 0 <= p <= 1 by steps of 0.1, and additionally fix any possible 
%   problems within matrices using weight_conversion.m with 'autofix'
%   option
%2: calls density_prop_mitsu.m to compute network density 
%3: calls efficiency_prop_mitsu.m to compute global efficiency
%4: calls clustering_coeff_prop_mitsu.m to compute CC
%5: calls distance_prop_mitsu.m to compute char path length


%   à voir: utliser la methode de supekar et al 2008 et choisir threshold 
%   matrices with function W = threshold_proportional(W,p) using 0.01 <= p <= 0.99
%   by steps of 0.01

clear all
addpath('D:\STROKE\MRIdata\BCT\2017_01_15_BCT');
BCTpath=('D:\STROKE\MRIdata\BCT\atlas_BNA\BNA_28_flipped_N32_retroicor_SBB4_prop_bin_window\ST');
NBSpath=('D:\STROKE\MRIdata\NBStoolbox\BNA_28_flipped_N32_retroicor\SBB4');

%% set!!!!!!!!!!!
session={'session01','session02','session03'};
TP={'ST01','ST02','ST03'};

%% retrieve CMs used for previous analyses
% session 01
for sess=1:numel(session)
    cd(fullfile(NBSpath,session{sess}))
    txtfiles=dir('P*');
    for t=1:numel(txtfiles)
        myfile=load(txtfiles(t).name);
        if sess==1
            CM3D28_norm_sbb4.ST01(:,:,t)=myfile;
        elseif sess==2
            CM3D28_norm_sbb4.ST02(:,:,t)=myfile;
        elseif sess==3
            CM3D28_norm_sbb4.ST03(:,:,t)=myfile;
        end
    end
end
cd(BCTpath)

%% thresholds matrices with function W = threshold_absolute(W, thr)
for tp=1:numel(TP)
    for s=1:15 
        myCM=CM3D28_norm_sbb4.(TP{tp})(:,:,s);
        myCM_abs_thres = threshold_absolute(myCM,0); % W_thr = threshold_absolute(W, thr);
        CM3D28_norm_sbb4.(TP{tp})(:,:,s)=myCM_abs_thres;% diag=0 and <0 elements removed
    end
end
save('CM3D28_norm_sbb4_posval.mat','CM3D28_norm_sbb4');%

%% thresholding matrices with function W = threshold_proportional(W, p); 50%
%  with 0 <= p <= 1 by steps of 0.1:
thr1=0:0.1:1;
thr2=0.1:0.1:1;
for tp=1:numel(TP)
    for s=1:15 
        for t=1:numel(thr2)
            myCM=CM3D28_norm_sbb4.(TP{tp})(:,:,s);
            lowbin = thr1(t);
            uppbin = thr2(t); 
            binname=['bin' num2str(t)];
            CM_win.(TP{tp}).(binname)(:,:,s)=threshold_bybin(myCM,lowbin);
            sumofweights(t) = sum(sum(CM_win.(TP{tp}).(binname)(:,:,s)))/2;
            nbofnodes(t) = numel(find(CM_win.(TP{tp}).(binname)(:,:,s)))/2;
        end
        figure
        plot(sumofweights)
        title (['subject ' num2str(s) ' TP' num2str(tp)])
        fname = ['winCM_subject' num2str(s) '_TP' num2str(tp)];
        saveas(gcf,fname)
        Ssumofweights(s,:) = sumofweights;
        Snbofnodes(s,:) = nbofnodes;
    end
    finalsumofweights.(TP{tp}) = Ssumofweights;
    finalnbofnodes.(TP{tp})= Snbofnodes;
end
save('CM3D28_norm_sbb4_posval_win.mat','CM3D28_norm_sbb4','CM_win','finalsumofweights','finalnbofnodes')

%% binarize matrices with function W = weight_conversion(W,'binarize');
for tp=1:numel(TP)
    for s=1:15
        for t=1:numel(thr2)
            lowbin = thr1(t);
            uppbin = thr2(t); 
            binname=['bin' num2str(t)];
            myCM=CM_win.(TP{tp}).(binname)(:,:,s);
            myCM_bin = weight_conversion(myCM,'binarize');
            CM_win_bin.(TP{tp}).(binname)(:,:,s) = myCM_bin;
        end
    end
end
save('CM3D28_norm_sbb4_posval_binwin.mat','CM_win_bin');%

%% additionally: fix any possible problems within matrices using weight_conversion.m
for tp=1:numel(TP)
    for s=1:15
        for t=1:numel(thr2)
            binname=['bin' num2str(t)];
            myCM_bin=CM_win_bin.(TP{tp}).(binname)(:,:,s);
            myCM_bin=weight_conversion(myCM_bin,'autofix');
            CM_win_bin.(TP{tp}).(binname)(:,:,s)=myCM_bin;
        end
    end
end
save('CM3D28_norm_sbb4_posval_binwin.mat','CM_win_bin')