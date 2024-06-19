function runMMIVA
% % Instructions to run this script:
% % Set current folder to be the MISA folder, e.g.:
% cd . % cd('~/MISA')
% % Add the folder containing this file to the path, e.g.:
% addpath('./scripts/tests')
% % Clear all and RUN:
% close all; clear all; clc
% % run this script
comps=30;
P=0;
addpath('/trdapps/linux-x86_64/matlab/toolboxes/spm12')
addpath(genpath('/trdapps/linux-x86_64/matlab/toolboxes/GroupICATv4.0b'))
addpath('./tests')
format compact
maxNumCompThreads(12);

MISAdir = '/data/users2/eswar/temp/MISA_analysis/MISA';
addpath(genpath(MISAdir));
UKBdataroot = '/data/users2/eswar/temp/MISA_analysis/UKB/derivatives';
outpath = '/data/users2/eswar/temp/MISA_analysis/UKB/results/hyb_/MIVA2step';%'/data/users2/rsilva/results/real/MISA_multimodal_eswar';
initpath=fullfile('/data/users2/eswar/temp/MISA_analysis/UKB/','results');
mask = '/data/users2/eswar/temp/MISA_analysis/UKB/amask_T1_3mm.nii';
slist = load('/data/users2/eswar/temp/MISA_analysis/UKB/scripts/goodID_UKBiobank_3497MISAfinal.txt');
nS = length(slist);
tab = readtable('/data/mialab/competition2019/UKBiobank/packNship/results/scores/new/ukb_unaffected_clean_select_log_with_sex.tab','Delimiter','\t','FileType','text');
tab = tab(ismember(tab.eid,slist),:);
tab.sex_f31_0_0 = double(strcmpi(tab.sex_f31_0_0,'Female'));
tab.genetic_sex_f22001_0_0 = double(strcmpi(tab.genetic_sex_f22001_0_0,'Female'));

% Load Data
load /data/users2/eswar/temp/MISA_analysis/UKB/scripts/ALFFdata_gmMaskNEW.mat % loads ALFFdata2
load /data/users2/eswar/temp/MISA_analysis/UKB/scripts/SMRIdata_gmMaskNEW.mat % loads SMRIdata2
load /data/users2/eswar/temp/MISA_analysis/UKB/scripts/DMRIdata_Sm6_wmMaskNew.mat
%/data/users2/eswar/temp/MISA_analysis/UKB/scripts/DMRIdata_anisoSm_wmMaskNew.mat


load('siteInfo_UKB.mat');
load('subID_2907.mat');
sitevecUKB_2907 = sitevecUKB_3497(lt5NanID);
tab = tab(lt5NanID,:);


% sMRI and ALFF mask:
Vol_sf = spm_vol('/data/users2/eswar/temp/MISA_analysis/UKB/gmMask_TPM_thrp2_fractc8_3mm.nii');
sf_3D = logical(spm_read_vols(Vol_sf));
sf_sz = size(sf_3D);

% Detect bad data
SMRI_drop = detect_bad_data(SMRIdata'); % Takes 10-15 seconds
ALFF_drop = detect_bad_data(ALFFdata'); % Takes 10-15 seconds

% dMRI mask:
Vol_d = spm_vol('/data/users2/eswar/temp/MISA_analysis/UKB/wmMask_TPM_thrp4_fractc8_3mm.nii');
d_3D = logical(spm_read_vols(Vol_d));
d_sz = size(d_3D);

% Detect bad data
DMRI_drop = detect_bad_data(DMRIdata'); % Takes 10-15 seconds
DMRIdata = DMRIdata(~DMRI_drop.bad_cols(:),:);
d_3D(d_3D(:)) = ~DMRI_drop.bad_cols(:);

% Apply variance normalization (per subject), and remove mean across subjects:
old_X = {SMRIdata(:,lt5NanID),ALFFdata(:,lt5NanID),DMRIdata(:,lt5NanID)};

[tmp, SMRI_vn_info] = do_variance_normalization(SMRIdata(:,lt5NanID)');
[tmp, SMRI_old_mean] = remove_mean(tmp);
SMRIdata = tmp';
[tmp, ALFF_vn_info] = do_variance_normalization(ALFFdata(:,lt5NanID)');
[tmp, ALFF_old_mean] = remove_mean(tmp);
ALFFdata = tmp';
[tmp, DMRI_vn_info] = do_variance_normalization(DMRIdata(:,lt5NanID)');
[tmp, DMRI_old_mean] = remove_mean(tmp);
DMRIdata = tmp';
vn_old_mean = {SMRI_vn_info.old_mean',ALFF_vn_info.old_mean',DMRI_vn_info.old_mean'};
vn_old_std = {SMRI_vn_info.old_std',ALFF_vn_info.old_std',DMRI_vn_info.old_std'};
old_mean = {SMRI_old_mean',ALFF_old_mean',DMRI_old_mean'};
% ALFFdata = remove_mean(do_variance_normalization(ALFFdata'))';
% DMRIdata = remove_mean(do_variance_normalization(DMRIdata'))';

%load('siteInfo_UKB.mat');
%load('subID_2907.mat');
%sitevecUKB_2907 = sitevecUKB_3497(lt5NanID);

Xs = [ones(length(sitevecUKB_2907),1) sitevecUKB_2907];
bh1 = Xs\SMRIdata';
SMRIdata = (SMRIdata' - Xs*bh1)';
bh2 = Xs\ALFFdata'; 
ALFFdata = (ALFFdata' - Xs*bh2)';
bh3 = Xs\DMRIdata';
DMRIdata = (DMRIdata' - Xs*bh3)';

% Define the number of datasets (here, the number of modalities)
X = {SMRIdata,ALFFdata,DMRIdata};
M = 1:length(X);
ut = utils;

%%%%% Preprocessing STARTS HERE
% comps = 10; P = 0;
rec_type = 'WT';

% [whtM, dwhtM] = doMMPCA(X);
[whtM, H] = doMMGPCA(X, comps, 'WT');
mmgpca1 = setup_MISA_MMIVA1(M,X,whtM,'kotz');

% Y = cellfun(@mtimes,whtM,X,'Un',0);

% figure
% imagesc(abs(corr(H',sum(cat(3,Y{:}),3)')),[0 1]);
% colorbar; colormap(hot); axis equal tight
% saveas(gcf,'test.png')

% t = H' - sum(cat(3,Y{:}),3)';
% sum(t(:).^2)
% sort(diag(corr(H',sum(cat(3,Y{:}),3)'))')
% figure
% plot(H',sum(cat(3,Y{:}),3)','.')
% axis equal tight
% saveas(gcf,'test.png')

% tt = H' ./ sum(cat(3,Y{:}),3)';
% figure
% h = histogram(tt(:));       
% h.Normalization = 'probability';
% h.BinWidth = 1e-3;
% saveas(gcf,'test.png')

gica1 = setup_MISA_MMIVA1(1,{H},{diag(1.813799364234218 ./ std(H,[],2))*eye(comps)},'kotz');

% Run GICA using Infomax
rng(202)
[W1,wht]=icatb_runica(H,'weights',gica1.W{1},'ncomps',comps,'sphering', 'off', 'verbose', 'off', 'posact', 'off', 'bias', 'on');
std_W1 = std(W1*H,[],2); % Ignoring wht because Infomax run with 'sphering' 'off' --> wht = eye(comps)
W1 = diag(1.813799364234218 ./ std_W1) * W1;

% RUN GICA using MISA: continuing from Infomax above...
%% Could use stochastic optimization, but not doing so because MISA does not implement bias weights (yet)...
% gica1.stochastic_opt('verbose', 'off', 'weights', gica1.W{1}, 'bias', 'off');%, 'block', 1100);
[wout,fval,exitflag,output] = run_MISA(gica1,{W1});
std_gica1_W1 = std(gica1.Y{1},[],2);
gica1.objective(ut.stackW({diag(1.813799364234218 ./ std_gica1_W1)*gica1.W{1}}));

% Compare Infomax and MISA: both are GICA components...
figure('color',[1 1 1],'visible','off')
imagesc(abs(corr((W1*H)',(gica1.W{1}*H)')),[0 1])
colorbar; colormap(hot); axis equal tight
saveas(gcf,fullfile(outpath,sprintf('C%03d_GICAcorr_Infomax_vs_MISA_P%03d_N2907.png',comps,P)))

% Combine MISA GICA with whitening matrices to initialize multimodal model
W = cellfun(@(w) gica1.W{1}*w,whtM,'Un',0);
W = cellfun(@(w,x) diag(1.813799364234218 ./ std(w*x,[],2))*w,W,X,'Un',0);
mmgica1 = setup_MISA_MMIVA1(M,X,W,'kotz');

%%%%% Preprocessing ENDS HERE

cs = comps;%min(comps,10);

%%%%% Multimodal IVA on reduced data STARTS HERE

% Y = cellfun(@(y) gica1.W{1}*y,Y,'Un',0);
Y = mmgica1.Y;
M = 1:length(Y);
nComp = comps;% size(Y{1},1);

% Identity initial W
for mm = M
    W0{mm} = [eye(nComp),zeros(nComp,size(Y{M(1)},1)-nComp)];
    % [u, s, v] = svd(randn(nComp,size(Y{M(1)},1)));
    % W0{mm} = u*s*v';
end

% Add noise based on P
seed = 100;
W0 = init_W_noise(P,W0,Y,seed);

mmiva1 = setup_MISA_MMIVA1(M,Y,W0,'kotz');

[wout,fval,exitflag,output] = run_MISA(mmiva1,W0(M));

% wout = ut.stackW(mmiva1.W);
% mmiva1.objective(wout);
min_r = diag(min(cat(3,abs(corr(mmiva1.Y{1}',mmiva1.Y{2}')),...
                        abs(corr(mmiva1.Y{2}',mmiva1.Y{3}')),...
                        abs(corr(mmiva1.Y{1}',mmiva1.Y{3}'))),[],3))';
sort(min_r) % Print top correlations

% Visualize cross-modal correlations:
figure('color',[1 1 1],'visible','off')
subplot(2,3,1)
imagesc(abs(corr(mmiva1.Y{1}',mmiva1.Y{2}')),[0 .5]), colormap(hot); colorbar
axis equal tight
subplot(2,3,2)
imagesc(abs(corr(mmiva1.Y{2}',mmiva1.Y{3}')),[0 .5]), colormap(hot); colorbar
axis equal tight
subplot(2,3,3)
imagesc(abs(corr(mmiva1.Y{1}',mmiva1.Y{3}')),[0 .5]), colormap(hot); colorbar
axis equal tight
subplot(2,3,4)
imagesc(abs(corr(mmiva1.Y{1}',mmiva1.Y{1}')),[0 .5]), colormap(hot); colorbar
axis equal tight
subplot(2,3,5)
imagesc(abs(corr(mmiva1.Y{2}',mmiva1.Y{2}')),[0 .5]), colormap(hot); colorbar
axis equal tight
subplot(2,3,6)
imagesc(abs(corr(mmiva1.Y{3}',mmiva1.Y{3}')),[0 .5]), colormap(hot); colorbar
axis equal tight
saveas(gcf,fullfile(outpath,sprintf('C%03d_mmiva1_corr_P%03d_N2907.png',comps,P)))

RRtab_v = abs(corr(tab{1:end,2:end},[mmiva1.Y{1}',mmiva1.Y{2}',mmiva1.Y{3}'],'Rows','pairwise','type','spearman'));
figure
imagesc(RRtab_v,[0 .5])
colorbar; colormap(hot); axis equal tight
saveas(gcf,sprintf('C%03d_corr_demog_overall_mmiva1_P%03d_N2907.png',nComp,P))
[a,ix_v] = sort(max(RRtab_v,[],2)','descend');
tab.Properties.VariableNames{ix_v(1:cs)+1} % Print top correlated demographics
% tab(1:2,ix_v(1:cs)+1)
% a(1:cs)

figure('color',[1 1 1],'visible','off')
for pick = 1:cs
    [~,ixm] = max(RRtab_v(:,1:nComp),[],2);
    subplot(1,3,1)
    plot(tab{:,ix_v(pick)+1},mmiva1.Y{1}(ixm(ix_v(pick)),:)','.')
    ylabel(sprintf('Component %d',ixm(ix_v(pick))))
    % xlabel(tab.Properties.VariableNames{ix_v(pick)+1})
    title(sprintf('sMRI Corr: %.2f',RRtab_v(ix_v(pick),ixm(ix_v(pick)))))
    axis square
    subplot(1,3,2)
    plot(tab{:,ix_v(pick)+1},mmiva1.Y{2}(ixm(ix_v(pick)),:)','.')
    % ylabel(sprintf('Component %d',ixm(ix_v(pick))))
    xlabel(tab.Properties.VariableNames{ix_v(pick)+1})
    title(sprintf('fMRI Corr: %.2f',RRtab_v(ix_v(pick),nComp+ixm(ix_v(pick)))))
    axis square
    subplot(1,3,3)
    plot(tab{:,ix_v(pick)+1},mmiva1.Y{3}(ixm(ix_v(pick)),:)','.')
    % ylabel(sprintf('Component %d',ixm(ix_v(pick))))
    % xlabel(tab.Properties.VariableNames{ix_v(pick)+1})
    title(sprintf('dMRI Corr: %.2f',RRtab_v(ix_v(pick),2*nComp+ixm(ix_v(pick)))))
    axis square
    saveas(gcf,fullfile(outpath,sprintf('C%03d_corr_demog_pick_top%03d_mmiva1_P%03d_N2907.png',nComp,pick,P)))
end

%%%%% Multimodal IVA on reduced data ENDS HERE



%%%%% Multimodal IVA on FULL data STARTS HERE

[wout,fval,exitflag,output] = run_MISA(mmgica1,mmgica1.W);

% wout = ut.stackW(mmgica1.W);
% mmgica1.objective(wout);
min_r = diag(min(cat(3,abs(corr(mmgica1.Y{1}',mmgica1.Y{2}')),...
                        abs(corr(mmgica1.Y{2}',mmgica1.Y{3}')),...
                        abs(corr(mmgica1.Y{1}',mmgica1.Y{3}'))),[],3))';
sort(min_r)
% Visualize cross-modal correlations:
figure('color',[1 1 1],'visible','off')
subplot(2,3,1)
imagesc(abs(corr(mmgica1.Y{1}',mmgica1.Y{2}')),[0 .5]), colormap(hot); colorbar
axis equal tight
subplot(2,3,2)
imagesc(abs(corr(mmgica1.Y{2}',mmgica1.Y{3}')),[0 .5]), colormap(hot); colorbar
axis equal tight
subplot(2,3,3)
imagesc(abs(corr(mmgica1.Y{1}',mmgica1.Y{3}')),[0 .5]), colormap(hot); colorbar
axis equal tight
subplot(2,3,4)
imagesc(abs(corr(mmgica1.Y{1}',mmgica1.Y{1}')),[0 .5]), colormap(hot); colorbar
axis equal tight
subplot(2,3,5)
imagesc(abs(corr(mmgica1.Y{2}',mmgica1.Y{2}')),[0 .5]), colormap(hot); colorbar
axis equal tight
subplot(2,3,6)
imagesc(abs(corr(mmgica1.Y{3}',mmgica1.Y{3}')),[0 .5]), colormap(hot); colorbar
axis equal tight
saveas(gcf,fullfile(outpath,sprintf('C%03d_mmgica1_corr_P%03d_N2907.png',comps,P)))

RR = max(cat(3,abs(corr(mmgica1.Y{1}',mmiva1.Y{1}')),abs(corr(mmgica1.Y{2}',mmiva1.Y{2}')),abs(corr(mmgica1.Y{3}',mmiva1.Y{3}'))),[],3);
ix_match = ut.munkres(1-abs(RR));
figure('color',[1 1 1],'visible','off')
imagesc(RR(:,ix_match),[0 1])
colorbar; colormap(hot); axis equal tight
saveas(gcf,fullfile(outpath,sprintf('C%03d_mmgica1_vs_mmiva1_P%03d_N2907.png',comps,P)))

RRtab_g = abs(corr(tab{1:end,2:end},[mmgica1.Y{1}',mmgica1.Y{2}',mmgica1.Y{3}'],'Rows','pairwise','type','spearman'));
figure('color',[1 1 1],'visible','off')
imagesc(RRtab_g,[0 .5])
colorbar; colormap(hot); axis equal tight
saveas(gcf,fullfile(outpath,sprintf('C%03d_corr_demog_overall_mmgica1_P%03d_N2907.png',nComp,P)))
[a,ix_g] = sort(max(RRtab_g,[],2)','descend');
tab.Properties.VariableNames{ix_g(1:cs)+1} % Print top correlated demographics
% tab(1,ix(1:cs)+1)
% a(1:cs)

% figure
% for pick = 1:cs
%     [~,ixm] = max(RRtab_g(:,1:nComp),[],2);
%     plot(tab{:,ix_g(pick)+1},mmgica1.Y{1}(ixm(ix_g(pick)),:)','.')
%     saveas(gcf,sprintf('C%03d_corr_demog_pick_top%03d_mmgica1_P%03d.png',nComp,pick,P))
% end
figure('color',[1 1 1],'visible','off')
for pick = 1:cs
    [~,ixm] = max(RRtab_g(:,1:nComp),[],2);
    subplot(1,3,1)
    plot(tab{:,ix_g(pick)+1},mmgica1.Y{1}(ixm(ix_g(pick)),:)','.')
    ylabel(sprintf('Component %d',ixm(ix_g(pick))))
    % xlabel(tab.Properties.VariableNames{ix_g(pick)+1})
    title(sprintf('sMRI Corr: %.2f',RRtab_g(ix_g(pick),ixm(ix_g(pick)))))
    axis square
    subplot(1,3,2)
    plot(tab{:,ix_g(pick)+1},mmgica1.Y{2}(ixm(ix_g(pick)),:)','.')
    % ylabel(sprintf('Component %d',ixm(ix_g(pick))))
    xlabel(tab.Properties.VariableNames{ix_g(pick)+1})
    title(sprintf('fMRI Corr: %.2f',RRtab_g(ix_g(pick),nComp+ixm(ix_g(pick)))))
    axis square
    subplot(1,3,3)
    plot(tab{:,ix_g(pick)+1},mmgica1.Y{3}(ixm(ix_g(pick)),:)','.')
    % ylabel(sprintf('Component %d',ixm(ix_g(pick))))
    % xlabel(tab.Properties.VariableNames{ix_g(pick)+1})
    title(sprintf('dMRI Corr: %.2f',RRtab_g(ix_g(pick),2*nComp+ixm(ix_g(pick)))))
    axis square
    saveas(gcf,fullfile(outpath,sprintf('C%03d_corr_demog_pick_top%03d_mmgica1_P%03d_N2907.png',nComp,pick,P)))
end

% Match reduced and full IVA results: some DEMOG. correlations go up, some go down
figure('color',[1 1 1],'visible','off')
for cc = 1:nComp
    subplot(2*nComp/10,5,cc)
    plot(RRtab_v(:,ix_match(cc)),RRtab_g(:,cc),'.r')
    hold on
    mi = min([RRtab_v(:,ix_match(cc));RRtab_g(:,cc)]);
    ma = max([RRtab_v(:,ix_match(cc));RRtab_g(:,cc)]);
    plot([mi ma],[mi ma],'--k')
    axis equal tight
    % set(gca,'XLim',[0 1],'YLim',[0 1])
end
saveas(gcf,fullfile(outpath,sprintf('C%03d_demog_corr_reduced_vs_full_matched_scatter_P%03d_N2907.png',nComp,P)))

%%%%% Multimodal IVA on FULL data ENDS HERE


old_W = W;

% Things to save:
% GPCA
step1 = mmgpca1;
% GICA
step2 = gica1;

% Reduced data approach
step3_reduced = mmiva1;
final_W = cellfun(@(w2,w1) w2 * w1, step3_reduced.W, old_W,'Un',0);
Y_adj = cellfun(@(w,x) w*x, final_W, old_X, 'Un', 0);
% Y_adj2 = cellfun(@(y,w,mu) y + repmat(w*mu,1,size(y,2)), step3_reduced.Y, final_W, old_mean, 'Un', 0);
% Y_adj2 = cellfun(@(y,w,mu,st) y.*repmat(st,size(y,1),1) + w*repmat(mu,size(w,2),1), Y_adj2, final_W, vn_old_mean, vn_old_std, 'Un', 0);

% plot((final_W{1}*old_X{1})',((step3_reduced.Y{1}+repmat(final_W{1}*old_mean{1},1,size(Dm_2,1))).*repmat(old_std_row,1,size(final_W{1},1))'+ final_W{1}*repmat(old_mean_row',size(final_W{1},2),1))','.')
% figure
% plot(Y_adj{1}',Y_adj2{1}','.')
% axis equal tight
% saveas(gcf,'test.png')

W = final_W;
A_tr = cellfun(@(w) w', W, 'Un', 0);
A_pinv = cellfun(@(w) ((w*w')\w)', W, 'Un', 0);
A_ls1 = cellfun(@(y,x) ((y*y')\(y*x'))', step3_reduced.Y, X, 'Un', 0);
A_ls2 = cellfun(@(y,x) ((y*y')\(y*x'))', Y_adj, X, 'Un', 0);
A_ls1_ox = cellfun(@(y,x) ((y*y')\(y*x'))', step3_reduced.Y, old_X, 'Un', 0);
A_ls2_ox = cellfun(@(y,x) ((y*y')\(y*x'))', Y_adj, old_X, 'Un', 0);
gpcasig = step1.Y;
gicasig = step2.Y;
icasig = step3_reduced.Y;

% Erase data from MISA objects
%tmpX = step1.X;
%tmpY = step1.Y;
%step1.X = {}; % --> these are variance normalized and mean removed, not original values
step1.Y = {};
step2.X = {};
step2.Y = {};
step3_reduced.X = {}; % --> these use the same Y as step1
step3_reduced.Y = {};

fname = sprintf('UKBMRI_MMIVAreduced_w0GICA_C%03d_results_goodSub2907_new_P%03d',step3_reduced.C(1),P);
save(fullfile(outpath,[fname '.mat']), ...
        'W','A_tr','A_pinv','A_ls1','A_ls2','A_ls1_ox','A_ls2_ox','gpcasig','gicasig','icasig','step1','step2','step3_reduced','vn_old_mean','vn_old_std','old_mean','outpath','-v7.3')


% Full data approach
step3_full = mmgica1;
final_W = step3_full.W;
Y_adj = cellfun(@(w,x) w*x, final_W, old_X, 'Un', 0);
% Y_adj2 = cellfun(@(y,w,mu) y + repmat(w*mu,1,size(y,2)), step3_full.Y, final_W, old_mean, 'Un', 0);
% Y_adj2 = cellfun(@(y,w,mu,st) y.*repmat(st,size(y,1),1) + w*repmat(mu,size(w,2),1), Y_adj2, final_W, vn_old_mean, vn_old_std, 'Un', 0);

% plot((final_W{1}*old_X{1})',((step3_full.Y{1}+repmat(final_W{1}*old_mean{1},1,size(Dm_2,1))).*repmat(old_std_row,1,size(final_W{1},1))'+ final_W{1}*repmat(old_mean_row',size(final_W{1},2),1))','.')
% figure
% plot(Y_adj{1}',Y_adj2{1}','.')
% axis equal tight
% saveas(gcf,'test.png')

W = final_W;
A_tr = cellfun(@(w) w', W, 'Un', 0);
A_pinv = cellfun(@(w) ((w*w')\w)', W, 'Un', 0);
A_ls1 = cellfun(@(y,x) ((y*y')\(y*x'))', step3_full.Y, X, 'Un', 0);
A_ls2 = cellfun(@(y,x) ((y*y')\(y*x'))', Y_adj, X, 'Un', 0);
A_ls1_ox = cellfun(@(y,x) ((y*y')\(y*x'))', step3_full.Y, old_X, 'Un', 0);
A_ls2_ox = cellfun(@(y,x) ((y*y')\(y*x'))', Y_adj, old_X, 'Un', 0);
gpcasig = step1.Y;
gicasig = step2.Y;
icasig = step3_full.Y;

% Erase data from MISA objects
%tmpX = step1.X;
%tmpY = step1.Y;
%step1.X = {}; % --> these are variance normalized and mean removed, not original values
% step1.Y = {};
% step2.X = {};
% step2.Y = {};
step3_full.X = {}; % --> these use the same X as step1
step3_full.Y = {};

fname = sprintf('UKBMRI_MMIVAfull_w0GICA_C%03d_results_goodSub2907_new_P%03d',step3_full.C(1),P);
save(fullfile(outpath,[fname '.mat']), ...
        'W','A_tr','A_pinv','A_ls1','A_ls2','A_ls1_ox','A_ls2_ox','gpcasig','gicasig','icasig','step1','step2','step3_full','vn_old_mean','vn_old_std','old_mean','outpath','-v7.3')

