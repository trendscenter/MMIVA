function runMMIVA_allHCSZ
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
BSNIPdataroot = '/data/users2/eswar/temp/MISA_analysis/HCSZ/BSNIP';
outpath = '/data/users2/eswar/temp/MISA_analysis/HCSZ/results/Combined/MIVA2step';%'/data/users2/rsilva/results/real/MISA_multimodal_eswar';
initpath=fullfile('/data/users2/eswar/temp/MISA_analysis/UKB/','results');
%mask = '/data/users2/eswar/temp/MISA_analysis/UKB/amask_T1_3mm.nii';
%slist = load('/data/users2/eswar/temp/MISA_analysis/UKB/scripts/goodID_UKBiobank_3497MISAfinal.txt');

%nS = length(slist);
%tab = readtable('/data/mialab/competition2019/UKBiobank/packNship/results/scores/new/ukb_unaffected_clean_select_log_with_sex.tab','Delimiter','\t','FileType','text');
%tab = tab(ismember(tab.eid,slist),:);
%tab.sex_f31_0_0 = double(strcmpi(tab.sex_f31_0_0,'Female'));
%tab.genetic_sex_f22001_0_0 = double(strcmpi(tab.genetic_sex_f22001_0_0,'Female'));

% Load Data

FB1=load('/data/users2/eswar/temp/MISA_analysis/HCSZ/scripts/SMRIdata_FBIRN_newMask.mat'); % loads SMRIdata2
FB2=load('/data/users2/eswar/temp/MISA_analysis/HCSZ/scripts/ALFFdata_FBIRN_newMask.mat'); % loads ALFFdata2
FB3=load('/data/users2/eswar/temp/MISA_analysis/HCSZ/scripts/DMRIdata_FBIRN_Sm6_newMask.mat');

CB1=load('/data/users2/eswar/temp/MISA_analysis/HCSZ/scripts/SMRIdata_COBRE_newMask.mat'); % loads SMRIdata2
CB2=load('/data/users2/eswar/temp/MISA_analysis/HCSZ/scripts/ALFFdata_COBRE_newMask.mat'); % loads ALFFdata2
CB3=load('/data/users2/eswar/temp/MISA_analysis/HCSZ/scripts/DMRIdata_COBRE_Sm6_newMask.mat');

MP1=load('/data/users2/eswar/temp/MISA_analysis/HCSZ/scripts/SMRIdata_MPRC_newMask.mat'); % loads SMRIdata2
MP2=load('/data/users2/eswar/temp/MISA_analysis/HCSZ/scripts/ALFFdata_MPRC_newMask.mat'); % loads ALFFdata2
MP3=load('/data/users2/eswar/temp/MISA_analysis/HCSZ/scripts/DMRIdata_MPRC_Sm6_newMask.mat');

BS1=load('/data/users2/eswar/temp/MISA_analysis/HCSZ/scripts/SMRIdata_BSNIP_newMask.mat'); % loads SMRIdata2
BS2=load('/data/users2/eswar/temp/MISA_analysis/HCSZ/scripts/ALFFdata_BSNIP_newMask.mat'); % loads ALFFdata2
BS3=load('/data/users2/eswar/temp/MISA_analysis/HCSZ/scripts/DMRIdata_BSNIP_Sm6_newMask.mat');


% sMRI and ALFF mask:
Vol_sf = spm_vol('/data/users2/eswar/temp/MISA_analysis/UKB/gmMask_TPM_thrp2_fractc8_3mm.nii');
sf_3D = logical(spm_read_vols(Vol_sf));
sf_sz = size(sf_3D);

% Detect bad data
SMRI_drop_FB = detect_bad_data(FB1.SMRIdata'); % Takes 10-15 seconds
ALFF_drop_FB = detect_bad_data(FB2.ALFFdata'); % Takes 10-15 seconds
SMRI_drop_CB = detect_bad_data(CB1.SMRIdata'); % Takes 10-15 seconds
ALFF_drop_CB = detect_bad_data(CB2.ALFFdata'); % Takes 10-15 seconds
SMRI_drop_MP = detect_bad_data(MP1.SMRIdata'); % Takes 10-15 seconds
ALFF_drop_MP = detect_bad_data(MP2.ALFFdata'); % Takes 10-15 seconds
SMRI_drop_BS = detect_bad_data(BS1.SMRIdata'); % Takes 10-15 seconds
ALFF_drop_BS = detect_bad_data(BS2.ALFFdata'); % Takes 10-15 seconds

% dMRI mask:
Vol_d = spm_vol('/data/users2/eswar/temp/MISA_analysis/UKB/wmMask_TPM_thrp4_fractc8_3mm.nii');
d_3D = logical(spm_read_vols(Vol_d));
d_sz = size(d_3D);

% Detect bad data
DMRI_drop_FB = detect_bad_data(FB3.DMRIdata'); % Takes 10-15 seconds
DMRI_drop_CB = detect_bad_data(CB3.DMRIdata'); % Takes 10-15 seconds
DMRI_drop_MP = detect_bad_data(MP3.DMRIdata'); % Takes 10-15 seconds
DMRI_drop_BS = detect_bad_data(BS3.DMRIdata'); % Takes 10-15 seconds
%DMRIdata = DMRIdata(~DMRI_drop.bad_cols(:),:);
%d_3D(d_3D(:)) = ~DMRI_drop.bad_cols(:);

%%
meanSMRI = mean([FB1.SMRIdata CB1.SMRIdata MP1.SMRIdata BS1.SMRIdata],2);
rFB1=corr(FB1.SMRIdata,meanSMRI); 
rCB1=corr(CB1.SMRIdata,meanSMRI); 
rMP1=corr(MP1.SMRIdata,meanSMRI); 
rBS1=corr(BS1.SMRIdata,meanSMRI); 

meanALFF = mean([FB2.ALFFdata CB2.ALFFdata MP2.ALFFdata BS2.ALFFdata],2);
rFB2=corr(FB2.ALFFdata,meanALFF);
rCB2=corr(CB2.ALFFdata,meanALFF);
rMP2=corr(MP2.ALFFdata,meanALFF);
rBS2=corr(BS2.ALFFdata,meanALFF); 

meanDMRI = mean([FB3.DMRIdata CB3.DMRIdata MP3.DMRIdata BS3.DMRIdata],2);
rFB3=corr(FB3.DMRIdata,meanDMRI);
rCB3=corr(CB3.DMRIdata,meanDMRI);
rMP3=corr(MP3.DMRIdata,meanDMRI);
rBS3=corr(BS3.DMRIdata,meanDMRI); 

rALL1 = [rFB1;rCB1;rMP1;rBS1];rStat1=[min(rALL1) mean(rALL1) median(rALL1) mode(rALL1) std(rALL1)];
rALL2 = [rFB2;rCB2;rMP2;rBS2];rStat2=[min(rALL2) mean(rALL2) median(rALL2) mode(rALL2) std(rALL2)];
rALL3 = [rFB3;rCB3;rMP3;rBS3];rStat3=[min(rALL3) mean(rALL3) median(rALL3) mode(rALL3) std(rALL3)];

rStatFB1 = [min(rFB1) mean(rFB1) median(rFB1) std(rFB1)];
rStatCB1 = [min(rCB1) mean(rCB1) median(rCB1) std(rCB1)];
rStatMP1 = [min(rMP1) mean(rMP1) median(rMP1) std(rMP1)];
rStatBS1 = [min(rBS1) mean(rBS1) median(rBS1) std(rBS1)];

rStatFB2 = [min(rFB2) mean(rFB2) median(rFB2) std(rFB2)];
rStatCB2 = [min(rCB2) mean(rCB2) median(rCB2) std(rCB2)];
rStatMP2 = [min(rMP2) mean(rMP2) median(rMP2) std(rMP2)];
rStatBS2 = [min(rBS2) mean(rBS2) median(rBS2) std(rBS2)];

rStatFB3 = [min(rFB3) mean(rFB3) median(rFB3) std(rFB3)];
rStatCB3 = [min(rCB3) mean(rCB3) median(rCB3) std(rCB3)];
rStatMP3 = [min(rMP3) mean(rMP3) median(rMP3) std(rMP3)];
rStatBS3 = [min(rBS3) mean(rBS3) median(rBS3) std(rBS3)];

meanSMRI_3D = to_vol(meanSMRI',sf_3D);
meanALFF_3D = to_vol(meanALFF',sf_3D);
meanDMRI_3D = to_vol(meanDMRI',d_3D);
addpath('/data/users2/eswar/software/utilities/')
%o1='meanSMRI_allHCSZ.nii';o2 ='meanALFF_allHCSZ.nii';o3='meanDMRI_allHCSZ.nii';
%create_4DNiftifile(fullfile(pwd,'meannii_newmask','SMRI',o1),meanSMRI_3D,Vol_d.mat);
%create_4DNiftifile(fullfile(pwd,'meannii_newmask','DMRI',o2),meanALFF_3D,Vol_d.mat);
%create_4DNiftifile(fullfile(pwd,'meannii_newmask','DMRI',o2),meanALFF_3D,Vol_d.mat);


nstd = 2;
outID_FB1 = find(rFB1 < rStatFB1(3)-nstd*rStatFB1(4));
outID_FB2 = find(rFB2 < rStatFB2(3)-nstd*rStatFB2(4));
outID_FB3 = find(rFB3 < rStatFB3(3)-nstd*rStatFB3(4));

outID_CB1 = find(rCB1 < rStatCB1(3)-nstd*rStatCB1(4));
outID_CB2 = find(rCB2 < rStatCB2(3)-nstd*rStatCB2(4));
outID_CB3 = find(rCB3 < rStatCB3(3)-nstd*rStatCB3(4));

outID_MP1 = find(rMP1 < rStatMP1(3)-nstd*rStatMP1(4));
outID_MP2 = find(rMP2 < rStatMP2(3)-nstd*rStatMP2(4));
outID_MP3 = find(rMP3 < rStatMP3(3)-nstd*rStatMP3(4));

outID_BS1 = find(rBS1 < rStatBS1(3)-nstd*rStatBS1(4));
outID_BS2 = find(rBS2 < rStatBS2(3)-nstd*rStatBS2(4));
outID_BS3 = find(rBS3 < rStatBS3(3)-nstd*rStatBS3(4));


outID_BSall  = unique([outID_BS1;outID_BS2;outID_BS3]);
outID_MPall  = unique([outID_MP1;outID_MP2;outID_MP3]);
outID_CBall  = unique([outID_CB1;outID_CB2;outID_CB3]);
outID_FBall  = unique([outID_FB1;outID_FB2;outID_FB3]);
keepID_FBall=setdiff(1:length(rFB1),outID_FBall);
keepID_CBall=setdiff(1:length(rCB1),outID_CBall);
keepID_BSall=setdiff(1:length(rBS1),outID_BSall); 
keepID_MPall=setdiff(1:length(rMP1),outID_MPall);

load('/data/users2/eswar/temp/MISA_analysis/HCSZ/stats/FBIRNdemo')
load('/data/users2/eswar/temp/MISA_analysis/HCSZ/stats/COBREdemo')
load('/data/users2/eswar/temp/MISA_analysis/HCSZ/stats/MPRCdemo')
load('/data/users2/eswar/temp/MISA_analysis/HCSZ/stats/BSNIPdemo')
[tf1,loc1] = ismember(FBDEMO.keepID,keepID_FBall);
FBfinalind = find(tf1==1);  FBfinalloc= loc1(loc1>0); %intersect(FBDEMO.keepID,keepID_FBall);
[tf2,loc2] = ismember(CBDEMO.keepID,keepID_CBall);
CBfinalind = find(tf2==1);CBfinalloc= loc2(loc2>0);%intersect(CBDEMO.keepID,keepID_CBall);
[tf3,loc3] = ismember(MPDEMO.keepID,keepID_MPall);
MPfinalind = find(tf3==1)';MPfinalloc= loc3(loc3>0);%intersect(MPDEMO.keepID,keepID_MPall);
[tf4,loc4] = ismember(BSDEMO.keepID,keepID_BSall);
BSfinalind = find(tf4==1);BSfinalloc= loc4(loc4>0);%intersect(BSDEMO.keepID,keepID_BSall);

SMRIdata = [FB1.SMRIdata(:,keepID_FBall(FBfinalloc))  CB1.SMRIdata(:,keepID_CBall(CBfinalloc)) MP1.SMRIdata(:,keepID_MPall(MPfinalloc)) BS1.SMRIdata(:,keepID_BSall(BSfinalloc)) ];
ALFFdata = [FB2.ALFFdata(:,keepID_FBall(FBfinalloc))  CB2.ALFFdata(:,keepID_CBall(CBfinalloc)) MP2.ALFFdata(:,keepID_MPall(MPfinalloc)) BS2.ALFFdata(:,keepID_BSall(BSfinalloc)) ];
DMRIdata = [FB3.DMRIdata(:,keepID_FBall(FBfinalloc))  CB3.DMRIdata(:,keepID_CBall(CBfinalloc)) MP3.DMRIdata(:,keepID_MPall(MPfinalloc)) BS3.DMRIdata(:,keepID_BSall(BSfinalloc)) ];

%SMRIdata = [FB1.SMRIdata(:,keepID_FBall)  CB1.SMRIdata(:,keepID_CBall) MP1.SMRIdata(:,keepID_MPall) BS1.SMRIdata(:,keepID_BSall) ];
%ALFFdata = [FB2.ALFFdata(:,keepID_FBall)  CB2.ALFFdata(:,keepID_CBall) MP2.ALFFdata(:,keepID_MPall) BS2.ALFFdata(:,keepID_BSall) ];
%DMRIdata = [FB3.DMRIdata(:,keepID_FBall)  CB3.DMRIdata(:,keepID_CBall) MP3.DMRIdata(:,keepID_MPall) BS3.DMRIdata(:,keepID_BSall) ];


%print(['nS = ', num2str(size(SMRIdata,2))])
%%

% Apply variance normalization (per subject), and remove mean across subjects:
old_X = {SMRIdata, ALFFdata,DMRIdata};

[tmp, SMRI_vn_info] = do_variance_normalization(SMRIdata');
[tmp, SMRI_old_mean] = remove_mean(tmp);
SMRIdata = tmp';
[tmp, ALFF_vn_info] = do_variance_normalization(ALFFdata');
[tmp, ALFF_old_mean] = remove_mean(tmp);
ALFFdata = tmp';
[tmp, DMRI_vn_info] = do_variance_normalization(DMRIdata');
[tmp, DMRI_old_mean] = remove_mean(tmp);
DMRIdata = tmp';
vn_old_mean = {SMRI_vn_info.old_mean',ALFF_vn_info.old_mean',DMRI_vn_info.old_mean'};
vn_old_std = {SMRI_vn_info.old_std',ALFF_vn_info.old_std',DMRI_vn_info.old_std'};
old_mean = {SMRI_old_mean',ALFF_old_mean',DMRI_old_mean'};
% ALFFdata = remove_mean(do_variance_normalization(ALFFdata'))';
% DMRIdata = remove_mean(do_variance_normalization(DMRIdata'))';

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
saveas(gcf,fullfile(outpath,sprintf('C%03d_GICAcorr_Infomax_vs_MISA_P%03d.png',comps,P)))

% Combine MISA GICA with whitening matrices to initialize multimodal model
W = cellfun(@(w) gica1.W{1}*w,whtM,'Un',0);
W = cellfun(@(w,x) diag(1.813799364234218 ./ std(w*x,[],2))*w,W,X,'Un',0);
mmgica1 = setup_MISA_MMIVA1(M,X,W,'kotz');

%%%%% Preprocessing ENDS HERE

cs = min(comps,10);

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
saveas(gcf,fullfile(outpath,sprintf('C%03d_mmiva1_corr_P%03d.png',comps,P)))

%RRtab_v = abs(corr(tab{1:end,2:end},[mmiva1.Y{1}',mmiva1.Y{2}',mmiva1.Y{3}'],'Rows','pairwise','type','spearman'));
%figure('color',[1 1 1],'visible','off')
%imagesc(RRtab_v,[0 .5])
%colorbar; colormap(hot); axis equal tight
%saveas(gcf,fullfile(outpath,sprintf('C%03d_corr_demog_overall_mmiva1_P%03d.png',nComp,P)))
%[a,ix_v] = sort(max(RRtab_v,[],2)','descend');
%tab.Properties.VariableNames{ix_v(1:cs)+1} % Print top correlated demographics
% tab(1:2,ix_v(1:cs)+1)
% a(1:cs)

%figure('color',[1 1 1],'visible','off')
%for pick = 1:cs
%    [~,ixm] = max(RRtab_v(:,1:nComp),[],2);
%    subplot(1,3,1)
%    plot(tab{:,ix_v(pick)+1},mmiva1.Y{1}(ixm(ix_v(pick)),:)','.')
%    ylabel(sprintf('Component %d',ixm(ix_v(pick))))
    % xlabel(tab.Properties.VariableNames{ix_v(pick)+1})
%    title(sprintf('sMRI Corr: %.2f',RRtab_v(ix_v(pick),ixm(ix_v(pick)))))
%    axis square
%    subplot(1,3,2)
%    plot(tab{:,ix_v(pick)+1},mmiva1.Y{2}(ixm(ix_v(pick)),:)','.')
    % ylabel(sprintf('Component %d',ixm(ix_v(pick))))
%    xlabel(tab.Properties.VariableNames{ix_v(pick)+1})
%    title(sprintf('fMRI Corr: %.2f',RRtab_v(ix_v(pick),nComp+ixm(ix_v(pick)))))
%    axis square
%    subplot(1,3,3)
%    plot(tab{:,ix_v(pick)+1},mmiva1.Y{3}(ixm(ix_v(pick)),:)','.')
    % ylabel(sprintf('Component %d',ixm(ix_v(pick))))
    % xlabel(tab.Properties.VariableNames{ix_v(pick)+1})
%    title(sprintf('dMRI Corr: %.2f',RRtab_v(ix_v(pick),2*nComp+ixm(ix_v(pick)))))
%    axis square
%    saveas(gcf,fullfile(outpath,sprintf('C%03d_corr_demog_pick_top%03d_mmiva1_P%03d.png',nComp,pick,P)))
%end

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
saveas(gcf,fullfile(outpath,sprintf('C%03d_mmgica1_corr_P%03d.png',comps,P)))

RR = max(cat(3,abs(corr(mmgica1.Y{1}',mmiva1.Y{1}')),abs(corr(mmgica1.Y{2}',mmiva1.Y{2}')),abs(corr(mmgica1.Y{3}',mmiva1.Y{3}'))),[],3);
ix_match = ut.munkres(1-abs(RR));
figure('color',[1 1 1],'visible','off')
imagesc(RR(:,ix_match),[0 1])
colorbar; colormap(hot); axis equal tight
saveas(gcf,fullfile(outpath,sprintf('C%03d_mmgica1_vs_mmiva1_P%03d.png',comps,P)))

%RRtab_g = abs(corr(tab{1:end,2:end},[mmgica1.Y{1}',mmgica1.Y{2}',mmgica1.Y{3}'],'Rows','pairwise','type','spearman'));
%figure('color',[1 1 1],'visible','off')
%imagesc(RRtab_g,[0 .5])
%colorbar; colormap(hot); axis equal tight
%saveas(gcf,fullfile(outpath,sprintf('C%03d_corr_demog_overall_mmgica1_P%03d.png',nComp,P)))
%[a,ix_g] = sort(max(RRtab_g,[],2)','descend');
%tab.Properties.VariableNames{ix_g(1:cs)+1} % Print top correlated demographics
% tab(1,ix(1:cs)+1)
% a(1:cs)

% figure
% for pick = 1:cs
%     [~,ixm] = max(RRtab_g(:,1:nComp),[],2);
%     plot(tab{:,ix_g(pick)+1},mmgica1.Y{1}(ixm(ix_g(pick)),:)','.')
%     saveas(gcf,sprintf('C%03d_corr_demog_pick_top%03d_mmgica1_P%03d.png',nComp,pick,P))
% end
%figure('color',[1 1 1],'visible','off')
%for pick = 1:cs
%    [~,ixm] = max(RRtab_g(:,1:nComp),[],2);
%    subplot(1,3,1)
%    plot(tab{:,ix_g(pick)+1},mmgica1.Y{1}(ixm(ix_g(pick)),:)','.')
%    ylabel(sprintf('Component %d',ixm(ix_g(pick))))
    % xlabel(tab.Properties.VariableNames{ix_g(pick)+1})
%    title(sprintf('sMRI Corr: %.2f',RRtab_g(ix_g(pick),ixm(ix_g(pick)))))
%    axis square
%    subplot(1,3,2)
%    plot(tab{:,ix_g(pick)+1},mmgica1.Y{2}(ixm(ix_g(pick)),:)','.')
    % ylabel(sprintf('Component %d',ixm(ix_g(pick))))
%    xlabel(tab.Properties.VariableNames{ix_g(pick)+1})
%    title(sprintf('fMRI Corr: %.2f',RRtab_g(ix_g(pick),nComp+ixm(ix_g(pick)))))
%    axis square
%    subplot(1,3,3)
%    plot(tab{:,ix_g(pick)+1},mmgica1.Y{3}(ixm(ix_g(pick)),:)','.')
    % ylabel(sprintf('Component %d',ixm(ix_g(pick))))
    % xlabel(tab.Properties.VariableNames{ix_g(pick)+1})
%    title(sprintf('dMRI Corr: %.2f',RRtab_g(ix_g(pick),2*nComp+ixm(ix_g(pick)))))
%    axis square
%    saveas(gcf,fullfile(outpath,sprintf('C%03d_corr_demog_pick_top%03d_mmgica1_P%03d.png',nComp,pick,P)))
%end

% Match reduced and full IVA results: some DEMOG. correlations go up, some go down
%figure('color',[1 1 1],'visible','off')
%for cc = 1:nComp
%    subplot(2,5,cc)
%    plot(RRtab_v(:,ix_match(cc)),RRtab_g(:,cc),'.r')
%    hold on
%    mi = min([RRtab_v(:,ix_match(cc));RRtab_g(:,cc)]);
%    ma = max([RRtab_v(:,ix_match(cc));RRtab_g(:,cc)]);
%    plot([mi ma],[mi ma],'--k')
%    axis equal tight
    % set(gca,'XLim',[0 1],'YLim',[0 1])
%end
%saveas(gcf,fullfile(outpath,sprintf('C%03d_demog_corr_reduced_vs_full_matched_scatter_P%03d.png',nComp,P)))

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

fname = sprintf('HCSZMRI_MMIVAreduced_w0GICA_C%03d_results_new_P%03d',step3_reduced.C(1),P);
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

fname = sprintf('HCSZMRI_MMIVAfull_w0GICA_C%03d_results_new_P%03d',step3_full.C(1),P);
save(fullfile(outpath,[fname '.mat']), ...
        'W','A_tr','A_pinv','A_ls1','A_ls2','A_ls1_ox','A_ls2_ox','gpcasig','gicasig','icasig','step1','step2','step3_full','vn_old_mean','vn_old_std','old_mean','outpath','-v7.3')

