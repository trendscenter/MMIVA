%clear all
%%
addpath('/home/eswar/Documents/MATLAB/utilities/')
addpath('/home/eswar/Documents/MATLAB/MANCOVAtb')
addpath('/home/eswar/Documents/MATLAB/spm12/');
addpath(genpath('/home/eswar/Documents/MATLAB/GroupICATv4.0b/'));
addpath(genpath('/home/eswar/Documents/MATLAB/MISA-master/'));
addpath('/home/eswar/Documents/MATLAB/export_fig')
addpath('/home/eswar/Documents/MATLAB/spm12/');

%%
FBroot = '/media/eswar/New Volume/Documents/MISA/FBIRN';
CBroot = '/media/eswar/New Volume/Documents/MISA/COBRE';
MProot = '/media/eswar/New Volume/Documents/MISA/MPRC';
BSroot = '/media/eswar/New Volume/Documents/MISA/BSNIP';

HCSZroot = '/media/eswar/New Volume/Documents/MISA/HCSZ';
load(fullfile(HCSZroot,'DEMOG_matchID.mat'));
load(fullfile(HCSZroot,'Data','keepID_HCSZ.mat'));

fulreg_HCSZ_UKBinit = load(fullfile(HCSZroot,'MMIVA2step_wUKBinit','C030_regsite_ukbinit','HCSZMRI_MMIVAfull_w0UKBinit2907_C030_results_new_P000_preregSite.mat'));
%fulreg_HCSZ_UKBinit = load(fullfile(HCSZroot,'MMIVA2step_wUKBinit','C030_regSite','HCSZMRI_MMIVAfull_w0UKBinit_C030_results_new_P000_preregSite.mat'));

%%
load(fullfile(FBroot,'FBIRNdemo.mat'));
load(fullfile(CBroot,'COBREdemo.mat'));
load(fullfile(MProot,'MPRCdemo.mat'));
load(fullfile(BSroot,'BSNIPdemo.mat'));

nFB = length(keepID_FBall(FBfinalloc));%length(FBDEMO.keepID) + length(FBDEMO.missID);
nCB = length(keepID_CBall(CBfinalloc));%length(CBDEMO.keepID) + length(CBDEMO.missID);
nMP = length(keepID_MPall(MPfinalloc));%length(MPDEMO.keepID) + length(MPDEMO.missID);
nBS = length(keepID_BSall(BSfinalloc));%length(BSDEMO.keepID) + length(BSDEMO.missID);


FBind = 1:nFB;
CBind = nFB+1:nFB+nCB;
MPind = nFB+nCB+1:nFB+nCB+nMP;
BSind = nFB+nCB+nMP+1:nFB+nCB+nMP+nBS;

%%
%fb_icasig_ukbinit_01 = fulreg_HCSZ_UKBinit.icasig{1}';%(:,FBind)';
%fb_icasig_ukbinit_01 = fulreg_HCSZ_UKBinit.icasig{1}';%(:,FBind)';

%fb_icasig_ukbinit_02 = fulreg_HCSZ_UKBinit.icasig{2}';%(:,FBind)';
%fb_icasig_ukbinit_02 = fulreg_HCSZ_UKBinit.icasig{2}';%(:,FBind)';

%fb_icasig_ukbinit_03 = fulreg_HCSZ_UKBinit.icasig{3}';%(:,FBind)';
%fb_icasig_ukbinit_03 = fulreg_HCSZ_UKBinit.icasig{3}';%(:,FBind)';


%cb_icasig_ukbinit_01 = fulreg_HCSZ_UKBinit.icasig{1}';%(:,CBind)';
%cb_icasig_ukbinit_01 = fulreg_HCSZ_UKBinit.icasig{1}';%(:,CBind)';

%cb_icasig_ukbinit_02 = fulreg_HCSZ_UKBinit.icasig{2}';%(:,CBind)';
%cb_icasig_ukbinit_02 = fulreg_HCSZ_UKBinit.icasig{2}';%(:,CBind)';

%cb_icasig_ukbinit_03 = fulreg_HCSZ_UKBinit.icasig{3}';%(:,CBind)';
%cb_icasig_ukbinit_03 = fulreg_HCSZ_UKBinit.icasig{3}';%(:,CBind)';

%mp_icasig_ukbinit_01 = fulreg_HCSZ_UKBinit.icasig{1}';%(:,MPind)';
%mp_icasig_ukbinit_01 = fulreg_HCSZ_UKBinit.icasig{1}';%(:,MPind)';

%mp_icasig_ukbinit_02 = fulreg_HCSZ_UKBinit.icasig{2}';%(:,MPind)';
%mp_icasig_ukbinit_02 = fulreg_HCSZ_UKBinit.icasig{2}';%(:,MPind)';

%mp_icasig_ukbinit_03 = fulreg_HCSZ_UKBinit.icasig{3}(:,MPind)';
%mp_icasig_ukbinit_03 = fulreg_HCSZ_UKBinit.icasig{3}(:,MPind)';

%bs_icasig_ukbinit_01 = fulreg_HCSZ_UKBinit.icasig{1}';%(:,BSind)';
%bs_icasig_ukbinit_01 = fulreg_HCSZ_UKBinit.icasig{1}';%(:,BSind)';

%bs_icasig_ukbinit_02 = fulreg_HCSZ_UKBinit.icasig{2}';%(:,BSind)';
%bs_icasig_ukbinit_02 = fulreg_HCSZ_UKBinit.icasig{2}';%(:,BSind)';

%bs_icasig_ukbinit_03 = fulreg_HCSZ_UKBinit.icasig{3}';%(:,BSind)';
%bs_icasig_ukbinit_03 = fulreg_HCSZ_UKBinit.icasig{3}';%(:,BSind)';

%%
%[tf1,loc1] = ismember(FBDEMO.keepID,keepID_FBall);
%FBfinalind = find(tf1==1);  FBfinalloc= loc1(loc1>0); %intersect(FBDEMO.keepID,keepID_FBall);
%[tf2,loc2] = ismember(CBDEMO.keepID,keepID_CBall);
%CBfinalind = find(tf2==1);CBfinalloc= loc2(loc2>0);%intersect(CBDEMO.keepID,keepID_CBall);
%[tf3,loc3] = ismember(MPDEMO.keepID,keepID_MPall);
%MPfinalind = find(tf3==1)';MPfinalloc= loc3(loc3>0);%intersect(MPDEMO.keepID,keepID_MPall);
%[tf4,loc4] = ismember(BSDEMO.keepID,keepID_BSall);
%BSfinalind = find(tf4==1);BSfinalloc= loc4(loc4>0);%intersect(BSDEMO.keepID,keepID_BSall);
%%
%[tf11,loc12] = ismember(keepID_FBall,FBDEMO.keepID);

%[tf21,loc22] = ismember(keepID_CBall,CBDEMO.keepID);

%[tf31,loc32] = ismember(keepID_MPall,MPDEMO.keepID);

%[tf41,loc42] = ismember(keepID_BSall,BSDEMO.keepID);

%all0finalind = find([tf11';tf21';tf31';tf41']==1);
%%
% %% FBIRN
% varnamesFB = FBDEMO.demo_col_names(1:4);
% % SMRI
% MODELHCSZ0_fb = make_design_matrix_HCSZ_fb(varnamesFB,FBDEMO.demo_extras(FBfinalind,1) ,{FBDEMO.demo_extravars{1}},FBDEMO.demo_extravars{1},FBfinalind);
% [DEMO0_FBIRN_UKBinit, MULT0_FBIRN_UKBinit, UNI0_FBIRN_UKBinit] = run_model_HCSZ(MODELHCSZ0_fb, 0.01,fb_icasig_ukbinit_01(FBfinalloc,:), [], [], 0);
% save FBIRN_UKBinit_preregSite_SMRI_Mancovan  *0_FBIRN_UKBinit
% 
% [DEMO0_FBIRN_UKBinit, MULT0_FBIRN_UKBinit, UNI0_FBIRN_UKBinit] = run_model_HCSZ(MODELHCSZ0_fb, 0.01,fb_icasig_ukbinit_01(FBfinalloc,:), [], [], 0);
% save FBIRN_UKBinit_preregSite_SMRI_Mancovan  *0_FBIRN_UKBinit
% 
% % ALFF
% MODELHCSZ1_fb = make_design_matrix_HCSZ_fb(varnamesFB,FBDEMO.demo_extras(FBfinalind,[2 4]) ,{FBDEMO.demo_extravars{2}, FBDEMO.demo_extravars{4}},{FBDEMO.demo_extravars(2),FBDEMO.demo_extravars(4)},FBfinalind);
% [DEMO1_FBIRN_UKBinit, MULT1_FBIRN_UKBinit, UNI1_FBIRN_UKBinit] = run_model_HCSZ(MODELHCSZ1_fb, 0.01,fb_icasig_ukbinit_02(FBfinalloc,:), [], [], 0);
% save FBIRN_UKBinit_preregSite_ALFF_Mancovan  *1_FBIRN_UKBinit
% 
% [DEMO1_FBIRN_UKBinit, MULT1_FBIRN_UKBinit, UNI1_FBIRN_UKBinit] = run_model_HCSZ(MODELHCSZ1_fb, 0.01,fb_icasig_ukbinit_02(FBfinalloc,:), [], [], 0);
% save FBIRN_UKBinit_preregSite_ALFF_Mancovan  *1_FBIRN_UKBinit
% 
% % DMRI
% MODELHCSZ2_fb = make_design_matrix_HCSZ_fb(varnamesFB,FBDEMO.demo_extras(FBfinalind,3) ,{FBDEMO.demo_extravars{3}},FBDEMO.demo_extravars{3},FBfinalind);
% [DEMO2_FBIRN_UKBinit, MULT2_FBIRN_UKBinit, UNI2_FBIRN_UKBinit] = run_model_HCSZ(MODELHCSZ2_fb, 0.01,fb_icasig_ukbinit_02(FBfinalloc,:), [], [], 0);
% save FBIRN_UKBinit_preregSite_DMRI_Mancovan  *2_FBIRN_UKBinit
% 
% [DEMO2_FBIRN_UKBinit, MULT2_FBIRN_UKBinit, UNI2_FBIRN_UKBinit] = run_model_HCSZ(MODELHCSZ2_fb, 0.01,fb_icasig_ukbinit_02(FBfinalloc,:), [], [], 0);
% save FBIRN_UKBinit_preregSite_DMRI_Mancovan  *2_FBIRN_UKBinit
% 
% 
% %%  COBRE
% 
% varnamesCB = CBDEMO.demo_col_names(1:3);
% 
% % SMRI
% MODELHCSZ0_cb = make_design_matrix_HCSZ_cb(varnamesCB,CBDEMO.demo_extras(CBfinalind,1) ,{CBDEMO.demo_extravars{1}},CBDEMO.demo_extravars{1},CBfinalind);
% 
% [DEMO0_COBRE_UKBinit, MULT0_COBRE_UKBinit, UNI0_COBRE_UKBinit] = run_model_HCSZ(MODELHCSZ0_cb, 0.01,cb_icasig_ukbinit_01(CBfinalloc,:), [], [], 0);
% save COBRE_UKBinit_preregSite_SMRI_Mancovan  *0_COBRE_UKBinit
% 
% [DEMO0_COBRE_UKBinit, MULT0_COBRE_UKBinit, UNI0_COBRE_UKBinit] = run_model_HCSZ(MODELHCSZ0_cb, 0.01,cb_icasig_ukbinit_01(CBfinalloc,:), [], [], 0);
% save COBRE_UKBinit_preregSite_SMRI_Mancovan  *0_COBRE_UKBinit
% 
% % ALFF
% MODELHCSZ1_cb = make_design_matrix_HCSZ_cb(varnamesCB,CBDEMO.demo_extras(CBfinalind,[2 4]) ,{CBDEMO.demo_extravars{2}, CBDEMO.demo_extravars{4}},{CBDEMO.demo_extravars(2),CBDEMO.demo_extravars(4)},CBfinalind);
% 
% [DEMO1_COBRE_UKBinit, MULT1_COBRE_UKBinit, UNI1_COBRE_UKBinit] = run_model_HCSZ(MODELHCSZ1_cb, 0.01,cb_icasig_ukbinit_02(CBfinalloc,:), [], [], 0);
% save COBRE_UKBinit_preregSite_ALFF_Mancovan  *1_COBRE_UKBinit
% 
% [DEMO1_COBRE_UKBinit, MULT1_COBRE_UKBinit, UNI1_COBRE_UKBinit] = run_model_HCSZ(MODELHCSZ1_cb, 0.01,cb_icasig_ukbinit_02(CBfinalloc,:), [], [], 0);
% save COBRE_UKBinit_preregSite_ALFF_Mancovan  *1_COBRE_UKBinit
% 
% % DMRI
% 
% MODELHCSZ2_cb = make_design_matrix_HCSZ_cb(varnamesCB,CBDEMO.demo_extras(CBfinalind,1) ,{CBDEMO.demo_extravars{3}},CBDEMO.demo_extravars{3},CBfinalind);
% 
% [DEMO2_COBRE_UKBinit, MULT2_COBRE_UKBinit, UNI2_COBRE_UKBinit] = run_model_HCSZ(MODELHCSZ2_cb, 0.01,cb_icasig_ukbinit_03(CBfinalloc,:), [], [], 0);
% save COBRE_UKBinit_preregSite_DMRI_Mancovan  *2_COBRE_UKBinit
% 
% [DEMO2_COBRE_UKBinit, MULT2_COBRE_UKBinit, UNI2_COBRE_UKBinit] = run_model_HCSZ(MODELHCSZ2_cb, 0.01,cb_icasig_ukbinit_03(CBfinalloc,:), [], [], 0);
% save COBRE_UKBinit_preregSite_DMRI_Mancovan  *2_COBRE_UKBinit
% 
% %% MPRC
% varnamesMP = MPDEMO.demo_col_names(1:4);
% 
% % SMRI
% MODELHCSZ0_mp = make_design_matrix_HCSZ_mp(varnamesMP,MPDEMO.demo_extras(MPfinalind,1) ,{MPDEMO.demo_extravars{1}},MPDEMO.demo_extravars{1},MPfinalind);
% 
% [DEMO0_MPRC_UKBinit, MULT0_MPRC_UKBinit, UNI0_MPRC_UKBinit] = run_model_HCSZ(MODELHCSZ0_mp, 0.01,mp_icasig_ukbinit_01(MPfinalloc,:), [], [], 0);
% save MPRC_UKBinit_preregSite_SMRI_Mancovan  *0_MPRC_UKBinit
% 
% [DEMO0_MPRC_UKBinit, MULT0_MPRC_UKBinit, UNI0_MPRC_UKBinit] = run_model_HCSZ(MODELHCSZ0_mp, 0.01,mp_icasig_ukbinit_01(MPfinalloc,:), [], [], 0);
% save MPRC_UKBinit_preregSite_SMRI_Mancovan  *0_MPRC_UKBinit
% 
% % ALFF
% MODELHCSZ1_mp = make_design_matrix_HCSZ_mp(varnamesMP,MPDEMO.demo_extras(MPfinalind,[2 4]) ,{MPDEMO.demo_extravars{2}, MPDEMO.demo_extravars{4}},{MPDEMO.demo_extravars(2),MPDEMO.demo_extravars(4)},MPfinalind);
% 
% [DEMO1_MPRC_UKBinit, MULT1_MPRC_UKBinit, UNI1_MPRC_UKBinit] = run_model_HCSZ(MODELHCSZ1_mp, 0.01,mp_icasig_ukbinit_02(MPfinalloc,:), [], [], 0);
% save MPRC_UKBinit_preregSite_ALFF_Mancovan  *1_MPRC_UKBinit
% 
% [DEMO1_MPRC_UKBinit, MULT1_MPRC_UKBinit, UNI1_MPRC_UKBinit] = run_model_HCSZ(MODELHCSZ1_mp, 0.01,mp_icasig_ukbinit_02(MPfinalloc,:), [], [], 0);
% save MPRC_UKBinit_preregSite_ALFF_Mancovan  *1_MPRC_UKBinit
% 
% % DMRI
% 
% MODELHCSZ2_mp = make_design_matrix_HCSZ_mp(varnamesMP,MPDEMO.demo_extras(MPfinalind,1) ,{MPDEMO.demo_extravars{3}},MPDEMO.demo_extravars{3},MPfinalind);
% 
% [DEMO2_MPRC_UKBinit, MULT2_MPRC_UKBinit, UNI2_MPRC_UKBinit] = run_model_HCSZ(MODELHCSZ2_mp, 0.01,mp_icasig_ukbinit_03(MPfinalloc,:), [], [], 0);
% save MPRC_UKBinit_preregSite_DMRI_Mancovan  *2_MPRC_UKBinit
% 
% [DEMO2_MPRC_UKBinit, MULT2_MPRC_UKBinit, UNI2_MPRC_UKBinit] = run_model_HCSZ(MODELHCSZ2_mp, 0.01,mp_icasig_ukbinit_03(MPfinalloc,:), [], [], 0);
% save MPRC_UKBinit_preregSite_DMRI_Mancovan  *2_MPRC_UKBinit
% 
% 
% %% BSNIP
% 
% varnamesBS = BSDEMO.demo_col_names(1:4);
% 
% %drpX=100;nBS1 = length(BSfinalind);
% %BSfinalind = BSfinalind(setdiff(1:nBS1,drpX));
% %BSfinalloc = BSfinalloc(setdiff(1:nBS1,drpX));
% % SMRI
% MODELHCSZ0_bs = make_design_matrix_HCSZ_bs(varnamesBS,BSDEMO.demo_extras(BSfinalind,1) ,{BSDEMO.demo_extravars{1}},BSDEMO.demo_extravars{1},BSfinalind);
% MODELHCSZ0_bs.labels{1} = {'Age'};MODELHCSZ0_bs.labels{5} = {'rSN_SMRI'};
% 
% MODELHCSZ0_bs.X(:,4) = MODELHCSZ0_bs.X(:,4)-1;
% MODELHCSZ0_bs.X(:,3) = MODELHCSZ0_bs.X(:,3)-1;
% [DEMO0_BSNIP_UKBinit, MULT0_BSNIP_UKBinit, UNI0_BSNIP_UKBinit] = run_model_HCSZ(MODELHCSZ0_bs, 0.01,bs_icasig_ukbinit_01(BSfinalloc,:), [], [], 0);
% save BSNIP_UKBinit_preregSite_SMRI_Mancovan  *0_BSNIP_UKBinit
% 
% [DEMO0_BSNIP_UKBinit, MULT0_BSNIP_UKBinit, UNI0_BSNIP_UKBinit] = run_model_HCSZ(MODELHCSZ0_bs, 0.01,bs_icasig_ukbinit_01(BSfinalloc,:), [], [], 0);
% save BSNIP_UKBinit_preregSite_SMRI_Mancovan  *0_BSNIP_UKBinit
% 
% % ALFF
% MODELHCSZ1_bs = make_design_matrix_HCSZ_bs(varnamesBS,BSDEMO.demo_extras(BSfinalind,[2 4]) ,{BSDEMO.demo_extravars{2}, BSDEMO.demo_extravars{4}},{BSDEMO.demo_extravars(2),BSDEMO.demo_extravars(4)},BSfinalind);
% MODELHCSZ1_bs.labels{1} = {'Age'};
% 
% MODELHCSZ1_bs.X(:,4) = MODELHCSZ0_bs.X(:,4)-1;
% MODELHCSZ1_bs.X(:,3) = MODELHCSZ0_bs.X(:,3)-1;
% 
% [DEMO1_BSNIP_UKBinit, MULT1_BSNIP_UKBinit, UNI1_BSNIP_UKBinit] = run_model_HCSZ(MODELHCSZ1_bs, 0.01,bs_icasig_ukbinit_02(BSfinalloc,:), [], [], 0);
% save BSNIP_UKBinit_preregSite_ALFF_Mancovan  *1_BSNIP_UKBinit
% 
% [DEMO1_BSNIP_UKBinit, MULT1_BSNIP_UKBinit, UNI1_BSNIP_UKBinit] = run_model_HCSZ(MODELHCSZ1_bs, 0.01,bs_icasig_ukbinit_02(BSfinalloc,:), [], [], 0);
% save BSNIP_UKBinit_preregSite_ALFF_Mancovan  *1_BSNIP_UKBinit
% 
% % DMRI
% 
% MODELHCSZ2_bs = make_design_matrix_HCSZ_bs(varnamesBS,BSDEMO.demo_extras(BSfinalind,1) ,{BSDEMO.demo_extravars{3}},BSDEMO.demo_extravars{3},BSfinalind);
% MODELHCSZ2_bs.labels{1} = {'Age'};MODELHCSZ2_bs.labels{5} = {'rSN_DMRI'};
% 
% MODELHCSZ2_bs.X(:,4) = MODELHCSZ0_bs.X(:,4)-1;
% MODELHCSZ2_bs.X(:,3) = MODELHCSZ0_bs.X(:,3)-1;
% 
% [DEMO2_BSNIP_UKBinit, MULT2_BSNIP_UKBinit, UNI2_BSNIP_UKBinit] = run_model_HCSZ(MODELHCSZ2_bs, 0.01,bs_icasig_ukbinit_03(BSfinalloc,:), [], [], 0);
% save BSNIP_UKBinit_preregSite_DMRI_Mancovan  *2_BSNIP_UKBinit
% 
% [DEMO2_BSNIP_UKBinit, MULT2_BSNIP_UKBinit, UNI2_BSNIP_UKBinit] = run_model_HCSZ(MODELHCSZ2_bs, 0.01,bs_icasig_ukbinit_03(BSfinalloc,:), [], [], 0);
% save BSNIP_UKBinit_preregSite_DMRI_Mancovan  *2_BSNIP_UKBinit
% %%
% load FBIRN_UKBinit_preregSite_SMRI_Mancovan
% load FBIRN_UKBinit_preregSite_SMRI_Mancovan.mat
% load FBIRN_UKBinit_preregSite_ALFF_Mancovan
% load FBIRN_UKBinit_preregSite_ALFF_Mancovan.mat
% load FBIRN_UKBinit_preregSite_DMRI_Mancovan
% load FBIRN_UKBinit_preregSite_DMRI_Mancovan.mat
% 
% load COBRE_UKBinit_preregSite_SMRI_Mancovan
% load COBRE_UKBinit_preregSite_SMRI_Mancovan.mat
% load COBRE_UKBinit_preregSite_ALFF_Mancovan
% load COBRE_UKBinit_preregSite_ALFF_Mancovan.mat
% load COBRE_UKBinit_preregSite_DMRI_Mancovan
% load COBRE_UKBinit_preregSite_DMRI_Mancovan.mat
% 
% load MPRC_UKBinit_preregSite_SMRI_Mancovan
% load MPRC_UKBinit_preregSite_SMRI_Mancovan.mat
% load MPRC_UKBinit_preregSite_ALFF_Mancovan
% load MPRC_UKBinit_preregSite_ALFF_Mancovan.mat
% load MPRC_UKBinit_preregSite_DMRI_Mancovan
% load MPRC_UKBinit_preregSite_DMRI_Mancovan.mat
% 
% load BSNIP_UKBinit_preregSite_SMRI_Mancovan
% load BSNIP_UKBinit_preregSite_SMRI_Mancovan.mat
% load BSNIP_UKBinit_preregSite_ALFF_Mancovan
% load BSNIP_UKBinit_preregSite_ALFF_Mancovan.mat
% load BSNIP_UKBinit_preregSite_DMRI_Mancovan
% load BSNIP_UKBinit_preregSite_DMRI_Mancovan.mat

%% sex
%SMRI
% [(myFDR(UNI0_FBIRN_UKBinit.p{1} ) <=0.05);(myFDR(UNI0_COBRE_UKBinit.p{1} ) <=0.05); (myFDR(UNI0_MPRC_UKBinit.p{2} ) <=0.05);(myFDR(UNI0_BSNIP_UKBinit.p{1}) <= 0.05)]
% %ALFF
% [(myFDR(UNI1_FBIRN_UKBinit.p{1} ) <=0.05);(myFDR(UNI1_COBRE_UKBinit.p{1} ) <=0.05); (myFDR(UNI1_MPRC_UKBinit.p{3} ) <=0.05);(myFDR(UNI1_BSNIP_UKBinit.p{1}) <= 0.05)]
% %DMRI
% [(myFDR(UNI2_FBIRN_UKBinit.p{1} ) <=0.05);(myFDR(UNI2_COBRE_UKBinit.p{1} ) <=0.05); (myFDR(UNI2_MPRC_UKBinit.p{2} ) <=0.05);(myFDR(UNI2_BSNIP_UKBinit.p{1}) <= 0.05)]
% 
% sum([(myFDR(UNI0_FBIRN_UKBinit.p{1} ) <=0.05);(myFDR(UNI0_COBRE_UKBinit.p{1} ) <=0.05); (myFDR(UNI0_MPRC_UKBinit.p{2} ) <=0.05);(myFDR(UNI0_BSNIP_UKBinit.p{1}) <= 0.05)])
% %ALFF
% sum([(myFDR(UNI1_FBIRN_UKBinit.p{1} ) <=0.05);(myFDR(UNI1_COBRE_UKBinit.p{1} ) <=0.05); (myFDR(UNI1_MPRC_UKBinit.p{3} ) <=0.05);(myFDR(UNI1_BSNIP_UKBinit.p{1}) <= 0.05)])
% %DMRI
% sum([(myFDR(UNI2_FBIRN_UKBinit.p{1} ) <=0.05);(myFDR(UNI2_COBRE_UKBinit.p{1} ) <=0.05); (myFDR(UNI2_MPRC_UKBinit.p{2} ) <=0.05);(myFDR(UNI2_BSNIP_UKBinit.p{1}) <= 0.05)]) 
% 
% %% age 
% sum([(myFDR(UNI0_FBIRN_UKBinit.p{4} ) <=0.05);(myFDR(UNI0_COBRE_UKBinit.p{3} ) <=0.05); (myFDR(UNI0_MPRC_UKBinit.p{3} ) <=0.05);(myFDR(UNI0_BSNIP_UKBinit.p{4}) <= 0.05)])
% %ALFF
% sum([(myFDR(UNI1_FBIRN_UKBinit.p{4} ) <=0.05);(myFDR(UNI1_COBRE_UKBinit.p{3} ) <=0.05); (myFDR(UNI1_MPRC_UKBinit.p{4} ) <=0.05);(myFDR(UNI1_BSNIP_UKBinit.p{4}) <= 0.05)])
% %DMRI
% sum([(myFDR(UNI2_FBIRN_UKBinit.p{4} ) <=0.05);(myFDR(UNI2_COBRE_UKBinit.p{2} ) <=0.05); (myFDR(UNI2_MPRC_UKBinit.p{3} ) <=0.05);(myFDR(UNI2_BSNIP_UKBinit.p{3}) <= 0.05)])



%%
% pVec = (UNI0_BSNIP_UKBinit.p{3}(:));
% pVecF = myFDR(pVec);
% pVecF = reshape(pVecF,36,20);
% sigI=find(pVecF<=0.05);
% [i1,i2]= ind2sub(size(pVecF),sigI);
% 
% sigContr_Diag_BSNIP_UKBinit = UNI0_BSNIP_UKBinit.stats{3}.Contrast(unique(i1));
% 
% %%
% pVec1 = (UNI0_BSNIP_UKBinit.p{3}(:));
% pVecF1 = myFDR(pVec1);
% pVecF1 = reshape(pVecF1,36,20);
% sigI1=find(pVecF1<=0.05);
% [i11,i21]= ind2sub(size(pVecF1),sigI1);
% 
% sigContr_Diag_BSNIP_UKBinit = UNI0_BSNIP_UKBinit.stats{3}.Contrast(unique(i11));
% %%
% sigIDb = find(myFDR(UNI1_BSNIP_UKBinit.p{3}(:))<=0.05);
% [i1b,i2b] = ind2sub(size(UNI1_BSNIP_UKBinit.p{3}),sigIDb);
% 
% [num2str(UNI1_BSNIP_UKBinit.p{3},'%1.4g') repmat(' ',36,1) strvcat(UNI1_BSNIP_UKBinit.stats{3}.Contrast')]

%% UKBinit

%% sex
%SMRI
% [(myFDR(UNI0_FBIRN_UKBinit.p{1} ) <=0.05);(myFDR(UNI0_COBRE_UKBinit.p{1} ) <=0.05); (myFDR(UNI0_MPRC_UKBinit.p{2} ) <=0.05);(myFDR(UNI0_BSNIP_UKBinit.p{1}) <= 0.05)]
% %ALFF
% [(myFDR(UNI1_FBIRN_UKBinit.p{1} ) <=0.05);(myFDR(UNI1_COBRE_UKBinit.p{1} ) <=0.05); (myFDR(UNI1_MPRC_UKBinit.p{3} ) <=0.05);(myFDR(UNI1_BSNIP_UKBinit.p{1}) <= 0.05)]
% %DMRI
% [(myFDR(UNI2_FBIRN_UKBinit.p{1} ) <=0.05);(myFDR(UNI2_COBRE_UKBinit.p{1} ) <=0.05); (myFDR(UNI2_MPRC_UKBinit.p{2} ) <=0.05);(myFDR(UNI2_BSNIP_UKBinit.p{1}) <= 0.05)]
% 
% %% age 
% sum([(myFDR(UNI0_FBIRN_UKBinit.p{4} ) <=0.05);(myFDR(UNI0_COBRE_UKBinit.p{3} ) <=0.05); (myFDR(UNI0_MPRC_UKBinit.p{3} ) <=0.05);(myFDR(UNI0_BSNIP_UKBinit.p{4}) <= 0.05)])
% %ALFF
% sum([(myFDR(UNI1_FBIRN_UKBinit.p{4} ) <=0.05);(myFDR(UNI1_COBRE_UKBinit.p{3} ) <=0.05); (myFDR(UNI1_MPRC_UKBinit.p{4} ) <=0.05);(myFDR(UNI1_BSNIP_UKBinit.p{4}) <= 0.05)])
% %DMRI
% sum([(myFDR(UNI2_FBIRN_UKBinit.p{2} ) <=0.05);(myFDR(UNI2_COBRE_UKBinit.p{2} ) <=0.05); (myFDR(UNI2_MPRC_UKBinit.p{3} ) <=0.05);(myFDR(UNI2_BSNIP_UKBinit.p{3}) <= 0.05)])

% %% diag
% sum([(myFDR(UNI0_FBIRN_UKBinit.p{2} ) <=0.05);(myFDR(UNI0_COBRE_UKBinit.p{2} ) <=0.05); (myFDR(UNI0_BSNIP_UKBinit.p{3}) <= 0.05)])
% %ALFF
% sum([(myFDR(UNI1_FBIRN_UKBinit.p{4} ) <=0.05);(myFDR(UNI1_COBRE_UKBinit.p{3} ) <=0.05); (myFDR(UNI1_MPRC_UKBinit.p{4} ) <=0.05);(myFDR(UNI1_BSNIP_UKBinit.p{4}) <= 0.05)])
% %DMRI
% sum([(myFDR(UNI2_FBIRN_UKBinit.p{2} ) <=0.05);(myFDR(UNI2_COBRE_UKBinit.p{2} ) <=0.05); (myFDR(UNI2_MPRC_UKBinit.p{3} ) <=0.05);(myFDR(UNI2_BSNIP_UKBinit.p{3}) <= 0.05)])
%%


%MODELHCSZ0_cb.X(:,3) = 2-MODELHCSZ0_cb.X(:,3); % 0-HC; 1: 1 SZ; 2 BP; 3 SZA
CBDEMO.demo_out(:,3) = 2 - CBDEMO.demo_out(:,3);% 0-HC; 1: 1 SZ; 2 BP; 3 SZA

BSDEMO.demo_out(:,4) = 3 - BSDEMO.demo_out(:,4);%   1-'BPP'}   2-{'BPR'}    0-{'NC'}   -1-{'SADBPP'}   -2: {'SADBPR'}  -3:{'SADDEPP'}   -4:{'SADDEPR'}   -5:{'SZP'}   -6: {'SZR'}
%MODELHCSZ0_bs.X(:,4) = 3 - MODELHCSZ0_bs.X(:,4); %   1-'BPP'}   2-{'BPR'}    0-{'NC'}   -1-{'SADBPP'}   -2: {'SADBPR'}  -3:{'SADDEPP'}   -4:{'SADDEPR'}   -5:{'SZP'}   -6: {'SZR'}
BSDEMO.demo_out( BSDEMO.demo_out(:,4)==1  ,4) = 12; % BPPtemp
BSDEMO.demo_out( BSDEMO.demo_out(:,4)==2  ,4) = 13; % BPRtemp

BSDEMO.demo_out( BSDEMO.demo_out(:,4)==-1  ,4) = 14; % SADBPPtemp
BSDEMO.demo_out( BSDEMO.demo_out(:,4)==-2  ,4) = 15; % SADBPRtemp
BSDEMO.demo_out( BSDEMO.demo_out(:,4)==-3  ,4) = 16; % SADEPPtemp
BSDEMO.demo_out( BSDEMO.demo_out(:,4)==-4  ,4) = 17; % SADEPRtemp
BSDEMO.demo_out( BSDEMO.demo_out(:,4)==-6  ,4) = 18; %  SZRtemp

BSDEMO.demo_out( BSDEMO.demo_out(:,4)==-5  ,4) = 1;  % SZP
BSDEMO.demo_out( BSDEMO.demo_out(:,4)==12  ,4) = 2;  % BPP
%                                                 3   SZA (in COBRE)  
BSDEMO.demo_out( BSDEMO.demo_out(:,4)==18  ,4) = 4;  % SZR
BSDEMO.demo_out( BSDEMO.demo_out(:,4)==13  ,4) = 5;  % BPR
BSDEMO.demo_out( BSDEMO.demo_out(:,4)==14  ,4) = 6;  % SADBPP
BSDEMO.demo_out( BSDEMO.demo_out(:,4)==15  ,4) = 7;  % SADBPR
BSDEMO.demo_out( BSDEMO.demo_out(:,4)==16  ,4) = 8;  % SADEPP
BSDEMO.demo_out( BSDEMO.demo_out(:,4)==17  ,4) = 9;  % SADEPR

%MODELHCSZ0_mp.X(:,2) = MODELHCSZ0_mp.X(:,2)+3;

% combined order: Age, sex, diag, rSN           %, site; 
% for site: COBRE : 0
%           FBIRN : 3     7     9    10    12    13    18;
%           MPRC  : 4   5    6
%           BSNIP : 1   2
MODELHCSZ0ns.X = [[FBDEMO.demo_out(FBfinalind,[1 2 3])  FBDEMO.demo_extras(FBfinalind,1)]; ...
                [CBDEMO.demo_out(CBfinalind,[1 2 3])  CBDEMO.demo_extras(CBfinalind,1)]; ...
                [MPDEMO.demo_out(MPfinalind,[3 4 1])  MPDEMO.demo_extras(MPfinalind,1)]; ...
                [BSDEMO.demo_out(BSfinalind,[1 2 4])  BSDEMO.demo_extras(BSfinalind,1)]; ...
                ];
                
MODELHCSZ0ns.group_ind = [2 3];MODELHCSZ0ns.cov_ind = [1 4];MODELHCSZ0ns.names = {'age','sex','diagnosis','rSN_SMRI'};
MODELHCSZ0ns.labels{1} =  {'age'};
MODELHCSZ0ns.labels{2} =  {'M' 'F'};
MODELHCSZ0ns.labels{3} = {'HC' 'SZ' 'BP' 'SZA' 'SZR' 'BPR' 'SADBPP' 'SADBPR' 'SADEPP' 'SADEPR'};
%                          0     1    2    3     4     5      6         7

%                      
MODELHCSZ0ns.labels{4} = {'rSN_SMRI'};
%MODELHCSZ0ns.labels{5} = {'CB0' 'BS1' 'BS2' 'FB1' 'MP1' 'MP2' 'MP3'  'FB2' 'FB3' 'FB4' 'FB5' 'FB6' 'FB7'};
%%
drpSubID =  mod(find(isnan(MODELHCSZ0ns.X)),length(MODELHCSZ0ns.X));
keepID = setdiff(1:length(MODELHCSZ0ns.X),drpSubID);
MODELHCSZ0ns.X = MODELHCSZ0ns.X(keepID,:);

%%
NMODELHCSZ0ns = MODELHCSZ0ns;
nL = length(MODELHCSZ0ns.labels{3});
ids =  cell(1,nL);
for ll = 1:nL 
    ids{ll} = find(NMODELHCSZ0ns.X(:,3)==ll-1);
end

for ll = 1:nL
   if MODELHCSZ0ns.labels{3}{ll}(end)=='R'
       NMODELHCSZ0ns.X(ids{ll},3) = 0; % HC
       NMODELHCSZ0ns.labels{3}{ll}='';
   end
end
NMODELHCSZ0ns.labels{3} = NMODELHCSZ0ns.labels{3}(~cellfun(@isempty,NMODELHCSZ0ns.labels{3}));
NMODELHCSZ0ns.X(ids{7},3) = 4;
NMODELHCSZ0ns.X(ids{9},3) = 5;     
    
%%


Y0_prereg_ukbinit = fulreg_HCSZ_UKBinit.icasig{1}';
%Y0_prereg_ukbinit = fulreg_HCSZ_UKBinit.icasig{1}';


Y0_prereg_ukbinit = Y0_prereg_ukbinit(keepID,:);
%Y0_prereg_ukbinit = Y0_prereg_ukbinit(keepID,:);
%%

%%
%[DEMO0_All_preregSite_UKBinit, MULT0_All_preregSite_UKBinit, UNI0_All_preregSite_UKBinit] = run_model_HCSZ(NMODELHCSZ0ns, 0.01,Y0_prereg_ukbinit, [], [], 0);
%save mancovaOuts_allHCSZcombined_preregSite_C30_SMRI_UKBinit *0_All_preregSite_UKBinit
[DEMO0_All_preregSite_UKBinit, MULT0_All_preregSite_UKBinit, UNI0_All_preregSite_UKBinit] = run_model_HCSZ(NMODELHCSZ0ns, 0.01,Y0_prereg_ukbinit, [], [], 0);
save mancovaOuts_allHCSZ_combinedRelatives_preregSite_C30_SMRI_UKBinit *0_All_preregSite_UKBinit  NMODELHCSZ0ns

%with interactions
%[DEMO0_AllwX_preregSite_UKBinit, MULT0_AllwX_preregSite_UKBinit, UNI0_AllwX_preregSite_UKBinit] = run_model_HCSZ_withX(MODELHCSZ0ns, 0.01,Y0_prereg_ukbinit, [], [], 0);
%save mancovaOuts_allHCSZcombined_wX_preregSite_C30_SMRI_UKBinit *0_AllwX_preregSite_UKBinit
[DEMO0_AllwX_preregSite_UKBinit, MULT0_AllwX_preregSite_UKBinit, UNI0_AllwX_preregSite_UKBinit] = run_model_HCSZ_withX(NMODELHCSZ0ns, 0.01,Y0_prereg_ukbinit, [], [], 0);
save mancovaOuts_allHCSZ_combinedRelatives_wX_preregSite_C30_SMRI_UKBinit *0_AllwX_preregSite_UKBinit NMODELHCSZ0ns

%%


NMODELHCSZ1ns =  NMODELHCSZ0ns;
temp = [FBDEMO.demo_extras(FBfinalind,2); CBDEMO.demo_extras(CBfinalind,2) ;MPDEMO.demo_extras(MPfinalind,2); BSDEMO.demo_extras(BSfinalind,2)];
temp(drpSubID) = []; 
NMODELHCSZ1ns.X(:,4) = temp;
NMODELHCSZ1ns.names{4} = 'rSN_ALFF';
NMODELHCSZ1ns.labels{4} = {'rSN_ALFF'};
temp = [FBDEMO.demo_extras(FBfinalind,4); CBDEMO.demo_extras(CBfinalind,4) ;MPDEMO.demo_extras(MPfinalind,4); BSDEMO.demo_extras(BSfinalind,4)];
temp(drpSubID) = [];

NMODELHCSZ1ns.X(:,5) = temp;%[FBDEMO.demo_extras(FBfinalind,4); CBDEMO.demo_extras(CBfinalind,4) ;MPDEMO.demo_extras(MPfinalind,4); BSDEMO.demo_extras(BSfinalind2,4)];
NMODELHCSZ1ns.names{5} = 'meanFD_FMRI';
NMODELHCSZ1ns.labels{5} = {'meanFD_FMRI'};
NMODELHCSZ1ns.cov_ind = [1 4 5];
%%
Y1_prereg_ukbinit = fulreg_HCSZ_UKBinit.icasig{2}';
%Y1_prereg_ukbinit = fulreg_HCSZ_UKBinit.icasig{2}';
Y1_prereg_ukbinit= Y1_prereg_ukbinit(keepID,:);
%Y1_prereg_ukbinit = Y1_prereg_ukbinit(keepID,:);
%%
%[DEMO1_All_preregSite_UKBinit, MULT1_All_preregSite_UKBinit, UNI1_All_preregSite_UKBinit] = run_model_HCSZ(NMODELHCSZ1ns, 0.01,Y1_prereg_ukbinit, [], [], 0);
%save mancovaOuts_allHCSZcombined_preregSite_C30_ALFF_UKBinit *1_All_preregSite_UKBinit
[DEMO1_All_preregSite_UKBinit, MULT1_All_preregSite_UKBinit, UNI1_All_preregSite_UKBinit] = run_model_HCSZ(NMODELHCSZ1ns, 0.01,Y1_prereg_ukbinit, [], [], 0);
save mancovaOuts_allHCSZ_combinedRelatives_preregSite_C30_ALFF_UKBinit *1_All_preregSite_UKBinit


%with interactions
%[DEMO1_AllwX_preregSite_UKBinit, MULT1_AllwX_preregSite_UKBinit, UNI1_AllwX_preregSite_UKBinit] = run_model_HCSZ_withX(NMODELHCSZ1ns, 0.01,Y1_prereg_ukbinit, [], [], 0);
%save mancovaOuts_allHCSZcombined_wX_preregSite_C30_ALFF_UKBinit *1_AllwX_preregSite_UKBinit
[DEMO1_AllwX_preregSite_UKBinit, MULT1_AllwX_preregSite_UKBinit, UNI1_AllwX_preregSite_UKBinit] = run_model_HCSZ_withX(NMODELHCSZ1ns, 0.01,Y1_prereg_ukbinit, [], [], 0);
save mancovaOuts_allHCSZ_combinedRelatives_wX_preregSite_C30_ALFF_UKBinit *1_AllwX_preregSite_UKBinit

%%
NMODELHCSZ2ns =  NMODELHCSZ0ns;
temp = [FBDEMO.demo_extras(FBfinalind,3); CBDEMO.demo_extras(CBfinalind,3) ;MPDEMO.demo_extras(MPfinalind,3); BSDEMO.demo_extras(BSfinalind,3)];
temp(drpSubID) = [];

NMODELHCSZ2ns.X(:,4) = temp;%[FBDEMO.demo_extras(FBfinalind,3); CBDEMO.demo_extras(CBfinalind,3) ;MPDEMO.demo_extras(MPfinalind,3); BSDEMO.demo_extras(BSfinalind2,3)];
NMODELHCSZ2ns.names{4} = 'rSN_DMRI';
NMODELHCSZ2ns.labels{4} = {'rSN_DMRI'};
%%
Y2_prereg_ukbinit = fulreg_HCSZ_UKBinit.icasig{3}';
%Y2_prereg_ukbinit = fulreg_HCSZ_UKBinit.icasig{3}';
Y2_prereg_ukbinit= Y2_prereg_ukbinit(keepID,:);
%Y2_prereg_ukbinit = Y2_prereg_ukbinit(keepID,:);

%%
%[DEMO2_All_preregSite_UKBinit, MULT2_All_preregSite_UKBinit, UNI2_All_preregSite_UKBinit] = run_model_HCSZ(NMODELHCSZ2ns, 0.01,Y2_prereg_ukbinit, [], [], 0);
%save mancovaOuts_allHCSZcombined_preregSite_C30_DMRI_UKBinit *2_All_preregSite_UKBinit
[DEMO2_All_preregSite_UKBinit, MULT2_All_preregSite_UKBinit, UNI2_All_preregSite_UKBinit] = run_model_HCSZ(NMODELHCSZ2ns, 0.01,Y2_prereg_ukbinit, [], [], 0);
save mancovaOuts_allHCSZ_combinedRelatives_preregSite_C30_DMRI_UKBinit *2_All_preregSite_UKBinit

%with interactions
%[DEMO2_AllwX_preregSite_UKBinit, MULT2_AllwX_preregSite_UKBinit, UNI2_AllwX_preregSite_UKBinit] = run_model_HCSZ_withX(NMODELHCSZ2ns, 0.01,Y2_prereg_ukbinit, [], [], 0);
%save mancovaOuts_allHCSZcombined_wX_preregSite_C30_DMRI_UKBinit *2_AllwX_preregSite_UKBinit
[DEMO2_AllwX_preregSite_UKBinit, MULT2_AllwX_preregSite_UKBinit, UNI2_AllwX_preregSite_UKBinit] = run_model_HCSZ_withX(NMODELHCSZ2ns, 0.01,Y2_prereg_ukbinit, [], [], 0);
save mancovaOuts_allHCSZ_combinedRelatives_wX_preregSite_C30_DMRI_UKBinit *2_AllwX_preregSite_UKBinit


%%
% SMRI
[sx00_comp,sx00_T,sx00_contr]= printouts(UNI0_AllwX_preregSite_UKBinit,1);
[d00_comp,d00_T,d00_contr]= printouts(UNI0_AllwX_preregSite_UKBinit,2);
[a00_comp,a00_T,a00_contr] = printouts(UNI0_AllwX_preregSite_UKBinit,3);
%[rsn00_comp,rsn00_T,rsn00_contr]=printouts(UNI0_AllwX_preregSite_UKBinit,5);
[sxd00_comp,sxd00_T,sxd00_contr]=printouts(UNI0_AllwX_preregSite_UKBinit,4);

% ALFF
%[sx11_comp,sx11_T,sx11_contr]=printouts(UNI1_AllwX_preregSite_UKBinit,1);
[d11_comp,d11_T,d11_contr]= printouts(UNI1_AllwX_preregSite_UKBinit,1);
%[st11_comp,st11_T,st11_contr]= printouts(UNI1_AllwX_preregSite_UKBinit,3);
[a11_comp,a11_T,a11_contr] = printouts(UNI1_AllwX_preregSite_UKBinit,2);
[rsn11_comp,rsn11_T,rsn11_contr]=printouts(UNI1_AllwX_preregSite_UKBinit,3);
[rfd11_comp,rfd11_T,rfd11_contr]=printouts(UNI1_AllwX_preregSite_UKBinit,4);
%[sxd11_comp,sxd11_T,sxd11_contr]=printouts(UNI1_AllwX_preregSite_UKBinit,6);


% DMRI
[sx22_comp,sx22_T,sx22_contr]=printouts(UNI2_AllwX_preregSite_UKBinit,1);
[d22_comp,d22_T,d22_contr]= printouts(UNI2_AllwX_preregSite_UKBinit,2);
%[st22_comp,st22_T,st22_contr]= printouts(UNI2_AllwX_preregSite_UKBinit,3);
[a22_comp,a22_T,a22_contr] = printouts(UNI2_AllwX_preregSite_UKBinit,3);
%[rsn22_comp,rsn22_T,rsn22_contr]=printouts(UNI2_AllwX_preregSite_UKBinit,5);
[sxd22_comp,sxd22_T,sxd22_contr]=printouts(UNI2_AllwX_preregSite_UKBinit,4);
%[axd22_comp,axd22_T,axd22_contr]=printouts(UNI2_AllwX_preregSite_UKBinit,5);
%%
aa{1} = UNI0_AllwX_preregSite_UKBinit;
aa{2} = UNI1_AllwX_preregSite_UKBinit;
aa{3} = UNI2_AllwX_preregSite_UKBinit;
%
dm{1} = DEMO0_All_preregSite_UKBinit;
dm{2} = DEMO1_All_preregSite_UKBinit;
dm{3} = DEMO2_All_preregSite_UKBinit;
%%
%plotuni(aa,dm,1,'sex',sx00_contr,sx00_comp,sx00_T)
%plotuni(aa,dm,2,'sex',sx11_contr,sx11_comp,sx11_T)
%plotuni(aa,dm,3,'sex',sx22_contr,sx22_comp,sx22_T)
outdir = fullfile(HCSZroot,'MMIVA2step_wUKBinit','C030_regsite_ukbinit','figs','stats');
plotuni(aa,dm,1,'age',a00_contr,a00_comp,a00_T,outdir)
plotuni(aa,dm,2,'age',a11_contr,a11_comp,a11_T,outdir)
plotuni(aa,dm,3,'age',a22_contr,a22_comp,a22_T,outdir)

plotuni(aa,dm,1,'diagnosis',d00_contr,d00_comp,d00_T,outdir)
plotuni(aa,dm,2,'diagnosis',d11_contr,d11_comp,d11_T,outdir)
plotuni(aa,dm,3,'diagnosis',d22_contr,d22_comp,d22_T,outdir)

plotuni(aa,dm,2,'rSN_ALFF',rsn11_contr,rsn11_comp,rsn11_T),outdir
%%

label = 'sex';
plotouts(aa,{'SMRI','ALFF','DMRI'},label,1)
export_fig(['UniResults_preregSite_wX_' label '.png'],gcf)
close(gcf)
%
label = 'age';
plotouts(aa,{'SMRI','ALFF','DMRI'},label,1)
export_fig(['UniResults_preregSite_wX_' label '.png'],gcf)
close(gcf)
%
label = 'diagnosis';
plotouts(aa,{'SMRI','ALFF','DMRI'},label,1)
export_fig(['UniResults_preregSite_wX_' label '.png'],gcf)
close(gcf)

%
label = 'sex_X_diagnosis';
plotouts(aa,{'SMRI','ALFF','DMRI'},label,1)
export_fig(['UniResults_preregSite_wX_' label '.png'],gcf)
close(gcf)

%%
atrm = 'age';
ageCol = find(contains(DEMO0_AllwXpreregSite_UKBinit.names,atrm)>0);
Age = DEMO0_AllwXpreregSite_UKBinit.full(:,ageCol);
cmAge = vals2colormap(Age,'jet',[min(Age) max(Age)]);
alocs0 = find(contains(UNI0_AllwXpreregSite_UKBinit.tests,atrm)>0);
alocs1 = find(contains(UNI1_AllwX_preregSite_UKBinit.tests,atrm)>0);
alocs2 = find(contains(UNI2_AllwX_preregSite_UKBinit.tests,atrm)>0);
% ag
agsig  = [(myFDR(UNI0_AllwXpreregSite_UKBinit.p{alocs0(1)} ) <=0.05);(myFDR(UNI1_AllwX_preregSite_UKBinit.p{alocs1(1)} ) <=0.05); (myFDR(UNI2_AllwX_preregSite_UKBinit.p{alocs2(1)} ) <=0.05)];
agsigID = find(sum(agsig)==3);
tmpp = [myFDR(UNI0_AllwXpreregSite_UKBinit.p{alocs0(1)} );myFDR(UNI1_AllwX_preregSite_UKBinit.p{alocs1(1)} );myFDR(UNI2_AllwX_preregSite_UKBinit.p{alocs2(1)} )];
agsig_p = tmpp(:,agsigID); 
%
%ageSig3_comp = intersect(intersect(a00_comp{1},a11_comp{1}),a22_comp{1})

for ii = 1:length(agsigID)
    figure('color',[1 1 1]);
    cc = agsigID(ii);
    scatter3(Y0_prereg_ukbinit(:,cc),Y1_prereg_ukbinit(:,cc),Y2_prereg_ukbinit(:,cc),[],cmAge);    
    view([30 10]);
    xlabel('SMRI');ylabel('ALFF');zlabel('DMRI');
    title(['Comp ' my_pad(cc,2)]);
    
    %set(cb,'FontName','Arial')
    set(gca,'FontName','Arial')
    export_fig('UKBicasig_full_prereg_UKBinit_C30_bysigfdrAge.pdf','-append',gcf);
    close all;
end

%% diagnosis

d0loc = find(cellfun(@isempty,d00_comp)==0);
d1loc = find(cellfun(@isempty,d11_comp)==0);
d2loc = find(cellfun(@isempty,d22_comp)==0);


%
%%
sxd0loc = find(cellfun(@isempty,sxd00_comp)==0);
sxd1loc = find(cellfun(@isempty,sxd11_comp)==0);
sxd2loc = find(cellfun(@isempty,sxd22_comp)==0);

%%

% -log10([p00_agefdr(ageSigInd);p11_agefdr(ageSigInd);p22_agefdr(ageSigInd)]+eps)

%   12.9625    7.2214    5.9038   15.6536    6.1938    5.4391   11.6836    3.8414    1.5718    3.5950    9.1811   15.6536    7.6948    5.4056   10.3801    3.9036
%    5.6656    4.7438    4.1856   15.6536    3.6144    3.3525    4.5007    2.6371    1.6982    4.0448    5.9298   15.6536    5.8948    3.4508    3.4508    1.6982
%    7.6058    5.2016    5.1901   15.6536    6.5546    3.7637   12.6361    3.3883    3.7976    5.0064    9.2074   15.6536    7.7077    5.2016    7.9205    4.6868
Tage_sigInd = [UNI0_AllwX_preregSite_UKBinit.t{3}(ageSigInd);UNI1_AllwX_preregSite_UKBinit.t{2}(ageSigInd);UNI2_AllwX_preregSite_UKBinit.t{3}(ageSigInd)]
%   -7.8486    5.6978   -5.0959  -10.1143   -5.2443   -4.8647    7.4074   -3.9864    2.4112   -3.8313    6.4946   16.8825   -5.9099   -4.8313    6.9432    4.0447
%   -5.1238    4.6592   -4.3172   -8.7336    3.9580   -3.7428   -4.5068   -3.2757    2.5364   -4.2167    5.3364   15.1580   -5.2673   -3.8202    3.8214    2.5507
%   -5.8749   -4.7800   -4.7473   -9.9024    5.4214   -3.9273    7.7486   -3.6953    3.9620   -4.6407    6.5674   17.3100   -5.9412   -4.7724    6.0544    4.4492


%% plotting


sloc = find(contains(DEMO0_AllwX_preregSite_UKBinit.names,'sex')>0);
x0id = find(DEMO0_AllwX_preregSite_UKBinit.full(:,sloc) == 0); % male
x1id = find(DEMO0_AllwX_preregSite_UKBinit.full(:,sloc) == 1); % male

figure('color',[1 1 1]);
for ii = 1:length(sxsigID)
    cc = sxsigID(ii);
    scatter3(Y0_prereg_ukbinit(x0id,cc),Y1_prereg_ukbinit(x0id,cc),Y2_prereg_ukbinit(x0id,cc),[],cmAge(x0id,:));
    hold on;
    scatter3(Y0_prereg_ukbinit(x1id,cc),Y1_prereg_ukbinit(x1id,cc),Y2_prereg_ukbinit(x1id,cc),[],cmAge(x1id,:),'x');
    view([30 10]);
    xlabel('SMRI');ylabel('ALFF');zlabel('DMRI');
    title(['Comp ' my_pad(cc,2)]);
    
    %set(cb,'FontName','Arial')
    set(gca,'FontName','Arial')
    %export_fig('UKBicasig_full_prereg_UKBinit_C30_bysigfdrAge.pdf','-append',gcf);
    %close all;
end

%% Sex X diagnosis
    

pp0 = UNI0_AllwXpreregSite_UKBinit.p{locs0(2)};
pp1 = UNI1_AllwXpreregSite_UKBinit.p{locs1(2)};
pp2 = UNI2_AllwXpreregSite_UKBinit.p{locs2(2)};
pp0f = reshape(myFDR(pp0(:)),size(pp0));
pp1f = reshape(myFDR(pp1(:)),size(pp1));
pp2f = reshape(myFDR(pp2(:)),size(pp2));

[(myFDR(UNI0_AllwXpreregSite_UKBinit.p{locs0(2)} ) <=0.05);(myFDR(UNI1_AllwXpreregSite_UKBinit.p{locs1(2)} ) <=0.05); (myFDR(UNI2_AllwXpreregSite_UKBinit.p{locs2(2)} ) <=0.05)]

%%

%% age 
ageVar = 'age';
tta1=find(strcmp(ageVar ,UNI0_AllwX_preregSite_UKBinit.tests));
tta2=find(strcmp(ageVar ,UNI1_AllwX_preregSite_UKBinit.tests));
tta3 = find(strcmp(ageVar ,UNI2_AllwX_preregSite_UKBinit.tests));

[myFDR(UNI0_AllwX_preregSite_UKBinit.p{tta1}) <= 0.05;myFDR(UNI1_AllwX_preregSite_UKBinit.p{tta2}) <= 0.05;myFDR(UNI2_AllwX_preregSite_UKBinit.p{tta3}) <= 0.05]
sigCompAGEfdr = find(sum([myFDR(UNI0_AllwX_preregSite_UKBinit.p{tta1}) <= 0.05;myFDR(UNI1_AllwX_preregSite_UKBinit.p{tta2}) <= 0.05;myFDR(UNI2_AllwX_preregSite_UKBinit.p{tta3}) <= 0.05])==3);
sigCompAGEbf = find(sum([(UNI0_AllwX_preregSite_UKBinit.p{tta1}) <= 0.05/20;(UNI1_AllwX_preregSite_UKBinit.p{tta2}) <= 0.05/20;(UNI2_AllwX_preregSite_UKBinit.p{tta3}) <= 0.05/20])==3);
%%
ageloc = find(contains(DEMO0_AllwX_preregSite_UKBinit.names,'age')>0);
sexloc = find(contains(DEMO0_AllwX_preregSite_UKBinit.names,'sex')>0);
if any(ageloc)
    ageloc = ageloc(1);
end
if any(sexloc)
    sexloc = sexloc(1);
end
sexCol = sexloc;
ageCol = ageloc;


Age = DEMO0_AllwX_preregSite_UKBinit.full(:,ageCol);
cmAge = vals2colormap(Age,'jet',[min(Age) max(Age)]);

sex = DEMO0_AllwX_preregSite_UKBinit.full(:,sexCol);
cmSex = vals2colormap(sex,'jet',[min(sex) max(sex)]);
x0id = find(sex==0);
x1id = find(sex==1);
%%
nC=30;
for cc= sigCompAGEfdr % [1 2 5 6 8 9 13 19]
    figure('color',[1 1 1]);
    %scatter3(Yful_SMRI(cc,:)',Yful_ALFF(cc,:)',Yful_DMRI(cc,:)',[],cmAge);
    scatter3(Y0_prereg_ukbinit(x0id,cc),Y1_prereg_ukbinit(x0id,cc),Y2_prereg_ukbinit(x0id,cc),[],cmAge(x0id,:));
    hold on;
    scatter3(Y0_prereg_ukbinit(x1id,cc),Y1_prereg_ukbinit(x1id,cc),Y2_prereg_ukbinit(x1id,cc),[],cmAge(x1id,:),'filled');
    view([30 10]);
    xlabel('SMRI');ylabel('ALFF');zlabel('DMRI');
    title(['Comp ' my_pad(cc,2)]);
    
    %set(cb,'FontName','Arial')
    set(gca,'FontName','Arial')
    export_fig('HCSZ_icasig_full_UKBinit_regsite_C30_bysigfdrAge.pdf','-append',gcf);
    close all;
end

%% Diagnosis

% SMRI
psigID = reshape(myFDR(UNI0_AllwXpreregSite_UKBinit.p{2}(:)) <=0.05,size(UNI0_AllwXpreregSite_UKBinit.p{2}));
for rr=1:45
    sigc = find(psigID(rr,:)==1);
    if any(sigc)
        disp([UNI0_AllwXpreregSite_UKBinit.stats{2}.Contrast{rr} ': ' num2str(sigc.*sign(UNI0_AllwXpreregSite_UKBinit.t{2}(sigc)) )]);
    end
end
%%
HCid = find(DEMO0_AllwX_preregSite_UKBinit.full(:,3)==0);
SZid = find(DEMO0_AllwX_preregSite_UKBinit.full(:,3)==1);
BPid = find(DEMO0_AllwX_preregSite_UKBinit.full(:,3)==2);
SZAid = find(DEMO0_AllwX_preregSite_UKBinit.full(:,3)==3);
SADBPPid = find(DEMO0_AllwX_preregSite_UKBinit.full(:,3)==4);
SADEPPid = find(DEMO0_AllwX_preregSite_UKBinit.full(:,3)==5);
%%


%%
% HCid = find(MODELHCSZ0ns.X(:,3)==0);SZid=find(MODELHCSZ0ns.X(:,3)==1);
% sigInd_diag_all_UKBinit_SMRI = find(myFDR(UNI0_AllpreregSite_UKBinit.p{1}(:))<=0.05);
% [i0_1g,i0_2g] = ind2sub(size(UNI0_AllpreregSite_UKBinit.p{1}),sigInd_diag_all_UKBinit_SMRI);
% 
% for ii=1:length(i0_1g)
%     if mod(ii,24)==1
%         figure('color',[1 1 1]);
%         jj=1;
%     end
%     contr =  UNI0_AllpreregSite_UKBinit.stats{1}.Contrast(i0_1g(ii));
%     contr = contr{1};
%     labst = strfind(contr,'(');labnd = strfind(contr,')');
%     lab1 = contr(labst(1)+1:labnd(1)-1);lab2 = contr(labst(2)+1:labnd(2)-1);
%     l1ind = find(strcmp(lab1,MODELHCSZ0ns.labels{3}));    
%     l2ind = find(strcmp(lab2,MODELHCSZ0ns.labels{3}));
%     
%     comp = i0_2g(ii);
%     
%     id1 = find(MODELHCSZ0ns.X(:,3)== l1ind-1);id2 = find(MODELHCSZ0ns.X(:,3)== l2ind-1);
%    figure('color',[1 1 1]); boxplot([Y0_ukbinit(id1,comp);Y0_ukbinit(id2,comp)],[zeros(length(id1),1);ones(length(id2),1)],'notch','on','labels' ,{lab1 lab2});title(['Comp ' my_pad(comp,2)])
%    export_fig('HCSZicasig_full_UKBinit_preregSite_C30_sigDiag.pdf','-append',gcf);close(gcf)
% end
% sigInd_diag_all_UKBinit_SMRI = find(myFDR(UNI0_AllpreregSite_UKBinit.p{1}(:))<=0.05);
% [i0_1u,i0_2u] = ind2sub(size(UNI0_AllpreregSite_UKBinit.p{1}),sigInd_diag_all_UKBinit_SMRI);
% 
% for ii=1:length(i0_1g)
%     if mod(ii,24)==1
%         figure('color',[1 1 1]);
%         jj=1;
%     end
%     contr =  UNI0_AllpreregSite_UKBinit.stats{1}.Contrast(i0_1g(ii));
%     contr = contr{1};
%     labst = strfind(contr,'(');labnd = strfind(contr,')');
%     lab1 = contr(labst(1)+1:labnd(1)-1);lab2 = contr(labst(2)+1:labnd(2)-1);
%     l1ind = find(strcmp(lab1,MODELHCSZ0ns.labels{3}));    
%     l2ind = find(strcmp(lab2,MODELHCSZ0ns.labels{3}));
%     
%     comp = i0_2g(ii);
%     
%     id1 = find(MODELHCSZ0ns.X(:,3)== l1ind-1);id2 = find(MODELHCSZ0ns.X(:,3)== l2ind-1);
%    figure('color',[1 1 1]); boxplot([Y0_ukbinit(id1,comp);Y0_ukbinit(id2,comp)],[zeros(length(id1),1);ones(length(id2),1)],'notch','on','labels' ,{lab1 lab2});title(['Comp ' my_pad(comp,2)])
%    set(gcs,'FontName','Arial')
%    export_fig('HCSZicasig_full_UKBinit_preregSite_C30_sigDiag.pdf','-append',gcf);close(gcf)
% end
% 
% %% with interaction
% 
% sigInd_diag_allwX_UKBinit_SMRI = find(myFDR(UNI0_AllwXpreregSite_UKBinit.p{2}(:))<=0.05);
% [i0x_1g,i0x_2g] = ind2sub(size(UNI0_AllwXpreregSite_UKBinit.p{2}),sigInd_diag_allwX_UKBinit_SMRI);
% for ii=1:length(i0x_1g)
%     if mod(ii,24)==1
%         figure('color',[1 1 1]);
%         jj=1;
%     end
%     contr =  UNI0_AllwXpreregSite_UKBinit.stats{2}.Contrast(i0x_1g(ii));
%     contr = contr{1};
%     labst = strfind(contr,'(');labnd = strfind(contr,')');
%     lab1 = contr(labst(1)+1:labnd(1)-1);lab2 = contr(labst(2)+1:labnd(2)-1);
%     l1ind = find(strcmp(lab1,MODELHCSZ0ns.labels{3}));    
%     l2ind = find(strcmp(lab2,MODELHCSZ0ns.labels{3}));
%     
%     comp = i0x_2g(ii);
%     
%     id1 = find(MODELHCSZ0ns.X(:,3)== l1ind-1);id2 = find(MODELHCSZ0ns.X(:,3)== l2ind-1);
%    figure('color',[1 1 1]); boxplot([Y0_ukbinit(id1,comp);Y0_ukbinit(id2,comp)],[zeros(length(id1),1);ones(length(id2),1)],'notch','on','labels' ,{lab1 lab2});title(['Comp ' my_pad(comp,2)])
%    export_fig('HCSZicasig_full_wX_UKBinit_preregSite_C30_sigDiag_SMRI.pdf','-append',gcf);close(gcf)
% end
% sigInd_diag_allwX_UKBinit_ALFF = find(myFDR(UNI1_AllwX_preregSite_UKBinit.p{2}(:))<=0.05);
% [i1x_1g,i1x_2g] = ind2sub(size(UNI1_AllwX_preregSite_UKBinit.p{2}),sigInd_diag_allwX_UKBinit_ALFF);
% 
% for ii=1:length(i1x_1g)
%     if mod(ii,24)==1
%         figure('color',[1 1 1]);
%         jj=1;
%     end
%     contr =  UNI1_AllwX_preregSite_UKBinit.stats{2}.Contrast(i1x_1g(ii));
%     contr = contr{1};
%     labst = strfind(contr,'(');labnd = strfind(contr,')');
%     lab1 = contr(labst(1)+1:labnd(1)-1);lab2 = contr(labst(2)+1:labnd(2)-1);
%     l1ind = find(strcmp(lab1,MODELHCSZ1ns.labels{3}));    
%     l2ind = find(strcmp(lab2,MODELHCSZ1ns.labels{3}));
%     
%     comp = i1x_2g(ii);
%     
%     id1 = find(MODELHCSZ1ns.X(:,3)== l1ind-1);id2 = find(MODELHCSZ1ns.X(:,3)== l2ind-1);
%    figure('color',[1 1 1]); boxplot([Y1_ukbinit(id1,comp);Y1_ukbinit(id2,comp)],[zeros(length(id1),1);ones(length(id2),1)],'notch','on','labels' ,{lab1 lab2});title(['Comp ' my_pad(comp,2)])
%    export_fig('HCSZicasig_full_wX_UKBinit_preregSite_C30_sigDiag_ALFF.pdf','-append',gcf);close(gcf)
% end
% 
% 
% 
% sigInd_diag_allwX_UKBinit_DMRI = find(myFDR(UNI2_AllwX_preregSite_UKBinit.p{2}(:))<=0.05);
% [i2x_1g,i2x_2g] = ind2sub(size(UNI2_AllwX_preregSite_UKBinit.p{2}),sigInd_diag_allwX_UKBinit_DMRI);
% for ii=1:length(i2x_1g)
%     if mod(ii,24)==1
%         figure('color',[1 1 1]);
%         jj=1;
%     end
%     contr =  UNI2_AllwX_preregSite_UKBinit.stats{2}.Contrast(i2x_1g(ii));
%     contr = contr{1};
%     labst = strfind(contr,'(');labnd = strfind(contr,')');
%     lab1 = contr(labst(1)+1:labnd(1)-1);lab2 = contr(labst(2)+1:labnd(2)-1);
%     l1ind = find(strcmp(lab1,MODELHCSZ2ns.labels{3}));    
%     l2ind = find(strcmp(lab2,MODELHCSZ2ns.labels{3}));
%     
%     comp = i2x_2g(ii);
%     
%     id1 = find(MODELHCSZ2ns.X(:,3)== l1ind-1);id2 = find(MODELHCSZ2ns.X(:,3)== l2ind-1);
%    figure('color',[1 1 1]); boxplot([Y2_ukbinit(id1,comp);Y2_ukbinit(id2,comp)],[zeros(length(id1),1);ones(length(id2),1)],'notch','on','labels' ,{lab1 lab2});title(['Comp ' my_pad(comp,2)])
%    export_fig('HCSZicasig_full_wX_UKBinit_preregSite_C30_sigDiag_DMRI.pdf','-append',gcf);close(gcf)
% end
% %%
% sigInd_diag_allwX_UKBinit_SMRI = find(myFDR(UNI0_AllwXpreregSite_UKBinit.p{2}(:))<=0.05);
% [i0x_1u,i0x_2u] = ind2sub(size(UNI0_AllwXpreregSite_UKBinit.p{2}),sigInd_diag_allwX_UKBinit_SMRI);
% for ii=1:length(i0x_1u)
%     if mod(ii,24)==1
%         figure('color',[1 1 1]);
%         jj=1;
%     end
%     contr =  UNI0_AllwXpreregSite_UKBinit.stats{2}.Contrast(i0x_1u(ii));
%     contr = contr{1};
%     labst = strfind(contr,'(');labnd = strfind(contr,')');
%     lab1 = contr(labst(1)+1:labnd(1)-1);lab2 = contr(labst(2)+1:labnd(2)-1);
%     l1ind = find(strcmp(lab1,MODELHCSZ0ns.labels{3}));    
%     l2ind = find(strcmp(lab2,MODELHCSZ0ns.labels{3}));
%     
%     comp = i0x_2u(ii);
%     
%     id1 = find(MODELHCSZ0ns.X(:,3)== l1ind-1);id2 = find(MODELHCSZ0ns.X(:,3)== l2ind-1);
%    figure('color',[1 1 1]); boxplot([Y0_ukbinit(id1,comp);Y0_ukbinit(id2,comp)],[zeros(length(id1),1);ones(length(id2),1)],'notch','on','labels' ,{lab1 lab2});title(['Comp ' my_pad(comp,2)])
%    export_fig('HCSZicasig_full_wX_UKBinit_preregSite_C30_sigDiag_SMRI.pdf','-append',gcf);close(gcf)
% end
% sigInd_diag_allwX_UKBinit_ALFF = find(myFDR(UNI1_AllwX_UKBinit.p{1}(:))<=0.05);
% [i1x_1u,i1x_2u] = ind2sub(size(UNI1_AllwX_UKBinit.p{1}),sigInd_diag_allwX_UKBinit_ALFF);
% 
% for ii=1:length(i1x_1u)
%     if mod(ii,24)==1
%         figure('color',[1 1 1]);
%         jj=1;
%     end
%     contr =  UNI1_AllwX_UKBinit.stats{1}.Contrast(i1x_1u(ii));
%     contr = contr{1};
%     labst = strfind(contr,'(');labnd = strfind(contr,')');
%     lab1 = contr(labst(1)+1:labnd(1)-1);lab2 = contr(labst(2)+1:labnd(2)-1);
%     l1ind = find(strcmp(lab1,MODELHCSZ1ns.labels{3}));    
%     l2ind = find(strcmp(lab2,MODELHCSZ1ns.labels{3}));
%     
%     comp = i1x_2u(ii);
%     
%     id1 = find(MODELHCSZ1ns.X(:,3)== l1ind-1);id2 = find(MODELHCSZ1ns.X(:,3)== l2ind-1);
%    figure('color',[1 1 1]); boxplot([Y1_ukbinit(id1,comp);Y1_ukbinit(id2,comp)],[zeros(length(id1),1);ones(length(id2),1)],'notch','on','labels' ,{lab1 lab2});title(['Comp ' my_pad(comp,2)])
%    export_fig('HCSZicasig_full_wX_UKBinit_preregSite_C30_sigDiag_ALFF.pdf','-append',gcf);close(gcf)
% end
% 
% 
% 
% sigInd_diag_allwX_UKBinit_DMRI = find(myFDR(UNI2_AllwX_UKBinit.p{2}(:))<=0.05);
% [i2x_1u,i2x_2u] = ind2sub(size(UNI2_AllwX_UKBinit.p{2}),sigInd_diag_allwX_UKBinit_DMRI);
% for ii=1:length(i2x_1u)
%     if mod(ii,24)==1
%         figure('color',[1 1 1]);
%         jj=1;
%     end
%     contr =  UNI2_AllwX_UKBinit.stats{2}.Contrast(i2x_1u(ii));
%     contr = contr{1};
%     labst = strfind(contr,'(');labnd = strfind(contr,')');
%     lab1 = contr(labst(1)+1:labnd(1)-1);lab2 = contr(labst(2)+1:labnd(2)-1);
%     l1ind = find(strcmp(lab1,MODELHCSZ2ns.labels{3}));    
%     l2ind = find(strcmp(lab2,MODELHCSZ2ns.labels{3}));
%     
%     comp = i2x_2u(ii);
%     
%     id1 = find(MODELHCSZ2ns.X(:,3)== l1ind-1);id2 = find(MODELHCSZ2ns.X(:,3)== l2ind-1);
%    figure('color',[1 1 1]); boxplot([Y2_ukbinit(id1,comp);Y2_ukbinit(id2,comp)],[zeros(length(id1),1);ones(length(id2),1)],'notch','on','labels' ,{lab1 lab2});title(['Comp ' my_pad(comp,2)])
%    export_fig('HCSZicasig_full_wX_UKBinit_preregSite_C30_sigDiag_DMRI.pdf','-append',gcf);close(gcf)
% end
% 
% 
% 
% %%
% sigInd_diag_all_UKBinit_wX_SMRI = find(myFDR(UNI0_AllwXpreregSite_UKBinit.p{1}(:))<=0.05);
% [i00_1g,i00_2g] = ind2sub(size(UNI0_AllwXpreregSite_UKBinit.p{1}),sigInd_diag_all_UKBinit_wX_SMRI);
% 
% sigInd_diag_all_UKBinit_wX_SMRI = find(myFDR(UNI0_AllwXpreregSite_UKBinit.p{1}(:))<=0.05);
% [i00_1u,i00_2u] = ind2sub(size(UNI0_AllwXpreregSite_UKBinit.p{1}),sigInd_diag_all_UKBinit_wX_SMRI);
% 
% 
% pAge_210_UKBinit =[myFDR(UNI2_AllwX_UKBinit.p{4});myFDR(UNI1_AllwX_UKBinit.p{3});myFDR(UNI0_AllwXpreregSite_UKBinit.p{4})];
% pAge_210_UKBinit =[myFDR(UNI2_AllwX_preregSite_UKBinit.p{4});myFDR(UNI1_AllwX_preregSite_UKBinit.p{4});myFDR(UNI0_AllwXpreregSite_UKBinit.p{4})];
% 
% pSex_210_UKBinit =[myFDR(UNI2_AllwX_preregSite_UKBinit.p{1});myFDR(UNI1_AllwX_preregSite_UKBinit.p{1});myFDR(UNI0_AllwXpreregSite_UKBinit.p{1})];
% 
% 
% sigIntInd2_sexXdiag_UKBinit = find(myFDR(UNI2_AllwX_preregSite_UKBinit.p{6}(:))<=0.05);
% [ix21g,ix22g] = ind2sub(size(UNI2_AllwX_preregSite_UKBinit.p{6}),sigIntInd2_sexXdiag_UKBinit);
% 
% %%
% Age = MODELHCSZ0ns.X(:,1);
% cmAge = vals2colormap(Age,'jet',[min(Age) max(Age)]);
% 
% nC=20;
% for cc=1:nC
% figure('color',[1 1 1]);scatter3(Y0_ukbinit(:,cc),Y1_ukbinit(:,cc),Y2_ukbinit(:,cc),[],cmAge);
% view([30 10]); 
% xlabel('SMRI');ylabel('ALFF');zlabel('DMRI');
% title(['Comp ' my_pad(cc,2)]);
% export_fig('HCSZicasig_full_UKBinit_preregSite_C30_byAge.pdf','-append',gcf);
% close all;
% end
% 
% for cc=1:nC
% figure('color',[1 1 1]);scatter3(Y0_ukbinit(:,cc),Y1_ukbinit(:,cc),Y2_ukbinit(:,cc),[],cmAge);
% view([30 10]); 
% xlabel('SMRI');ylabel('ALFF');zlabel('DMRI');
% title(['Comp ' my_pad(cc,2)]);
% export_fig('HCSZicasig_full_UKBinit_preregSite_C30_byAge.pdf','-append',gcf);
% close all;
% end
% 
% %
% 
% site = MODELHCSZ0ns.X(:,5);
% 
% % combined order: Age, sex, diag, rSN , site; 
% % for site: COBRE : 0
% %           FBIRN : 3     7     9    10    12    13    18;
% %           MPRC  : 4   5    6
% %           BSNIP : 1   2
% 
% site(site==4) = 16;site(site==5) = 24; site(site==6) = 30;
% site(site==1) = 48; site(site==2) = 56; 
% site(site==3) = 68;site(site==7) = 74;site(site==9) = 80;site(site==10) = 88;site(site==12) = 94;site(site==13) = 100;site(site==18) = 108;
% 
% % site color code: bottom to top
% % {'COBRE'  'MPRC1' 'MPRC2' 'MPRC3' 'BSNIP1' 'BSNIP2' 'FBIRN3' 'FBIRN7'
% % 'FBIRN9' 'FBIRN10' 'FBIRN12' 'FBIRN13' 'FBIRN18'}
% 
% cmSite = vals2colormap(site,'jet',[min(site) max(site)]);
% 
% for cc=1:nC
% figure('color',[1 1 1]);scatter3(Y0_ukbinit(:,cc),Y1_ukbinit(:,cc),Y2_ukbinit(:,cc),[],cmSite);
% view([30 10]); 
% xlabel('SMRI');ylabel('ALFF');zlabel('DMRI');
% title(['Comp ' my_pad(cc,2)]);
% export_fig('HCSZicasig_full_UKBinit_preregSite_C30_bySite.pdf','-append',gcf);
% close all;
% end
% 
% for cc=1:nC
% figure('color',[1 1 1]);scatter3(Y0_ukbinit(:,cc),Y1_ukbinit(:,cc),Y2_ukbinit(:,cc),[],cmSite);
% view([30 10]); 
% xlabel('SMRI');ylabel('ALFF');zlabel('DMRI');
% title(['Comp ' my_pad(cc,2)]);
% export_fig('HCSZicasig_full_UKBinit_preregSite_C30_bySite.pdf','-append',gcf);
% close all;
% end
% 
% diag = MODELHCSZ0ns.X(:,3);
% % {'HC' 'SZ' 'BP' 'SZA' 'SZR' 'BPR' 'SADBPP' 'SADBPR' 'SADEPP' 'SADEPR'};
% diag(diag==1) = 120;diag(diag==3) = 112; diag(diag==4) = 106; % SZ SZA SZR
% diag(diag==2) = 32;diag(diag==5) = 40; % BPP BPR
% diag(diag==6) = 56; diag(diag==7) = 64; % 'SADBPP' 'SADBPR'
% diag(diag==8) = 88; diag(diag==9) = 94; % SADEPP' 'SADEPR
% %  legend({'HC','BPP', 'BPR','SADBPP', 'SADBPR','SADEPP',
% %  'SADEPR','SZR','SZA','SZ'}) % DIAG color order bottom to top
% cmDiag = vals2colormap(diag,'jet',[min(diag) max(diag)]);
% 
% for cc=1:nC
% figure('color',[1 1 1]);scatter3(Y0_ukbinit(:,cc),Y1_ukbinit(:,cc),Y2_ukbinit(:,cc),[],cmDiag);
% view([30 10]); 
% xlabel('SMRI');ylabel('ALFF');zlabel('DMRI');
% title(['Comp ' my_pad(cc,2)]);
% export_fig('HCSZicasig_full_UKBinit_preregSite_C30_byDiag.pdf','-append',gcf);
% close all;
% end
% 
% for cc=1:nC
% figure('color',[1 1 1]);scatter3(Y0_ukbinit(:,cc),Y1_ukbinit(:,cc),Y2_ukbinit(:,cc),[],cmDiag);
% view([30 10]); 
% xlabel('SMRI');ylabel('ALFF');zlabel('DMRI');
% title(['Comp ' my_pad(cc,2)]);
% export_fig('HCSZicasig_full_UKBinit_preregSite_C30_byDiag.pdf','-append',gcf);
% close all;
% end
% 
% %%
% Age1u=find(strcmp('age' ,MULT0_AllwXpreregSite_UKBinit.final_terms));
% Age2u=find(strcmp('age' ,MULT1_All_UKBinit.final_terms));
% Age3u = find(strcmp('age' ,MULT2_AllwX_UKBinit.final_terms));
% [myFDR(UNI0_AllwXpreregSite_UKBinit.p{Age1u}) <= 0.05;myFDR(UNI1_AllwX_UKBinit.p{Age2u}) <= 0.05;myFDR(UNI2_AllwX_UKBinit.p{Age3u}) <= 0.05]
% sigCompAge_UKBinitwX = find(sum([myFDR(UNI0_AllwXpreregSite_UKBinit.p{Age1u}) <= 0.05;myFDR(UNI1_AllwX_UKBinit.p{Age2u}) <= 0.05;myFDR(UNI2_AllwX_UKBinit.p{Age3u}) <= 0.05])==3);
% %%
% Sex1u=find(strcmp('sex' ,MULT0_AllwXpreregSite_UKBinit.final_terms));
% Sex2u=find(strcmp('sex' ,MULT1_All_UKBinit.final_terms));
% Sex3u = find(strcmp('sex' ,MULT2_AllwX_UKBinit.final_terms));
% [myFDR(UNI0_AllwXpreregSite_UKBinit.p{Sex1u}) <= 0.05;myFDR(UNI2_AllwX_UKBinit.p{Sex3u}) <= 0.05]
% sigCompSex_UKBinitwX = find(sum([myFDR(UNI0_AllwXpreregSite_UKBinit.p{Sex1u}) <= 0.05; myFDR(UNI2_AllwX_UKBinit.p{Sex3u}) <= 0.05])==2);
% %%
% Site1u=find(strcmp('site' ,MULT0_AllwXpreregSite_UKBinit.final_terms));
% Site2u=find(strcmp('site' ,MULT1_All_UKBinit.final_terms));
% Site3u = find(strcmp('site' ,MULT2_AllwX_UKBinit.final_terms));
% %[myFDR(UNI0_AllwX_UKBinit.p{Site1u}) <= 0.05;myFDR(UNI1_AllwX_UKBinit.p{Site2u}) <= 0.05;myFDR(UNI2_AllwX_UKBinit.p{Site3u}) <= 0.05]
% %sigCompSite_UKBinitwX = find(sum([myFDR(UNI0_AllwX_UKBinit.p{Site1u}) <= 0.05; myFDR(UNI1_AllwX_UKBinit.p{Site2u}) <= 0.05;myFDR(UNI2_AllwX_UKBinit.p{Site3u}) <= 0.05])==2);
% 
% 
% %%
% figure;plot(1:length(diag),diag,cmDiag)
%%
Tage_all = [UNI0_AllwX_preregSite_UKBinit.t{3};UNI1_AllwX_preregSite_UKBinit.t{2};UNI2_AllwX_preregSite_UKBinit.t{3}];
Tage_all_SgnSum =  sum(sign(Tage_all));
signMaps = ones(1,30);
signMaps(Tage_all_SgnSum==3) = -1;
%%
fulreg_HCSZ_UKBinitSgnCor = fulreg_HCSZ_UKBinit;
for mm=1:3,fulreg_HCSZ_UKBinitSgnCor.A_ls1{mm} = bsxfun(@times,fulreg_HCSZ_UKBinitSgnCor.icasig{mm},signMaps');end
for mm=1:3,fulreg_HCSZ_UKBinitSgnCor.A_ls1{mm} = bsxfun(@times,fulreg_HCSZ_UKBinit.A_ls1{mm},signMaps);end
for mm=1:3,fulreg_HCSZ_UKBinitSgnCor.A_ls1_ox{mm} = bsxfun(@times,fulreg_HCSZ_UKBinit.A_ls1_ox{mm},signMaps);end
for mm=1:3,fulreg_HCSZ_UKBinitSgnCor.A_ls2{mm} = bsxfun(@times,fulreg_HCSZ_UKBinit.A_ls2{mm},signMaps);end
for mm=1:3,fulreg_HCSZ_UKBinitSgnCor.A_ls2_ox{mm} = bsxfun(@times,fulreg_HCSZ_UKBinit.A_ls2_ox{mm},signMaps);end
for mm=1:3,fulreg_HCSZ_UKBinitSgnCor.A_pinv{mm} = bsxfun(@times,fulreg_HCSZ_UKBinit.A_pinv{mm},signMaps);end
for mm=1:3,fulreg_HCSZ_UKBinitSgnCor.A_tr{mm} = bsxfun(@times,fulreg_HCSZ_UKBinit.A_tr{mm},signMaps);end
for mm=1:3,fulreg_HCSZ_UKBinitSgnCor.W{mm} = bsxfun(@times,fulreg_HCSZ_UKBinit.W{mm},signMaps');end
save HCSZMRI_MMIVAfull_w0GICA_C020_results_new_P000_preregSITE_signCorrected fulreg_HCSZ_UKBinitSgnCor
%%
Y0_prereg_ukbinit_sgnCor = bsxfun(@times,Y0_prereg_ukbinit,signMaps);
Y1_prereg_ukbinit_sgnCor = bsxfun(@times,Y1_prereg_ukbinit,signMaps);
Y2_prereg_ukbinit_sgnCor = bsxfun(@times,Y2_prereg_ukbinit,signMaps);

%%
nC=30;
for cc= sigCompAGEfdr % [1 2 5 6 8 9 13 19]
    figure('color',[1 1 1]);
    %scatter3(Yful_SMRI(cc,:)',Yful_ALFF(cc,:)',Yful_DMRI(cc,:)',[],cmAge);
    scatter3(Y0_prereg_ukbinit_sgnCor(x0id,cc),Y1_prereg_ukbinit_sgnCor(x0id,cc),Y2_prereg_ukbinit_sgnCor(x0id,cc),[],cmAge(x0id,:));
    hold on;
    scatter3(Y0_prereg_ukbinit_sgnCor(x1id,cc),Y1_prereg_ukbinit_sgnCor(x1id,cc),Y2_prereg_ukbinit_sgnCor(x1id,cc),[],cmAge(x1id,:),'filled');
    view([30 10]);
    xlabel('SMRI');ylabel('ALFF');zlabel('DMRI');
    title(['Comp ' my_pad(cc,2)]);
    
    %set(cb,'FontName','Arial')
    set(gca,'FontName','Arial')
    export_fig('HCSZ_icasig_full_UKBinit_regsite_C30_bysigfdrAge_sgnCor.pdf','-append',gcf);
    close all;
end
%%
for cc = 20
    figure('color',[1 1 1]); 
    %scatter3(Y0_prereg_ukbinit_sgnCor(HCid,cc),Y1_prereg_ukbinit_sgnCor(HCid,cc),Y2_prereg_ukbinit_sgnCor(HCid,cc),[],'b');%cmAge(HCid,:));
    %hold on;
    scatter3(Y0_prereg_ukbinit_sgnCor(SZid,cc),Y1_prereg_ukbinit_sgnCor(SZid,cc),Y2_prereg_ukbinit_sgnCor(SZid,cc),[],'r');%cmAge(x1id,:),'x');
    hold on;
    scatter3(Y0_prereg_ukbinit_sgnCor(BPid,cc),Y1_prereg_ukbinit_sgnCor(BPid,cc),Y2_prereg_ukbinit_sgnCor(BPid,cc),[],'g');
    view([30 10]);
    xlabel('SMRI');ylabel('ALFF');zlabel('DMRI');
    title(['Comp ' my_pad(cc,2)]);
    
    %set(cb,'FontName','Arial')
    set(gca,'FontName','Arial')
    %export_fig('UKBicasig_full_prereg_UKBinit_C30_bysigfdrAge.pdf','-append',gcf);
    
    %close all;
end
%%
odir ='/media/eswar/New Volume/Documents/MISA/HCSZ/MMIVA2step_wUKBinit/C030_regsite_ukbinit/figs/stats';%'/media/eswar/New Volume/Documents/MISA/HCSZ/MMIVA2step_wUKBinit/C030_regsite_ukbinit/figs/stats_png'; 
plotuni(aa,dm,1,'age',a00_contr,a00_comp,a00_T,signMaps,odir)
plotuni(aa,dm,2,'age',a11_contr,a11_comp,a11_T,signMaps,odir)
plotuni(aa,dm,3,'age',a22_contr,a22_comp,a22_T,signMaps,odir)

plotuni(aa,dm,1,'diagnosis',d00_contr,d00_comp,d00_T,signMaps,odir)
plotuni(aa,dm,2,'diagnosis',d11_contr,d11_comp,d11_T,signMaps,odir)
plotuni(aa,dm,3,'diagnosis',d22_contr,d22_comp,d22_T,signMaps,odir)

plotuni(aa,dm,2,'rSN_ALFF',rsn11_contr,rsn11_comp,rsn11_T,signMaps,odir)
%%

fulreg_HCSZ_UKBinit_A_ls1_sgnCor = fulreg_HCSZ_UKBinit.A_ls1;
fulreg_HCSZ_UKBinit_A_ls1_sgnCor{1} = bsxfun(@times,fulreg_HCSZ_UKBinit_A_ls1_sgnCor{1},signMaps);
fulreg_HCSZ_UKBinit_A_ls1_sgnCor{2} = bsxfun(@times,fulreg_HCSZ_UKBinit_A_ls1_sgnCor{2},signMaps);
fulreg_HCSZ_UKBinit_A_ls1_sgnCor{3} = bsxfun(@times,fulreg_HCSZ_UKBinit_A_ls1_sgnCor{3},signMaps);
%%
hdr = spm_vol(fullfile(pwd,'MMIVA2step_wUKBinit/C030_regsite_ukbinit/outs_A_ls1_DMRI_ful.nii'));
dta1 = spm_read_vols(hdr);
szdta = size(dta1);
dta1 = reshape(dta1,prod(szdta(1:3)),szdta(4));
dta1 =  bsxfun(@times,dta1,signMaps);
dta1 = reshape(dta1,szdta);
outname = fullfile(pwd,'MMIVA2step_wUKBinit/C030_regsite_ukbinit/outs_A_ls1_DMRI_ful_signCorr.nii');
create_4DNiftifile(outname,dta1,hdr(1).mat);
%%


