addpath('/data/users2/eswar/software/utilities/')
%addpath('/home/eswar/Documents/MATLAB/MANCOVAtb')
%addpath('/home/eswar/Documents/MATLAB/spm12/');
%addpath(genpath('/home/eswar/Documents/MATLAB/GroupICATv4.0b/'));
addpath(genpath('/data/users2/eswar/temp/MISA_analysis/MISA-master/'));
addpath('/data/users2/eswar/software/export_fig')
addpath('/home/eswar/Documents/MATLAB/spm12/');
scrpth = '/data/users2/eswar/temp/MISA_analysis/UKB/postMMIVAresults/scripts';
origscrdir = '/data/users2/eswar/temp/MISA_analysis/UKB/scripts';
demogdir = '/data/users2/eswar/temp/MISA_analysis/UKB/demographics';
ddir1 = '/data/users2/eswar/temp/MISA_analysis/UKB/postMMIVAresults/C030_UKB_regSite';
id2907 = load(fullfile(demogdir,'subID_2907.mat'));
load(fullfile(demogdir, 'DEMO_2907_curated.mat'))
%%
%demog = readtable('/home/eswar/Documents/projects/MISA/UKB/demographics/ukb_unaffected_clean_select_log_with_sex.dat');
%slist = load('/home/eswar/Documents/projects/MISA/UKB/demographics/goodID_UKBiobank_3497MISAfinal.txt');
%nS = length(slist);
%tab = readtable('/home/eswar/Documents/projects/MISA/UKB/demographics/ukb_unaffected_clean_select_log_with_sex.dat','Delimiter','\t','FileType','text');
%tab = tab(ismember(tab.eid,slist),:);
%tab.sex_f31_0_0 = double(strcmpi(tab.sex_f31_0_0,'Female'));
%tab.genetic_sex_f22001_0_0 = double(strcmpi(tab.genetic_sex_f22001_0_0,'Female'));


%%
rSNall = load(fullfile(origscrdir,'spatNorm_3497UKB.mat'));
motUKB = load(fullfile(origscrdir,'motpar_summary.mat'));

%%
fid=fopen(origscrdir,'goodID_UKBiobank_3500.txt'),'r');
g3500=textscan(fid,'%s');fclose(fid);
fid=fopen(origscrdir),'goodID_UKBiobank_3497MISAfinal.txt'),'r');
g3497=textscan(fid,'%s');fclose(fid);
g3500=g3500{1};
g3497=g3497{1};
[tfu locu] = ismember(g3497,g3500);
%% FD
for ii=1:length(locu)
    FDaf = compute_FD(motUKB.motpars_all_nfilt(:,:,locu(ii)),2);
    FDa = compute_FD(motUKB.motpars_all(:,:,locu(ii)),2);
    if ii==1
        FDall_UKB3497 = zeros(size(motUKB.motpars_all_nfilt,1)-1,length(locu));
        FDfall_UKB3497 = zeros(size(motUKB.motpars_all_nfilt,1)-1,length(locu));
        mxMot_UKB3497 = zeros(length(locu),1);
    end
    FDall_UKB3497(:,ii) = FDa;FDfall_UKB3497(:,ii) = FDaf;
    mxMot_UKB3497(ii) = max(max(motUKB.motpars_all_nfilt(:,:,locu(ii))));
end
meanFD_notchfiltmp_UKB3497 = mean(FDfall_UKB3497)';
meanFD_mp_UKB3497 = mean(FDall_UKB3497)';

meanFD_mp_UKB2907 = meanFD_mp_UKB3497(id2907.lt5NanID);
meanFD_nfiltmp_UKB2907 = meanFD_notchfiltmp_UKB3497(id2907.lt5NanID);

rSN_ALFF_UKB2907 = rSNall.rALFF(id2907.lt5NanID);
rSN_SMRI_UKB2907 = rSNall.rSMRI(id2907.lt5NanID);
rSN_DMRI_UKB2907 = rSNall.rDMRI(id2907.lt5NanID);




%%
%SMRI
%[num2str((1:54)') repmat('  ', 54,1) strvcat(DEMO.demo_col_names)]
ful = load(fullfile(ddir1,'UKBMRI_MMIVAfull_w0GICA_C030_results_goodSub2907_new_P000.mat'));% loads ful

%%

keepHeaderCol = DEMO_2907.demo_col_names(1:54);
%keepHeaderCol = DEMO.demo_col_names([6 21 24 25 26 52:54]);

%%
PEind =[4:18 28:31 38:46];
Xpe = DEMO0_ful.full(:,PEind);
PElabs = DEMO0_ful.names(PEind);
[COEFFpe, SCOREpe, LATENTpe, TSQUAREDpe, EXPLAINEDpe] = pca(Xpe);
figure;plot(COEFFpe(:,1),COEFFpe(:,2),'*'); %plotting first and second prin comps
gname %select outs pca

figure;biplot(COEFFpe(:,1:2),'Scores',SCOREpe(:,1:2),'Varlabels',PElabs);
PEpcainds = [9  11  19  6  13 2  17 4 ];
PElabs(PEpcainds)
figure;imagesc(corr(Xpe(:,PEpcainds)));
figure;pareto(EXPLAINEDpe); % shows first 8 PC's capture most variance
Xpe_pcBEST =  SCOREpe(:,1:8);

%%
ageORyearIND = [find(contains(DEMO0_ful.names,'age')>0) find(contains(DEMO0_ful.names,'years')>0)];
% 26    {'age_first_had_sexual_intercourse_f2139_2_0'             }
% 27    {'age_started_wearing_glasses_or_contact_lenses_f2217_2_0'}
% 37    {'age_when_attended_assessment_centre_f21003_2_0'         }
% 47    {'years_since_first_sexual_intercourse_2_0'               }
% 48    {'years_since_started_wearing_glasses_2_0'                }

% Correlation between these
%         26         27      37          47       48  
% 26|    1.0000   -0.0597    0.1773   -0.0761    0.1221
% 27|   -0.0597    1.0000    0.1257    0.1306   -0.8777
% 37|    0.1773    0.1257    1.0000    0.7369    0.2963
% 47|   -0.0761    0.1306    0.7369    1.0000    0.1853
% 48|    0.1221   -0.8777    0.2963    0.1853    1.0000
% DROP  columns 47,48 as they are highly correlated with other age
% variables

ageDRPind = [47 48];
%%
nonPEind_age = setdiff(1:length(rDEMO),[PEind ageDRPind]);

drpFIind = [36 53];
drpLOGind = [2 49 50 51 ];

finalGOODind = setdiff(1:length(rDEMO),[PEind ageDRPind drpFIind drpLOGind]);
%finalGOODind = finalGOODind(1:end-1);


%%
Yful_SMRI = ful.icasig{1};


[DEMO0s_ful, MULT0s_ful, UNI0s_ful] = run_model_UKB_wX(MODELUKB0s_ful, 0.01,Yful_SMRI', [], [], 0);
save UKB_MMIVA_C30_preregSite_SMRI_MancovanOuts_wX_FINAL  *0s_ful 

Yful_ALFF = ful.icasig{2};
C30_outs11 = load(fullfile(fileparts(ddir),'UKB_MMIVA_C30_ALFF_MancovanOuts_FINAL'));

MODELUKB1s_ful = C30_outs11.MODELUKB1_ful;

[DEMO1s_ful, MULT1s_ful, UNI1s_ful] = run_model_UKB_wX(MODELUKB1s_ful, 0.01,Yful_ALFF', [], [], 0);
save UKB_MMIVA_C30_preregSite_ALFF_MancovanOuts_wX_FINAL  *1s_ful


Yful_DMRI = ful.icasig{3};
C30_outs22 = load(fullfile(fileparts(ddir),'UKB_MMIVA_C30_DMRI_MancovanOuts_FINAL'));

MODELUKB2s_ful = C30_outs22.MODELUKB2_ful;

[DEMO2s_ful, MULT2s_ful, UNI2s_ful] = run_model_UKB_wX(MODELUKB2s_ful, 0.01,Yful_DMRI', [], [], 0);
save UKB_MMIVA_C30_preregSite_DMRI_MancovanOuts_wX_FINAL  *2s_ful

%%
%XX = DEMO0_ful.full;
% XX=UNI0_ful.stats{2}.X;
% colcnt=cell(size(XX,2),1);
% for col =1:size(XX,2),
%     leng=length(unique(XX(:,col)));
%     if leng<17, 
%         uv=unique(XX(:,col));
%         for uu=1:length(uv) 
%            cnt(uu)=sum(XX(:,col)==uv(uu));
%         end
%     end
%     if exist('cnt','var'),colcnt{col}=cnt;clear cnt;end
% end
% Xdrp = []; 
% for cc=1:length(colcnt) 
%     if ~isempty(colcnt{cc})
%         uv=unique(XX(:,cc));
%         cnt=colcnt{cc};
%         lt5=find(cnt<6);
%         if ~isempty(lt5), 
%              for ll=1:length(lt5)
%                       Xdrp = [Xdrp; find(XX(:,cc)==uv(lt5(ll)))];
%              end
%         end
%      end
% end
% Xdrp = unique(Xdrp);
%%
%Xkeep = setdiff(1:
sigTerms = intersect(MULT0s_ful.final_terms,intersect(MULT1s_ful.final_terms,MULT2s_ful.final_terms))';

%%
% SMRI
[sx00_comp,sx00_T,sx00_contr]= printouts(UNI0s_ful,find(strcmp(UNI0s_ful.tests,sigTerms{4})>0));
[tta00_comp,tta00_T,tta00_contr]= printouts(UNI0s_ful,find(strcmp(UNI0s_ful.tests,sigTerms{6})>0));
[a00_comp,a00_T,a00_contr] = printouts(UNI0s_ful,find(strcmp(UNI0s_ful.tests,sigTerms{1})>0));
[fi00_comp,fi00_T,fi00_contr]= printouts(UNI0s_ful,find(strcmp(UNI0s_ful.tests,sigTerms{2})>0));
[rsn00_comp,rsn00_T,rsn00_contr]=printouts(UNI0s_ful,find(strcmp(UNI0s_ful.tests,sigTerms{3})>0));
[sxa01_comp,sxa01_T,sxa01_contr]=printouts(UNI0s_ful,find(strcmp(UNI0s_ful.tests,sigTerms{5})>0));

% ALFF
[sx11_comp,sx11_T,sx11_contr]= printouts(UNI1s_ful,find(strcmp(UNI1s_ful.tests,sigTerms{4})>0));
[tta11_comp,tta11_T,tta11_contr]= printouts(UNI1s_ful,find(strcmp(UNI1s_ful.tests,sigTerms{6})>0));
[a11_comp,a11_T,a11_contr] = printouts(UNI1s_ful,find(strcmp(UNI1s_ful.tests,sigTerms{1})>0));
[fi11_comp,fi11_T,fi11_contr] = printouts(UNI1s_ful,find(strcmp(UNI1s_ful.tests,sigTerms{2})>0));
[rsn11_comp,rsn11_T,rsn11_contr]=printouts(UNI1s_ful,find(strcmp(UNI1s_ful.tests,sigTerms{3})>0));
[sxa12_comp,sxa12_T,sxa12_contr]=printouts(UNI1s_ful,find(strcmp(UNI1s_ful.tests,sigTerms{5})>0));

%
% DMRI
[sx22_comp,sx22_T,sx22_contr]= printouts(UNI2s_ful,find(strcmp(UNI2s_ful.tests,sigTerms{4})>0));
[tta22_comp,tta22_T,tta22_contr]= printouts(UNI2s_ful,find(strcmp(UNI2s_ful.tests,sigTerms{6})>0));
[a22_comp,a22_T,a22_contr] = printouts(UNI2s_ful,find(strcmp(UNI2s_ful.tests,sigTerms{1})>0));
[fi22_comp,fi22_T,fi22_contr] = printouts(UNI2s_ful,find(strcmp(UNI2s_ful.tests,sigTerms{2})>0));
[rsn22_comp,rsn22_T,rsn22_contr]=printouts(UNI2s_ful,find(strcmp(UNI2s_ful.tests,sigTerms{3})>0));
[sxa23_comp,sxa23_T,sxa23_contr]=printouts(UNI2s_ful,find(strcmp(UNI2s_ful.tests,sigTerms{5})>0));
%%
all3ageSigInd = intersect(a11_comp{1},intersect(a00_comp{1},a22_comp{1}));

T0_ag_all3Sig = a00_T{1}(ismember(a00_comp{1},all3ageSigInd ));
T1_ag_all3Sig = a11_T{1}(ismember(a11_comp{1},all3ageSigInd));
T2_ag_all3Sig = a22_T{1}(ismember(a22_comp{1},all3ageSigInd));
signComps = ones(1,30);
signComps([1 3 10 12 16  17  19 26])=-1;% making sure age slope is -ve
disp('---- age -----');
[all3ageSigInd;T0_ag_all3Sig.*signComps(all3ageSigInd);T1_ag_all3Sig.*signComps(all3ageSigInd);T2_ag_all3Sig.*signComps(all3ageSigInd)]
%%
sxAllSigID = intersect(sx22_comp{1},intersect(sx00_comp{1},sx11_comp{1}));
T0_sx_all3Sig = sx00_T{1}(ismember(sx00_comp{1},sxAllSigID));
T1_sx_all3Sig = sx11_T{1}(ismember(sx11_comp{1},sxAllSigID));
T2_sx_all3Sig = sx22_T{1}(ismember(sx22_comp{1},sxAllSigID));
disp('---- sex -----');
[sxAllSigID;T0_sx_all3Sig;T1_sx_all3Sig;T2_sx_all3Sig]
%%
ttaAllSigID = intersect(tta22_comp{1},intersect(tta00_comp{1},tta11_comp{1}));
T0_tta_all3Sig = tta00_T{1}(ismember(tta00_comp{1},ttaAllSigID));
T1_tta_all3Sig = tta11_T{1}(ismember(tta11_comp{1},ttaAllSigID));
T2_tta_all3Sig = tta22_T{1}(ismember(tta22_comp{1},ttaAllSigID));
disp('---- tta -----');
[ttaAllSigID;T0_tta_all3Sig;T1_tta_all3Sig;T2_tta_all3Sig]
%%
fiAllSigID = intersect(fi22_comp{1},intersect(fi00_comp{1},fi11_comp{1}));
T0_fi_all3Sig = fi00_T{1}(ismember(fi00_comp{1},fiAllSigID));
T1_fi_all3Sig = fi11_T{1}(ismember(fi11_comp{1},fiAllSigID));
T2_fi_all3Sig = fi22_T{1}(ismember(fi22_comp{1},fiAllSigID));
disp('---- fi -----');
[fiAllSigID;T0_fi_all3Sig;T1_fi_all3Sig;T2_fi_all3Sig]
%%
aa{1} = UNI0s_ful;
aa{2} = UNI1s_ful;
aa{3} = UNI2s_ful;

dm{1} = DEMO0s_ful;
dm{2} = DEMO1s_ful;
dm{3} = DEMO2s_ful;

%%
odir = fullfile(ddir1,'figs/stats');
plotuni(aa,dm,1,sigTerms{1},a00_contr,a00_comp,a00_T,signComps,odir)
plotuni(aa,dm,2,sigTerms{1},a11_contr,a11_comp,a11_T,signComps,odir)
plotuni(aa,dm,3,sigTerms{1},a22_contr,a22_comp,a22_T,signComps,odir)

plotuni(aa,dm,1,sigTerms{2},fi00_contr,fi00_comp,fi00_T,signComps,odir)
plotuni(aa,dm,2,sigTerms{2},fi11_contr,fi11_comp,fi11_T,signComps,odir)
plotuni(aa,dm,3,sigTerms{2},fi22_contr,fi22_comp,fi22_T,signComps,odir)

plotuni(aa,dm,1,sigTerms{4},sx00_contr,sx00_comp,sx00_T,signComps,odir)
plotuni(aa,dm,2,sigTerms{4},sx11_contr,sx11_comp,sx11_T,signComps,odir)
plotuni(aa,dm,3,sigTerms{4},sx22_contr,sx22_comp,sx22_T,signComps,odir)

plotuni(aa,dm,1,sigTerms{6},tta00_contr,tta00_comp,tta00_T,signComps,odir)
plotuni(aa,dm,2,sigTerms{6},tta11_contr,tta11_comp,tta11_T,signComps,odir)
plotuni(aa,dm,3,sigTerms{6},tta22_contr,tta22_comp,tta22_T,signComps,odir)

%%


%%
dofi = {'sex_f31_0_0' 'age_when_attended_assessment_centre_f21003_2_0' 'fluid_intelligence_score_f20016_2_0' 'log_pm_score_2_0' 'time_to_answer_f4288_2_0' 'rSpatNorm'  'meanFD'  };
[strvcat( MULT0_ful.final_terms) repmat('  ',length(MULT0_ful.t),1) num2str(MULT0_ful.t)]
[strvcat( MULT1_ful.final_terms) repmat('  ',length(MULT1_ful.t),1) num2str(MULT1_ful.t)]
[strvcat( MULT2_ful.final_terms) repmat('  ',length(MULT2_ful.t),1) num2str(MULT2_ful.t)]
%%
ageVar = 'age_when_attended_assessment_centre_f21003_2_0';
tta1=find(strcmp(ageVar ,MULT0_ful.final_terms));
tta2=find(strcmp(ageVar ,MULT1_ful.final_terms));
tta3 = find(strcmp(ageVar ,MULT2_ful.final_terms));
[myFDR(UNI0_ful.p{tta1}) <= 0.05;myFDR(UNI1_ful.p{tta2}) <= 0.05;myFDR(UNI2_ful.p{tta3}) <= 0.05]
sigCompAGEfdr = find(sum([myFDR(UNI0_ful.p{tta1}) <= 0.05;myFDR(UNI1_ful.p{tta2}) <= 0.05;myFDR(UNI2_ful.p{tta3}) <= 0.05])==3);
sigCompAGEbf = find(sum([(UNI0_ful.p{tta1}) <= 0.05/20;(UNI1_ful.p{tta2}) <= 0.05/20;(UNI2_ful.p{tta3}) <= 0.05/20])==3);

%%
sexVar = 'sex_f31_0_0';
ttx1=find(strcmp(sexVar ,MULT0s_ful.final_terms));
ttx2=find(strcmp(sexVar ,MULT1s_ful.final_terms));
ttx3 = find(strcmp(sexVar ,MULT2s_ful.final_terms));
[myFDR(UNI0s_ful.p{ttx1}) <= 0.05;myFDR(UNI1s_ful.p{ttx2}) <= 0.05;myFDR(UNI2s_ful.p{ttx3}) <= 0.05]
sigCompSEXfdr = find(sum([myFDR(UNI0s_ful.p{ttx1}) <= 0.05;myFDR(UNI1s_ful.p{ttx2}) <= 0.05;myFDR(UNI2s_ful.p{ttx3}) <= 0.05])==3);
sigCompSEXbf = find(sum([(UNI0s_ful.p{ttx1}) <= 0.05/20;(UNI1s_ful.p{ttx2}) <= 0.05/20;(UNI2s_ful.p{ttx3}) <= 0.05/20])==3);
%%
tt1=find(strcmp('time_to_answer_f4288_2_0' ,MULT0s_ful.final_terms));
tt2=find(strcmp('time_to_answer_f4288_2_0' ,MULT1s_ful.final_terms));
tt3 = find(strcmp('time_to_answer_f4288_2_0' ,MULT2s_ful.final_terms));
[myFDR(UNI0s_ful.p{tt1}) <= 0.05;myFDR(UNI1s_ful.p{tt2}) <= 0.05;myFDR(UNI2s_ful.p{tt3}) <= 0.05]
sigCompTTA = find(sum([myFDR(UNI0s_ful.p{tt1}) <= 0.05/20;myFDR(UNI1s_ful.p{tt2}) <= 0.05/20;myFDR(UNI2s_ful.p{tt3}) <= 0.05/20])==3);
%%
isM_hdr  = spm_vol(fullfile(ddir,'outs_A_ls1_SMRI_ful.nii'));
iaM_hdr  = spm_vol(fullfile(ddir,'outs_A_ls1_ALFF_ful.nii'));
idM_hdr  = spm_vol(fullfile(ddir,'outs_A_ls1_DMRI_ful.nii'));

dta = spm_read_vols(isM_hdr);
dtas = bsxfun(@times,reshape(dta,prod(isM_hdr(1).dim),30),signComps);
dtas = reshape(dtas,[isM_hdr(1).dim ,30]);
outName = 'outs_A_ls1_SMRI_ful_sgnCorr.nii';
create_4DNiftifile(outName,dtas,isM_hdr(1).mat)

dta = spm_read_vols(iaM_hdr);
dtas = bsxfun(@times,reshape(dta,prod(isM_hdr(1).dim),30),signComps);
dtas = reshape(dtas,[isM_hdr(1).dim ,30]);
outName = 'outs_A_ls1_ALFF_ful_sgnCorr.nii';
create_4DNiftifile(outName,dtas,isM_hdr(1).mat)

dta = spm_read_vols(idM_hdr);
dtas = bsxfun(@times,reshape(dta,prod(isM_hdr(1).dim),30),signComps);
dtas = reshape(dtas,[isM_hdr(1).dim ,30]);
outName = 'outs_A_ls1_DMRI_ful_sgnCorr.nii';
create_4DNiftifile(outName,dtas,isM_hdr(1).mat)
clear dta dtas 
%%
ageCol=find(strcmp('age_when_attended_assessment_centre_f21003_2_0',MODELUKB0s_ful.names));
Age = MODELUKB0s_ful.X(:,ageCol);
cmAge = vals2colormap(Age,'jet',[min(Age) max(Age)]);

sexCol = find(strcmp('sex_f31_0_0',MODELUKB0s_ful.names));
sex = MODELUKB0s_ful.X(:,sexCol);
cmSex = vals2colormap(sex,'jet',[min(sex) max(sex)]);
x0id = find(sex==0);
x1id = find(sex==1);
%%
nC=30;
for cc= all3ageSigInd % [ 2 3 11 15  18  21 23]
    figure('color',[1 1 1]);
    %scatter3(Yful_SMRI(cc,:)',Yful_ALFF(cc,:)',Yful_DMRI(cc,:)',[],cmAge);
    scatter3(Yful_SMRI(cc,x0id)',Yful_ALFF(cc,x0id)',Yful_DMRI(cc,x0id)',[],cmAge(x0id,:));
    hold on;
    scatter3(Yful_SMRI(cc,x1id)',Yful_ALFF(cc,x1id)',Yful_DMRI(cc,x1id)',[],cmAge(x1id,:),'x');
    view([30 10]);
    xlabel('SMRI');ylabel('ALFF');zlabel('DMRI');
    title(['Comp ' my_pad(cc,2)]);
    
    %set(cb,'FontName','Arial')
    set(gca,'FontName','Arial')
    export_fig(fullfile(ddir,'figs','UKBicasig_full_GICAinit_C30_preregSite_bysigfdrAge.pdf'),'-append',gcf);
    %close all;
end
%%
POS =[104         450        1572         413];% [805   186   659   689];
POSS = [99          31        1207         820];
for cc=sigCompSEXfdr
    figure('color',[1 1 1],'position',POS);
    subplot(1,3,1);plot(Age(sex==0),Yful_SMRI(cc,sex==0)','bo',Age(sex==1),Yful_SMRI(cc,sex==1)','rx');lsline;title(['SMRI: Comp ' my_pad(cc,2)]);xlabel('Age');ylabel('a.u.');legend({'M','F'})
    subplot(1,3,2);plot(Age(sex==0),Yful_ALFF(cc,sex==0)','bo',Age(sex==1),Yful_ALFF(cc,sex==1)','rx');lsline;title(['ALFF: Comp ' my_pad(cc,2)])
    subplot(1,3,3);plot(Age(sex==0),Yful_DMRI(cc,sex==0)','bo',Age(sex==1),Yful_DMRI(cc,sex==1)','rx');lsline;title(['DMRI: Comp ' my_pad(cc,2)])
    
    set(gca,'FontName','Arial')
    export_fig('UKBicasig_full_GICAinit_C30_preregSite_bysigfdrSexN.pdf','-append',gcf);
end
%%
for cc=9%1:20
    figure('color',[1 1 1]);
    % subplot(3,1,1);
    scatter3(Yful_SMRI(cc,:)',MODELUKB0_ful.X(:,37),MODELUKB0_ful.X(:,32),[],cmTTA);xlabel('SMRI');ylabel('Age');zlabel('TTA');
    view([0 2]);
    figure('color',[1 1 1]);
    %subplot(3,1,2);
    scatter3(Yful_ALFF(cc,:)',MODELUKB0_ful.X(:,37),MODELUKB0_ful.X(:,32),[],cmTTA);xlabel('ALFF');ylabel('Age');zlabel('TTA');
    %subplot(3,1,3);
    figure('color',[1 1 1]);
    scatter3(Yful_DMRI(cc,:)',MODELUKB0_ful.X(:,37),MODELUKB0_ful.X(:,32),[],cmTTA);xlabel('DMRI');ylabel('Age');zlabel('TTA');
end  
%%
ttaCol=find(strcmp('time_to_answer_f4288_2_0' ,MODELUKB0s_ful.names));
tta = MODELUKB0s_ful.X(:,ttaCol);
cmTTA = vals2colormap(tta,'jet',[min(tta) max(tta)]);
for cc=ttaAllSigID
figure('color',[1 1 1],'position',POSS);
subplot(2,3,1:3)
scatter3(Yful_SMRI(cc,:)',Yful_ALFF(cc,:)',-1*Yful_DMRI(cc,:)',[],cmTTA);
%hold on;
%scatter3(Yful_SMRI(cc,sex==1)',Yful_ALFF(cc,sex==1)',Yful_DMRI(cc,sex==1)',[],cmAge(sex==1,:),'x');
view([30 10]); 
xlabel('SMRI');ylabel('ALFF');zlabel('DMRI');
title(['Comp ' my_pad(cc,2)]);

%set(cb,'FontName','Arial')
set(gca,'FontName','Arial')

subplot(2,3,4);plot(tta,Yful_SMRI(cc,:)','ko');title('SMRI');xlabel('TTA (ms)');ylabel('source [a.u.]')
subplot(2,3,5);plot(tta,Yful_ALFF(cc,:)','ko');title('ALFF');xlabel('TTA (ms)');ylabel('source [a.u.]')
subplot(2,3,6);plot(tta,Yful_DMRI(cc,:)','ko');title('DMRI');xlabel('TTA (ms)');ylabel('source [a.u.]')

export_fig('UKBicasig_full_GICAinit_C30_preregSite_bysigfdrTTA_wSub.pdf','-append',gcf);

close all;
end

%%
fi1=find(strcmp('fluid_intelligence_score_f20016_2_0' ,MULT0s_ful.final_terms));
fi2=find(strcmp('fluid_intelligence_score_f20016_2_0' ,MULT1s_ful.final_terms));
fi3 = find(strcmp('fluid_intelligence_score_f20016_2_0' ,MULT2s_ful.final_terms));
[myFDR(UNI0s_ful.p{fi1}) <= 0.05;myFDR(UNI1s_ful.p{fi2}) <= 0.05;myFDR(UNI2s_ful.p{fi3}) <= 0.05]
sigCompFI = find(sum([myFDR(UNI0s_ful.p{fi1}) <= 0.05;myFDR(UNI1s_ful.p{fi2}) <= 0.05;myFDR(UNI2s_ful.p{fi3}) <= 0.05])==3);
sigCompFIb = find(sum([(UNI0s_ful.p{fi1} <= 0.05/20);(UNI1s_ful.p{fi2} <= 0.05/20);(UNI2s_ful.p{fi3} <= 0.05/20)])==3);
%%
fintCol =  find(strcmp('fluid_intelligence_score_f20016_2_0' ,MODELUKB0s_ful.names));
fi = MODELUKB0s_ful.X(:,fintCol);
cmFI = vals2colormap(fi,'jet',[min(fi) max(fi)]);
for cc=sigCompFI
figure('color',[1 1 1],'position',POSS);
subplot(2,3,1:3)
scatter3(Yful_SMRI(cc,:)',Yful_ALFF(cc,:)',Yful_DMRI(cc,:)',[],cmFI);
%hold on;
%scatter3(Yful_SMRI(cc,sex==1)',Yful_ALFF(cc,sex==1)',Yful_DMRI(cc,sex==1)',[],cmAge(sex==1,:),'x');
view([30 10]); 
xlabel('SMRI');ylabel('ALFF');zlabel('DMRI');
title(['Comp ' my_pad(cc,2)]);
    subplot(2,3,4);plot(fi,Yful_SMRI(cc,:)','bo');lsline;title(['SMRI: Comp ' my_pad(cc,2)]);
    subplot(2,3,5);plot(fi,Yful_ALFF(cc,:)','bo');lsline;title(['ALFF: Comp ' my_pad(cc,2)])
    subplot(2,3,6);plot(fi,Yful_DMRI(cc,:)','bo');lsline;title(['DMRI: Comp ' my_pad(cc,2)])
    xlabel('Fluid Intelligence Score');ylabel('source intensity a.u.')
    set(gca,'FontName','Arial')

%set(cb,'FontName','Arial')
set(gca,'FontName','Arial')
export_fig('UKBicasig_full_GICAinit_C30_preregsite_bysigfdrFIn_wSub.pdf','-append',gcf);
close all;
end
%%
figure('color',[1 1 1],'position',POS);
for cc=sigCompFI
    figure('color',[1 1 1],'position',POS);
    subplot(3,1,1);plot(fi,Yful_SMRI(cc,:)','bo');lsline;title(['SMRI: Comp ' my_pad(cc,2)]);
    subplot(3,1,2);plot(fi,Yful_ALFF(cc,:)','bo');lsline;title(['ALFF: Comp ' my_pad(cc,2)])
    subplot(3,1,3);plot(fi,Yful_DMRI(cc,:)','bo');lsline;title(['DMRI: Comp ' my_pad(cc,2)])
    xlabel('Fluid Intelligence Score');ylabel('source intensity a.u.')
    set(gca,'FontName','Arial')
    export_fig('UKBicasig_full_GICAinit_C30_preregsite_bysigfdrSN.pdf','-append',gcf);
    close all;
end

%%
axs0 = find(strcmp('sex_f31_0_0_X_age_when_attended_assessment_centre_f21003_2_0',MULT0s_ful.final_terms));
axs1 = find(strcmp('sex_f31_0_0_X_age_when_attended_assessment_centre_f21003_2_0',MULT1s_ful.final_terms));
axs2 = find(strcmp('sex_f31_0_0_X_age_when_attended_assessment_centre_f21003_2_0',MULT2s_ful.final_terms));
[myFDR(UNI0s_ful.p{axs0}) <= 0.05;myFDR(UNI1s_ful.p{axs1}) <= 0.05;myFDR(UNI2s_ful.p{axs2}) <= 0.05]
sigCompAxS = find(sum([myFDR(UNI0s_ful.p{axs0}) <= 0.05/20;myFDR(UNI1s_ful.p{axs1}) <= 0.05/20;myFDR(UNI2s_ful.p{axs2}) <= 0.05/20])==3)

%%
sn1=find(strcmp('rSpatNorm' ,MULT0s_ful.final_terms));
sn2=find(strcmp('rSpatNorm' ,MULT1s_ful.final_terms));
sn3 = find(strcmp('rSpatNorm' ,MULT2s_ful.final_terms));
[myFDR(UNI0s_ful.p{sn1}) <= 0.05;myFDR(UNI1s_ful.p{sn2}) <= 0.05;myFDR(UNI2s_ful.p{sn3}) <= 0.05]
sigCompSN = find(sum([myFDR(UNI0s_ful.p{sn1}) <= 0.05;myFDR(UNI1s_ful.p{sn2}) <= 0.05;myFDR(UNI2s_ful.p{sn3}) <= 0.05])==3);

snSCol =  find(strcmp('rSpatNorm' ,MODELUKB0s_ful.names));
snS = MODELUKB0s_ful.X(:,snSCol);
snACol =  find(strcmp('rSpatNorm' ,MODELUKB1s_ful.names));
snA = MODELUKB1s_ful.X(:,snACol);
snDCol = find(strcmp('rSpatNorm' ,MODELUKB2s_ful.names));
snD = MODELUKB2s_ful.X(:,snDCol);


%%


%cmSN = vals2colormap(fi,'hot',[min(fi) max(fi)]);
for cc=sigCompSN
    figure('color',[1 1 1],'position',POS);
    subplot(3,1,1);plot(snS,Yful_SMRI(cc,:)','bo');lsline;title(['SMRI: Comp ' my_pad(cc,2)]);
    subplot(3,1,2);plot(snA,Yful_ALFF(cc,:)','bo');lsline;title(['ALFF: Comp ' my_pad(cc,2)])
    subplot(3,1,3);plot(snD,Yful_DMRI(cc,:)','bo');lsline;title(['DMRI: Comp ' my_pad(cc,2)])
    xlabel('rSN');ylabel('a.u.')
    set(gca,'FontName','Arial')
    export_fig('UKBicasig_full_GICAinit_C30_preregsite_bysigfdrSN.pdf','-append',gcf);
    close all;
end

%%


sn1=find(strcmp('rSpatNorm' ,MULT0_ful.final_terms));
sn2=find(strcmp('rSpatNorm' ,MULT1_ful.final_terms));
sn3 = find(strcmp('rSpatNorm' ,MULT2_ful.final_terms));
[myFDR(UNI0_ful.p{sn1}) <= 0.05;myFDR(UNI1_ful.p{sn2}) <= 0.05;myFDR(UNI2_ful.p{sn3}) <= 0.05]
sigCompSN = find(sum([myFDR(UNI0_ful.p{sn1}) <= 0.05;myFDR(UNI1_ful.p{sn2}) <= 0.05;myFDR(UNI2_ful.p{sn3}) <= 0.05])==3);

rSN0 = MODELUKB0_ful.X(:,55);
%cmAge = vals2colormap(Age,'jet',[min(Age) max(Age)]);
rSN1 = MODELUKB1_ful.X(:,55);
rSN2 = MODELUKB2_ful.X(:,55);
mID = find(MODELUKB1_ful.X(:,54) == 0);
fID = find(MODELUKB1_ful.X(:,54) == 1);
for cc = sigCompSN
    figure('color',[1 1 1]);
    subplot(3,1,1);
    plot(rSN0(mID),Yful_SMRI(cc,mID)','b*',rSN0(fID),Yful_SMRI(cc,fID)','r*');lsline;title(['SMRI: Comp' my_pad( cc,2)]);
    subplot(3,1,2);
    plot(rSN1(mID),Yful_ALFF(cc,mID)','b*',rSN1(fID),Yful_ALFF(cc,fID)','r*');lsline;title(['ALFF: Comp' my_pad( cc,2)]);
    subplot(3,1,3);
    plot(rSN2(mID),Yful_DMRI(cc,mID)','b*',rSN2(fID),Yful_DMRI(cc,fID)','r*');lsline;title(['DMRI: Comp' my_pad( cc,2)]);
    xlabel('r(SpatialNormToMeanMap)');ylabel('a.u.')
    export_fig('UKBicasig_full_GICAinit_C30_bysigfdrSNbySex.pdf','-append',gcf);
    close(gcf)
end
%%

for cc= all3ageSigInd % [1 2 5 6 8 9 13 19]
    figure('color',[1 1 1],'position',POSS);
    subplot(2,3,1:3)
    scatter3(Yful_SMRI(cc,:)'*signComps(cc),Yful_ALFF(cc,:)'*signComps(cc),Yful_DMRI(cc,:)'*signComps(cc),[],cmAge);
    %scatter3(Yful_SMRI(cc,x0id)',Yful_ALFF(cc,x0id)',Yful_DMRI(cc,x0id)',[],cmAge(x0id,:));
    %hold on;
    %scatter3(Yful_SMRI(cc,sex==1)',Yful_ALFF(cc,sex==1)',Yful_DMRI(cc,sex==1)','x',cmAge(x1id,:));
    view([30 10]);
    xlabel('SMRI');ylabel('ALFF');zlabel('DMRI');
    title(['Comp ' my_pad(cc,2)]);
    
    subplot(2,3,4);ppp=plot(Age(x0id),Yful_SMRI(cc,x0id)'*signComps(cc),'bx',Age(x1id),Yful_SMRI(cc,x1id)'*signComps(cc),'ro');legend({'M','F'},'AutoUpdate','off');lsline;title('SMRI');xlabel('Age');ylabel('source [a.u.]');
    subplot(2,3,5);plot(Age(x0id),Yful_ALFF(cc,x0id)'*signComps(cc),'bx',Age(x1id),Yful_ALFF(cc,x1id)'*signComps(cc),'ro');LG1=lsline;title('ALFF');xlabel('Age');ylabel('source [a.u.]');
    subplot(2,3,6);plot(Age(x0id),Yful_DMRI(cc,x0id)'*signComps(cc),'bx',Age(x1id),Yful_DMRI(cc,x1id)'*signComps(cc),'ro');LG2=lsline;title('DMRI');xlabel('Age');ylabel('source [a.u.]');
    
    %set(cb,'FontName','Arial')
    set(gca,'FontName','Arial')
    export_fig(['UKBicasig_full_GICAinit_C30_bysigfdrAge_PlotbySx_COMP' my_pad(cc,2) '.pdf'],gcf);
    close(gcf);
end

%% Compare model order 20 and 30
ful20 = load('C020/UKBMRI_MMIVAfull_w0GICA_C020_results_goodSub3497_new_P000.mat');
A_r1=corr(ful.icasig{1}',ful20.icasig{1}');
A_r2=corr(ful.icasig{2}',ful20.icasig{2}');
A_r3=corr(ful.icasig{3}',ful20.icasig{3}');
[mxv0,mxi0] = max(abs(A_r1),[],2);
[mxv1,mxi1] = max(abs(A_r2),[],2);
[mxv2,mxi2] = max(abs(A_r3),[],2);
disp([max(abs(A_r1),[],2) max(abs(A_r2),[],2) max(abs(A_r3),[],2)]')
disp([(1:30)' mxi0 mxi1 mxi2]')

A_r11=corr(ful.A_ls1{1},ful20.A_ls1{1});
A_r12=corr(ful.A_ls1{2},ful20.A_ls1{2});
A_r13=corr(ful.A_ls1{3},ful20.A_ls1{3});

[mxv00,mxi00] = max(abs(A_r11),[],2);
[mxv11,mxi11] = max(abs(A_r12),[],2);
[mxv22,mxi22] = max(abs(A_r13),[],2);
disp([mxv00 mxv11 mxv22]')
disp([(1:30)' mxi00 mxi11 mxi22]')

%% Compare HCSZ and UKB

ful30HCSZ = load('/media/eswar/New Volume/Documents/MISA/HCSZ/MMIVA2step_wGICApre/C030_regSite/HCSZMRI_MMIVAfull_w0GICA_C020_results_new_P000_preregSITE_signCorrected.mat');
%AA_r1=corr(ful.icasig{1}',ful30HCSZ.fulreg_HCSZ_GICAinitSgnCor.icasig{1}');
%AA_r2=corr(ful.icasig{2}',ful30HCSZ.fulreg_HCSZ_GICAinitSgnCor.icasig{2}');
%AA_r3=corr(ful.icasig{3}',ful30HCSZ.fulreg_HCSZ_GICAinitSgnCor.icasig{3}');
%[mxvA0,mxiA0] = max(abs(AA_r1),[],2);
%[mxvA1,mxiA1] = max(abs(AA_r2),[],2);
%[mxvA2,mxiA2] = max(abs(AA_r3),[],2);
%disp([max(abs(AA_r1),[],2) max(abs(AA_r2),[],2) max(abs(AA_r3),[],2)]')
%disp([(1:30)' mxiA0 mxiA1 mxiA2]')

AA_r11=corr(ful.A_ls1{1},ful30HCSZ.fulreg_HCSZ_GICAinitSgnCor.A_ls1{1});
AA_r12=corr(ful.A_ls1{2},ful30HCSZ.fulreg_HCSZ_GICAinitSgnCor.A_ls1{2});
AA_r13=corr(ful.A_ls1{3},ful30HCSZ.fulreg_HCSZ_GICAinitSgnCor.A_ls1{3});

[mxvA00,mxiA00] = max(abs(AA_r11),[],2);
[mxvA11,mxiA11] = max(abs(AA_r12),[],2);
[mxvA22,mxiA22] = max(abs(AA_r13),[],2);
disp([(1:30)' mxvA00 mxvA11 mxvA22]')
disp([(1:30)' mxiA00 mxiA11 mxiA22]')
%%
figure('color',[1 1 1]);
subplot(1,3,2);imagesc(corr(Yful_DMRI'),[-1 1]);axis square;title('DMRI')
subplot(1,3,2);imagesc(corr(Yful_DMRI'),[-1 1]);axis square;title('DMRI')
subplot(1,3,3);imagesc(corr(Yful_ALFF'),[-1 1]);axis square;title('ALFF')

figure('color',[1 1 1]);
subplot(1,3,1);imagesc(corr(Yful_SMRI',Yful_ALFF'),[-1 1]);axis square;title('SMRI-ALFF')
subplot(1,3,2);imagesc(corr(Yful_SMRI',Yful_DMRI'),[-1 1]);axis square;title('SMRI-DMRI')
subplot(1,3,3);imagesc(corr(Yful_ALFF',Yful_DMRI'),[-1 1]);axis square;title('ALFF-DMRI')

figure('color',[1 1 1]);
plot([abs(diag(corr(Yful_SMRI',Yful_ALFF')))';abs(diag(corr(Yful_SMRI',Yful_DMRI')))';abs(diag(corr(Yful_ALFF',Yful_DMRI')))']')
%%
figure;plot(Yful_SMRI([8],:)',Yful_SMRI([19],:)','*')
figure;histogram(Yful_SMRI(8,:),40)

%%
mask = '/home/eswar/Documents/projects/MISA/Brain_Masks/gmMask_TPM_thrp2_fractc8_3mmbin.nii';
maskd = '/home/eswar/Documents/projects/MISA/Brain_Masks/wmMask_TPM_thrp4_fractc8_3mmbin.nii';


smask = spm_read_vols(spm_vol(mask));
dmask = spm_read_vols(spm_vol(maskd));

comp12_SMRI  = spm_read_vols(spm_vol('comp12_SMRI_Clust_thr0.04_N400.nii'));
comp12_DMRI  = spm_read_vols(spm_vol('comp12_DMRI_Clust_thr0.02_N40.nii'));
comp12_ALFF  = spm_read_vols(spm_vol('comp12_ALFF_Clust_thr0.02_N40_mask.nii'));

comp17_SMRI  = spm_read_vols(spm_vol('comp17_SMRI_Clust_thr0.04_N400_mask.nii'));
comp17_DMRI  = spm_read_vols(spm_vol('comp17_DMRI_Clust_thr0.02_N40_mask.nii'));
comp17_ALFF  = spm_read_vols(spm_vol('comp17_ALFF_Clust_thr0.04_N1000_mask.nii'));
%%

smriData = load(fullfile(fileparts(origddir),'data','SMRIdata_gmMaskNEW.mat'));
dmriData = load(fullfile(fileparts(origddir),'data','DMRIdata_Sm6_wmMaskNew.mat'));
alffData = load(fullfile(fileparts(origddir),'data','ALFFdata_gmMaskNEW.mat'));

SMRIdata = smriData.SMRIdata(:,id2907.lt5NanID);
DMRIdata = dmriData.DMRIdata(:,id2907.lt5NanID);
ALFFdata = alffData.ALFFdata(:,id2907.lt5NanID);
%%

c12_smriMask = comp12_SMRI(smask>0);
c12_alffMask = comp12_ALFF(smask>0);
c12_dmriMask = comp12_DMRI(dmask>0);

unique(c12_smriMask)
c12_y00_l05 = prctile(Yful_SMRI(12,:)',5);
c12_y00_h95 = prctile(Yful_SMRI(12,:)',95);

c12_y11_l05 = prctile(Yful_ALFF(12,:)',5);
c12_y11_h95 = prctile(Yful_ALFF(12,:)',95);

c12_y22_l05 = prctile(Yful_DMRI(12,:)',5);
c12_y22_h95 = prctile(Yful_DMRI(12,:)',95);

y00_l05id = find(Yful_SMRI(12,:) <= c12_y00_l05);
y00_h95id = find(Yful_SMRI(12,:) >= c12_y00_h95);

y11_l05id = find(Yful_ALFF(12,:) <= c12_y11_l05);
y11_h95id = find(Yful_ALFF(12,:) >= c12_y11_h95);

y22_l05id = find(Yful_DMRI(12,:) <= c12_y22_l05);
y22_h95id = find(Yful_DMRI(12,:) >= c12_y22_h95);
%%
uq12_00 = unique(c12_smriMask);
for ii=2:length(uq12_00)
    c12_mnLow_SMRI(ii-1,:) = mean(SMRIdata(c12_smriMask==uq12_00(ii),y00_l05id));
    c12_mnHigh_SMRI(ii-1,:) = mean(SMRIdata(c12_smriMask==uq12_00(ii),y00_h95id));
end

uq12_11 = unique(c12_alffMask);
for ii=2:length(uq12_11)
    c12_mnLow_ALFF(ii-1,:) = mean(ALFFdata(c12_alffMask==uq12_11(ii),y11_l05id));
    c12_mnHigh_ALFF(ii-1,:) = mean(ALFFdata(c12_alffMask==uq12_11(ii),y11_h95id));
end

uq12_22 = unique(c12_dmriMask);
for ii=2:length(uq12_22)
    c12_mnLow_DMRI(ii-1,:) = mean(DMRIdata(c12_dmriMask==uq12_22(ii),y22_l05id));
    c12_mnHigh_DMRI(ii-1,:) = mean(DMRIdata(c12_dmriMask==uq12_22(ii),y22_h95id));
end
%%


c12_mnRAW3dLow_SMRI = mean(SMRIdata(:,y00_l05id),2);
c12_mnRAW3dHigh_SMRI = mean(SMRIdata(:,y00_h95id),2);

c12_mnRAW3dLow_ALFF = mean(ALFFdata(:,y11_l05id),2);
c12_mnRAW3dHigh_ALFF = mean(ALFFdata(:,y11_h95id),2);

c12_mnRAW3dLow_DMRI = mean(DMRIdata(:,y22_l05id),2);
c12_mnRAW3dHigh_DMRI = mean(DMRIdata(:,y22_h95id),2);

hdr = spm_vol(mask);

tmp00 = zeros([numel(smask) 3]);
tmp00(smask>0,1) = c12_mnRAW3dHigh_SMRI;
tmp00(smask>0,2) = c12_mnRAW3dLow_SMRI;
tmp00(smask>0,3) = c12_mnRAW3dHigh_SMRI - c12_mnRAW3dLow_SMRI;
tmp00 = reshape(tmp00,[size(smask) 3]);
create_4DNiftifile('c12_SMRI_qH_qL_qdiff.nii',tmp00,hdr.mat);


tmp00 = zeros([numel(smask) 3]);
tmp00(smask>0,1) = c12_mnRAW3dHigh_ALFF;
tmp00(smask>0,2) = c12_mnRAW3dLow_ALFF;
tmp00(smask>0,3) = c12_mnRAW3dHigh_ALFF - c12_mnRAW3dLow_ALFF;
tmp00 = reshape(tmp00,[size(smask) 3]);
create_4DNiftifile('c12_ALFF_qH_qL_qdiff.nii',tmp00,hdr.mat);

tmp00 = zeros([numel(dmask) 3]);
tmp00(dmask>0,1) = c12_mnRAW3dHigh_DMRI;
tmp00(dmask>0,2) = c12_mnRAW3dLow_DMRI;
tmp00(dmask>0,3) = c12_mnRAW3dHigh_DMRI - c12_mnRAW3dLow_DMRI;
tmp00 = reshape(tmp00,[size(dmask) 3]);
create_4DNiftifile('c12_DMRI_qH_qL_qdiff.nii',tmp00,hdr.mat);


%%
ttaCol=find(strcmp('time_to_answer_f4288_2_0' ,DEMO_2907.demo_col_names));
tta = DEMO_2907.demo_out(:,ttaCol);
rat0 = ([sum(tta(y00_h95id) < 1000)/length(y00_h95id) sum(tta(y00_h95id) > 1000)/length(y00_h95id) sum(tta(y00_l05id) < 1000)/length(y00_l05id) sum(tta(y00_l05id) > 1000)/length(y00_l05id)]);
POSSS = [116         149        1625         776];
figure('color',[1 1 1],'position',POSSS);
for uu = 2:length(uq12_00)
    subplot(3,5,uu-1);
    plot(tta(y00_h95id),c12_mnHigh_SMRI(uu-1,:),'ro',tta(y00_l05id),c12_mnLow_SMRI(uu-1,:)','bo')
    if uu == 2
        legend({'high','low'});
        xlabel('tta')
        ylabel('mean GM density')
    end
    title(['Comp12 SMRI: cluster ' num2str(uu-1)])
end
%
rat1 = ([sum(tta(y11_h95id) < 1000)/length(y11_h95id) sum(tta(y11_h95id) > 1000)/length(y11_h95id) sum(tta(y11_l05id) < 1000)/length(y11_l05id) sum(tta(y11_l05id) > 1000)/length(y11_l05id)]);

% figure('color',[1 1 1])
for uu=2:length(uq12_11)
    subplot(3,5,5+uu-1);
    plot(tta(y11_h95id),c12_mnHigh_SMRI(uu-1,:),'ro',tta(y11_l05id),c12_mnLow_SMRI(uu-1,:)','bo');
    if uu == 2
        legend({'high','low'});
        xlabel('tta')
        ylabel('mean ALFF')
    end
    title(['Comp12 ALFF: cluster ' num2str(uu-1)])

end
%  
rat2 = ([sum(tta(y22_h95id) < 1000)/length(y22_h95id) sum(tta(y22_h95id) > 1000)/length(y22_h95id) sum(tta(y22_l05id) < 1000)/length(y22_l05id) sum(tta(y22_l05id) > 1000)/length(y22_l05id)]); 

%figure('color',[1 1 1])
for uu=2:length(uq12_22)
    subplot(3,5,10+uu-1);
    plot(tta(y22_h95id),c12_mnHigh_SMRI(uu-1,:),'ro',tta(y22_l05id),c12_mnLow_SMRI(uu-1,:)','bo')
    if uu == 2
        legend({'high','low'});
        xlabel('tta')
        ylabel('mean FA')
    end
    title(['Comp12 FA: cluster ' num2str(uu-1)])

end
export_fig('comp12_tta_rawActivation_icasig_quartile5_vs_quartile95.png');
%%

c17_smriMask = comp17_SMRI(smask>0);
c17_alffMask = comp17_ALFF(smask>0);
c17_dmriMask = comp17_DMRI(dmask>0);

unique(c17_smriMask)
c17_y00_l05 = prctile(Yful_SMRI(17,:)',5);
c17_y00_h95 = prctile(Yful_SMRI(17,:)',95);

c17_y11_l05 = prctile(Yful_ALFF(17,:)',5);
c17_y11_h95 = prctile(Yful_ALFF(17,:)',95);

c17_y22_l05 = prctile(Yful_DMRI(17,:)',5);
c17_y22_h95 = prctile(Yful_DMRI(17,:)',95);

y00_l05id_c17 = find(Yful_SMRI(17,:) <= c17_y00_l05);
y00_h95id_c17 = find(Yful_SMRI(17,:) >= c17_y00_h95);

y11_l05id_c17 = find(Yful_ALFF(17,:) <= c17_y11_l05);
y11_h95id_c17 = find(Yful_ALFF(17,:) >= c17_y11_h95);

y22_l05id_c17 = find(Yful_DMRI(17,:) <= c17_y22_l05);
y22_h95id_c17 = find(Yful_DMRI(17,:) >= c17_y22_h95);
%%
uq17_00 = unique(c17_smriMask);
for ii=2:length(uq17_00)
    c17_mnLow_SMRI(ii-1,:) = mean(SMRIdata(c17_smriMask==uq17_00(ii),y00_l05id_c17));
    c17_mnHigh_SMRI(ii-1,:) = mean(SMRIdata(c17_smriMask==uq17_00(ii),y00_h95id_c17));
end

uq17_11 = unique(c17_alffMask);
for ii=2:length(uq17_11)
    c17_mnLow_ALFF(ii-1,:) = mean(ALFFdata(c17_alffMask==uq17_11(ii),y11_l05id_c17));
    c17_mnHigh_ALFF(ii-1,:) = mean(ALFFdata(c17_alffMask==uq17_11(ii),y11_h95id_c17));
end

uq17_22 = unique(c17_dmriMask);
for ii=2:length(uq17_22)
    c17_mnLow_DMRI(ii-1,:) = mean(DMRIdata(c17_dmriMask==uq17_22(ii),y22_l05id_c17));
    c17_mnHigh_DMRI(ii-1,:) = mean(DMRIdata(c17_dmriMask==uq17_22(ii),y22_h95id_c17));
end
%%
c17_mnRAW3dLow_SMRI = mean(SMRIdata(:,y00_l05id_c17),2);
c17_mnRAW3dHigh_SMRI = mean(SMRIdata(:,y00_h95id_c17),2);

c17_mnRAW3dLow_ALFF = mean(ALFFdata(:,y11_l05id_c17),2);
c17_mnRAW3dHigh_ALFF = mean(ALFFdata(:,y11_h95id_c17),2);

c17_mnRAW3dLow_DMRI = mean(DMRIdata(:,y22_l05id_c17),2);
c17_mnRAW3dHigh_DMRI = mean(DMRIdata(:,y22_h95id_c17),2);

tmp00 = zeros([numel(smask) 3]);
tmp00(smask>0,1) = c17_mnRAW3dHigh_SMRI;
tmp00(smask>0,2) = c17_mnRAW3dLow_SMRI;
tmp00(smask>0,3) = c17_mnRAW3dHigh_SMRI - c17_mnRAW3dLow_SMRI;
tmp00 = reshape(tmp00,[size(smask) 3]);
create_4DNiftifile('c17_SMRI_qH_qL_qdiff.nii',tmp00,hdr.mat);


tmp00 = zeros([numel(smask) 3]);
tmp00(smask>0,1) = c17_mnRAW3dHigh_ALFF;
tmp00(smask>0,2) = c17_mnRAW3dLow_ALFF;
tmp00(smask>0,3) = c17_mnRAW3dHigh_ALFF - c17_mnRAW3dLow_ALFF;
tmp00 = reshape(tmp00,[size(smask) 3]);
create_4DNiftifile('c17_ALFF_qH_qL_qdiff.nii',tmp00,hdr.mat);

tmp00 = zeros([numel(dmask) 3]);
tmp00(dmask>0,1) = c17_mnRAW3dHigh_DMRI;
tmp00(dmask>0,2) = c17_mnRAW3dLow_DMRI;
tmp00(dmask>0,3) = c17_mnRAW3dHigh_DMRI - c17_mnRAW3dLow_DMRI;
tmp00 = reshape(tmp00,[size(dmask) 3]);
create_4DNiftifile('c17_DMRI_qH_qL_qdiff.nii',tmp00,hdr.mat);

%%

rat0_c17 = ([sum(tta(y00_h95id_c17) < 1000)/length(y00_h95id_c17) sum(tta(y00_h95id_c17) > 1000)/length(y00_h95id_c17) sum(tta(y00_l05id_c17) < 1000)/length(y00_l05id_c17) sum(tta(y00_l05id_c17) > 1000)/length(y00_l05id_c17)]);

figure('color',[1 1 1],'position',POSSS);
for uu = 2:length(uq17_00)
    subplot(3,5,uu-1);
    plot(tta(y00_h95id_c17),c17_mnHigh_SMRI(uu-1,:),'ro',tta(y00_l05id_c17),c17_mnLow_SMRI(uu-1,:)','bo')
    if uu == 2
        legend({'high','low'});
        xlabel('tta')
        ylabel('mean GM density')
    end
    title(['Comp17 SMRI: cluster ' num2str(uu-1)])
end
%
rat1_c17 = ([sum(tta(y11_h95id_c17) < 1000)/length(y11_h95id_c17) sum(tta(y11_h95id_c17) > 1000)/length(y11_h95id_c17) sum(tta(y11_l05id_c17) < 1000)/length(y11_l05id_c17) sum(tta(y11_l05id_c17) > 1000)/length(y11_l05id_c17)]);

% figure('color',[1 1 1])
for uu=2:length(uq17_11)
    subplot(3,5,5+uu-1);
    plot(tta(y11_h95id_c17),c17_mnHigh_SMRI(uu-1,:),'ro',tta(y11_l05id_c17),c17_mnLow_SMRI(uu-1,:)','bo');
    if uu == 2
        legend({'high','low'});
        xlabel('tta')
        ylabel('mean ALFF')
    end
    title(['Comp17 ALFF: cluster ' num2str(uu-1)])

end
%  
rat2_c17 = ([sum(tta(y22_h95id_c17) < 1000)/length(y22_h95id_c17) sum(tta(y22_h95id_c17) > 1000)/length(y22_h95id_c17) sum(tta(y22_l05id_c17) < 1000)/length(y22_l05id_c17) sum(tta(y22_l05id_c17) > 1000)/length(y22_l05id_c17)]); 

%figure('color',[1 1 1])
for uu=2:length(uq17_22)
    subplot(3,5,10+uu-1);
    plot(tta(y22_h95id_c17),c17_mnHigh_SMRI(uu-1,:),'ro',tta(y22_l05id_c17),c17_mnLow_SMRI(uu-1,:)','bo')
    if uu == 2
        legend({'high','low'});
        xlabel('tta')
        ylabel('mean FA')
    end
    title(['Comp17 FA: cluster ' num2str(uu-1)])

end
export_fig('comp17_tta_rawActivation_icasig_quartile5_vs_quartile95.png');
%%  with site

sigTermsN = intersect(MULT0s_ful.final_terms,intersect(MULT1s_ful.final_terms,MULT2s_ful.final_terms))';

%%
% SMRI
[sx00_comp,sx00_T,sx00_contr]= printouts(UNI0s_ful,find(strcmp(UNI0s_ful.tests,sigTermsN{4})>0));
[tta00_comp,tta00_T,tta00_contr]= printouts(UNI0s_ful,find(strcmp(UNI0s_ful.tests,sigTermsN{7})>0));
[a00_comp,a00_T,a00_contr] = printouts(UNI0s_ful,find(strcmp(UNI0s_ful.tests,sigTermsN{1})>0));
[fi00_comp,fi00_T,fi00_contr]= printouts(UNI0s_ful,find(strcmp(UNI0s_ful.tests,sigTermsN{2})>0));
[rsn00_comp,rsn00_T,rsn00_contr]=printouts(UNI0s_ful,find(strcmp(UNI0s_ful.tests,sigTermsN{3})>0));
[sxa01_comp,sxa01_T,sxa01_contr]=printouts(UNI0s_ful,find(strcmp(UNI0s_ful.tests,sigTermsN{5})>0));
[st00_comp,st00_T,st00_contr]= printouts(UNI0s_ful,find(strcmp(UNI0s_ful.tests,sigTermsN{6})>0));

% ALFF
[sx11_comp,sx11_T,sx11_contr]= printouts(UNI1s_ful,find(strcmp(UNI1s_ful.tests,sigTermsN{4})>0));
[tta11_comp,tta11_T,tta11_contr]= printouts(UNI1s_ful,find(strcmp(UNI1s_ful.tests,sigTermsN{7})>0));
[a11_comp,a11_T,a11_contr] = printouts(UNI1s_ful,find(strcmp(UNI1s_ful.tests,sigTermsN{1})>0));
[fi11_comp,fi11_T,fi11_contr] = printouts(UNI1s_ful,find(strcmp(UNI1s_ful.tests,sigTermsN{2})>0));
[rsn11_comp,rsn11_T,rsn11_contr]=printouts(UNI1s_ful,find(strcmp(UNI1s_ful.tests,sigTermsN{3})>0));
[sxa12_comp,sxa12_T,sxa12_contr]=printouts(UNI1s_ful,find(strcmp(UNI1s_ful.tests,sigTermsN{5})>0));
[st11_comp,st11_T,st11_contr]= printouts(UNI1s_ful,find(strcmp(UNI1s_ful.tests,sigTermsN{6})>0));

%
% DMRI
[sx22_comp,sx22_T,sx22_contr]= printouts(UNI2s_ful,find(strcmp(UNI2s_ful.tests,sigTermsN{4})>0));
[tta22_comp,tta22_T,tta22_contr]= printouts(UNI2s_ful,find(strcmp(UNI2s_ful.tests,sigTermsN{7})>0));
[a22_comp,a22_T,a22_contr] = printouts(UNI2s_ful,find(strcmp(UNI2s_ful.tests,sigTermsN{1})>0));
[fi22_comp,fi22_T,fi22_contr] = printouts(UNI2s_ful,find(strcmp(UNI2s_ful.tests,sigTermsN{2})>0));
[rsn22_comp,rsn22_T,rsn22_contr]=printouts(UNI2s_ful,find(strcmp(UNI2s_ful.tests,sigTermsN{3})>0));
[sxa23_comp,sxa23_T,sxa23_contr]=printouts(UNI2s_ful,find(strcmp(UNI2s_ful.tests,sigTermsN{5})>0));
[st22_comp,st22_T,st22_contr]= printouts(UNI2s_ful,find(strcmp(UNI2s_ful.tests,sigTermsN{6})>0));

%%
odir = '/media/eswar/New Volume/Documents/MISA/UKB/MMIVA2step_wGICApre/C030_UKB/figs/stats/wSite';

aa{1} = UNI0s_ful;
aa{2} = UNI1s_ful;
aa{3} = UNI2s_ful;

dm{1} = DEMO0s_ful;
dm{2} = DEMO1s_ful;
dm{3} = DEMO2s_ful;

%%
plotuni(aa,dm,1,sigTermsN{1},a00_contr,a00_comp,a00_T,signComps,odir)
plotuni(aa,dm,2,sigTermsN{1},a11_contr,a11_comp,a11_T,signComps,odir)
plotuni(aa,dm,3,sigTermsN{1},a22_contr,a22_comp,a22_T,signComps,odir)

plotuni(aa,dm,1,sigTermsN{2},fi00_contr,fi00_comp,fi00_T,signComps,odir)
plotuni(aa,dm,2,sigTermsN{2},fi11_contr,fi11_comp,fi11_T,signComps,odir)
plotuni(aa,dm,3,sigTermsN{2},fi22_contr,fi22_comp,fi22_T,signComps,odir)

plotuni(aa,dm,1,sigTermsN{4},sx00_contr,sx00_comp,sx00_T,signComps,odir)
plotuni(aa,dm,2,sigTermsN{4},sx11_contr,sx11_comp,sx11_T,signComps,odir)
plotuni(aa,dm,3,sigTermsN{4},sx22_contr,sx22_comp,sx22_T,signComps,odir)

plotuni(aa,dm,1,sigTermsN{6},st00_contr,st00_comp,st00_T,signComps,odir)
plotuni(aa,dm,2,sigTermsN{6},st11_contr,st11_comp,st11_T,signComps,odir)
plotuni(aa,dm,3,sigTermsN{6},st22_contr,st22_comp,st22_T,signComps,odir)

plotuni(aa,dm,1,sigTermsN{7},tta00_contr,tta00_comp,tta00_T,signComps,odir)
plotuni(aa,dm,2,sigTermsN{7},tta11_contr,tta11_comp,tta11_T,signComps,odir)
plotuni(aa,dm,3,sigTermsN{7},tta22_contr,tta22_comp,tta22_T,signComps,odir)

