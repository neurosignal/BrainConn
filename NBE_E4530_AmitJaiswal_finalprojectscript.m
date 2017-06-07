%% Script for source level connectivity analysis under course NBEHBC (E4530)
% Author: Amit K. Jaiswal (PhD Student NBE group)
% Email: amit.jaiswal@aalto.fi
% Date: 23/04/2017

% Details:
% The script is written for source level functional connectivity analysis.
% To get the source level data, LCMV beamformer source estimation method 
% has been used.


%% Section 1: Source localization part (Using Elekta analysis software)
% Step1: MaxFilter the data (in Elekta MaxFilter 2.0)
% Step2: MRI/MEG Coregistration (in Elekta Mrilab)
% Step3: MRI Segmentation (Elekta Seglab)
% Step4: BEM modeling (Elekta Xfit)
% Step5: Source localization using LCMV (Elekta Beamformer 2.0)
% Step6: Virtual electrode signal calculation (Elekta Beamformer 2.0)
% The virtual electrode signal file was then saved for further network analysis 


%% Section 2: Source level functional connectivity calculation and visualization 
% Set path
restoredefaultpath
clc; clear; close all;
addpath('/net/bonsai/home/amit/Documents/MATLAB/fieldtrip-master')
ft_defaults
cd /net/bonsai/home/amit/Documents/MATLAB/

for finger=2:5 % All the 4 dataset (D2-D5) were calculated for 4 different datasets but only D5 is presented in the report
    megfile=['/neuro/data/beamformer_data/Beamformer_share/SEF_laser/MH/MH_raw/' 'MH_D' num2str(finger) '_raw.fif'];
    VEfile= ['/net/bonsai/home/amit/Documents/0Aalto_PhD_work/BrainConn/D' num2str(finger) '_raw_sss_eventbf_bpno_0.01-0.08_VE.fif'];
    mrifile='/neuro/data/beamformer_data/Beamformer_share/SEF_laser/MH/mri_coreg/sets/hall_michael_1_01-amit-170424.fif';
    segmrifile='/neuro/data/beamformer_data/Beamformer_share/SEF_laser/MH/mri_coreg/sets/hall_michael_1_01-amit-170424-seg.mat';

    stimchan='STI101';
    megchan= 'VIRT*';
    stimval=1;
    trialwin=[-0.5 0.5];
    bpfreq=[2 95];
    mri_seg='no';
    
% Browse continuous raw data          
    cfg = [];
    cfg.dataset  = VEfile;
    data_c = ft_preprocessing(cfg);
    cfg.viewmode = 'vertical';
    cfg.channel  = megchan;
    cfg.preproc.demean = 'yes';
    cfg.blocksize      = 5;
    cfg.ylim           = [-5e-8 5e-8];
    cfg.hpfilter       = 'yes';
    cfg.hpfreq         = 2;
    cfg.demean         ='yes';
    ft_databrowser(cfg, data_c);

% Define trials
    cfg                     = [];                   
    cfg.dataset             = VEfile;
    cfg.channel             = megchan;
    cfg.trialdef.eventtype  = stimchan;
    cfg.trialdef.eventvalue = 1;                     
    cfg.trialdef.prestim    = abs(trialwin(1));                     
    cfg.trialdef.poststim   = trialwin(2); 
    cfg.trialfun            = 'EN_phanton_trialfun';
    cfg.trig_min_gap        = 0.1; % in second
    cfg.minlength = cfg.trialdef.prestim + cfg.trialdef.poststim;
    cfg = ft_definetrial(cfg);
    % Leave the first and last trigger
    cfg.trl(1:2,:)  =[];
    cfg.trl(end-1:end,:)=[];
    
% Trial preprocessing  
    cfg.demean     = 'yes';
    cfg.bpfilter   = 'yes';
    cfg.bpfiltord  = 2;
    cfg.bpfilttype = 'but';
    cfg.bpfreq     = bpfreq;
    cfg.bsfilter   = 'yes';
    cfg.bsfreq     = [49 51];
    data = ft_preprocessing(cfg);
    
% Interactive data browser  
    cfg.blocksize  = cfg.trialdef.prestim + cfg.trialdef.poststim;
    cfg.continuous = 'no';
    cfg.viewmode   = 'vertical'; 
    ft_databrowser(cfg, data);
    
% Connectivity (correlation based)
for tr=1:size(data.trial,2)
    for chr=1:length(data.label)
        for chc=1:length(data.label)
            CorMat(chr,chc,tr)=corr(data.trial{1,tr}(chr,:)',data.trial{1,tr}(chc,:)');
        end
    end
end
Cor=squeeze(mean(CorMat,3));

% Corelation plot
    addpath('~/Documents/MATLAB/fieldtrip-master/external/bct/') % Add path to BCT toolbox
    p=0.8;
    thr_Cor = threshold_proportional(Cor, p);
    thr_Cor=triu(thr_Cor);
    figure('Color',[1 1 1]);
    imagesc(thr_Cor); axis xy
    caxis([0 1])
    tick=1:1:length(data.label);
    xticklabel=char(data.label);
    set(gca,'XTick',tick,'XTickLabel',xticklabel)
    set(gca,'YTick',tick,'YTickLabel',xticklabel)
    title('Source level connectivity', 'FontSize', 20); % set title
    colorbar('FontSize', 15)
    
% Preprocess and average virtual channel signals
    cfg=[];
    cfg.preproc.hpfilter='yes';
    cfg.preproc.hpfreq=5;
    cfg.preproc.lpfilter='yes';
    cfg.preproc.lpfreq=100;
    tlkvc=ft_timelockanalysis(cfg, data);

    xx=sqrt(length(tlkvc.label));
    yy=ceil(xx); zz=fix(xx);
    figure('Color', 'w')
    for ii=1:length(tlkvc.label)
        ax(ii)=subplot_tight(yy,zz,ii,0.04); 
        plot(tlkvc.time, tlkvc.avg(ii,:), 'Parent', ax(ii)); 
        xlim([tlkvc.time(1) tlkvc.time(end)]);
        title(tlkvc.label(ii), 'FontSize', 12)
    end
    linkaxes(ax,'x');
            
% Plot connections over 3D brain volume
if isequal(mri_seg, 'yes')
    mri = ft_read_mri(mrifile);
% Apply volume realignment
    cfg=[];
    cfg.method    = 'interactive';
    cfg.coordsys  = 'neuromag';
    cfg.parameter = 'anatomy';
    cfg.viewresult= 'yes' ;
    [mri] = ft_volumerealign(cfg, mri);
% Apply segmentation
    cfg          = [];
    cfg.output   = 'brain';
    segmentedmri = ft_volumesegment(cfg, mri);
    segmentedmri.transform = mri.transform;
    segmentedmri.anatomy   = mri.anatomy;
else
    load(segmrifile)
% Plot segmented MRI
    cfg              = [];
    cfg.funparameter = 'brain';
    ft_sourceplot(cfg, segmentedmri);
end

% compute the subject's headmodel/volume conductor model
    cfg                = [];
    cfg.method         = 'singleshell';
    headmodel          = ft_prepare_headmodel(cfg, segmentedmri);
% Extract fiducial points location and label
    headshape=ft_read_headshape(megfile);
    headshape=ft_convert_units(headshape, 'm');
    headshape.pos=[];
    headshape.label=[];
% Plot connectivity over coregistered 3D brain
    figure
    ft_plot_vol(headmodel, 'edgecolor', 'w', 'facecolor','none', 'facealpha', 0.3); material dull; camlight; % plot volumetric heahmodel
    ft_plot_headshape(headshape, 'fidcolor', 'r', 'fidmarker', '*', 'fidlabel', 'yes', 'FontSize', 20)
    hold on
    nodes=importdata(['/net/bonsai/home/amit/Documents/0Aalto_PhD_work/BrainConn/D' num2str(finger) '_raw_sss_eventbf_bpno_0.01-0.08.pts'])/1000;
    plot3(nodes(:,1), nodes(:,2), nodes(:,3), 'bo','MarkerSize', 8,'linewidth', 3); axis off
    text(nodes(:,1), nodes(:,2), nodes(:,3), data.label, 'FontSize', 20)
    rotate3d
    
% connect lines
    [ch_pr_r,ch_pr_c,ch_pr_v]=find(thr_Cor);

    minVal = min(ch_pr_v);  maxVal = max(ch_pr_v);
    norm_ch_pr_v = (ch_pr_v - minVal) / ( maxVal - minVal );

    cnt=0;
for chr=1:length(data.label)
        for chc=1:length(data.label)
            if ~(thr_Cor(chr,chc)==0)
                cnt=cnt+1;
                line([nodes(chr,1),nodes(chc,1)],[nodes(chr,2),nodes(chc,2)],[nodes(chr,3),nodes(chc,3)], 'linewidth', 1, 'color',[0.5 norm_ch_pr_v(cnt) 0.5])
            else
            end
        end
end 

% Auto rotate the view    
    for ii=1:500
        view(ii, -20)
        pause(0.01)
    end
    close all
end


%% Section 3: Calculate Network metrics of above functional network

CIJ= thr_Cor; % or load the save adjacency matrix

% Adjacency matrix plot
figure('color', [1 1 1])
imagesc(CIJ)
set(gca,'XTick',tick,'XTickLabel',xticklabel, 'FontSize', 20)
set(gca,'YTick',tick,'YTickLabel',xticklabel, 'FontSize', 20)
xlabel('Node number', 'FontSize', 20)
ylabel('Node number', 'FontSize', 20)
title('Adjacency matrix plot', 'FontSize', 20)

% Joint degree distribution matrix
[J,J_od,J_id,J_bl] = jdegree(CIJ); 
figure('color', [1 1 1])
imagesc(J);
set(gca,'XTick',tick,'XTickLabel',xticklabel, 'FontSize', 20)
set(gca,'YTick',tick,'YTickLabel',xticklabel, 'FontSize', 20)
xlabel('Nodes --->', 'FontSize', 20)
ylabel('Nodes --->', 'FontSize', 20)
title('Joint degree distribution matrix', 'FontSize', 20)

% Calculate inward and  outward links (Degree)
[id,od,deg2] = degrees_dir(CIJ); 
figure('color', [1 1 1])
subplot(2,1,1); stem(id, 'Linewidth', 2.5);
set(gca,'XTick',tick,'XTickLabel',xticklabel, 'FontSize', 20)
xlabel('Node number', 'FontSize', 20)
ylabel('Number of inward links', 'FontSize', 20)
title('Inward Node degree ', 'FontSize', 20)
subplot(2,1,2); stem(od,'Linewidth', 2.5);
set(gca,'XTick',tick,'XTickLabel',xticklabel, 'FontSize', 20)
xlabel('Node number', 'FontSize', 20)
ylabel('Number of outward links', 'FontSize', 20)
title('Outward Node degree ', 'FontSize', 20)

% Clustering & Transitivity
C=clustering_coef_bd(CIJ); % Clustering coeff
T=transitivity_bd(CIJ); % Transitivity
figure('color', [1 1 1])
subplot(2,1,1); stem(C, 'Linewidth', 2.5);
set(gca,'XTick',tick,'XTickLabel',xticklabel, 'FontSize', 20)
xlabel('Node number', 'FontSize', 20);
ylabel('Clustering coeficient', 'FontSize', 20);
title('Clustering Coeff', 'FontSize', 20);
subplot(2,1,2)
stem(T, 'Linewidth', 2.5);
ylabel('Transitivity', 'FontSize', 20);
title('Transitivity', 'FontSize', 20);

% Optimal overlapping community structure
M = link_communities(CIJ,'complete');
figure('color', [1 1 1])
imagesc(M);
xlabel('Node number', 'FontSize', 20);
ylabel('Communities', 'FontSize', 20)
title('Nodal community-affiliation matrix', 'FontSize', 20);

% Modularity
[Ci,Q2]=modularity_und(CIJ);
figure('color', [1 1 1])
stem(Ci, 'Linewidth', 2.5);
set(gca,'XTick',tick,'XTickLabel',xticklabel, 'FontSize', 20)
xlabel('Node number', 'FontSize', 20);
ylabel('Communities', 'FontSize', 20)
title('Optimal community structure', 'FontSize', 20);

% Rich club coefficients
[R,Nk,Ek] = rich_club_bd(CIJ,25); 
figure('color', [1 1 1])
subplot(3,1,1)
stem(R, 'Linewidth', 1.5);
xlabel('Node number')
ylabel('R for levels 1-25')
title('Rich-club coefficient')
subplot(3,1,2)
stem(Nk, 'Linewidth', 1.5);
xlabel('Node number')
ylabel('Nk for degree>25')
title('Number of nodes with degree>25')
subplot(3,1,3)
stem(Ek, 'Linewidth', 1.5);
xlabel('Node number')
ylabel('Ek for degree>25')
title('Number of edges remaining in subgraph with degree>25')

% Assortativity coefficient (assortativity_bin)
r0 = assortativity_bin(CIJ,0);
r1 = assortativity_bin(CIJ,1);
r2 = assortativity_bin(CIJ,2);
r3 = assortativity_bin(CIJ,3);
r4 = assortativity_bin(CIJ,4);
figure('Color',[1 1 1]);
stem([r0; r1; r2; r3; r4], 'linewidth', 2)
xlim([0.5 5.5])
tick=1:5;
xticklabel={'flag=0','flag=1','flag=2','flag=3','flag=4'};
set(gca,'XTick',tick,'XTickLabel',xticklabel, 'FontSize', 20)
ylabel('Assortativity coefficient', 'FontSize', 20)
title('Assortativity coefficient (assortativity__bin)', 'FontSize', 20)

% Distance matrix for shortest distance between nodes
D=distance_bin(CIJ); 
figure('Color',[1 1 1]);
imagesc(D);
set(gca,'XTick',tick,'XTickLabel',xticklabel, 'FontSize', 20)
set(gca,'YTick',tick,'YTickLabel',xticklabel, 'FontSize', 20)
xlabel('Node number', 'FontSize', 20)
ylabel('Node number', 'FontSize', 20)
title('Shortest path length between node-pair' , 'FontSize', 20)

% Node betweenness centrality
BC=betweenness_bin(CIJ);
figure('Color',[1 1 1]);
stem(BC, 'Linewidth', 2.5);
set(gca,'XTick',tick,'XTickLabel',xticklabel, 'FontSize', 20)
xlabel('node number', 'FontSize', 20);
ylabel('Betweenness centrality', 'FontSize', 20);
title('Node betweenness centrality', 'FontSize', 20);

% Edge betweenness centrality
[EBC,BC_EDGE]=edge_betweenness_bin(CIJ);
figure('Color',[1 1 1]);
subplot(1,2,1)
stem(BC_EDGE, 'Linewidth', 2.5);
set(gca,'XTick',tick,'XTickLabel',xticklabel, 'FontSize', 20)
xlabel('node number', 'FontSize', 20);
ylabel('Edge betweenness centrality', 'FontSize', 20);
title('Edge betweenness centrality', 'FontSize', 20);
subplot(1,2,2)
imagesc(EBC);
set(gca,'XTick',tick,'XTickLabel',xticklabel, 'FontSize', 20)
set(gca,'YTick',tick,'YTickLabel',xticklabel, 'FontSize', 20)
xlabel('node number', 'FontSize', 20);
ylabel('node number', 'FontSize', 20);
title('Edge betweenness centrality matrix', 'FontSize', 20);

% Global and local diffusion efficiency
[GEdiff,Ediff] = diffusion_efficiency(CIJ);
figure('Color',[1 1 1]);
subplot(1,2,1); imagesc(Ediff); caxis([0 0.05])
xlabel('Node number')
ylabel('Node number')
title('Pair-wise diffusion efficiency');
subplot(1,2,2); stem(GEdiff, 'Linewidth', 1.5);
ylabel('Global diffusion efficiency')
title('Global diffusion efficiency');

% Global and local efficiency
Eloc = efficiency_bin(CIJ,1);
Eglob = efficiency_bin(CIJ); 
figure('Color',[1 1 1]);
subplot(2,1,1); stem(Eloc, 'Linewidth', 2.5);
set(gca,'XTick',tick,'XTickLabel',xticklabel, 'FontSize', 20)
xlabel('Node number', 'FontSize', 20)
ylabel('Local efficiency', 'FontSize', 20)
title('Local efficiency', 'FontSize', 20);
subplot(2,1,2); stem(Eglob, 'Linewidth', 2.5);
ylabel('Global efficiency', 'FontSize', 20)
title('Global efficiency', 'FontSize', 20);

%% %%%%%%%%%%%%%%%% END %%%%%%%%%%%%%%%%%%%%%%%