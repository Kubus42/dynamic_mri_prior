clear; close all; clc,
load('BrainWeb_area_fMRI_data_50_randomSpokes18.mat'); 
load('anatomical_T1_prior_109x91_alpha_ap125.mat');
addpath(genpath('../../'));

%% Choose the number of samplings 
samp = 1:50;

%% 
for i = samp;
    sampling_vis(f{i});
    
    set(gca, 'visible', 'off')
    set(gcf, 'color', 'w');
    export_fig(['images/smp/smp_', num2str(i)]);
    
end

%% 

for i = samp;
    recon = abs(K{i,2}(f{i}));
    
    figure; imagesc(recon); colormap gray; axis image; caxis([0,1]);
    set(gca, 'visible', 'off')
    set(gcf, 'color', 'w');
    export_fig(['images/ift/ift_', num2str(i)]);
    
end

%% Subgradient
figure;
imagesc(p_aux); colormap gray; axis image; 
set(gca, 'visible', 'off')
set(gcf, 'color', 'w');
export_fig('subgradient');

%% Prior
sampling_vis(data_ap);

set(gca, 'visible', 'off')
set(gcf, 'color', 'w');
export_fig('smp_prior');

%% 
figure; imagesc(img_ap); colormap gray; axis image; caxis([0,1]);

set(gca, 'visible', 'off')
set(gcf, 'color', 'w');
export_fig('prior');


%% Reconstructions

rec = 1:50;

for i = rec;
    recon = abs(u_clean{i});
    
    figure; imagesc(recon); colormap gray; axis image; caxis([0,1]);
    set(gca, 'visible', 'off');
    set(gcf, 'color', 'w');
    export_fig(['images/gt/gt_', num2str(i)]);
    
    close all
    
end



