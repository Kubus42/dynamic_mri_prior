% RUN IT
clear; close all; clc;
% Add the necessary folders
addpath(genpath('.'));

vec = @(x) x(:);

%% Load the phantoms 
load('prior.mat'); % Anatomical prior
load('dynamic_phantom.mat'); % T2 dynamic data

% The phantom has a resolution of 109x91 px and 50 time frames
T  = numel(u_clean);
sz = size(u0_clean); 

%% Create sampling operators and an artificial data set for the prior

% Create a sampling S_0 for the anatomical prior (here: full k-space)
S_0  = sampling_geom(u0_clean, 'full');

% Setup the operator K_0 and its adjoint
K_0{1,1} = @(x) sampling_op(x,S_0);
K_0{1,2} = @(y) sampling_adj(y,S_0);

% Create data f_0 with complex valued Gaussian noise 
% of the anatomical prior
energy = 0.05; % noise level
f0 = fourier_data(K_0{1,1},u0_clean,energy);

% Visualize the sampling 
sampling_vis(f0);

%% Do a TV reconstruction of the prior u0

param_u0  = paramfile_tv; % new parameter file
alpha_u0 = 50; % regularization parameter
[u0, hist_u0] = l2_tv(f0,alpha_u0,K_0,'param',param_u0,'show',false);

% Show the reconstruction and the convergence criteria
figure; 
subplot(131); imagesc(abs(u0)); title('Prior recon (abs)'); axis image; colormap gray;
subplot(132); plot(hist_u0.en); title('Primal energy');
subplot(133); plot(hist_u0.pdres); title('PD residual');

% Compute the gradient with edge parameter 'eta'
eta = 0.05; 
grad_u0 = gradient_direction(u0,'eta',eta);

% Compare the subgradients visually 
figure; 
imagesc(sqrt(sum(abs(grad_u0).^2,3))); 
colormap gray; 


%% Create samplings S and data f for the dynamic phantom u

for t = 1:T
    % Samplings
    nSpokes     = 5;
    S{t}        = sampling_geom(u_clean{t}, 'spokes5','nSpokes',...
        nSpokes,'number',(t-1)*nSpokes + 1);
    K{t,1}      = @(x) sampling_op(x,S{t});
    K{t,2}      = @(y) sampling_adj(y,S{t});
    % Data 
    f{t} = fourier_data(K{t,1},u_clean{t},energy);
end

% Show two consecutive samplings, and their overlap
figure; 
subplot(131); sampling_vis(S{1},1); title('Sampling 1');
subplot(132); sampling_vis(S{2},1); title('Sampling 2');
subplot(133); sampling_vis(S{1} | S{2},1); title('Overlap');

%% One data set
set{1} = 1:60; 
w{1} = 0.1 * ones(numel(set{1}),1);
gamma{1} = 25 * ones(numel(set{1}),1);
alpha{1} = 50 * ones(numel(set{1}),1);

%% Divide the data sets into overlapping bits 
% set{1} = 1:10; 
% set{2} = 11:20;
% set{3} = 21:30; 
% set{4} = 31:40; 
% set{5} = 41:50;
% set{6} = 51:60;
% 
% % Choose the reconstruction parameters 
% % Weights:
% w{1} = 0.1 * ones(numel(set{1}),1);
% w{2} = 0.1 * ones(numel(set{2}),1);
% w{3} = 0.1 * ones(numel(set{3}),1);
% w{4} = 0.1 * ones(numel(set{4}),1);
% w{5} = 0.1 * ones(numel(set{5}),1);
% w{6} = 0.1 * ones(numel(set{6}),1);
% 
% % Temporal regularity: 
% gamma{1} = 25 * ones(numel(set{1}),1);
% gamma{2} = 25 * ones(numel(set{2}),1);
% gamma{3} = 25 * ones(numel(set{3}),1);
% gamma{4} = 25 * ones(numel(set{4}),1);
% gamma{5} = 25 * ones(numel(set{5}),1);
% gamma{6} = 25 * ones(numel(set{6}),1);
% 
% % Data fidelity:
% alpha{1} = 50 * ones(numel(set{1}),1);
% alpha{2} = 50 * ones(numel(set{2}),1);
% alpha{3} = 50 * ones(numel(set{3}),1);
% alpha{4} = 50 * ones(numel(set{4}),1);
% alpha{5} = 50 * ones(numel(set{5}),1);
% alpha{6} = 50 * ones(numel(set{6}),1);

%% Run the reconstruction with dTV
x0 = zeros(size(u0_clean));
for i = 1:numel(set)    
    [u{i},hist_u{i}] = jr_fmri_dtv(f(set{i}), alpha{i}, gamma{i}, w{i}, K(set{i},:), grad_u0,'param.pdtol',1e-03, 'show',false,'param.niter',50,'x0',x0);
    x0 = u{i}{10};
end

% Check the convergence (e.g. for the first set{1})
figure; 
subplot(121); plot(hist_u{1}.en); title('Primal energy');
subplot(122); plot(hist_u{1}.pdres); title('PD residual');

%% Save the results

% Store the results in an array 
res = zeros([sz,T]);
k = 1;
for i= 1:1%numel(set)
    for j = 1:numel(u{i})
        res(:,:,k) = abs(u{i}{j});
        k = k + 1;
    end
end

jfScrollImage(res);

% Note that I did not eliminate the overlap here! 

save_name = fullfile('results','name');
save(save_name', 'u', 'hist_u', 'res', 'S');













