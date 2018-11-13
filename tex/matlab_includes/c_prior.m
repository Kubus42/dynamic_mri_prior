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
% sampling_vis(f0);

%% Do a TV reconstruction of the prior u0

param_u0  = paramfile_tv; % new parameter file
alpha_u0  = 75; % regularization parameter
[u0, hist_u0] = l2_tv(f0,alpha_u0,K_0,'param',param_u0,'show',false);

% Compute the true subgradient from the optimality condition
p_true = K_0{1,2}(alpha_u0 * (f0 - K_0{1,1}(u0)));

% Create an artificial subgradient with edge parameter 'eta'
eta = 0.025; 
p_art = artificial_subgradient(u0,'eta',eta);
