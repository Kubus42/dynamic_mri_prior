% RUN IT
% A description of how this code works can be found in the 'tex'-folder.
clear; close all; clc;
% Add the necessary folders
addpath(genpath('.'));

vec = @(x) x(:);

%% Load the phantoms and the tv-prior
load('prior.mat'); % Anatomical prior
load('dynamic_phantom.mat'); % T2 dynamic data

T  = numel(u_clean);
sz = size(u0_clean); 

load('tv_prior.mat');
% load('data5.mat');

%% Setup the samplings

energy = 0.05; % noise level

for t = 1:T
    % Samplings
    nSpokes     = 10;
    S{t}        = sampling_geom(u_clean{t}, 'spokes5','nSpokes',...
        nSpokes,'number',(t-1)*nSpokes + 1);
    K{t,1}      = @(x) sampling_op(x,S{t});
    K{t,2}      = @(y) sampling_adj(y,S{t});
    % Data 
    f{t} = fourier_data(K{t,1},u_clean{t},energy);
end

  
%% Run temporal + TV

a = 500; b = 500;

set{1} = 1:T;
alpha{1} = a * ones(numel(set{1}),1);
gamma{1} = b * ones(numel(set{1}),1);

for i = 1:numel(set)    
    [u_tt{i},hist_u_tt{i}] = jr_temp_tv(f(set{i}), alpha{i}, gamma{i}, K(set{i},:),'show',false,'param.niter',1e5);
end


%% Visualize and save the results
% This is still to do, but I guess for now you will find a way to do it
% yourself. But I will prepare it. For now we just save the results.

% Store the results in an array 
res = zeros([sz,T]);
k = 1;
for i= 1:numel(set)
    for j = 1:numel(u_tt{i})
        res(:,:,k) = abs(u_tt{i}{j});
       % res(:,:,k) = (res(:,:,k) - min(vec(res(:,:,k))) ) / ...
        %    ( max(vec(res(:,:,k))) - min(vec(res(:,:,k))) );
        k = k + 1;
    end
end

% jfScrollImage(res);

% Note that I did not eliminate the overlap here! 
name = ['result',num2str(nSpokes),'ga_tmp_tv_alpha',num2str(a),'_gamma',num2str(b),'.mat'];
save(name, 'u_tt', 'hist_u_tt', 'res', 'S');
