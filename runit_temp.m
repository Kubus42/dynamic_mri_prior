% RUN IT
% A description of how this code works can be found in the 'tex'-folder.
clear; close all; clc;
% Add the necessary folders
addpath(genpath('.'));

vec = @(x) x(:);

%% Load the phantoms 
load('prior.mat'); % Anatomical prior
load('dynamic_phantom.mat'); % T2 dynamic data

T  = numel(u_clean);
sz = size(u0_clean); 

load('tv_prior.mat');
load('data5.mat');

%% Run only temporal
set{1} = 1:T;
alpha{1} = 50 * ones(numel(set{1}),1);
gamma{1} = 15 * ones(numel(set{1}),1);

for i = 1:numel(set)    
    [u_temp{i},hist_u_temp{i}] = jr_temp(f(set{i}), alpha{i}, gamma{i}, K(set{i},:),'show',false,'param.niter',25000);
end


%% Visualize and save the results
% This is still to do, but I guess for now you will find a way to do it
% yourself. But I will prepare it. For now we just save the results.

% Store the results in an array 
res = zeros([sz,T]);
k = 1;
for i= 1:numel(set)
    for j = 1:numel(u_temp{i})
        res(:,:,k) = abs(u_temp{i}{j});
       % res(:,:,k) = (res(:,:,k) - min(vec(res(:,:,k))) ) / ...
        %    ( max(vec(res(:,:,k))) - min(vec(res(:,:,k))) );
        k = k + 1;
    end
end

% jfScrollImage(res);

% Note that I did not eliminate the overlap here! 

save('result5ga_tmp_alpha50_gamma15.mat', 'u_temp', 'hist_u_temp', 'res', 'S');
