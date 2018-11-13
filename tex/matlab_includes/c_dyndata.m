%% Create samplings S and data f for the dynamic phantom u
for t = 1:T
    % Samplings
    nSpokes     = 18;
    S{t}        = sampling_geom(u_clean{t}, 'spokes5','nSpokes',...
        nSpokes,'number',(t-1)*nSpokes + 1);
    K{t,1}      = @(x) sampling_op(x,S{t});
    K{t,2}      = @(y) sampling_adj(y,S{t});
    % Data 
    f{t} = fourier_data(K{t,1},u_clean{t},energy);
end

% Show two consecutive samplings, and their overlap
figure; 
subplot(131); sampling_vis(S{1},1); 
subplot(132); sampling_vis(S{2},1);
subplot(133); sampling_vis(S{1} | S{2},1);