%% Run the reconstruction
for i = 1:numel(set)
    [u{i},hist_u{i}] = jr_fmri(f(set{i}), alpha{i}, gamma{i}, w{i}, K(set{i},:), p_art,'param.pdtol',5e-02, 'show',false);
end

% Check the convergence (e.g. for the first set{i})
figure; 
subplot(121); plot(hist_u{1}.en); title('Primal energy');
subplot(122); plot(hist_u{1}.pdres); title('PD residual');
