function [x, hist] = jr_temp(data, alpha, gamma, op, varargin) 

% Computes a complex-valued minimizer x = [x1, ..., xn] of 
%
%   sum_{i=1}^n alpha_i/2 * |B_i x_i - g_i|_C^2 
%                    + sum_{i=0}^n gamma_i/2 * |x_{i+1} - x_i|^2
%
% where |.|_C denotes a the complex 2-norm.
%
% Input: 
% data       == cell array containing the MR Fourier data sets 'g_i'
% alpha      == array containing the regularization parameters
% gamma      == temporal regularization parameter
% op         == cell array: forward operators and their adjoints
%
% varargin   == param  : parameter file
%               show   : shows iteration history: 'true'/'false'
%               x0     : starting point for temporal regularization 
%               gamma0 : weighting for temporal starting point x0
% 
% Output: 
% x          == The computed minimizer 'x' 
% hist       == Struct containing iteration history:
%               - en         == Primal energy
%               - pdres      == Primal-dual residual
%               - it         == Number of iterations
%
% Questions to julian.rasch@wwu.de 

%% Overload varargin
% Defaults 
param           = struct;
param.norm_type = 'iso';
param.pdtol     = 1e-4;
param.niter     = 1e5;
param.int       = 50;
param.tau       = 0.02; 
param.sigma     = 0.02;

show            = true;
gamma0          = 0;


if isempty(varargin) == 0  
    for i=1:2:length(varargin) % overwrites default parameter
	eval([varargin{i},'=varargin{',int2str(i+1),'};']);
    end
end

%% Definitions 
tic; 

% Proximal operators
p_prox = @(x, u1, u2, gamma_1, gamma_2, tau) (x + tau * gamma_1 * u1 + tau * gamma_2 * u2) / ( tau * (gamma_1 + gamma_2) + 1); 
d_prox = @(y,alpha,sigma) (alpha * y) / (alpha + sigma);

% Step sizes 
tau   = param.tau;
sigma = param.sigma; 

% Helper 
vec = @(x) x(:);

%% Check the sizes 
n  = numel(data);

%%  Initialize
for i = 1:n
    % Primal variables
    x{i} = op{1,2}(data{i}); x_old{i} = x{i};
    
    % Dual variables
    y{i} = op{1,1}(x{i}); y_old{i} = y{i};

    % Store the calculations to save computation time
    % Evaluated primals
    B_x{i} = op{1,1}(x{i}); B_x_old{i} = B_x{i};

    % Evaluated duals
    Bt_y{i} = op{1,2}(y{i}); Bt_y_old{i} = Bt_y{i};

end

% History 
hist       = struct;
hist.en    = zeros(param.niter/param.int,1);
hist.pdres = zeros(param.niter/param.int,1);
hist.it    = zeros(param.niter/param.int,1);

%% Do the work
stop_pdres = false;
it = 1;

while (it <= param.niter && ~stop_pdres)
    
    for i = 1:n
        % Compute the updates
        % Dual update
        B_x{i} = op{i,1}(x{i});
        B_x_bar = 2 * B_x{i} - B_x_old{i}; 
        dual_aux = y{i} + sigma .* (B_x_bar - data{i});
        y{i} = d_prox(dual_aux,alpha(i),sigma); % L2-data term

        % Primal update
        Bt_y{i} = op{i,2}(y{i});
        primal_aux = x{i} - tau .* Bt_y{i};
        
        if (i == 1 && exist('x0','var'))
            x{i} = (primal_aux + tau * gamma0 * x0 + tau * gamma(i) * x{i+1}) / (tau * (gamma0 + gamma(i)) + 1); 
        elseif (i == 1 && ~exist('x0','var'))
            x{i} = (primal_aux + tau * gamma(i) * x{i+1}) / (tau * gamma(i) + 1);
        elseif (i == n)
            x{i} = (primal_aux + tau * gamma(i-1) * x{i-1}) / (tau * gamma(i-1) + 1);
        else
            x{i} = p_prox(primal_aux,x{i-1},x{i+1},gamma(i-1),gamma(i),tau);
        end

    end    
    
    % Compute convergence criteria
    if (mod(it, param.int) == 0)
        pdres = 0;
        en    = 0; 
        % Primal-dual residual
        for j = 1:n
            % primal
            pr = norm(vec( (x_old{j} - x{j})./tau - Bt_y_old{j} + Bt_y{j} ),1) / numel(x{j});
            % dual
            dr = norm(vec( (y_old{j} - y{j})./sigma - B_x_old{j} + B_x{j} ),1) / numel(y{j});
            % residual
            pdres = pdres + pr + dr;  

        end
        % Primal energy
        for j = 1:n
            if (j == 1 && exist('x0','var'))
                en = en + alpha(j)/2 * norm(B_x{j}(:) - data{j}(:),2)^2 ...
                    + gamma0/2 * norm(vec(x0 - x{j}))^2;
            elseif (j == 1 && ~exist('x0','var'))
                en = en + alpha(j)/2 * norm(B_x{j}(:) - data{j}(:),2)^2;
            elseif (j==n)
                en = en + alpha(j)/2 * norm(B_x{j}(:) - data{j}(:),2)^2;
            else
                en = en + alpha(j)/2 * norm(B_x{j}(:) - data{j}(:),2)^2 ... 
                    + gamma(j)/2 * norm(vec(x{j} - x{j+1}))^2; 
            end
        end
            
        % Save in history 
        hist.en(it/param.int)    = en;
        hist.pdres(it/param.int) = pdres;
        hist.it(it/param.int)    = it;
        
        % Check whether residual is small enough
        if (pdres < param.pdtol )
            stop_pdres = true;
        end
        % Some output
        fprintf('It.: %6.6d. Primal energy: %6.6d. PD residual: %6.6d > %8.2E. \n', it, en, pdres, param.pdtol);
        
    end
    
    for i=1:n
        % SAVE THE ITERATES
        % Primals
        x_old{i} = x{i};

        % Duals
        y_old{i} = y{i};

        % Evaluated primals
        B_x_old{i} = B_x{i};
     
        % Evaluated duals
        Bt_y_old{i} = Bt_y{i};

   
    end   
    

    % Update iteration number
    it = it + 1;

end

% Shorten history
hist.en = nonzeros(hist.en);
hist.pdres = nonzeros(hist.pdres);
hist.it = nonzeros(hist.it);

% Show history 
if (show)
    figure; 
    subplot(121); plot(hist.en); title('Primal energy'); 
    subplot(122); plot(hist.pdres); title('PD residual');
end

fprintf('Computation time %6.6d minutes.\n',toc/60);

