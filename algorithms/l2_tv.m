function [x, hist] = l2_tv(data, alpha, op, varargin) 

% Computes the solution of 
%
%   min_{x in C} alpha/2 * |Bx - g|_C^2 + |\grad x|_1  
%
% where |.|_C denotes a the complex 2-norm.
%
% Input: 
% data       == MR Fourier data 'g'
% alpha      == regularization parameter
% op         == cell array containing the forward operator and its adjoint
%
% varargin   == param: Parameter file (create with paramfile.m)
% 
% Output: 
% result     == The computed minimizer 'x' 
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
param.niter     = 5e4;
param.int       = 100;
param.tau       = 0.05; 
param.sigma     = 0.05;

show           = false;

if isempty(varargin) == 0  
    for i=1:2:length(varargin) % overwrites default parameter
	eval([varargin{i},'=varargin{',int2str(i+1),'};']);
    end
end

%% Definitions 
tic; 

% Proximal operators
p_prox  = @(x) x;
d_prox1 = @(y1, sigma) alpha./(alpha+sigma) .* (y1-sigma .* data); % alpha/2 * \|Bv-g\|_2^2
d_prox2 = @(y2, reg) prox_dual_l1(y2, reg, param.norm_type); % TV

% Step sizes 
tau   = param.tau;
sigma = param.sigma; 

% Helper 
vec = @(x) x(:);

%%  Initialize
% Primal variables
x = op{1,2}(data); x_old = x;

% Dual variables
y1 = op{1,1}(x); y1_old = y1;
y2 = grad(x);    y2_old = y2;

% Store the calculations to save computation time
% Evaluated primals
B_x = op{1,1}(x); B_x_old = B_x;
G_x = grad(x);    G_x_old = G_x;

% Evaluated duals
Bt_y1 = op{1,2}(y1); Bt_y1_old = Bt_y1;
Gt_y2 = -div(y2);    Gt_y2_old = Gt_y2;

% History 
hist       = struct;
hist.en    = zeros(param.niter/param.int,1);
hist.pdres = zeros(param.niter/param.int,1);
hist.it    = zeros(param.niter/param.int,1);

%% Do the work
stop_pdres = false;
it = 1;

while (it <= param.niter && ~stop_pdres)
   
    % Compute the updates
    % Dual update
    B_x = op{1,1}(x);
    G_x = grad(x);
    
    B_x_bar = 2 * B_x - B_x_old; 
    G_x_bar = 2 * G_x - G_x_old;
    
    y1 = d_prox1(y1 + sigma .* B_x_bar,sigma); % L2-data term
    y2 = d_prox2(y2 + sigma .* G_x_bar,1); % TV
    
    % Primal update
    Bt_y1 = op{1,2}(y1);
    Gt_y2 = -div(y2);

    x = p_prox(x - tau .* (Bt_y1 + Gt_y2)); 
    
    % Compute convergence criteria
    if (mod(it, param.int) == 0)
        
        % Primal-dual residual
        % primal
        pr = norm(vec((x_old - x)./tau - (Bt_y1_old + Gt_y2_old) + (Bt_y1 + Gt_y2)),1);
        % dual
        dr1 = norm(vec((y1_old - y1)./sigma - B_x_old + B_x),1);
        dr2 = norm(vec((y2_old - y2)./sigma - G_x_old + G_x),1);  
        % residual
        pdres = pr/numel(x) + dr1/numel(y1) + dr2/numel(y2);  
     
        % Check whether residual is small enough
        if (pdres < param.pdtol )
            stop_pdres = true;
        end
        
        % Primal energy
        en = alpha/2 * norm(B_x(:) - data(:),2)^2 + tv(x,param.norm_type);
        
        % Some output
        fprintf('It.: %6.6d. Primal energy: %6.6d. PD residual: %6.6d > %8.2E. \n', it, en, pdres, param.pdtol);
        
        % Save in history 
        hist.en(it/param.int)    = en;
        hist.pdres(it/param.int) = pdres;
        hist.it(it/param.int)    = it;
    
    end
   
    % SAVE THE ITERATES
    % Primals
    x_old = x;
    
    % Duals
    y1_old = y1;
    y2_old = y2;
    
    % Evaluated primals
    B_x_old = B_x;
    G_x_old = G_x;
    
    % Evaluated duals
    Bt_y1_old = Bt_y1;
    Gt_y2_old = Gt_y2;
   
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

