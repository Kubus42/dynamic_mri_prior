function [x, hist] = jr_temp_tv(data, alpha, gamma, op, varargin) 

% Computes a complex-valued minimizer x = [x1, ..., xn] of 
%
%   sum_{i=1}^n alpha_i/2 * |B_i x_i - g_i|_C^2 + |\grad x_i |_1
%   + gamma_i/2 * |x_{i+1} - x_i|^2
%
% where |.|_C denotes a the complex 2-norm.
%
% Input: 
% data       == cell array containing the MR Fourier data sets 'g_i'
% alpha      == array containing the regularization parameters
% gamma      == temporal regularization parameter
% op         == cell array containing the forward operators and their adjoints
%
% varargin   == param: Parameter file
%               show: Shows iteration history: 'true'/'false'
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
gamma_0         = gamma(1);
x0              = zeros(size(op{1,2}(data{1})));

show            = true;

if isempty(varargin) == 0  
    for i=1:2:length(varargin) % overwrites default parameter
	eval([varargin{i},'=varargin{',int2str(i+1),'};']);
    end
end

%% Definitions 
tic; 

% Proximal operators
p_prox  = @(x, u1, u2, gamma_1, gamma_2, tau) (x + tau * gamma_1 * u1 + tau * gamma_2 * u2) / ( tau * (gamma_1 + gamma_2) + 1); 
d_prox1 = @(y, data, alpha, sigma) alpha./(alpha+sigma) .* (y - sigma .* data); % alpha * 1/2 * \|Bv-g\|_2^2
d_prox2 = @(y, reg) prox_dual_l1(y, reg, param.norm_type); % TV

% Functions 
TV_iso = @(x) tv(x,param.norm_type);

% Step sizes 
tau   = param.tau;
sigma = param.sigma; 

% Helper 
vec = @(x) x(:);

%% Compute the strong convexity parameter mu of the primal part
% Set up the Hessian
m_diag = gamma + [gamma_0; gamma(1:end-1)];
M = diag(m_diag,0) + diag(-gamma(1:end-1),-1) + diag(-gamma(1:end-1),1);
mu = min(eig(M));

% Check for zeros
if mu <= 1e-4
    mu = 0;
    disp('No strong convexity detected: Using standard version.');
else
    disp(['Strong convexity detected: Using accelerated version with mu = ',num2str(mu),'.']);
end

%% Check the sizes 
n  = numel(data);

%%  Initialize
for i = 1:n
    % Primal variables
    x{i} = op{1,2}(data{i}); x_old{i} = x{i};
    
    % Dual variables
    y1{i} = op{1,1}(x{i}); y1_old{i} = y1{i};
    y2{i} = grad(x{i});    y2_old{i} = y2{i};

    % Store the calculations to save computation time
    % Evaluated primals
    B_x{i} = op{1,1}(x{i}); B_x_old{i} = B_x{i};
    G_x{i} = grad(x{i});    G_x_old{i} = G_x{i};

    % Evaluated duals
    Bt_y1{i} = op{1,2}(y1{i}); Bt_y1_old{i} = Bt_y1{i};
    Gt_y2{i} = -div(y2{i});    Gt_y2_old{i} = Gt_y2{i};

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
        G_x{i} = grad(x{i});
                

        B_x_bar = 2 * B_x{i} - B_x_old{i}; 
        G_x_bar = 2 * G_x{i} - G_x_old{i};
                
        dual_1 = y1{i} + sigma .* B_x_bar; 
        dual_2 = y2{i} + sigma .* G_x_bar; 
        
        y1{i} = d_prox1(dual_1,data{i},alpha(i),sigma); % L2-data term
        y2{i} = d_prox2(dual_2,1); % TV

        % Primal update
        Bt_y1{i} = op{i,2}(y1{i});
        Gt_y2{i} = -div(y2{i});
        
        primal = x{i} - tau .* (Bt_y1{i} + Gt_y2{i});
        
        if (i == 1 && n == 1)
            x{i} = primal;
        elseif (i == 1 && n > 1) % if x0 is zero, x1 might be not so good
                x{i} = p_prox(primal,x0,x{i+1},gamma_0,gamma(i),tau);
        elseif (i == n && n > 1 ) % set gamma_{n+1} = 0
            x{i} = p_prox(primal,x{i-1},x{i-1},gamma(i-1),0,tau);
        else
            x{i} = p_prox(primal,x{i-1},x{i+1},gamma(i-1),gamma(i),tau);
        end
    end    
    
    % Compute convergence criteria
    if (mod(it, param.int) == 0)
        pdres = 0;
        en    = 0; 
        % Primal-dual residual
        for j = 1:n
            % primal
            pr = norm(vec((x_old{j} - x{j})./tau - (Bt_y1_old{j} + Gt_y2_old{j}) + (Bt_y1{j} + Gt_y2{j})),1) / numel(x{j});
            % dual
            dr1 = norm(vec((y1_old{j} - y1{j})./sigma - B_x_old{j} + B_x{j}),1) / numel(y1{j});
            dr2 = norm(vec((y2_old{j} - y2{j})./sigma - G_x_old{j} + G_x{j}),1) / numel(y2{j}); 
            % residual
            pdres = pdres +  pr + dr1 + dr2;  

        end
        % Primal energy
        for j = 1:n
            if j == n
                en = en + alpha(j)/2 * norm(B_x{j}(:) - data{j}(:),2)^2 ... 
                    + TV_iso(x{j});
            else
                en = en + alpha(j)/2 * norm(B_x{j}(:) - data{j}(:),2)^2 ... 
                    + TV_iso(x{j}) ...
                    + gamma(i)/2 * norm(vec(x{j} - x{j+1}))^2; 
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
        y1_old{i} = y1{i};
        y2_old{i} = y2{i};
        
        % Evaluated primals
        B_x_old{i} = B_x{i};
        G_x_old{i} = G_x{i};
                
        % Evaluated duals
        Bt_y1_old{i} = Bt_y1{i};
        Gt_y2_old{i} = Gt_y2{i};
           
    end   
    
    % Update the step sizes 
    theta = 1/sqrt(1+2*mu*tau);
    tau   = tau * theta;
    sigma = sigma / theta;

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

