function nn = tv(f, varargin)
%
% Computes | ||grad f|| |_1, 
%
% where ||.|| can be an l1-norm or l2-norm, i.e.
% anisotropic (default) or isotropic.
%
% Input: 
% f          ==   data
% varargin   ==   'norm_type': 'aniso' / 'iso' 

grad_f = grad(f);
sz     = size(grad_f);

% Helper 
vec = @(x) x(:);

% Anisotropic 
if (isempty(varargin) == 1 || strcmp(varargin{1}, 'aniso'))
    nn = norm(vec(grad_f),1); 
% Isotropic 
elseif (strcmp(varargin{1}, 'iso') && numel(sz(sz~=1)) > 1 )
    nn   = reshape(grad_f,[prod(sz(1:end-1)),sz(end)]);
    nn   = sum(vec(sqrt(sum(abs(nn).^2,2))));
else 
    error('Check input!');
end
end
