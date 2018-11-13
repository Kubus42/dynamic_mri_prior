function img = sampling_adj(data,geom)

% MRI adjoint/inverse operator for sampling_op.m.
%
% Input: 
% data   ==   2-dimensional Fourier data set 
% geom   ==   2-dimensional sampling geometry
%
% Output  
% img    ==   Zero-padded inverse FT of the data set 'data'
%
% Questions to julian.rasch@gmx.de

% Check the input 
sz = size(data);
if ( sz(2) == 1 || numel(data) ~= numel(geom) ) 
    error('Input does not match or has the wrong dimension.');
end

% Apply geometry
data = data .* geom;

% Do the adjoint/inverse Fourier transform
% img = sqrt(prod(sz)) * real(ifft2(data));
img = sqrt(prod(sz)) * ifft2(data);

end
