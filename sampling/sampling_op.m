function data = sampling_op(img, geom)

% MRI forward operator for geometries 'geom' provided by
% createSamplingGeometry.m. The binary map 'geom' defines the frequencies
% that are sampled.
%
% Input: 
% img    ==   2-dimensional image to be sampled
% geom   ==   2-dimensional sampling geometry
%
% Output: 
% data   ==   Fourier data set sampled at frequencies in 'geom'
% 
% Use sampling_adj.m for the adjoint.
% 
% Questions to julian.rasch@wwu.de

% Check the input 
sz = size(img);
if ( sz(2) == 1 || numel(img) ~= numel(geom) ) 
    error('Input does not match or has the wrong dimension.');
end

% Do the Fast Fourier Transform
data = fft2(img)/sqrt(prod(sz));

% Apply geometry
data = data .* geom;

end

