function sampling_vis(data,varargin) 

% Visualizes the Fourier frequencies in 'data'.
%
% Input: 
% sampling == Binary map for the Fourier sampling
% data     == Fourier data
%
% Questions to julian.rasch@wwu.de

% Move zero frequency to the middle
data = fftshift(data);

% Display
if (isempty(varargin) == 1)
figure; 
end
imagesc(log(abs(data))); colormap hot; axis image; 
title('k-space data');

end