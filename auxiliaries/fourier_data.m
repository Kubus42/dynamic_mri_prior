function data = fourier_data(B, v, energy)
% 
% Creates Fourier data of the input image 'v', corrupted by Gaussian 
% noise with an energy of 'energy'. The operator 'B' is the Fourier 
% operator including possible undersampling. 
%
% Input: 
% B      == Fourier/MR operator including undersampling
% v      == Input image 
% energy == Energy of the noise
%
% Questions to: julian.rasch@wwu.de

% Project
data_c = B(v);
mask   = (data_c ~= 0);

% Calculate normalization 
data_en = norm(data_c(:));
sz = size(data_c);

noise_re = randn(sz) .* mask;
noise_im = randn(sz) .* mask;

% noise_en = norm(noise_re(:))/2 + norm(noise_im(:))/2;

noise_en = norm(noise_re(:) + noise_im(:))/2;

% Standard deviation split into real and imaginary part
sigma = energy * data_en/noise_en;
noise = sigma/2 * noise_re + 1i * sigma/2 * noise_im;

noise_aux = norm(noise(:));
% fprintf('Ratio: %d\n', noise_aux/data_en);

% Add noise to the data
data = data_c + noise;

end