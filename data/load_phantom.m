function [t1, t2] = load_phantom(varargin)

% This function loads a T1 and a T2 phantom from a BrainWeb data set (URL).
% The original resolution is 181x217x181 px. Choose the slice you want to
% consider with 'slice', and change the resolution via the parameter 's'.
%
% Input: 
%
% sz        ==   size of the phantom in px
% slice     ==   slice of the phantom to consider (choose one) 
%                    slice_x 
%                    slice_y
%                    slice_z    
% s         ==   factor for the change of resolution: 
%                   s<1 is downsampling, s>1 is upsampling 
% show      ==   show the phantoms in a figure
% print_res ==   print the results
% 
%
% 

% Defaults
sz        = [181,217,181];
slice_x   = 1:sz(1);
slice_y   = 1:sz(2);
slice_z   = 1:sz(3);
s         = 1; 
show      = true; 
print_res = true; 

if isempty(varargin) == 0  
    for i=1:2:length(varargin) % overwrites default parameter
	eval([varargin{i},'=varargin{',int2str(i+1),'};']);
    end
end

% Open the 3D phantoms
% t1 contrast
fid1 = fopen('t1_phantom.rawb','r');
I1=fread(fid1,prod(sz));
t1_tmp = reshape(I1, sz);
fclose(fid1);

% t2 contrast
fid2 = fopen('t2_phantom.rawb','r');
I2=fread(fid2,prod(sz));
t2_tmp = reshape(I2, sz);
fclose(fid2);

% If not specified, select default transversal slice 
if (numel(slice_x) > 1 && numel(slice_y) > 1 && numel(slice_z) > 1)
    slice_z = 64;
end

% Select the slice 
t1 = squeeze(t1_tmp(slice_x,slice_y,slice_z)); 
t2 = squeeze(t2_tmp(slice_x,slice_y,slice_z));

% Flip direction
t1 = t1';
t2 = t2';
t1 = t1(end:-1:1,:);
t2 = t2(end:-1:1,:);

% Assign name
if (numel(slice_x) == 1 && numel(slice_y) > 1 && numel(slice_z) > 1)
    s_name = 'sagittal';
elseif (numel(slice_x) > 1 && numel(slice_y) == 1 && numel(slice_z) > 1)
    s_name = 'coronal';
elseif (numel(slice_x) > 1 && numel(slice_y) > 1 && numel(slice_z) == 1)
    s_name = 'transversal';
end

% Scale the images to [0,1]
t1 = (t1 - min(t1(:))) / (max(t1(:)) - min(t1(:)));
t2 = (t2 - min(t2(:))) / (max(t2(:)) - min(t2(:)));

% Threshold the images to be piecewise constant
thresh_t1 = multithresh(t1,4);
thresh_t2 = multithresh(t2,4);

val_t1 = [0, 0.2, 0.45, 0.6, 1];
val_t2 = [0, 0.3, 0.4, 0.6, 1];

t1 = imquantize(t1,thresh_t1, val_t1);
t2 = imquantize(t2,thresh_t2,val_t2);

% Resize the images
t1 = imresize(t1,s,'nearest');
t2 = imresize(t2,s,'nearest');

% Show them 
if (show)
    figure;
    subplot(121);
    imagesc(t1); axis image; colormap gray; colorbar; title(['phantom T1 contrast ',s_name]);
    subplot(122);
    imagesc(t2); axis image; colormap gray; colorbar; title(['phantom T2 contrast ',s_name]);
end

% Print them 
if (print_res && exist('results','dir'))
    print_name = fullfile('results', ['phantom_t1_', s_name]);
    writeImage(print_name,shrinkImage(t1,0,1),gray(256));    
    print_name = fullfile('results', ['phantom_t2_', s_name]);
    writeImage(print_name,shrinkImage(t2,0,1),gray(256));   
elseif (print_res && ~exist('results','dir'))
    mkdir('results')
    print_name = fullfile('results', ['phantom_t1_', s_name]);
    writeImage(print_name,shrinkImage(t1,0,1),gray(256));    
    print_name = fullfile('results', ['phantom_t2_', s_name]);
    writeImage(print_name,shrinkImage(t2,0,1),gray(256));
end
    
    
    