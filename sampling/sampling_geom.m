function [geom , numberOfSamplingPoints ]= sampling_geom(img, type, varargin)

% Creates a sampling geometry for sampling_op.m and sampling_adj.m.
%
% Provides a binary map with the same size as 'img', which defines a sampling
% geometry for the k-space, i.e. a pixel has value 1, thus is a sampling point,
% if it is contained in the specified geometry 'type', 0 else. 
%
% Types:
% full    == Full sampling of the k-space
% half    == Every second line of the k-space
% lines1  == Uniformly random lines of the k-space
% lines2  == Gaussian random lines of the k-space
% center1 == Center block of the k-space
% center2 == Long center block 
% spokes1 == Sampling on radial uniform spokes
% spokes2 == Sampling on radial logarithmic spokes
% spokes3 == Sampling on radial sqrt spokes
% spokes4 == Sampling on random uniform spokes
% spokes5 == Sampling on spokes chosen via the golden ratio
% spiral1 == Sampling on a regular spiral 
% spiral2 == Sampling on a parabolic spiral
% spiral3 == Sampling on a ?? spiral
% spiral4 == Sampling on a regular spiral + dense center
% random1 == Fully random sampling 
% random2 == Random + dense center sampling
%
% The following parameters can be overloaded via ('parameterName', value):
%
% nSpokes         == Number of spokes 
% nSamplesOnSpoke == Number of samples per spoke
% nTwists         == Winding number of the spirals 
% K               == Number of random spokes for 'lines1/2'
% sigma           == Standard deviation for 'lines2'
% C = [C1,C2]     == Size of the center block for 'center1/2'
% flp             == Flip lines to vertical direction '1' or not '0'
% show            == Show the sampling '1' or not '0'
% nRandSamples    == Number of random sampling points for 'random1/2'
% radCenter       == Size of the dense center for 'random2'
% number          == Starting spoke for the golden ratio 
% 
% Example:
% geom = createSamplingGeometry(img, 'spiral1', 'nTwists', 55);
%
% Questions to: julian.rasch@wwu.de


% =========================================================================
% EVALUATE INPUT

N = size(img);

% The following parameters can be overloaded with varargin:
nTwists         = 10; % Number of twists for the spiral
show            = 0; % Show the sampling
nSpokes         = 20; % Number of radial spokes
nSamplesOnSpoke = N(1); % Number of samples on spoke
K               = floor(N(1)/5); % Number of random 'lines'
sigma           = 0.25; % Standard deviation of the Gaussian for 'lines2'
flp             = 0; % Flip lines to vertical lines 
C               = [floor(N(1)/3),floor(N(2)/3)]; % Size of the center block
nRandSamples    = floor(prod(N)/2); % Number of random samples for 'random1/2'
radCenter       = floor(N(1)/20); % Size of the dense center for 'random2'
number          = 1; % Starting number for the golden ratio

% Overload
for i=1:2:length(varargin) % overwrites default parameter
eval([varargin{i},'=varargin{',int2str(i+1),'};']);
end

% =========================================================================
switch type
% Full sampling
    case 'full'
    geom = ones(N);
    fprintf('Full k-space sampling.\n');

% Every second line
    case 'half'
    geom = ones(N);
    if (flp ==1)
        for i = 2:2:N(2)
            geom(:,i) = 0;
        end
    else
        for i = 2:2:N(1)
            geom(i,:) = 0;
        end
    end
    fprintf('Half k-space sampling.\n');

% K uniformly distributed random lines 
    case 'lines1'
    geom = zeros(N);
    if (flp == 1) % Flipped sampling
        n = N(2);
        randmask = zeros([1,n]);
        ii = 1;
        while (numel(find(randmask ~= 0)) <= K)
            aux = randi([1 n], 1, 1);
            if (isempty(find(randmask == aux,1)))
                randmask(ii) = aux;
                ii = ii + 1;
            end
        end

        geom(:, randmask(1:K)) = 1;
    else % Normal sampling
        n = N(1);
        randmask = zeros([n,1]);
        ii = 1;
        while (numel(find(randmask ~= 0)) <= K)
            aux = randi([1 n], 1, 1);
            if (isempty(find(randmask == aux,1)))
                randmask(ii) = aux;
                ii = ii + 1;
            end
        end

        geom(randmask(1:K), :) = 1;
    end
    fprintf(['Uniformly random ',num2str(K),'-line k-space sampling.\n']);
    
% K normally distributed random lines 
    case 'lines2'
    geom = zeros(N);
    mind = zeros([K,1]);
    
    if (flp == 1)
        % Draw from a Gaussian distribution, but only within [-1,1] 
        ii = 1; 
        while (numel(find(mind == 0)) > 0 )
            aux = sigma * randn;
            if (abs(aux)<=1)
                % Scale to [1,N] and take the nearest integer
                aux = uint16((aux + 1) * floor((N(1)-1)/2) + 1); 
                if (isempty(find(mind == aux,1)))
                    mind(ii) = aux; 
                    ii = ii + 1;
                end
            end
        end
        geom(mind,:) = 1; 
    else
        % Draw from a Gaussian distribution, but only within [-1,1] 
        ii = 1; 
        while (numel(find(mind == 0)) > 0 )
            aux = sigma * randn;
            if (abs(aux)<=1)
                % Scale to [1,N] and take the nearest integer
                aux = uint16((aux + 1) * floor((N(2)-1)/2) + 1); 
                if (isempty(find(mind == aux,1)))
                    mind(ii) = aux; 
                    ii = ii + 1;
                end
            end
        end
        geom(:,mind) = 1; 
    end
    fprintf(['Gauss random ',num2str(K),'-line k-space sampling.\n']);
    
    
% Center block
    case 'center1'
    geom = zeros(N);
    geom(1:C(1),1:C(2)) = 1;
    geom = circshift(geom, [floor(N(1)/2)-floor(C(1)/2),ceil(N(2)/2)-floor(C(2)/2)]);
    fprintf(['Center block of ',num2str(C(1)),'x',num2str(C(2)),' px.\n']);
    
% Long center block
    case 'center2'
    geom = zeros(N);
    if (flp == 1)
        geom(1:C(1),:) = 1;
        geom = circshift(geom, [floor(N(1)/2)-floor(C(1)/2),0]);
    else 
        geom(:,1:C(2)) = 1;
        geom = circshift(geom, [0,floor(N(2)/2)-floor(C(2)/2)]);
    end
    fprintf(['Center block of ',num2str(C(1)),'x',num2str(C(2)),' px.\n']);

% Radial standard spokes
    case 'spokes1'
    r = -1:1/floor(nSamplesOnSpoke/2):1; r(end) = [];
    phi = 0:pi/ceil(nSpokes):pi; phi(end) = [];
    
    f = @(x) x;
    X = floor(N(1)/2) * f(r)' * cos(phi);
    Y = floor(N(2)/2) * f(r)' * sin(phi);

    X = uint16(X+floor(N(1)/2)) - 1;
    Y = uint16(Y+floor(N(2)/2)) - 1;

    geom = zeros(N(1),N(2));
    s = size(X);
    for j = 1:s(1)
        for i = 1:s(2)
           x = X(j,i) + 1;
           y = Y(j,i) + 1;
           geom(x,y) = 1;
        end
    end
    
% Radial logarithmic spokes    
    case 'spokes2'
    r = -1:1/floor(nSamplesOnSpoke/2):1; r(end) = [];
    phi = 0:pi/ceil(nSpokes):pi; phi(end) = [];

    f = @(x) sign(x) .* abs(x).^2;

    X = floor(N(1)/2) * f(r)' * cos(phi);
    Y = floor(N(2)/2) * f(r)' * sin(phi);

    X = uint16(X+floor(N(1)/2)) - 1;
    Y = uint16(Y+floor(N(2)/2)) - 1;

    geom = zeros(N(1),N(2));
    s = size(X);
    for j = 1:s(1)
        for i = 1:s(2)
           x = X(j,i) + 1;
           y = Y(j,i) + 1;
           geom(x,y) = 1;
        end
    end
    
% Radial sqrt spokes    
    case 'spokes3'
    r = -1:1/floor(nSamplesOnSpoke/2):1; r(end) = [];
    phi = 0:pi/ceil(nSpokes):pi; phi(end) = [];

    f = @(x) sign(x) .* abs(x).^(0.2);

    X = floor(N(1)/2) * f(r)' * cos(phi);
    Y = floor(N(2)/2) * f(r)' * sin(phi);

    X = uint16(X+floor(N(1)/2)) - 1;
    Y = uint16(Y+floor(N(2)/2)) - 1;

    geom = zeros(N(1),N(2));
    s = size(X);
    for j = 1:s(1)
        for i = 1:s(2)
           x = X(j,i) + 1;
           y = Y(j,i) + 1;
           geom(x,y) = 1;
        end
    end
    geom(64,64) = 1;
    
% Radial random spokes
    case 'spokes4'
    r = -1:1/floor(nSamplesOnSpoke/2):1; r(end) = [];
    p1 = 0;
    p2 = pi;
    phi = sort((p2-p1).*rand(nSpokes+1,1) + p1)';   
    
    % phiphi = 0:pi/ceil(nSpokes):pi; phi(end) = [];
    
    f = @(x) x;
    X = floor(N(1)/2) * f(r)' * cos(phi);
    Y = floor(N(2)/2) * f(r)' * sin(phi);

    X = uint16(X+floor(N(1)/2)) - 1;
    Y = uint16(Y+floor(N(2)/2)) - 1;

    geom = zeros(N(1),N(2));
    s = size(X);
    for j = 1:s(1)
        for i = 1:s(2)
           x = X(j,i) + 1;
           y = Y(j,i) + 1;
           geom(x,y) = 1;
        end
    end
    
% Radial spokes according to the Golden Angle 
    case 'spokes5'
        
    r = -1:1/floor(nSamplesOnSpoke/2):1; r(end) = [];
    ga = (2*pi) / (1 + sqrt(5));
    
    phi = number:(number+nSpokes-1); % Determine the amount of spokes
    phi = mod(phi * ga, 2*pi);

    f = @(x) x;
    X = floor(N(1)/2) * f(r)' * cos(phi);
    Y = floor(N(2)/2) * f(r)' * sin(phi);

    X = uint16(X+floor(N(1)/2)) - 1;
    Y = uint16(Y+floor(N(2)/2)) - 1;

    geom = zeros(N(1),N(2));
    s = size(X);
    for j = 1:s(1)
        for i = 1:s(2)
           x = X(j,i) + 1;
           y = Y(j,i) + 1;
           geom(x,y) = 1;
        end
    end     
    
% Standard spiral    
    case 'spiral1'
    f = @(x) x;
    x= -1*pi*nTwists : 0.005 : pi*nTwists;
    r= 0:1/(length(x)-1):1;
    
    X=sin(x).*f(r);  Y=cos(x).*f(r);     
    % Normalize X,Y
    X = X/max(X);   Y = Y/max(Y);
    
    X = uint16(N(1)/2 *(X+1)) - 1;
    Y = uint16(N(2)/2 *(Y+1)) - 1;

    geom = zeros(N(1),N(2));

    for i = 1:length(X)
       x = X(i)+1;
       y = Y(i)+1;
       geom(x,y) = 1;
    end

% Parabolic spiral
    case 'spiral2'
    f = @(x) x.^2;
    x= -1*pi*nTwists : 0.05 : pi*nTwists;
    r= 0:1/(length(x)-1):1;
    
    X=sin(x).*f(r);  Y=cos(x).*f(r);     
    % Normalize X,Y
    X = X/max(X);   Y = Y/max(Y);
    
    X = uint16(N(1)/2 *(X+1)) - 1;
    Y = uint16(N(2)/2 *(Y+1)) - 1;

    geom = zeros(N);

    for i = 1:length(X)
       x = X(i)+1;
       y = Y(i)+1;
       geom(x,y) = 1;
    end
    
% Another spiral
    case 'spiral3'
    f = @(x) x.^(1/4);
    x= -1*pi*nTwists : 0.005 : pi*nTwists;
    r= 0:1/(length(x)-1):1;
    
    X=sin(x).*f(r);  Y=cos(x).*f(r);     
    % Normalize X,Y
    X = X/max(X);   Y = Y/max(Y);
    
    X = uint16(N(1)/2 *(X+1)) - 1;
    Y = uint16(N(2)/2 *(Y+1)) - 1;

    geom = zeros(N);

    for i = 1:length(X)
       x = X(i)+1;
       y = Y(i)+1;
       geom(x,y) = 1;
    end
    % Add some center points
    geom(64, 64) = 1;
    
% Regular spiral + dense center  
    case 'spiral4'
    f = @(x) x;
    x= -1*pi*nTwists : 0.005 : pi*nTwists;
    r= 0:1/(length(x)-1):1;
    
    X=sin(x).*f(r);  Y=cos(x).*f(r);     
    % Normalize X,Y
    X = X/max(X);   Y = Y/max(Y);
    
    X = uint16(N(1)/2 *(X+1)) - 1;
    Y = uint16(N(2)/2 *(Y+1)) - 1;

    geom = zeros(N(1),N(2));

    for i = 1:length(X)
       x = X(i)+1;
       y = Y(i)+1;
       geom(x,y) = 1;
    end
    
    % Add the center
    [X, Y] = meshgrid(-1: 2/(N(2)-1) :1, -1: 2/(N(1)-1) : 1);
    unitGrid = zeros([N(1), N(2), 2]);
    unitGrid(:,:,1) = X;
    unitGrid(:,:,2) = Y;
    clear X Y;

    h = 2/(N(2)-1);
    radCenter = radCenter * h; 

    Z = unitGrid(:,:,1).^2 + unitGrid(:,:,2).^2;
    mask = (Z<radCenter^2);
    geom(mask==1) = 1;
    
% Random sampling
    case 'random1'
        geom = zeros([prod(N),1]);
        mind = randsample(prod(N),nRandSamples);
        geom(mind) = 1;
        geom = reshape(geom,N);
        fprintf('Random sampling.\n');
        
% Random sampling with dense center 
    case 'random2'
        geom = zeros([prod(N),1]);
        mind = randsample(prod(N),nRandSamples);
        geom(mind) = 1;
        geom = reshape(geom,N);
        
        % Add the center
        [X, Y] = meshgrid(-1: 2/(N(2)-1) :1, -1: 2/(N(1)-1) : 1);
        unitGrid = zeros([N(1), N(2), 2]);
        unitGrid(:,:,1) = X;
        unitGrid(:,:,2) = Y;
        clear X Y;
        
        h = 2/(N(2)-1);
        radCenter = radCenter * h; 
        
        Z = unitGrid(:,:,1).^2 + unitGrid(:,:,2).^2;
        mask = (Z<radCenter^2);
        geom(mask==1) = 1;
        
        fprintf('Random sampling.\n');
        
    otherwise
    error('Please enter valid type.')
end

geom = fftshift(geom);
numberOfSamplingPoints = numel(geom(geom~=0));
% fprintf(['Number of samples: ', num2str(numberOfSamplingPoints),'.\n']);

% Display geometry
if ( show == 1 )
    figure; imagesc(ifftshift(geom)); axis image; colormap gray;
    title('Sampling geometry');
end

