function u3d = make_3d(u) 

% Input: nx1 (1xn) cell containing n images 
%        or 
%        dx1 cell containing d nx1 cells with n images

% Distinguish input
if (iscell(u{1}) == 0)
    n = max(size(u));
    sz_im = size(u{1,1}); 
    u3d = zeros([sz_im,n]);
    for i = 1:n
        u3d(:,:,i) = u{i};
    end

else 
    d = max(size(u));
    sz_im = size(u{1}{1});
    % Determine number of frames
    T = zeros(d,1);
    for i = 1:d
        T(i) = max(size(u{i}));
    end 
    u3d = zeros([sz_im,sum(T)]);
    
    % Fill 
    t = 1;
    for i = 1:d 
        for j = 1:T(i)
            u3d(:,:,t) = u{i}{j};
            t = t + 1;
        end
    end
    
    
end