%SHRINKIMAGE is a function.
%
%   [shrinked, num] = shrinkImage(img,a,b) shrinks the values of the input 
% image 'img' to the range [a b].
%
% Input:    
%   img [matrix]              
%         scalar valued image
%
%   a   [float; optional]              
%         lower bound for projection. DEFAULT = 0
%
%   b   [float; optional]              
%         upper bound for projection. DEFAULT = 1
%
% Output:
%   img [matrix]              
%         scalar valued image with values projected to [a b]
%
%   num [float]
%         ratio of shrinked values
%
%   just run:
%       SHRINKIMAGE
%
%   to see a minimal example.


%   Author:     Julian Rasch
%   E-Mail:     julian.rasch@uni-muenster.de
%   Institute:  Westfälische Wilhelms Universität (WWU) Münster
%   Date:       07-Oct-2016

function [img, num] = shrinkImage(img,a,b)
    
    if nargin == 0
        runExample;
        return;
    end
    
    if nargin < 2
        a = 0;
        b = 1;
    end
    
    l = find(img<a);
    u = find(img>b);
    
    num = (numel(l) + numel(u)) / numel(img);
    img(l) = a; 
    img(u) = b;
    
    
    function runExample()
        %defining an example starts here
        in = sprintf('This minimal example of ''shrinkImage'' just shows you the ''shrinkImage'' help text:\n');
        disp(shrinkImage(in));
        help(mfilename);
        %defining an example ends here
    end
end