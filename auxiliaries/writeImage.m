%WRITEIMAGE is a function.
%
%
%   just run:
%       WRITEIMAGE
%
%   to see a minimal example.


%   Author:     Julian Rasch
%   E-Mail:     julian.rasch@uni-muenster.de
%   Institute:  Westfälische Wilhelms Universität (WWU) Münster
%   Date:       19-Sep-2016

function writeImage(filename,img, cmap)
    
    if nargin == 0
        runExample;
        return;
    end
    
    if nargin < 3
        cmap = gray(256);
    end
    
    % img = (img + min(0,min(img(:))) ) / (max(img(:)) + min(0,min(img(:))) );
    
    imwrite(255*img, cmap, [filename '.png']) 
    
    function runExample()
        %defining an example starts here
        in = sprintf('This minimal example of ''writeImage'' just show you the ''writeImage'' help text:\n');
        disp(writeImage(in));
        help(mfilename);
        %defining an example ends here
    end
end