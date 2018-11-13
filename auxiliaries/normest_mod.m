function [e,cnt] = normest_mod(S,St,sz,tol)
%NORMEST Estimate the matrix 2-norm.
%
%   MODIFICATION 
%   Uses a function handle instead!!
%
%
% 
%   NORMEST(S) is an estimate of the 2-norm of the matrix S.
%   NORMEST(S,tol) uses relative error tol instead of 1.e-6.
%   [nrm,cnt] = NORMEST(..) also gives the number of iterations used.
%
%   This function is intended primarily for sparse matrices,
%   although it works correctly and may be useful for large, full
%   matrices as well.  Use NORMEST when your problem is large
%   enough that NORM takes too long to compute and an approximate
%   norm is acceptable.
%
%   Class support for input S:
%      float: double, single
%
%   See also NORM, COND, RCOND, CONDEST.

%   Copyright 1984-2014 The MathWorks, Inc. 

if nargin < 4
    tol = 1.e-6;
end
maxiter = 200; % set max number of iterations to avoid infinite loop 
% x = full(sum(abs(S),1))';
cnt = 0;
% e = norm(x);
% if e == 0
%     return;
% end
% x = x/e;

x = randn(sz);
e = norm(x); 
x = x/e;

e0 = 0;
while abs(e-e0) > tol*e
   e0 = e;
   Sx = S(x);
   x = St(Sx);
   normx = norm(x(:));
   e = normx/norm(Sx(:));
   x = x/normx;
   cnt = cnt+1;
   if cnt > maxiter
      warning(message('MATLAB:normest:notconverge', maxiter, sprintf('%g',tol)));
      break;
   end
end
