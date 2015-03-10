function [x,succ] = newt(g,h,x,tol,al,fl)
%
% function [x,succ] = newt(g,h,x,tol,al,fl)
%
% Upon success, this function returns in x a value x_n such that
%      | x_n - x_{n-1} | < tol(1+|x_n|), and
%      |g(x_n)| < tol.
% n is the number of iterations required.
% In this case, set succ = 1.
% However, if for some iteration |g(x_k)| > 0.5 |g(x_{k-1})|
% then succ = 0 and thus "failure" is declared.

xo = al;
afo = 5*abs(fl); % times 5 because input x is not Newton's
succ = 0;

while succ == 0
  fx = feval(g,x);
  afx = abs(fx);
  % check successful termination
  if ( abs(x-xo) < tol*(1+abs(x)) ) * ( afx < tol )
    succ = 1;
    return
  % check insufficient progress
  elseif afx >= 0.5 * afo
    succ = 0;
    return
  else
    % Newton step
    xo = x;
    afo = afx;
    fxp = feval(h,x);
    x = x - fx / fxp;
  end
end
