% Weak line search with gradient descent
% Try different initial guesses 
% Use two line search strategies
% Uses functions phi.m and funv.m

clear all
% exact solution
xe = [.695884386117764,-1.34794219305888]';
fe = phi(xe);
nmax = 2000; tol = 1.e-6; alphamin = 1.e-6; sigma = 1.e-4;

% first initial guess
 x0 = [.75,-1.25]';
% second initial guess
 x0 = [0,0.3]';
linesearch = input(' which line search? 0 = backtrachking, 1 = quad : ');

fprintf ('k      alpha(k-1)      ||x_k - x*||        f_k - f(x*)      g_kp_k  \n')

% apply grad descent
x = x0; fx = phi(x); alpha = 0; 
fc = 1; % number of calls to evaluate f = phi
for k=1:nmax  
  gx = funv(x); % the gradient at x
  p = -gx; pg = p'*gx; alphamax = 1/sqrt(-pg);
  fprintf ('%d     %e     %e     %e     %e \n',k-1,alpha,norm(x-xe),fx-fe,gx'*p)
  alpha = alphamax;
  
  if linesearch == 1,
    goon = 1;
    while (alpha > alphamin) & goon
      % quadratic interpolant minimization  
      xn = x + alpha * p;
      fxn = phi(xn);
      fc = fc + 1;
      if (fxn > fx + sigma * alpha * pg)
        mu = -0.5 * pg * alpha / (fxn - fx - alpha * pg );
        mu = max(mu,.01); % don't trust quadratic interpolation from far away
        alpha = mu * alpha;
      else
         x = xn; fx = fxn; goon = 0;
      end
    end 
    
  else
      
    goon = 1;
    while (alpha > alphamin) & goon
      % backtracking
      xn = x + alpha*p;
      fxn = phi(xn);
      fc = fc + 1;
      if (fxn > fx + sigma * alpha * pg)
        alpha = alpha/2;
      else
        x = xn; fx = fxn; goon = 0;
      end
    end
  end
  if norm(p) < tol*(1+norm(x)) break, end
end
fprintf('number of function calls upon convergence = %d \n',fc)

