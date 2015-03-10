% Example 8.8 : minimum norm solution

format long g

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% defined a well-conditioned 5x3 problem
A = [1,0,1; 2,3,5; 5,3,-2;3,5,4; -1,6,3];
b = [4,-2,5,-2,1]';

prob = input(' enter 1,2 or 3 : ');

if (prob == 1),    
display('----  solution of well-conditioned least squares problem ---')

x = A \ b
res = norm(A*x - b)
normx = norm(x)

elseif (prob == 2)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%define a rank deficient 5x4 problem
AA = [A,sum(A,2)];

display('----  qr solution of rank-deficient least squares problem ---')

xx = AA \ b
resres = norm(AA*xx - b)
normxx = norm(xx)

display('----  svd solution of rank-deficient least squares problem ---')

[U,S,V] = svd(AA);
z = U'*b;
sig = diag(S)
y = z(1:3) ./ sig(1:3);
y = [y;0];
xs = V*y
ress = norm(AA*xs - b)
normxs = norm(xs)

elseif (prob == 3)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%define an almost rank deficient 5x4 problem
AA = [A,sum(A,2)];  
AAA = AA + 1.e-3*randn(size(AA));

display('----  qr solution of almost-deficient least squares problem ---')

xxx = AAA \ b
resresres = norm(AAA*xxx - b)
normxxx = norm(xxx)

display('----  svd solution of almost-deficient least squares problem ---')

[U,S,V] = svd(AAA);
z = U'*b;
sigp = diag(S)
y = z(1:3) ./ sigp(1:3);
y = [y;0];
xsp = V*y
ressp = norm(AAA*xsp - b)
normxsp = norm(xsp)

end