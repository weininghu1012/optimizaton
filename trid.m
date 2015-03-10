function x = trid (md,ld,ud,b)
%
%  Solve Ax = b for a tridiagonal A
%  md = diagonal vecotor of A (length n)
%  ld = lower diagonal of A (length (n-1)
%  ud = upper diagonal of A (length (n-1)
%  b = right hand side vector

% LU decomposition (overwriting entries)
n = length(b);
for k=1:n-1
    ld(k)=ld(k)/md(k);            % compute multiplier
    md(k+1)=md(k+1)-ld(k)*ud(k);  % update pivot
    b(k+1)=b(k+1)-ld(k)*b(k);     % update r-h-s
end

% backward substitution
x=b;
x(n)=x(n)/md(n);
for k=n-1:-1:1
    x(k)=(x(k)-ud(k)*x(k+1))/md(k);
end