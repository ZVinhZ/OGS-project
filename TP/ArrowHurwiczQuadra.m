% TP1
function [ u0,miu0,it ] = ArrowHurwiczQuadra( A,b,C,d,u0,ksi,rho,itMax,tol )
%u1 and u2 positif for inequality constraints
%   Detailed explanation goes here
m=size(C,1);
miu0=zeros(m,1);
it=1; 
err=1;
while it<itMax && err>tol
    u1=max(0,u0-ksi*(A*u0-b+C'*miu0));
    err=norm(u1-u0);
    
    miu0=max(0,miu0+rho*(C*u1-d));
    u0=u1;
    it=it+1;
end
end

