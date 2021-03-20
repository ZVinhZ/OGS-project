% TP1
function [ u0,miu0,it ] = UzawaQuadra( A,b,C,d,u0,rho,itMax,tol )
%Uzawa algorithm for the quadratic case with affine inequality constraints
%   u0 : initial point
m=size(C,1);
miu0=zeros(m,1);
it=1; 
err=1;
while it<itMax && err>tol
    u1=A\(b-C'*miu0);
    err=norm(u1-u0);
    miu0=max(0,miu0+rho*(C*u1-d));
    u0=u1;
    it=it+1;
end

end

