% TP2
function [ u0,lambda0,it ] = UzawaQuadraEquality( A,b,C,d,u0,rho,itMax,tol )
%Uzawa algorithm for the quadratic case with affine equality constraints
%   u0 : initial point
m=size(C,1);
lambda0=zeros(m,1);
it=1; 
err=1;
while it<itMax && err>tol
    u1=A\(b-C'*lambda0);
    err=norm(u1-u0);
    lambda0=lambda0+rho*(C*u1-d);
    u0=u1;
    it=it+1;
end

end


