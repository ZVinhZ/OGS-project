% TP2
function [ u0,lambda0,it ] = decompRessource( A,b,C,u0, v0,rho,rho_U,tol,itMax )
%decomposition by ressource algorithm
%   d is the second member of ctr
%   v0 is the initial allocation
it=1;
err=1;
n=length(b); % number of variables
m=size(C,1); % number of constraints

lambda0=zeros(m,n);
u1=zeros(n,1);
while (err>tol && it<itMax)
    for i=1:n
        [u1(i),lambda0(:,i), ~]=UzawaQuadraEquality(A(i,i),b(i),C(:,i),v0(:,i),u0(i),rho_U,itMax,0.1*tol);
    end
    % coordination
    v0=v0+rho*(lambda0-repmat(mean(lambda0,2),1,n));
    it=it+1;
    err=norm(u1-u0); 
    u0=u1;
end

end
