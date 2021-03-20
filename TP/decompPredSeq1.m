function [ u0,lambda0,it ] = decompPredSeq1(A,b,C,d,i0,u0,v0, gamma,rho_U,tol,itMax )

%   i0 : imposed production
%   v0 : initial allocation for i0
%   gamma >0
it=1;
err=1;
n=length(b); % number of variables

u1=zeros(n,1);
C_temp=C;
C_temp(:,i0)=[];
while (err>tol && it<itMax)
    % i0
    [u1(i0),lambda0, ~]=UzawaQuadraEquality(A(i0,i0),b(i0),C(:,i0),v0,u0(i0),rho_U,itMax,0.1*tol);
    % coordi
    p1=lambda0;
    for i=1:n
        if i~=i0
            u1(i)=A(i,i)\(b(i)-p1'*C(:,i));
        end
    end
    u1_temp=u1;
    u1_temp(i0)=[];
    v0=(1-gamma)*v0+gamma*(d-C_temp*u1_temp);
    it=it+1;
    err=norm(u1-u0); 
    u0=u1;
end
end

