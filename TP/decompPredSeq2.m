function [ u0,lambda0,it ] = decompPredSeq2(A,b,C,d,i0,u0,p0,beta,rho_U,tol,itMax )
%UNTITLED2 Summary of this function goes here
%   i0 : imposed production
%   p0 : initial price
%   beta >0
it=1;
err=1;
n=length(b); % number of variables

u1=zeros(n,1);
C_temp=C;
C_temp(:,i0)=[];
while (err>tol && it<itMax)
    for i=1:n
        if i~=i0
            u1(i)=A(i,i)\(b(i)-p0'*C(:,i));
        end
    end
    u1_temp=u1;
    u1_temp(i0)=[];
    v0=d-C_temp*u1_temp;
    % i0
    [u1(i0),lambda0, ~]=UzawaQuadraEquality(A(i0,i0),b(i0),C(:,i0),v0,u0(i0),rho_U,itMax,0.1*tol);
    % else
    p0=(1-beta)*p0+beta*lambda0;
    it=it+1;
    err=norm(u1-u0); 
    u0=u1;
end
end

