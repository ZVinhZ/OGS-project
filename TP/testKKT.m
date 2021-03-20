% TP2
function [ KKT1,KKT2 ] = testKKT( A,b,C,d,u1,lambda1,tol,param )
%   param : 1 if decomp ressource, 2 if decomp prediction

if param==1 
    grad=ones(4,1);
    for i=1:4
        grad(i)=A(i,i)*u1(i)-b(i)+C(:,i)'*lambda1(:,i);
    end
    KKT1= all(abs(grad)<=10*tol);
    KKT2= all(abs(C*u1-d)<=1);
elseif param==2
    KKT1=A*u1-b+C'*lambda1;
    KKT2=C*u1-d;
end

end

