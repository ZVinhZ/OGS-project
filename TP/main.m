clear all
clc
% 1.
A=[5 0;0 3];
b=[-3;4];
C=[-1 2; 2 -1];
d=[1;1];
[u1,u2] = meshgrid(-4:0.1:1,-4:0.1:1);

J=5*u1.^2+3*u2.^2+3*u1-4*u2;
surf(u1,u2,J,'FaceAlpha',0.5,'EdgeColor','none');
xlabel('u1')
ylabel('u2')
zlabel('J')
%%
% 2.

rhoMax=2*min(svds(A))/9;
rho=0.1;
itMax=100;
u0=[0;0];
itVec=zeros(6,1);
tolVec=zeros(6,1);

for i=3:8
    tolVec(i-2)=10^-(i);
    [u,miu,it]=UzawaQuadra( A,b,C,d,u0,rho,itMax,tolVec(i-2) ); 
    itVec(i-2)=it;
end

plot(tolVec,itVec);
%%
% 3.
J_star=dot(A*u,u)/2-dot(b,u); % -1.993
grad=A*u-b+C'*miu; %=0
inequ=C*u-d; %<=0
cond_ineq=miu.*(C*u-d); %=0
miu; %>=0
%%
% 4. 
tol=1e-8;
rho=0.1;
while all(abs(grad)<=0.1) && all(inequ<=0.1) && all(abs(cond_ineq) <=0.1) && all(miu>=0)
    rho=rho+0.1;
    [u,miu,it]=UzawaQuadra( A,b,C,d,u0,rho,itMax,tol ); 
    grad=A*u-b+C'*miu; %=0
    inequ=C*u-d; %<=0
    cond_ineq=miu.*(C*u-d); %=0
end
% rho=1.29
%%
rho=0.1;
tol=10^-(3);
[u,miu,it]=UzawaQuadra( A,b,C,d,u0,rho,itMax,tol ); 
grad=A*u-b+C'*miu; %=0
inequ=C*u-d; %<=0
cond_ineq=miu.*(C*u-d); %=0
while all(abs(grad)<=0.1) && all(inequ<=0.1) && all(abs(cond_ineq) <=0.1) && all(miu>=0)
    rho=rho+0.1;
    [u,miu,it]=UzawaQuadra( A,b,C,d,u0,rho,itMax,tol ); 
    grad=A*u-b+C'*miu; %=0
    inequ=C*u-d; %<=0
    cond_ineq=miu.*(C*u-d); %=0
end
% rho=1.29
%%
% exo2
% 1.
n = 6;
D = sparse(1:n,1:n,2*ones(1,n),n,n);
E = sparse(2:n,1:n-1,-1*ones(1,n-1),n,n);
A = E+D+E';
b = -1*ones(n,1);
C = [3 1 0 -1 0 0; -1 2 1 0 0 0 ; 0 0 0 1 -1 1];
d = [0;1;0];
rhoMax = 2*min(svds(A))/max(svds(C));
u0=zeros(6,1);
rho = 0.1;
[u,miu,it]=UzawaQuadra( A,b,C,d,u0,rho,itMax,tol );
J_star=dot(A*u,u)/2-dot(b,u);

%%
% 2.
grad=A*u-b+C'*miu; %=0
inequ=C*u-d; %<=0
cond_ineq=miu.*(C*u-d); %=0
miu; %>=0

%%
% exo supplementaire
A=-2*eye(2,2);
b=[-12;-9];
C=[1 2 ; 1 1; -2 1];
d=[12;7;12];
u0=[1;1];
ksi=0.1;
rho=0.01;
itMax=1000;
tol=10^(-5);
% run with A H
[u,miu,it]=ArrowHurwiczQuadra( A,b,C,d,u0,ksi,rho,itMax,tol );
% it = 4, on a la vraie sol

%%
% run with uzawa
C=[1 2 ; 1 1; -2 1; -1 0; 0 -1];
d=[12;7;12;0;0];
rhoMax=4/64;
rho=0.00001;
[u,miu,it]=UzawaQuadra( A,b,C,d,u0,rho,itMax,tol );
% difficile d'obtenir la vraie sol