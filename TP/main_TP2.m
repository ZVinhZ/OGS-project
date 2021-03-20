clear all
clc

% resolution with uzawa
A=eye(4,4)*4;  
b=[1;-2;3;-1]; 
C=[1 1 -1 -1; 1 -1 1 -1];
d=[1;5];
u0=[0;0;0;0];
rho=0.5;
tol=10^(-4);
itMax=1000;
[u,lambda, it]=UzawaQuadraEquality(A,b,C,d,u0,rho,itMax,tol);

%%
%u0=[3; 0; 2; 0]; % pas d'influence
%u0=[3; -2; 0; 0];
u0=[1; 1;1; 1];
%v0=[3 0 -2 0; 3 0 2 0];
v0=[4 0 -3 0; 4 0 1 0]; % doit verifier CG
rho=0.01;
rho_U=1;
% rho_U=0.1; % ca ne cvg pas, it1=290 ?
tol=10^(-3);
[ u1,lambda1,it1]=decompRessource( A,b,C,u0, v0,rho,rho_U,tol,itMax );
% rho est lie au epsilon  
%%
% 2. KKT
[grad,equa]=testKKT(A,b,C,d,u1,lambda1,tol,1);
% grad est de meme precision que tol de uzawa
% equa 0.1
%%
% 3. Encadrement pour rho
rho=0.01;
grad=1;
equa=1;
while grad+equa ==2
    rho=rho+0.01;
    [ u1,lambda1,it1]=decompRessource( A,b,C,u0, rho,rho_U,tol,itMax );
    [grad,equa]=testKKT(A,b,C,d,u1,lambda1,tol,1 );
end
% rho=0.22 pour tol=0.001
%%
% 4.
rho=0.01;
list_tol=zeros(1,11);
list_it=zeros(1,11);
for k=3:13
    list_tol(k-2)=10^(-k);
    [ u1,lambda1,list_it(k-2)]=decompRessource( A,b,C,u0, rho,rho_U,list_tol(k-2),itMax );
end
plot(list_tol,list_it)
%%
% exo2
% Seq1
i0=3;
%v0=[2;2];
%v0=[-2;2]; % pas de lien pour la cvgence
v0=[1;1];
u0=[3; 0; 2; 0]; % pas de lien pour la cvgence
gamma=0.1;
[ u1,lambda1,it1 ] = decompPredSeq1(A,b,C,d,i0,u0,v0, gamma,rho_U,tol,itMax );
% eps fixe, tout le reste depend de ca
% l'algo qui fait appeal a un autre algo, error pour l'autre doit etre plus
% petit
% cvg independant du point depart (min locaux ? --> Dichonomie)
% stabilite : petite perturbation, algo marche tjr?
% consistance etc
% nb it
[grad,equa]=testKKT(A,b,C,d,u1,lambda1,tol,2 );
%%
% gamma optimal pour seq1
gamma=0.01;
K=14;
list_grad=ones(4,K);
list_equa=ones(2,K);
list_it=ones(1,K);
for i=1:K
    [ u1,lambda1,list_it(i) ] = decompPredSeq1(A,b,C,d,i0,u0,v0, gamma,rho_U,tol,itMax );
    [list_grad(:,i),list_equa(:,i)]=testKKT(A,b,C,d,u1,lambda1,tol,2 );
    gamma=gamma+0.01;
end
gamma=0.12;
% critere : min it
%%
beta=0.1;
p0=[1;1];
[ u2,lambda2,it2 ] = decompPredPara(A,b,C,d,i0,u0,v0, p0,gamma,beta,rho_U,tol,itMax );
[grad,equa]=testKKT(A,b,C,d,u2,lambda2,tol,2 );
% gain en temps de calcul, resultat le plus mauvais
%%
% beta optimal pour Para
beta=0.01;
K=17;
list_grad=ones(4,K);
list_equa=ones(2,K);
list_it=ones(1,K);
for i=1:K
    [ u2,lambda2,list_it(i) ] = decompPredPara(A,b,C,d,i0,u0,v0, p0,gamma,beta,rho_U,tol,itMax );
    [list_grad(:,i),list_equa(:,i)]=testKKT(A,b,C,d,u2,lambda2,tol,2 );
    beta=beta+0.01;
end
vec_grad=sum(abs(list_grad));
vec_equa=sum(abs(list_equa));

%%
[ u3,lambda3,it3 ] = decompPredSeq2(A,b,C,d,i0,u0, p0,beta,rho_U,tol,itMax );
[grad,equa]=testKKT(A,b,C,d,u3,lambda3,tol,2 );