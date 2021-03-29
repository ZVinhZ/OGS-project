%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Exercice 1 (Kefan S., Hiba S., Victor K. Vinh N.  %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N = 4;

A = eye(N);
b = ((-1).^(0:(N-1)))';
C = eye(N) + 2*diag(ones(N-1,1),1);
d = zeros(N,1);

%% Présentation des solutions admissibles pour N=2
U1 = -5:.05:5;
U2 = -5:.05:5;
U1_ad = 0;
U2_ad = 0;
k = 1;
for u1 = U1
  for u2 = U2
    if (u1 + 2*u2 <= 0) && (u2 <= 0)
      U1_ad(k) = u1;
      U2_ad(k) = u2;
      k = k+1;
    endif
  end
end
%{
figure(1)
plot(U1_ad,U2_ad,'.')
axis([-5,5, -5,5])
hold on 
plot([0,0],[-5,5], ':k','linewidth',.1)
plot([-5,5],[0,0], ':k','linewidth',.1)
hold off
%}
%% Résolution avec la décomposition par prix (Uzawa) :
mu0 = zeros(N,1);
rho = .01;
tol = 10^-4;

[ u_k, lam_k, mu_k, k ] = Uzawa(A,b,0,0,C,d, 0,mu0, rho,tol);
u_k
k
conditions_KKT(u_k,0,mu_k, A,b,0,0,C,d, .001, true);

%% Résolution avec la décomposition par quantités :
v0  = d/N*ones(1,N);
rho = .4;
tol = 10^-6;
rho_uzawa = .1;
tol_uzawa = 10^-4;

[ U, V, Lam, Mu, k ] = decomp_ressources( A,b,0,0,C,d, v0, rho,tol, rho_uzawa,tol_uzawa );
U(:,k)
k

% on regarde la convergence et les KKT
seuil = 10^-2;

figure(2)
for u = U(:,k)'
  plot([1,k],[u,u], '--k')
  axis([1,k, min(U(:,k))-3,max(U(:,k))+3])
  hold on
end

for i = 1:k
  if conditions_KKT(U(:,i),0,Mu(:,2,i), A,b,0,0,C,d, seuil, false)
    plot(i*ones(N,1),U(:,i), 'og','MarkerFaceColor','g')
  else
    plot(i*ones(N,1),U(:,i), 'or','MarkerFaceColor','r')
  end
  pause(.01)
end
hold off

%{
%% Résolution avec la décomposition par prédiction (parallèle) :
i0    = 1;
v0    = ones(N,1);
mu0   = ones(N,1);
beta  = .5;
gamma = .5;
tol   = 10^-5;
rho_usawa = .1;
tol_usawa = 10^-4;

[ U, V, Lam, Mu, k ] = decomp_pred_para( A,b,0,0,C,d, i0, v0,0,mu0, beta,gamma, tol, rho_usawa,tol_usawa );
U(:,k)

% on regarde la convergence et les KKT
seuil = 10^-2;

figure(3)
for u = U(:,k)'
  plot([1,k],[u,u], '--k')
  axis([1,k, min(U(:,k))-3,max(U(:,k))+3])
  hold on
end

for i = 1:k
  if conditions_KKT(U(:,i),0,Mu(:,i), A,b,0,0,C,d, seuil, false)
    plot(i*ones(N,1),U(:,i), 'og','MarkerFaceColor','g')
  else
    plot(i*ones(N,1),U(:,i), 'or','MarkerFaceColor','r')
  end
  pause(.1)
end
hold off
%}

