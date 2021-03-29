%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Exercice 3.1 (Kefan S., Hiba S., Victor K. Vinh N.  %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Génération de données
N = 4;

rand ("seed", 12); % Graine pour la génération de nombres aléatoires

e = rand(N,1);       % Gain moyen de chaque action
Q = diag(rand(N,1)); % Risque encouru pour chaque action
Re = rand/2;         % Rendement minimal autorisé


%% Problème 1 : imposer un rendement minimum

% Écriture du problème quadratique :
A = Q;
b = zeros(N,1);
C_e = ones(1,N);
d_e = 1;
C_i = [-e' ; -eye(N)];
d_i = [-Re ; zeros(N,1)];


%== Résolution avec la décomposition par prix (Uzawa) :
lam0 = 0;
mu0 = zeros(N+1,1);
rho = .01;
tol = 10^-5;

[ u_k, lam_k, mu_k, k ] = Uzawa(A,b,C_e,d_e,C_i,d_i, lam0,mu0, rho,tol);
u_k
k
conditions_KKT(u_k,lam_k,mu_k, A,b,C_e,d_e,C_i,d_i, .001, true);


%== Résolution avec la décomposition par quantités :
v0  = [d_e;d_i]/N*ones(1,N);
rho = .4;
tol = 10^-6;
rho_uzawa = .1;
tol_uzawa = 10^-4;

[ U, V, Lam, Mu, k ] = decomp_ressources( A,b,C_e,d_e,C_i,d_i, v0, rho,tol, rho_uzawa,tol_uzawa );
U(:,k)
k
conditions_KKT(U(:,k),Lam(:,k),Mu(:,k), A,b,C_e,d_e,C_i,d_i, .001, true);


%== Résolution avec la décomposition par prédiction :
i0 = 1;
v0 = zeros(N+2,1);
lam0 = 0;
mu0 = zeros(N+1,1);
beta = .11;
gamma = 1;
tol = 10^-4;
rho_uzawa = .1;
tol_uzawa = 10^-4;

[ U, V, Lam, Mu, k ] = decomp_pred_para( A,b,C_e,d_e,C_i,d_i, i0, v0,lam0,mu0, beta,gamma, tol, rho_uzawa,tol_uzawa );
U(:,k)
k
conditions_KKT(U(:,k),Lam(:,k),Mu(:,k), A,b,C_e,d_e,C_i,d_i, .001, true);


