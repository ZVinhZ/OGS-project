function [conditions_verfiees] = conditions_KKT(u,lam,mu, A,b,C_e,d_e,C_i,d_i, seuil, affichage)
%% Fonction qui permet de vérifier si les conditions KKT sont satisfaites à un certain seuil
%% pour un problème quadratique de type :
%%   minimiser  J(u) := 1/2<A*u,u> - <b,u>
%%      sous     C_e*u == d_e        (e pour contraintes d'égalité)
%%               C_i*u <= d_i        (i pour contraintes d'inégalité)

conditions_verfiees = all(abs(A*u - b + C_e'*lam + C_i'*mu) < seuil) ...
                      && all(abs(C_e*u - d_e) < seuil) ...                                                      % pour les contraintes d'égalité
                      && all(C_i*u - d_i < seuil) && all(mu > -seuil) && all(abs(mu .* (C_i*u - d_i)) < seuil); % pour les contraintes d'inégalité

if affichage
  if conditions_verfiees
    disp(['Les conditions KKT sont vérifiées au seuil d''erreur de ',num2str(seuil)])
  else
    disp(['!!! ATTENTION !!! Les conditions KKT NE sont PAS vérifiées au seuil d''erreur de ',num2str(seuil)])    
  end
end
  
end