function [a, n, R2, mu_qm, Sigma_qm] = uncertaintyVieilleDos(p, rb)
% TODO MODIFICARE COMMENTI PER NUOVA FUNZIONE



%
% Original file Incertezze.m V 1.02
% Modified to output Covariance Matrix for Monte Carlo simulations
%
% syntax: [a, n, R2, mu_qm, Sigma_qm] = uncertaintyVieille(p, rb)
%
% rb = a * p^n -> log(rb) = log(a*p^n) = log(a) + n * log(p)
% Y = log(rb); X = log(p); q = log(a); m = n
% Y = m*X + q
%
X = log(p);                                                         % conversione valori di pressione in scala logaritmica
N = length(p);                                                      % lunghezza vettore delle pressioni
Y = log(rb);                                                        % conversione valori di rb in scala logaritmica

delta = N * sum(X.^2)-(sum(X))^2;                                   % denominatore per il calcolo di q e m
m = ( N*sum(X.*Y) - sum(X)*sum(Y) ) / delta;                        % coefficiente angolare della retta Y = mX + q con m = n
n = m;                                                              % dal passaggio sopra
q = ( sum(X.^2)*sum(Y) - sum(X)*sum(X.*Y)) / delta;                 % intercetta della retta Y = mX + q con m = n

Y_eval = m .* X + q;                                                % valori valutati tramite equazione della retta appena calcolata
sig_sq = ( sum( ( Y - Y_eval ).^2) ) / (N-2);                       % Varianza (scarto quadratico medio al quadrato)


% Multivariate state calc
Var_q  = sig_sq * ( sum( X.^2 ) / delta );                          % Varianza di q (intercetta)
Var_m  = sig_sq * ( N / delta );                                    % Varianza di m (pendenza)
Cov_qm = sig_sq * ( -sum( X ) / delta );                            % Covarianza tra q e m

mu_qm = [q, m];                                                     % Vettore delle medie in spazio logaritmico
Sigma_qm = [Var_q, Cov_qm;                                          % Matrice di Covarianza [2x2]
            Cov_qm, Var_m];

a = exp(q);                                                         % calcolo del valore nominale "a"

M_Y = mean(Y);                                                      % calcolo della media di Y
dev_reg = sum ( (Y_eval - M_Y).^2 );                                % calcolo della devianza di regressione
dev_tot = sum ( (Y - M_Y).^2 );                                     % calcolo della devianza totale
R2 = dev_reg / dev_tot;                                             % calcolo del coefficiente R^2

end