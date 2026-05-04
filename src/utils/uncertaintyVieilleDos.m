function [a, n, R2, qmMu, qmSigma] = uncertaintyVieilleDos(p, rb)
% uncertaintyVieilleDos - Computes Vieille's law parameters and covariance matrix via logarithmic linear regression
%
% This function determines the parameters 'a' and 'n' for Vieille's ballistic law (rb = a * p^n)
% by performing a least-squares linear regression on the logarithmic data: log(rb) = log(a) + n * log(p).
% In addition to the nominal parameters and the goodness of fit (R^2), it calculates the mean vector
% and covariance matrix of the log-space parameters for use in Monte Carlo uncertainty simulations.
%
% SYNTAX:
%  [a, n, R2, qmMu, qmSigma] = uncertaintyVieilleDos(p, rb)
%
% INPUT:
%  p        - Array of experimental pressure values.
%  rb       - Array of experimental burn rate values corresponding to 'p'.
%
% OUTPUT:
%  a        - Nominal pre-exponential coefficient of Vieille's law.
%  n        - Pressure exponent of Vieille's law.
%  R2       - Coefficient of determination (R-squared) for the regression fit.
%  qmMu     - 1x2 mean vector [q, m] in logarithmic space, where q = log(a) and m = n.
%  qmSigma  - 2x2 covariance matrix of the log-space parameters [q, m].

% The comments were left in italian as the original function uncertaintyVielle.m


        X = log(p);                                                         % conversione valori di pressione in scala logaritmica
        N = length(p);                                                      % lunghezza vettore delle pressioni
        Y = log(rb);                                                        % conversione valori di rb in scala logaritmica

        delta = N * sum(X.^2)-(sum(X))^2;                                   % denominatore per il calcolo di q e m
        m = ( N*sum(X.*Y) - sum(X)*sum(Y) ) / delta;                        % coefficiente angolare della retta Y = mX + q con m = n
        n = m;                                                              % dal passaggio sopra
        q = ( sum(X.^2)*sum(Y) - sum(X)*sum(X.*Y)) / delta;                 % intercetta della retta Y = mX + q con m = n

        Y_eval = m .* X + q;                                                % valori valutati tramite equazione della retta appena calcolata
        sig_sq = ( sum( ( Y - Y_eval ).^2) ) / (N-2);                       % Varianza (scarto quadratico medio al quadrato)


        % Calcolo multivariabile
        Var_q  = sig_sq * ( sum( X.^2 ) / delta );                          % Varianza di q (intercetta)
        Var_m  = sig_sq * ( N / delta );                                    % Varianza di m (pendenza)
        Cov_qm = sig_sq * ( -sum( X ) / delta );                            % Covarianza tra q e m

        qmMu = [q, m];                                                      % Vettore delle medie in spazio logaritmico
        qmSigma = [Var_q, Cov_qm;                                           % Matrice di Covarianza [2x2]
                    Cov_qm, Var_m];

        a = exp(q);                                                         % calcolo del valore nominale "a"

        M_Y = mean(Y);                                                      % calcolo della media di Y
        dev_reg = sum ( (Y_eval - M_Y).^2 );                                % calcolo della devianza di regressione
        dev_tot = sum ( (Y - M_Y).^2 );                                     % calcolo della devianza totale
        R2 = dev_reg / dev_tot;                                             % calcolo del coefficiente R^2

end