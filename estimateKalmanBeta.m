function [QR, distAlpha, kfBeta, Pt,Yhat, Alpha, kfAlpha, Lik] = estimateKalmanBeta(QR0, varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
%
%               ETF STATISTICAL ARBITARGE STRATEGY
%                           VERSION 2.2.2
%                   JONATHAN KINLAY JAN 15, 2015
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
% USAGE:
%
% [QR, distAlpha, kfBeta,Yhat, Alpha, kfAlpha, Lik]] =
%            estimateKalmanBeta_ver_2_1_1(Y, X, QR0)
%
% SYSTEM DYNAMICS:
%
% The system evolves according to the following difference equations,
% where quantities are further defined below:
%
% b(t) = b(t-1) + w     Beta, the unobserved state variable
%                       follows a random walk
% Y(t) = b(t)X(t) + v   The observed processes of stock prices
%                        Y(t) and X(t)
% where w ~ N(0,Q) meaning w is gaussian noise with variance Q
%       v ~ N(0,R) meaning v is gaussian noise with variance R
%
% The signal noise ratio Q/R is critical, generally around 1E-6/1E-7
%
%
% =========================================================================
%INPUTS:
%
% QR0(1x2)                  - initial estimates of Q and R
%
% Prices (nobs x 2)          - Matrix of prices
%  OR
% [Y(nobs,1) X(nobs,2)]       - two, single column vectors of prices
%
% =========================================================================
%OUTPUTS:
%
% QR(1)                     - updated estimate of variance Q
% QR(2)                     - updated estimate of variance R
%
% distaAlpha(1)             - estimated mean of Alpha
% distAlpha (2)             - estimated std of Alpha
%
% kfBeta(nobs,1)            - vector of beta(t) estimates
%
% Pt                        - see eqn below
%
% Yhat(nobs x 1)            - vector of estimated Y(t) prices
%
% Alpha(nobs x 1)           - vector of Alpha(t)=Y(t) - Yhat(t)
%
% kfAlpha(nobs x 1)         - vector of standardized alpha(t)
%
% Lik (scalar)              - estimated likelihood function value

    if nargin < 3
        Prices = cell2mat(varargin(1));
        Y = Prices(:, 1);
        X = Prices(:, 2);
    else
        Y = cell2mat(varargin(1));
        X = cell2mat(varargin(2));
    end

    cleanPrices = [Y, X];
    cleanPrices = cleanPrices(~any(isnan(cleanPrices), 2), :);
    nobs = length(Y);
    kfBeta = ones(nobs, 1);
    Yhat = NaN(nobs, 1);

    if ~isnan(Y(1)) && ~isnan(X(1))
        kfBeta(1) = Y(1) / X(1);
    end
    Yhat(1) = Y(1);
    Pt = 0;

%%%%%%%%%%%%%%%%%  GA OPTIMIZATION %%%%%%%%%%%%%%
%     kfBetaLik = @(QR)kfBetaLikelihood(Y, X, QR);
%     options = gaoptimset('Generations', 60, 'TolFun', 1E-10);
%     [QR, Lik] = ga(kfBetaLik, 2, [], [], [], [], [0, 0], [0.1, 0.1], [], options);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%  PARTICLE SWARM OPTIMIZATION %%%%%%%%%%%

    kfBetaLik = @(QR) kfBetaLikelihood(cleanPrices(:, 1), cleanPrices(:, 2), QR);
    options = optimoptions('particleswarm', 'SwarmSize', 100, 'HybridFcn', @fmincon, 'TolFun', 1E-10, 'Display', 'off');
    [QR, Lik] = particleswarm(kfBetaLik, 2, [0, 0], [1, 1], options);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    Q = QR(1);
    R = QR(2);
    
%%%%%%%%%%%%  APPLY THE KALMAN FILTER %%%%%%%%%%%%
    for t = 2 : nobs
        if isnan(X(t)) || isnan(Y(t))
            kfBeta(t) = kfBeta(t-1);
        else
            Yhat(t) = kfBeta(t - 1) * X(t);
            vt = Y(t)- Yhat(t);
            Ft = Pt * X(t) ^ 2 + R;
            kfBeta(t) = kfBeta(t-1) + Pt * X(t) * vt / Ft;
            Pt = Pt - (Pt*X(t))^2 / Ft + Q;
        end
    end

%%%% Calculate Estimated Prices and Alphas %%%%%%
    Alpha = Y - Yhat;
    distAlpha(1) = nanmean(Alpha);
    distAlpha(2) = nanstd(Alpha);
    kfAlpha = (Alpha - distAlpha(1)) / distAlpha(2);
end

function Lik = kfBetaLikelihood(Y, X, QR)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
%
%               ETF STATISTICAL ARBITARGE STRATEGY
%                           VERSION 2.1.2
%                   JONATHAN KINLAY JAN 12, 2015
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
% ESTIMATES THE LIKELIHOOD FOR THE KALMAN BETA MODEL
%
%
% SYSTEM DYNAMICS:
%
% The system evolves according to the following difference equations,
% where quantities are further defined below:
%
% b(t) = b(t-1) + w     Beta, the unobserved state variable
%                       follows a random walk
% Y(t) = b(t)X(t) + v   The observed processes of stock prices
%                        Y(t) and X(t)
% where w ~ N(0,Q) meaning w is gaussian noise with variance Q
%       v ~ N(0,R) meaning v is gaussian noise with variance R
%
% The signal noise ratio Q/R is critical, generally around 1E-7

    nobs = length(Y);
    Q = QR(1);
    R = QR(2);
    Bt = Y(1) / X(1);
    Pt = 0;
    Lik = -0.5 * nobs * log(2 * pi);

    for t = 2 : nobs
        vt = Y(t) - Bt * X(t);
        Ft = Pt * X(t) ^ 2 + R ;
        Bt = Bt + Pt * X(t) * vt / Ft;
        Pt = Pt - (Pt * X(t)) ^ 2 / Ft + Q;
        Lik = Lik - 0.5 * log(Ft) - 0.5 * vt ^ 2 / Ft ;
    end
    
    Lik = -Lik;
    
end