function [Y,Ynoise,P,F,Pnoise,Fnoise] = simulateFactorModel(beta,corM,factorCorr,noiseVar,span)

N = span;
delta = 1/N;


[d,rmax] = size(beta);

% Simulate factor model: dP = beta*dF + dZ

%First, simulate rmax factor log-prices F

B = [zeros(1,rmax);cumsum(sqrt(delta).*randn(N,rmax))]; %Standard rmax-dimensional Brownian Motion

U = chol(factorCorr,'lower');
W = (U*B')'; %Correlated r-dimensional Brownian Motion

% Simulate stochastic factor volatilities using CIR processes. Parameters
% adapted from Ait-Sahalia et al. (2020)

kappa = 1;
theta = 0.0001;
eta = .01;



% kappa = 2;
% theta = 0.001;
% eta = .045;

feller = 2*kappa*theta >= eta^2;
if ~feller
    das
end
sigma2 = zeros(N,rmax);
sigma0 = theta;
for t = 1:N
    if t == 1
        sigma2(t,:) = sigma0 + kappa*(theta-sigma0).*delta + eta*sqrt(sigma0.*delta).*randn(1,rmax);
    else
        sigma2(t,:) = sigma2(t-1,:)  + kappa*(theta-sigma2(t-1,:)).*delta + eta*sqrt(sigma2(t-1,:).*delta).*randn(1,rmax);
    end
end
% plot(sigma2)

% var(sigma2(end,:))
% theta*exp(-kappa) + theta*(1-exp(-kappa))
% (theta*(eta^2/kappa))*(exp(-kappa)-exp(-2*kappa))+((theta*eta^2)/(2*kappa))*((1-exp(-kappa))^2)

df = sqrt(sigma2).*diff(W);
f = cumsum(df);

U = chol(corM,'lower');

B = [zeros(1,d);cumsum(sqrt(delta).*randn(N,d))]; %Standard d-dimensional Brownian Motion

W = (U*B')'; %Block-Correlated d-dimensional Brownian Motion


kappa = 1;
theta = 0.00015;
eta = .01;

% kappa = 2;
% theta = 0.001;
% eta = .045;


feller = 2*kappa*theta >= eta^2;
if ~feller
    das
end
sigma2 = zeros(N,d);
sigma0 = theta;
for t = 1:N
    if t == 1
        sigma2(t,:) = sigma0  + kappa*(theta-sigma0).*delta + eta*sqrt(sigma0.*delta).*randn(1,d);
    else
        sigma2(t,:) = sigma2(t-1,:)  + kappa*(theta-sigma2(t-1,:)).*delta + eta*sqrt(sigma2(t-1,:).*delta).*randn(1,d);
    end
end

dz = sqrt(sigma2).*diff(W);
z = cumsum(dz);

p = f*beta' + z;

P = exp(p);
F = exp(f);

Y = [P,F];

% nvar = 0.001*sqrt(sum(diff(p).^4))
% sqrt(nvar)

Pnoise = exp(p + sqrt(noiseVar).*randn(N,d));
Fnoise = exp(f + sqrt(noiseVar).*randn(N,rmax));
Ynoise = [Pnoise, Fnoise];

% blockInfo.dim = dim;
% blockInfo.rho = rho;


