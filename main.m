clc;
clear;
cd('C:\Users\simon\Dropbox\PhD\CovarianceModelling\Replication code') %Change this
addpath(genpath('Functions'))

%Simulate data and estimate the daily covariance matrix
nObs = 50; %The number of days in the sample
N = 40; %The number of stocks
totalFactors = 10; %The number of factors in the data

%Simulation settings
targetObs = 15000;
freq = 1/2; %30 seconds sampling
blockCorr.LB = 0.03; blockCorr.UB = 0.40;
nSimulations = 1;
noiseVar = 7e-08;
nBlocks = 10; %Number of blocks in the data

tic;
%Simulate prices and generate the full covariance matrix (i.e., for individual stocks and factors) 
% The covariance matrix is in vech form here, and the inverse vech operator
% must be used to get the covariance matrix for each day
COV = simulationData(N,totalFactors,nBlocks,blockCorr,noiseVar,targetObs,nObs,nSimulations,freq);
toc;

%% Forecasting analysis

nFac = 10; %The number of factors to be used in the model. nFac = 1: SPY, nFac = 10: All SPDR ETFs 

% Calculate factor variances and corelations + the idiosyncratic covariance matrix 
rvF = zeros(nObs,nFac); %Factor variances
corrF = zeros(nObs,nFac*(nFac-1)/2); %Factor correlations
selN = ~triu(true(nFac));
beta = zeros(nObs,N*nFac); %Betas
S = cell(nObs,N); %Idiosyncratic covariance matrix
for t = 1:nObs
    if mod(t,250)==0
        disp(t);
    end
    C = ivech(COV{t});
    rvcov = C(N+(1:nFac),N+(1:nFac));
    rvF(t,:) = diag(rvcov);
    [~,tmp] = cov2corr(rvcov);
    corrF(t,:) = tmp(selN);
    b = C(1:N,N+(1:nFac))*(rvcov\eye(nFac));
    beta(t,:) = reshape(b,nFac*N,1);
    S{t} = C(1:N,1:N) - b*rvcov*b';
end

inSampleObs = nObs-1; %Number of in-sample observations
nForecast = 1; %Number of out-of-sample forecasts
horizon = 1; %Forecasting horizon


%Cluster estimation and construction of variances and correlations within
%each block

clusterMetric = 'complete'; %Complete linkage for hierachical clustering
metric = 'standard'; % Dismilliarity measure 1-abs(rho)
[nGroups, sectorClust] = deal([]); %No pre-specification of the number of clusters in the idiosyncratic covariance matrix
cutoff = 0.9999; %Threshold to define clusters, i.e., tau in the paper

lags = {[1;5;22];[1;5;22];[1;5;22];[1;5;22];[1;5;22]}; %Lags for the HAR model, (1): betas, (2) factor variances, 
% (3) factor correlations, (4) idiosyncratic variances, (5) idiosyncratic correlations

%Estimation of the block structure and of the corresponding idiosyncratic
%variances and correlations
[yFV,xFV,yCorr,xCorr,order,dim,Ng] = getClusterReady(S,N,1:inSampleObs,metric,cutoff,lags,horizon,sectorClust,clusterMetric, nGroups);

% Prepare for HAR estimation of factor variances, factor correlations, and betas 


[yBeta,xBeta] = HARMATRIX_H(beta,0,lags{1},horizon,'average');

[yFacRV, xFacRV] = deal(cell(nFac,1));
for m = 1:nFac
    [yFacRV{m},xFacRV{m}] = HARMATRIX_H(log(rvF(:,m)),1,lags{2},horizon,'average');
end

[yFacCorr,xFacCorr] = HARMATRIX_H(corrF,0,lags{3},horizon,'average');

%Forecast factor variances
sim = 10000;
index = unidrnd(inSampleObs-max(cell2mat(lags))-horizon+1,sim,nForecast);
iota = ones(1,sim);
RVF = zeros(nFac,1);
for m = 1:nFac
    y = yFacRV{m}(1:(end-horizon),:);
    x = xFacRV{m}(1:(end-horizon),:);
    xf = xFacRV{m}(end,:);
    RVF(m) = alternativeForecastingV2(y,x,xf,'average',index,iota);
end

%Forecast betas

aY = yBeta(1:(end-horizon),:);
aX = xBeta(1:(end-horizon),:,:);
aXF = xBeta(end,:,:);
mC = mean(aY);
y = reshape(aY-mC,size(aY,1)*nFac*N,1);
x = reshape(permute(aX,[1 3 2]),size(aY,1)*nFac*N,size(lags{1},1));
x = x - reshape(repmat(mC,size(aY,1),1),size(aY,1)*nFac*N,1);
if size(lags{1},1)>1
    xf = squeeze(aXF)'-mC';
else
    xf = squeeze(aXF)-mC';
end
    
BETA = reshape(xf*(x\y) + mC',N,nFac);

%Forecast factor correlations

if nFac>1
    aY = yFacCorr(1:(end-horizon),:);
    aX = xFacCorr(1:(end-horizon),:,:);
    aXF = xFacCorr(end,:,:);

    mC = mean(aY);
    y = reshape(aY-mC,size(aY,1)*(nFac*(nFac-1)/2),1);
    x = reshape(permute(aX,[1 3 2]),size(aY,1)*(nFac*(nFac-1)/2),size(lags{3},1));
    x = x - reshape(repmat(mC,size(aY,1),1),size(aY,1)*(nFac*(nFac-1)/2),1);
    if nFac*(nFac-1)/2>1
        if size(lags{3},1)>1
            xf = squeeze(aXF)'-mC';
        else
            xf = squeeze(aXF)-mC';
        end
    else
        xf = squeeze(aXF)-mC;
    end
    facCorrMat = corr_ivech(xf*(x\y) +  mC');
else
    facCorrMat = 1;
end

%Combined forecast of the factor covariance matrix
SIGMAF = facCorrMat.*sqrt(RVF*RVF');

%Forecast idiosyncratic variances

FV = zeros(N,1);
for m = 1:N
    y = yFV{m}(1:(end-horizon),:);
    x = xFV{m}(1:(end-horizon),:);
    xf = xFV{m}(end,:);
    FV(m) = alternativeForecastingV2(y,x,xf,'average',index,iota);
end

%Forecast idiosyncratic block correlations

RHO = cell(Ng,1);
for k = 1:Ng
    if dim(k)> 0
        aY = yCorr{k}(1:(end-horizon),:);
        aX = xCorr{k}(1:(end-horizon),:,:);
        aXF = xCorr{k}(end,:,:);
        mC = mean(aY);
        y = reshape(aY-mC,size(aY,1)*dim(k),1);
        x = reshape(permute(aX,[1 3 2]),size(aY,1)*dim(k),size(lags{5},1));
        x = x - reshape(repmat(mC,size(aY,1),1),size(aY,1)*dim(k),1);
        if dim(k)>1
            if size(lags{5},1)>1
                xf = squeeze(aXF)'-mC';
            else
                xf = squeeze(aXF)-mC';
            end
        else
            xf = squeeze(aXF)-mC;
        end
        tmp = corr_ivech(xf*(x\y) +  mC');
    else
        tmp = 1;
    end

    RHO{k} = tmp.*sqrt(FV(order(:,2)==k)*FV(order(:,2)==k)');
end

% systematic part
market = BETA*SIGMAF*BETA';
tmp = blkdiag(RHO{:});
    
%idiosyncratic part
[~,invOrder] = sort(order(:,1));    
idio = tmp(invOrder,invOrder);

%Forecast of the full N x N covariance matric
SIGMA = market + idio;
