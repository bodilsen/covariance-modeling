function [Y,X,XT] = HARMATRIX_H(data,constant,lags,horizon,varargin)
%lags is a vector with the chosen lag order. Typically in a HAR regression
%we have lags = [1;5;22], corresponding to the daily, weekly and monthly
%components

% Y is the dependent variable in a HAR regression
% X is the design matrix 
% XT is the design matrix to construct out-of-sample forecast

if isempty(horizon)
    horizon = 1;
end

[obs,dim] = size(data);

k = max(lags);
if  nargin<5
    Y = data(k+horizon:end,:);
else
    temp = movmean(data,[horizon-1 0]);
    Y = temp(k+horizon:end,:);
end

Xtemp = zeros(obs,size(lags,1),dim);
for i = 1:size(lags,1)
    for d = 1:dim
        Xtemp(:,i,d) = movmean(data(:,d),[lags(i)-1 0]);
    end
end

X = Xtemp(k:end-horizon,:,:);
XT = Xtemp(end,:,:);
if constant
    X   =   [ones(size(X,1),1,dim) X];
    XT  =   [ones(1,1,dim) XT];
end


end