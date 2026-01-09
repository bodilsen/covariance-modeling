function res = alternativeForecastingV2(Y,X,XF,foreMethod,index,iota)
% rng(1);
% Compared to standard alternativeforecasting I generate indices a priori
T = length(Y);
b = X\Y;
yhat = X*b;
u = Y-yhat;
if strcmpi(foreMethod,'naive')
    res = exp(XF*b);
elseif strcmpi(foreMethod,'gaussian')
    sig2 = sum(u.^2)/(T-size(b,1));
    res = exp(XF*b+0.5*sig2);
else
    %     index = unidrnd(T,sim,1);
    %     sum(index)
    sim = size(iota,2);
    res = exp(XF*b)*Mean(exp(DemeanV2(u(index),iota,sim)),iota,sim);
    %     res = exp(XF*b)*mean(exp(demean(u(index))));
end
end

