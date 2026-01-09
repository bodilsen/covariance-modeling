function [EstOmega2,nNoise] = subNoiseRVEst(Z,useF,All_or_Chg)

[nobs, w_est, w_N]= deal(length(Z),0,0);

for i=1:useF
    if i==1
        sray   = unique([1 i:useF:nobs nobs]);
    else
        sray   = [1 i:useF:nobs nobs];
    end
    rets   = diff(1*log(Z(sray)));
    
    N     = ifelse(All_or_Chg,length(rets),sum(1-(rets==0)));
    
    w_est = w_est + sum(rets.^2)/(2*N);
    w_N   = w_N   + N;
end
[EstOmega2,nNoise] = deal(w_est/useF, w_N/useF);
end
