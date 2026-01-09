function [SubRv] = SubRV_ef(Z,useF)

[nobs, rv]= deal(length(Z),0);

for i=1:useF
    if i==1
        sray   = unique([1 i:useF:nobs nobs]);
    else
       sray   = [1 i:useF:nobs nobs];
    end
    rets   = diff(1*log(Z(sray)));
    rv = rv + sum(rets.^2);
end
SubRv = rv/useF;
end
