function COV = MultiRK(useDim,prices,times,async,freq)

%         tic
%         disp(d/days);
RV = zeros(useDim,1);  LfreqRV = RV;
% nNonZero = RV;
% 
% for m = 1:useDim
%     ret = diff(log(prices{m}));
%     %         sray30min = 1:60*30:23400;
%     %         sray30min = [sray30min 23400];
%     LfreqRet  =  diff(log(realized_price_filter(prices{m},times{m},'seconds','CalendarTime',60*30)));
% 
%     nNonZero(m,1) = sum(~ret==0);
%     RV(m,1)       = ret'*ret;
%     LfreqRV(m,1)  = LfreqRet'*LfreqRet;
% end
% 
% NtoS     = bsxfun(@rdivide,bsxfun(@rdivide,RV,2*nNonZero),LfreqRV);

useOmega2 = zeros(useDim,1);
for m = 1:useDim
    [EstOmega2,nNoise]  = subNoiseRVEst(prices{m}, 15, 0);
    LfreqRV(m)    = SubRV_ef(prices{m},10*60);
    NtoS       = EstOmega2/LfreqRV(m);
    [~, MKern] = MultivarRKernel(1/2,1,[times{m},prices{m}] ,NtoS,'parzen',1);
    useOmega2(m)  = exp(log(EstOmega2)-(MKern./10000)/(2*nNoise*EstOmega2));
end
NtoS = useOmega2./ LfreqRV;

tim = (1:60*freq:23400)';
if max(tim)<23400
    tim = [tim;23400];
end
if async
    RF = tim;
else
    RF = (1:23400)';
end
for i = 1:useDim
    if async
        fp =  realized_price_filter(prices{i},times{i},'seconds','fixed',tim);
        % mean(diff(fp)==0)
        RF = [RF,fp];
    else
        RF = [RF,prices{i}];
    end
end


[H, RK] = MultivarRKernel(3/5,1,RF,mean(NtoS),'parzen',0);

COV = vech(RK);

