function COV = simulationData(nAssets,nFactors,nBlocks,blockCorr,noiseVar,targetObs,nDays,nSimulations,freq)
span = 23400;

for s = 1:nSimulations

    % Simulate idiosyncratic component
    dim = ones(nBlocks,1);
    alpha1 = .9+(1.5-.9).*rand();
    alpha2 = .9+(1.5-.9).*rand();
    while sum(dim) ~= nAssets
        indi = ceil(nBlocks.*betarnd(alpha1,alpha2));
        dim(indi) = dim(indi) + 1;
    end
    % disp(dim)
    ndim = dim.*(dim-1)./2;
    trueMembership = ones(nAssets,1);
    if nBlocks >= 2
        for b = 2:nBlocks
            trueMembership(sum(dim(1:b-1))+1: sum(dim(1:b))) = b;
        end
    end

    RHO = cell(nBlocks,1);
    for bl = 1:nBlocks
        b = -1;
        a = eye(dim(bl));
        if dim(bl)>1
            while (min(b)<blockCorr.LB || max(b)>blockCorr.UB) == 1
                pd = makedist('InverseGaussian',0.08,(-0.0009634*dim(bl)+0.6463)/(dim(bl)+3.461));
                t = truncate(pd,0,0.25-0.00025*dim(bl));
                rho = random(t,ndim(bl),1);
                a = logCorr(corr_ivech(rho));
                b = corr_vech(a);
            end
        end
        RHO{bl} = a;
    end

    corM = blkdiag(RHO{:});

    %Simulate factor loadings - assumed to be constant here
    rbar = nFactors-1;
    beta = [0.25+ (1.75-0.25).*rand(nAssets,1), -.55 + (0.25-(-.55)).*rand(nAssets,rbar)];

    corr_struct = -.8 + (0.8-(-.8)).*rand(nFactors*(nFactors-1)/2,1);
    factorCorr = logCorr(corr_ivech(corr_struct));

    useDim = nAssets + nFactors;

    COV = cell(nDays,1);
    async = 1;

    lambda = 2*(span/targetObs-1)*rand(1,useDim);
    for t = 1:nDays
        rng(10+t)
        [Y,Ynoise,~,~,~,~] = simulateFactorModel(beta,corM,factorCorr,noiseVar,span);
        [prices, pricesNoise,pricesAsync,timesAsync, times,OpenClose]  = deal(cell(useDim,1));

        cng = -0.3 + (0.3 + 0.3).*rand(1,useDim);
        result = cumsum(poissrnd(repmat(max(0,lambda+cng), span, 1))+1);
        for n = 1:useDim
            prices{n} = 100*Y(:,n);
            pricesNoise{n} = 100*Ynoise(:,n);
            times{n} = (1:length(Y))';
            up = pricesNoise{n};
            useInd = result(result(:,n)<=span,n);
            pricesAsync{n} = round(up(useInd),2);
            timesAsync{n} = useInd;
            OpenClose{n} = [prices{n}(1),prices{n}(end)];
        end
        COV{t} = MultiRK(useDim,pricesAsync,timesAsync,async,freq);
    end

end