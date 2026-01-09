function [group, hej, Z, leafOrder, dissimilarity] = clusterEstimation(S,N,range,metric,clusterMetric,cutoff,nGroups)
sel = ~triu(true(N));

if strcmp(metric,'standard')
    f = zeros(N);
    hej = zeros(N);
    for t = range
        f = f + S{t}./length(range);
        [~,corrMat] = cov2corr(S{t});
        hej = hej + (abs(corrMat)>=0.15)/length(range);
    end
    [~,corrMat] = cov2corr(f);
    dissimilarity  = 1-abs((corrMat(sel)))';
else
    if strcmp(metric,'cos')
        dist = @(x) (0.5*pi - abs(0.5*pi-acos(x)));
    elseif strcmp(metric,'sqrtabs')
        dist = @(x) sqrt(1-abs(x));
    elseif strcmp(metric,'sqrtDist')
        dist = @(x) sqrt(0.5.*(1-x));
    else
        error("Wrong input")
    end
    f = zeros(N);
    hej = zeros(N);
    sel = ~triu(true(N));
    for t = range
        [~,corrMat] = cov2corr(S{t});
        f = f + dist(corrMat)./length(range);
        hej = hej + (abs(corrMat)>=0.2)/length(range);
    end
    dissimilarity  = f(sel)';
end

Z = linkage(dissimilarity,clusterMetric);
leafOrder = optimalleaforder(Z, dissimilarity);
group = [(1:N)', cluster(Z,'cutoff',cutoff,'criterion','distance')];
if ~isempty(nGroups)
    group = [(1:N)', cluster(Z,'maxclust',nGroups,'criterion','distance')];
end

