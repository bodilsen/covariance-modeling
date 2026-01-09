function [yFV,xFV,yCorr,xCorr,order,dim,Ng] = getClusterReady(S,N,range,metric,cutoff,lags,horizon,sectorClust,clusterMetric, nGroups)

T = size(S,1);
if isempty(sectorClust)
    group = clusterEstimation(S,N,range,metric,clusterMetric,cutoff,nGroups);
else
    group = sectorClust;
end
    
order = sortrows(group,2);
Ngroup = tabulate(group(:,2));
Ngroup = Ngroup(:,2);
dim = Ngroup.*(Ngroup-1)./2;
Ng = size(dim,1);
corrM = cell(Ng,1);
for k = 1:Ng
    if dim(k)==0
        corrM{k} = zeros(T,1);
    else
        corrM{k} = zeros(T,dim(k));
    end
end
fv = zeros(T,N);
for t = 1:T
    for k = 1:Ng
        if dim(k)>0
            tmp = S{t}(group(:,2)==k,group(:,2)==k);
            [~,tmp] = cov2corr(tmp);
            selK = ~triu(true(Ngroup(k)));
            corrM{k}(t,:) = tmp(selK);
        else
            corrM{k}(t,:) = 1;
        end
    end
    tmp = diag(S{t});
    fv(t,:) = tmp(order(:,1));    
end


[yFV, xFV] = deal(cell(N,1));
[yCorr, xCorr] = deal(cell(Ng,1));
for m = 1:N
    [yFV{m},xFV{m}] = HARMATRIX_H(log(fv(:,m)),1,lags{4},horizon,'avearge');
end

for k = 1:Ng
    if dim(k)> 0
        [yCorr{k},xCorr{k}] = HARMATRIX_H(corrM{k},0,lags{5},horizon,'average');
    end
end