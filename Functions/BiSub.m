function [RQsubs] = BiSub(Z, frac, eta)
%
%  ttime: array of time stamps of prices
%      X: array of prices
%   frac: fraction that gives the number of subsamples to apply
%    eta: the estimate of \omega^2 to use
%
    ttime  = Z(:,1)';
	X      = 100*log(Z(:,2)');
	npr    = length(X);			  		% number of prices
	nsub   = floor((npr+1)^frac);       % number of subsamples
	nsamp  = ceil((npr)/nsub)-1;        % number of prices in each subsample
	scaler = zeros(1,nsub);
	x2     = 0;

	for j=1:nsub

		xindx  = Subsample(npr,j,nsub);  % subsample times  
		xindx  = xindx(1:nsamp);         % chop off to make sure everything is same size
		Xsub   = X(xindx);			     % subsample prices
		xsub   = diff(Xsub);             % subsample returns, remove the 0

		x2sub  = (xsub.*xsub)-(2.0*eta); % squared return - bias correction
		x2     = x2 +x2sub;
			
		scaler(j) = (max(ttime)-min(ttime)+1)/(max(ttime(xindx))-min(ttime(xindx)));  % period of vol calculation. Corrects for truncation of the day
    end

	x2   = (scaler(1:length(x2))./nsub).*x2;
	x2_2 = x2(1:end-2);
	x22  = x2(3:end);
	
	RQsubs =  (length(x2)^2/length(x22))*sum(x22.*x2_2);  % subsampled BPV
end

function [ray] = Subsample(n, istart, nsub)

    ray = istart:nsub:n-1; 
end
