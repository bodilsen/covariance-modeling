function [H MKer] = MultivarRKernel(rate, mJit, TimePrices, NtoS, KVers, topF)
% Computes kernel weights for Multivariate Jittered Realized Kernel
%
% USAGE:
%
% INPUTS:
%   -rate        (1x1) the rate power for n^(rate)
%   -mJit        (1x1) #observatioins to jitter in the begining and end of sample
%   -TimePrices  (nx(d+1)) time stamps followed by d-dimensional matrix of Prices
%   -NtoS        (1x1) Noise to signal ratio omega^2/sqrt(Quaticity)
%   -KVers	     (string) which kernel to use
%				  "cubic","TH2","TH16","parzen","5th","6th","7th","8th","bartlett"
% OUTPUTS:
%   -H            bandwidth
%   -MKer         Jittered Realized Kernel
%
% COMMENTS:
%   ....
%
% EXAMPLES:
%   ....
%
%  See also ....
%
% Copyright: Asger Lunde
% alunde@asb.dk
% Revision: 1    Date: 7/1/2010

	Mker=0;
	tm = TimePrices(:,1);
	Pr = TimePrices(:,2:end);
	
	timrg = max(tm)-min(tm); % length of periods in seconds
	dimPr = size(Pr,2);      % get the number of price series

	% Jitter first m and last m observations

    jitPr = [ mean(Pr(1:mJit,:),1); Pr(mJit+1:end-mJit,:); mean(Pr(end-mJit+1:end,:),1)];

	rets  = 100*diff(log(jitPr));
    nRt   = size(rets,1);

	% Get scaling factor
	if (mJit>1)
		jittm  = [ tm(1:mJit), fliplr(tm(end-mJit+1:end))];
		sc     = timrg/(timrg-[mJit-1:-1:1]/mJit*sum(diff(jittm),2));
    else
        sc = 1;
    end;

	% Get H and kernel weights 
	H   = UseH_MKern(rate,NtoS,size(rets,1),KVers,topF);
	H   = min([H nRt-1]);
    k_h = MKerWeight(H,KVers,topF);

    %fprintf('Bandwidth is %6i \n', H);
    %disp(k_h);
	% Get gamma_0,....,gamma_H, and compute Mkernel
	for h=1:H+1
        if h>nRt || H==1
    		Mker = rets'*rets;
        else
            gamma = rets'*[ rets(h:nRt,:); zeros(h-1,dimPr)];
            Mker = Mker + k_h(h)*(gamma+gamma');
        end;
    end;

	% return average of subsampled kernels
	MKer = Mker*sc;
end