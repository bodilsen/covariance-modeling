function [prPro,spt,trueVals,vEnd] = SimFactorSV1F(d, NSecs, parm, useVol)
% Simulate d-variate Factor 1-Factor Stochastic Volatility Model%
%
% USAGE:
%   [DATA] = SimFactorSV1F(DIMENSION,NSECS,PARM)
%
% INPUTS:
%   DIMENSION        - Dimension of price process
%   NSECS            - coarseness of discretization
%   PARM             - Process parameters.
%
% OUTPUTS:
%   DATA      - NSecs by d of simulated prices 
%
% COMMENTS:
%   .... 
%
% EXAMPLES:
%   ....
%
%  See also ....

	prPro = zeros(NSecs,d);
	Delta  = 1/NSecs;
	mu     = parm{1}(1);
	beta0  = parm{1}(2);
	beta1  = parm{1}(3);
	alpha  = parm{1}(4);
	levEf  = parm{1}(5);
	rho    = parm{3};
	p0 	   = 100*log(parm{4});
	spt    = zeros(NSecs,d);	
	ns 	   = NSecs;

	W  = randn(NSecs+2,1);
	B  = randn(NSecs+2,d);
	Wp = bsxfun(@times,rho,B) + bsxfun(@times,sqrt(1-rho.^2),W);

    for i=1:d
        switch useVol
        case 'ConstVol'
            spt(:,i)     = Delta;	
            prPro(:,i)   = cumsum([ p0(i); mu*Delta+zeros(NSecs-1,1)] + bsxfun(@times,sqrt(spt(1:ns,i)),Wp(1:ns,i)));
            vEnd         = Delta;
        case 'OuVol'
            [lPr, vol]   = SimOuSVF1_v01(ns, ns, parm{1}, parm{2}(1,i), p0(i), [B(:,i); Wp(:,i)]);
            spt(:,i)     = exp(2*(beta0+beta1*vol))*Delta;	
            prPro(:,i)   = lPr;
            vEnd         = vol(end);
        otherwise
            disp('Do not know that vol!');
        end
    end;

   	TrueRV = sum(diff(prPro).^2);
	TrueQr = NSecs*sum(diff(prPro).^4)/3;

 	TrueIV = ifelse(strcmp(useVol,'ConstVol'),ones(1,d),sum(spt));
 	TrueQu = ifelse(strcmp(useVol,'ConstVol'),ones(1,d),NSecs*sum(spt.^2));

   	TrueRC = sum(prod(diff(prPro),2),1);

    trueVals = [ TrueRV TrueQr TrueIV TrueQu rho.^2 TrueRC];
end
