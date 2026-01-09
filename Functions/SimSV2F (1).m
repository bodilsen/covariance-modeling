function [prPro,spt,trueVals,vEnd] = SimSV2F(NSecs, parm, useVol)
% Simulate Stochastic Volatility Model%
%
% USAGE:
%   [DATA] = SimFactorSV1F(DIMENSION,NSECS,PARM)
%
% INPUTS:
%   NSECS     - coarseness of discretization
%   PARM      - Process parameters.
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

	prPro = zeros(NSecs,1);
	Delta  = 1/NSecs;
	mu     = parm{1}(1);
	beta0  = parm{1}(2);
	beta1  = parm{1}(3);
	beta2  = parm{1}(4);
	alph1  = parm{1}(5);
	alph2  = parm{1}(6);
	phi    = parm{1}(7);
	rho1   = parm{1}(8);
	rho2   = parm{1}(9);
	p0 	   = 100*log(parm{3});
	spt    = zeros(NSecs,2);	
	ns 	   = NSecs;

	W  = randn(NSecs+2,1);
	B  = randn(NSecs+2,2);
	Wp = bsxfun(@times,rho1,B(:,1)) + bsxfun(@times,rho2,B(:,2)) + bsxfun(@times,sqrt(1-rho1.^2-rho2.^2),W);

    [lPr, vol1, vol2]   = SimOuSVF2_v01(ns, ns, parm{1}, parm{2}(1:2), p0, B(:,1), B(:,2), Wp);
    spt          = exp(beta0+beta1*vol1+beta2*vol2)*Delta;	
    prPro        = lPr;
    vEnd         = [vol1(end) vol2(end)];

   	[TrueRV,TrueQr] = deal(sum(diff(prPro).^2),NSecs*sum(diff(prPro).^4)/3);

 	[TrueIV,TrueQu] = deal(sum(spt),NSecs*sum(spt.^2));

    trueVals = [ TrueRV TrueQr TrueIV TrueQu];
end
