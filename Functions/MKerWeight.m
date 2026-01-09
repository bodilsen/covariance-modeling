function [v] = MKerWeight(H, Vers, topF)
% Computes kernel weights for Multivariate Jittered Realized Kernel
%
% USAGE:
%
% INPUTS:
%	-H       (1x1) number of autocovariances to be weighted
%	-KVers	 (string) which kernel to use
%				"cubic","TH2","TH16","parzen","5th","6th","7th","8th","bartlett"
% OUTPUTS:
%   -v       (1xH+1) kernel weigths including gamma_0
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

    if topF>1, error('Only Flat-top 0 or 1 allowed'); end;
    v  = [];
	
	if (topF<1)
		if (H>0)
			hs = [1:1:H]/(H+1);
            switch lower(Vers)
                case 'th2',         v  = (1-cos(pi.*(1-hs).^2))/2;
                case 'parzen',      v  = 1-6*(hs.^2)+6*(hs.^3);
										if (H>2) 
                                            v(floor(H/2)+1:end)=2*(1-hs(floor(H/2)+1:end)).^3;	
                                        end;
                otherwise
                    error('Invalid kernel choice');
            end;
        end;
		v = [0.5 v];
    else
		if H>topF
			hs = [1:1:H-topF]/(H-topF+1);
            switch lower(Vers)
                case 'th2',         v  = (1-cos(pi.*(1-hs).^2))/2;
    			case 'parzen',      v  = 1-6*(hs.^2)+6*(hs.^3);
										if (H-topF>3) 
                                            v(floor((H-topF)/2)+1:end)=2*(1-hs(floor((H-topF)/2)+1:end)).^3;
                                        end;
                otherwise
                    error('Invalid kernel choice');
            end;
        end;
		if (H==0)        v = 0.5;
        elseif (H<=topF) v = [0.5 ones(1,H-1)];
        else             v = [0.5 ones(1,topF) v];
        end;
    end;
end
