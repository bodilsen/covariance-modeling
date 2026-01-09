function [H] = UseH_MKern(rate, NtoS, N, Vers, topF)
%
% Copyright: Asger Lunde
% alunde@asb.dk
% Revision: 1    Date: 7/1/2010

    if rate==1/2
        if topF==0, error('Invalid rate/Flat-top combination'); end;
        switch lower(Vers)
            case 'th2',     H = ceil(5.74*(NtoS*N)^rate);
            case 'parzen',  H = ceil(4.77*(NtoS*N)^rate);
            otherwise
                error('Invalid kernel choice');
        end
    elseif rate==3/5
        if topF==1, error('Invalid rate/Flat-top combination'); end;
        switch lower(Vers)
            case 'th2',     H = ceil(4.467*NtoS^(2/5)*N^rate);
            case 'parzen',  H = ceil(3.513*NtoS^(2/5)*N^rate);
            otherwise
                error('Invalid kernel choice');
        end
    else
        error('Invalid rate choice');
    end
end
