
function [Xlag]=funcLag2(X,p)

    [Traw,N]=size(X);
    lgp = length(p);
    mp  = max(p);
    Xlag=zeros(Traw,lgp);

    plc = 1;
    for ii=p
        Xlag(mp+1:Traw,plc)=X(mp+1-ii:Traw-ii,1);
        plc = plc+1;
    end
end


