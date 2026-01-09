function rtn = MkPrices(TimePr,sray) 

time  = TimePr(:,1);
pr    = TimePr(:,2);

npr = size(TimePr,1);
nsr = length(sray);
rtn = [sray' zeros(nsr,1)];

j=1;
for i=1:nsr
    if time(j) <= rtn(i,1)
        while time(j) < rtn(i,1) && j<npr, j=j+1; end;
        if time(j)~=rtn(i,1), j=j-1; end;
        rtn(i,2) = pr(j);
    end;
    if j==npr, break; end; 
    if time(end)< rtn(i,1), break; end; 
end;
rtn = rtn(:,2);


