function rtn = CorM(x)

invol = diag(sqrt(diag(x)).^(-1));
rtn = invol*x*invol;
rtn = (rtn+rtn')/2;

