function [C,x,count] = logCorr(R,threshold)

N = size(R,1);
if nargin < 2
    threshold = 1e-6;
end

x = -0.5*ones(N,1);
diagIdx = 1:N+1:N^2;

count = 0;
while true
    R(diagIdx) = x;

    % Symmetric eigendecomposition (exact for symmetric R)
    [V, lam] = eig(R, 'vector');   % lam is N×1 eigenvalues

    % diag(expm(R)) computed exactly without expm:
    d = (V.^2) * exp(lam);

    xnew = x - log(d);
    count = count + 1;

    if norm(xnew - x, inf) < threshold
        x = xnew;
        break
    end
    x = xnew;
end


% Build C = expm(R) via eig, avoiding expm:
eLam = exp(lam);
C = (V .* eLam.') * V.';   % equals V*diag(eLam)*V'
end
