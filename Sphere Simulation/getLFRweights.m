%% Compute the Local Frechet regression weights for p-variate covariate

% Inputs:  1) xin    - a nxp matrix of covariates
%          2) x0     - a px1 vector of arbitrary covariate point 
%          3) h      - a scalar bandwidth used for product kernel

function w = getLFRweights(xin, x0, h)
      n = size(xin, 1);
    aux = K((xin - repmat(x0', n, 1)), h*ones(1, p));
    mu0 = mean(aux);
    mu1 = mean(aux .* (xin - repmat(x0', n, 1)));
    mu2 = 0;
    for i=1:n
        mu2 = mu2 + aux(i) .* ((xin(i,:)' - x0) * (xin(i,:) - x0'))./n;
    end
    sL = zeros(n,1);
    
    for i= 1:n
      sL(i) = aux(i) .* (1 - mu1 * ((mu2)^(-1)) * (xin(i,:)' - x0));
    end
    
    s = sum(sL);
    w =sL./s;

end

