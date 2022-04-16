%% Kernel function to compute weights for LFR when p>1. This function is to used in 'getLFRweights.m' function.

% Inputs:  1) x : an nxp matrix, n being the number of observations, p
%                 being number of covariates.
%          2) h : a vector of length p of bandwidths.


function k = K(xtemp, h)
    k=1;
    p = size(xtemp,2);
    for i=1:p
        k = k .* normpdf(xtemp(:,i)./h(i));
    end    
end
