%% Gaussian kernel function to compute weights for LFR when p>1. This function is to used in 'getLFRweights.m' function. Also computes the first and second derivatives if needed.

% Inputs:  1) x : an nxp matrix, n being the number of observations, p
%                 being number of covariates.
%          2) h : a vector of length p of bandwidths.
%          3) kern : a string indicating the kernel choice ('gauss' or
%          'epan') are the only options currently accepted. If missing,
%          'gauss' is used
%
% Outputs: 1) k: the value of the scaled kernel evaluated at xtemp
%          
%   Also, ONLY IF p = 1 AND nargin > 1,
%
%          2) kder1: the derivative of the scaled kernel evaluated at xtemp
%          3) kder2: the second derivative of the scaled kernel evaluated
%          at xtemp
% 

function [k, kder1, kder2] = K(xtemp, h, kern)

    if(nargin < 3)
        kern = 'gauss';
    elseif(~strcmp(kern, 'gauss') && ~ strcmp(kern, 'epan'))
        error("Invalid kernel choice");
    end

    p = size(xtemp,2);

    if(nargout > 1 && p > 1)
        error('This function does not compute derivatives when p > 1')
    end
    
    k=1; 
    if(nargout > 1)
        kder1=1; 
    end
    if(nargout > 2)
        kder2=1;
    end
    
    for i=1:p
        
        switch kern
            
            case 'gauss'
                
                xh = xtemp(:, i)./h(i);
                kcur = normpdf(xh)./h(i);
                k = k .* kcur;
                
                if(nargout > 1)
                    kder1 = kder1 .* (-xh.*kcur)./h(i);
                end
                if(nargout > 2)
                    kder2 = kder2 .* (xh.^2 - 1).*kcur./(h(i)^2); 
                end
                
            otherwise
                
                xh = xtemp(:, i)./h(i);
                k = k .* 0.75 .* (1 - xh.^2)./h(i) .* (abs(xh) <= 1);
                
                if(nargout > 1)
                    kder1 = kder1 .*(-1.5).* xh .* (abs(xh) <= 1) ./ (h(i)^2);
                end
                if(nargout > 2)
                    kder2 = (-1.5)*(abs(xh) <= 1)/(h(i)^3);
                end
        end
    end    
end
