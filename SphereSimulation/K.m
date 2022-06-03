%% Gaussian kernel function to compute weights for LFR when p>1. This function is to used in 'getLFRweights.m' function. Also computes the first and second derivatives if needed.

% Inputs:  1) x : an nxp matrix, n being the number of observations, p
%                 being number of covariates.
%          2) h : a vector of length p of bandwidths.
%          3) kern : a string indicating the kernel choice ('gauss' or
%          'epan') are the only options currently accepted. If missing,
%          'epan' is used


function [k, kder1, kder2] = K(xtemp, h, kern)

    if(nargin < 3)
        kern = 'epan';
    elseif(~strcmp(kern, 'gauss') && ~ strcmp(kern, 'epan'))
        error("Invalid kernel choice");
    end

    k=1; 
    if(nargout > 1)
        kder1=1; 
    end
    if(nargout > 2)
        kder2=1;
    end
    
    p = size(xtemp,2);
    for i=1:p
        
        switch kern
            
            case 'gauss'
                
                kcur = normpdf(xtemp(:, i)./h(i))/h(i);
                k = k .* kcur;
                
                if(nargout > 1)
                    kcur1 = ((-xtemp(:, i)/h(i)) .* kcur);
                    kder1 = kder1 .* kcur1 / (h(i)^2);
                end
                if(nargout > 2)
                    kder2 = kder2 .* -(kcur + (xtemp(:, i)/h(i)).*kcur1) /(h(i)^3); 
                end
                
            otherwise
                
                k = k .* 0.75 .* (1 - (xtemp(:, i)/h(i)).^2)/h(i) .* (abs(xtemp(:, i)) <= h(i));
                
                if(nargout > 1)
                    kder1 = kder1 .*(-1.5).*(xtemp(:, i)/h(i))/(h(i)^2) .* (abs(xtemp(:, i)) <= h(i));
                end
                if(nargout > 2)
                    kder2 = (-1.5)*(abs(xtemp(:, i)) <= h(i))/(h(i)^3);
                end
        end
    end    
end
