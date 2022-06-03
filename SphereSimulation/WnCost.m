%% Cost function and (approximate) gradient and Hessian for estimating the index parameter in the FSI model

% eta: ploar coordinates in radians for current iterate.
% Y  : 3x n response belongs to metric space: surface of a unit sphere.
% x  : n x p covariate, n=#of observations, p=#of covariates.
% h  : bandwidth.
% useAltCost : Binary - if true, will use normalized local linear estimates
%              instead of local Frechet in computing the cost function 
%              (this is much faster, especially for large sample sizes), 
%              but can be less accurate



function [WnVal, WnGrad, WnHess] = WnCost(eta, Y, x, h, useAltCost)
    
    if(nargin < 4 || isempty(useAltCost))
        
        useAltCost = 0; % Use original cost function as default
        
    end

    n = size(x,1); % find the sample size
    
    theta = polar2cart(eta);   % compute the cartesian coordinates from polar ones
    ztemp = x*theta;      % project covariates
    
    if(~useAltCost)
                
        Y_lf = get_sphere_fit_LF(Y, ztemp, h);  % frechet mean using 'trustregion'
        WnVal = mean((arrayfun(@(g) acos(Y(:, g)'*Y_lf(:, g)), 1:n)).^2); % prediction error on sphere
    
        if(nargout == 2)
            
            [~, WnGrad] = getAltCost(eta, Y, x, h);
            
        end
        
        if(nargout == 3)
            
            [~, WnGrad, WnHess] = getAltCost(eta, Y, x, h);
            
        end
        
    else
        
        if(nargout == 1)
            
            WnVal = getAltCost(eta, Y, x, h);
            
        end
        
        if(nargout == 2)
            
            [WnVal, WnGrad] = getAltCost(eta, Y, x, h);
            
        end
        
        if(nargout == 3)
            
            [WnVal, WnGrad, WnHess] = getAltCost(eta, Y, x, h);
            
        end
        
    end

end

