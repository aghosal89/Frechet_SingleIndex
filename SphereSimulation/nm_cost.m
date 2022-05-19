%% Cost function for Nelder Mead optimization

% eta: ploar coordinates in radians for current iterate.
% Y  : 3x n response belongs to metric space: surface of a unit sphere.
% x  : n x p covariate, n=#of observations, p=#of covariates.
% h  : bandwidth.

function WnVal = nm_cost(eta, Y, x, h)
    
    n = size(x,1); % find the sample size
    
    theta = polar2cartGen(eta);   % compute the cartesian coordinates from polar ones then 
    ztemp = x*theta;      % compute single index
    Y_lf = get_sphere_fit_LF(Y, ztemp, h);  % frechet mean using 'trustregion'
    WnVal = mean((arrayfun(@(g) acos(Y(:, g)'*Y_lf(:, g)), 1:n)).^2); % prediction error on sphere
        
end
