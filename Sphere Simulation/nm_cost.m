%% Cost function for Nelder Mead optimization

% x  : n x p covariate, n=#of observations, p=#of covariates.
% Y  : response belongs to metric space: surface of a unit sphere.
% h  : bandwidth.
% eta: ploar coordinates in radian.

function dist_pe = nm_cost(Y, eta, h, x)
    
    n = size(x,1); % find the sample size
    p = size(x,2); % find out dimension of covariates 
    
    theta = polar2cart(eta, 1);   % compute the cartesian coordinates from polar ones then 
    ztemp = x*theta/sqrt(p);      % compute single index
    Y_fr = get_sphere_fit_FSI(Y, ztemp, h);  % frechet mean using 'trustregion'
    dist_pe = mean((arrayfun(@(g) acos(Y(:, g)'*Y_fr(:,g)),1:n)).^2); % prediction error on sphere
        
end
