%% Get sphere fits for simulations for general p covariates

function lfr_fit = get_sphere_fit_pcov(Y, x, x0, h)
    
    n = size(x,1); 
    lfr_fit = zeros(3, 1);
    M = spherefactory(3);
    ops.verbosity = 0;
    fr.M = M;
    
    w = getLFRweights(x, x0, h);
    
    % compute the nadaraya watson smoothed mean
    y0 = sum(cell2mat(arrayfun(@(k) w(k)*Y(:, k), 1:n, 'UniformOutput', false))')'; 
    % compute the NW smoother with norm = 1, so that it is on the unit sphere surface 
    y0 = y0./norm(y0); % use it as initial guess for trustregion algorithm below~

    %Compute cost and Euclidean gradient

    if(length(find(w)) < 2) % if number of nonzero elements in the array in < 2, then return the j-th column 
        % as NaN values
        
      lfr_fit = repmat(NaN, 1, 3);
    
    else
        
      % run the main algorithm to obtain frechet mean ~
      fr.cost = @(y) get_cost(w, Y, y, M);
      fr.egrad = @(y) get_egrad(w, Y, y, M);
      fr.ehess = @(y, u) get_ehess(w, Y, y, M, u);

      lfr_fit = trustregions(fr, y0, ops);

      %[x, fgradx, cost, info, options]= trustregions(fr, y0, ops);
      
    end

end