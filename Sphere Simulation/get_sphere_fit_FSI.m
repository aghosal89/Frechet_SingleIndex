%% Get sphere fits for simulations

% Inputs: Y = response data 3xn 
%         z = nx1 index vector
%         h = bandwidth

% Output: fr_fit = a 3xn matrix, each row is the fitted value for i-th
%         observation, fitted using an iterative algorithm, is the Frechet Fit. 

function fr_fit = get_sphere_fit_FSI(Y, z, h)

    n = length(z); 
    fr_fit = zeros(3, length(z));
    M = spherefactory(3);
    ops.verbosity = 0;
    fr.M = M; 
    
    % Get weights
    zdiff = repmat(z, 1, length(z)) - repmat(z', length(z), 1); % creating the matrix of pairwise differences 
    %Kmat = 0.75*(1-(zdiff/h).^2).*(abs(zdiff)<=h);
    Kmat = ((2*pi)^(-.5))*exp((-0.5)*(zdiff/h).^2);  % getting Gaussian kernel 
    % computing the parameters for computing Nadaraya-Watson smoother
    mu0 = mean(Kmat);  
    mu1 = mean(Kmat.*zdiff); mu2 = mean(Kmat.*(zdiff.^2)); 
    sig2 = mu0.*mu2 - mu1.^2; 
    
    % computing the weight matrix
    w = Kmat.*(repmat(mu2, length(z), 1) - repmat(mu1, length(z), 1).*zdiff)./repmat(sig2, length(z), 1);
    w = w - diag(diag(w),0); % main diagonal is zero
    
  for j = 1:n
    
    % compute the nadaraya watson smoothed mean
    y0 = sum(cell2mat(arrayfun(@(k) w(k, j)*Y(:, k), 1:n, 'UniformOutput', false))')'; 
    % compute the NW smoother with norm = 1, so that it is on the unit sphere surface 
    y0 = y0/norm(y0); % use it as initial guess for trustregion algorithm below~

    %Compute cost and Euclidean gradient

    if(length(find(w(:, j))) < 2) % if number of nonzero elements in the array in < 2, then return the j-th column 
        % as NaN values
        
      fr_fit(:, jj) = repmat(NaN, 1, 3);
    
    else
        
      % run the main algorithm to obtain frechet mean ~
      fr.cost = @(y) get_cost(w(:, j), Y, y, M);
      fr.egrad = @(y) get_egrad(w(:, j), Y, y, M);
      fr.ehess = @(y, u) get_ehess(w(:, j), Y, y, M, u);

      fr_fit(:, j) = trustregions(fr, y0, ops);

      %[x, fgradx, cost, info, options]= trustregions(fr, y0, ops);
      
    end
    
  end

end

