%% Gets conditional Frechet mean estimates using local Frechet regression for spherical response data, generalized for vector predictors

% This implementation uses a product Gaussian kernel with a single bandwidth for
% simplicity
%
% Inputs: Y = response data 3xn 
%         x = nxp vector of predictors corresponding to observations in Y
%         h = 1xH bandwidth vector
%         xout = mxp vector of predictors to be used for computing fitted
%                values (if missing, set xout = x so m = n)

% Output: lfr_fit = 3xmxH array of fitted values using local Frechet
%                  regression for the predictor values in xout and
%                  bandwidths in h

function lfr_fit = get_sphere_fit_LFpcov(Y, x, h, xout)
 
%% Check Inputs

    if(nargin < 3 || any([isempty(Y) isempty(x) isempty(h)]))
        error('Must provide Y, x and at least one bandwidth')
    end
    
    if(size(Y, 2) ~= size(x, 1))
        error('Dimensions of x and Y do not match')
    end

    n = size(Y, 2); p = size(x, 2); H = length(h);

    if(nargin < 4 || isempty(xout))
        xout = x;
    end
    
    m = size(xout, 1);

%% Set up for estimation

    lfr_fit = zeros(3, 1, H);
    M = spherefactory(3);
    ops.verbosity = 0;
    fr.M = M;

for l = 1:H
    
    
    hh = h(l);
    
  
  for j = 1:m
    
    x0 = xout(j, :)';
    % Get Weights
    w = getLFRweights(x, x0, hh);
    Kvec = K((x - repmat(x0', n, 1)), h*ones(1, p));
    Kvec(j) = 0;

    y0 = sum(cell2mat(arrayfun(@(k) Kvec(k)*Y(:, k), 1:n, 'UniformOutput', false))')'; y0 = y0/norm(y0); % Initial guess (leave-one-out NW)

    if(length(find(w)) < 2) % if there are less than two points, the fitted value cannot be computed

      lfr_fit(:, j, l) = NaN(1, 3); 

    else

      % Compute cost and Euclidean gradient

      fr.cost = @(y) get_cost(w, Y, y, M);
      fr.egrad = @(y) get_egrad(w, Y, y, M);
      fr.ehess = @(y, u) get_ehess(w, Y, y, M, u);

      lfr_fit(:, j, l) = trustregions(fr, y0, ops);

    end

  end

end

end