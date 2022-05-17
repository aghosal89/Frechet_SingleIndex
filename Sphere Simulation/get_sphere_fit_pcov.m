%% Gets conditional Frechet mean estimates using local Frechet regression for spherical response data, generalized for vector predictors

% This implementation uses a product Gaussian kernel with a single bandwidth for
% simplicity
%
% Inputs: Y = response data 3xn 
%         x = nxp vector of predictors corresponding to observations in Y
%         h = bandwidth
%         xout = mxp vector of predictors to be used for computing fitted
%                values (if missing, set xout = x so m = n)

% Output: lfr_fit = 3xm matrix of fitted values using local Frechet
%                  regression for the predictor values in xout

function lfr_fit = get_sphere_fit_pcov(Y, x, h, xout)
 
%% Check Inputs

    if(nargin < 3 || any([isempty(Y) isempty(x) isempty(h)))
        error('Must provide Y, x and a bandwidth')
    end
    
    if(length(h) > 1)
        warning('Only one bandwidth supported - taking first element')
        h = h(1);
    end
    
    if(size(Y, 2) ~= size(x, 1))
        error('Dimensions of x and Y do not match')
    end

    n = size(Y, 2); p = size(x, 2);

    if(nargin < 4)
        xout = x;
    end
    
    m = size(xout, 1);

%% Set up for estimation

    lfr_fit = zeros(3, 1);
    M = spherefactory(3);
    ops.verbosity = 0;
    fr.M = M;
    
    % Get Weights
    w = getLFRweights(x, xout, h);
    Kmat = K((x - repmat(xout', n, 1)), h*ones(1, p));
    KmatLO = Kmat - diag(diag(Kmat));
  
  for j = 1:m

    y0 = sum(cell2mat(arrayfun(@(k) KmatLO(k, j)*Y(:, k), 1:n, 'UniformOutput', false))')'; y0 = y0/norm(y0); % Initial guess (leave-one-out NW)

    if(length(find(w(:, j))) < 2) % if there are less than two points, the fitted value cannot be computed

      lfr_fit(:, j) = repmat(NaN, 1, 3); 

    else

      % Compute cost and Euclidean gradient

      fr.cost = @(y) get_cost(w(:, j)', Y, y, M);
      fr.egrad = @(y) get_egrad(w(:, j)', Y, y, M);
      fr.ehess = @(y, u) get_ehess(w(:, j)', Y, y, M, u);

      lfr_fit(:, j) = trustregions(fr, y0, ops);

    end

  end

end