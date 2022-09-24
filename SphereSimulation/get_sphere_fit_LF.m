%% Gets conditional Frechet mean estimates using local Frechet regression for spherical response data and scalar predictor

% Inputs: Y = response data 3xn 
%         x = nx1 vector of predictors corresponding to observations in Y
%         h = bandwidth
%         xout = mx1 vector of predictors to be used for computing fitted
%                values (if missing, set xout = x so m = n)

% Output: lf_fit = 3xm matrix of fitted values using local Frechet
%                  regression for the predictor values in xout

function lf_fit = get_sphere_fit_LF(Y, x, h, xout)

%% Check Inputs

    if(nargin < 3 || any([isempty(Y) isempty(x) isempty(h)]))
        error('Must provide Y, x and a bandwidth')
    end
    
    if(length(h) > 1)
        warning('Only one bandwidth supported - taking first element')
        h = h(1);
    end
    
    if(size(Y, 2) ~= length(x))
        error('Dimensions of x and Y do not match')
    end

    n = size(Y, 2);

    if(nargin < 4)
        xout = x;
    end
    
    m = length(xout);

%% Set up for estimation

  lf_fit = zeros(3, m);
  M = spherefactory(3);
  ops.verbosity = 0;
  fr.M = M;

  % Get weights
  
  for j = 1:m
      
    x0 = xout(j, :)';
    % Get Weights
    w = getLFRweights(x, x0, h);
    Kvec = K((x - repmat(x0', n, 1)), h);
    Kvec(j) = 0;

    y0 = sum(cell2mat(arrayfun(@(k) Kvec(k)*Y(:, k), 1:n, 'UniformOutput', false))')'; y0 = y0/norm(y0); % Initial guess (leave-one-out NW)

    if(length(find(w)) < 2) % if there are less than two points, the fitted value cannot be computed

      lf_fit(:, j) = NaN(1, 3); 

    else

      % Compute cost and Euclidean gradient

      fr.cost = @(y) get_cost(w, Y, y, M);
      fr.egrad = @(y) get_egrad(w, Y, y, M);
      fr.ehess = @(y, u) get_ehess(w, Y, y, M, u);

      lf_fit(:, j) = trustregions(fr, y0, ops);

    end

  end

end
