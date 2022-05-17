%% Fits the FSI model for spherical response data

% Inputs: Y = response data 3xn 
%         x = nxp matrix of predictors corresponding to observations in Y
%         h = 1xH vector of bandwidths
%         xout = mxp matrix of predictors to be used for computing fitted
%                values (if missing, set xout = x so m = n)
%         theta_init = pxq matrix of q starting values for optimizing over the
%                      index parameter (if missing, 10 random starting points 
%                      will be drawn)
%         optRet = flag for returning optimizer info (1 = yes [default], 0
%         = no)

% Output: fsi_fit = a struct object with the following fields:
%           thetaHat = pxH vector of the index parameter estimates, each
%                      corresponding to one of the bandwidths in h
%           Yout = 3xmxH array of fitted values using thetaHat estimates and local Frechet
%                  regression for the predictor values in xout for all
%                  bandwidths in h
%           optInfo = if optRet = 1, then this is a cell array with four
%                     elements.  The first contains the matrix of starting
%                     values, the second contains the array of optimum values returned 
%                     by the algorithm for various h, the third contains the exit flags for each
%                     starting value, and the fourth contains the additional
%                     information for each optimization problem

function fr_fit = get_sphere_fit_FSI(Y, x, xout, h, theta_init, optRet)


%% Check Inputs

    if(nargin < 3 || any([isempty(Y) isempty(x) isempty(h)]))
        error('Must provide Y, x and at least one bandwidth')
    end
    
    if(size(Y, 2) ~= size(x, 1))
        error('Dimensions of x and Y do not match')
    end

    n = size(Y, 2); p = size(x, 2); H = length(h);

    if(nargin < 4)
        xout = x;
    end
    
    m = size(xout, 2);
       
    if(nargin < 5 || isempty(theta_init) || size(theta_init, 1) ~= p)
        warning('Initial values are missing are invalid - randomly generating 10 initial values')
        eta_init = pi*rand(p - 1, 10) - pi/2;
        theta_init = cell2mat(arrayfun(@(j) polar2cartGen(eta_init(:, j)), 1:10, 'UniformOutput', false));
    end
    
    q = size(theta_init, 2);
    
    if(nargin < 6 || isempty(optRet) || length(optRet) ~= 1 || ~ismember(optRet, [0 1]))
        warning('Input optRet is missing or invalid - setting to 1')
        optRet = 1;
    end

%% Create objects for storing results
    thetaHat = zeros(p, H); Yout = zeros(3, m, H);
    if(optRet)
        optInfo = cell(1, 4);
        optInfo{1} = theta_init;
        optInfo{2} = zeros(p, q, H); % one p-vector at convergence per starting value per bandwidth
        optInfo{3} = zeros(q, H); % one flag per staring value per bandwidth
        optInfo{3} = cell(q, H); % one set of output info per starting value per bandwidth
    else
        optInfo = [];
    end

    

%% Execute Estimation

    for l = 1:H

        % Define cost function for optimization
        hh = h(l);
        cost_fun = @(eta) nm_cost(eta, Y, x, hh);
        WnVec = zeros(1, q); % for storing the criterion value
        thetaMat = zeros(p, q);

        for j = 1:q

            etaCur = cart2polarGen(theta_init(:, j));
            [etaOpt, WnVal, fl, op] = fminsearchbnd(cost_fun, etaCur, -pi/2*ones(1, p - 1), pi/2*ones(1, p - 1));
        
            % Store current optimizer and optimum
            thetaMat(:, j) = polar2cartGen(etaOpt);
            WnVec(j) = WnVal;

            if(optRet)
        
                optInfo{2}(:, j, l) = thetaMat(:, j);
                optInfo{3}(j, l) = fl;
                optInfo{4}{j, l} = op;

            end

        end

        [~, ind] = min(WnVec);
        thetaHat(:, l) = thetaMat(:, ind);
        Yout(:, :, l) = get_sphere_fit_LF(Y, x*thetaHat(:, l), hh, xout*thetaHat(:, l));

    end

    fsi_fit = struct('thetaHat', thetaHat, 'Yout', Yout, 'optInfo', optInfo);

end