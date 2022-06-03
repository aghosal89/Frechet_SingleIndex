%% Fits the FSI model for spherical response data

% Inputs: Y = response data 3xn 
%         x = nxp matrix of predictors corresponding to observations in Y
%         h = 1xH vector of bandwidths
%         xout = mxp matrix of predictors to be used for computing fitted
%                values (if missing, set xout = x so m = n)
%         theta_init = pxq matrix of q starting values for optimizing over the
%                      index parameter OR an integer indicating the number 
%                      of random starting points to draw (if missing, 10 random 
%                      starting points will be drawn)
%         opts = a structure array with the following fields:
%
%           Algorithm = can be either 'NM' (for Nelder-Mead optimization
%                       using fminsearchbnd) or 'TR' (for trust region
%                       optimization using fmincon) - default is 'TR'
%           Ret = binary flag for returning optimizer info (1 = yes [default], 0 = no)
%           maxIter = integer indicating maximum number of iterations to be
%                     executed (default is 100)
%           useAltCost = binary flag indicating whether or not an
%                        alternative cost function should be used based on 
%                        normalized local linear estimates (1 = yes, 0 = no
%                        [default])
%           useAltGrad = binary flag indicating whether or not the gradient
%                        for the alternative cost function should be used
%                        (1 = yes, 0 = no [default]).  This is only
%                        utilized if Algorithm = 'TR'.  Note that if
%                        Algorithm = 'TR' and useAltCost = 1, then the
%                        alternative gradient is the exact gradient for the
%                        alternative cost.
%           useAltHess = binary flag indicating whether or not the Hessian
%                        for the alternative cost function should be used
%                        (1 = yes, 0 = no [default]).  This is only
%                        utilized if Algorithm = 'TR'.  Note that if
%                        Algorithm = 'TR' and useAltCost = 1, then the
%                        alternative Hessian is the exact Hessian for the
%                        alternative cost.

% Output: fsi_fit = a structure array with the following fields:
%           thetaHat = pxH vector of the index parameter estimates, each
%                      corresponding to one of the bandwidths in h
%           Yout = 3xmxH array of fitted values using thetaHat estimates and local Frechet
%                  regression for the predictor values in xout for all
%                  bandwidths in h
%           optInfo = if opts.Ret = 1, then this is a cell array with four
%                     elements.  The first contains the matrix of starting
%                     values, the second contains the array of optimum values returned 
%                     by the algorithm for various h, the third contains the exit flags for each
%                     starting value, and the fourth contains the additional
%                     information for each optimization problem

function fsi_fit = get_sphere_fit_FSI(Y, x, h, xout, theta_init, opts)


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
    
    if(nargin < 5 || isempty(theta_init) || (length(theta_init) ~= 1 && size(theta_init, 1) ~= p))
        warning('Initial values are missing are invalid - randomly generating 10 initial values')
        eta_init = pi*rand(p - 1, 10) - pi/2;
        theta_init = cell2mat(arrayfun(@(j) polar2cart(eta_init(:, j)), 1:10, 'UniformOutput', false));
    elseif(length(theta_init) == 1)
        q = length(theta_init);
        clear theta_init; % remove, then replace with actual theta values
        eta_init = pi*rand(p - 1, q) - pi/2;
        theta_init = cell2mat(arrayfun(@(j) polar2cart(eta_init(:, j)), 1:q, 'UniformOutput', false));
    end
    
    q = size(theta_init, 2);
    optsDef = struct('Algorithm', 'NM', 'Ret', 1, 'maxIter', 100, 'useAltCost', 0, 'useAltGrad', 0, 'useAltHess', 0);
    
    
    if(nargin < 6 || isempty(opts) || ~isstruct(opts))
        warning('Invalid input for opts - using defaults')
        opts = optsDef;
    else
        if(isfield(opts, 'Algorithm'))
            optsDef.Algorithm = opts.Algorithm;
        end
        if(isfield(opts, 'Ret'))
            optsDef.Ret = opts.Ret;
        end
        if(isfield(opts, 'maxIter'))
            optsDef.maxIter = opts.maxIter;
        end
        if(isfield(opts, 'useAltCost'))
            optsDef.useAltCost = opts.useAltCost;
        end
        if(isfield(opts, 'useAltGrad'))
            optsDef.useAltGrad = opts.useAltGrad;
        end
        if(isfield(opts, 'useAltHess'))
            optsDef.useAltHess = opts.useAltHess;
        end
        opts = optsDef; clear optsDef; % reset and clear
    end
    
    if(strcmp(opts.Algorithm, 'TR')) % Define options input for fmincon
        options = optimoptions('fmincon', 'Display', 'off', 'Algorithm', 'trust-region-reflective', 'MaxIterations', opts.maxIter);
        if(opts.useAltGrad)
            options = optimoptions(options, 'SpecifyObjectiveGradient', true);
        end
        if(opts.useAltHess)
            options = optimoptions(options, 'HessianFcn', 'objective');
        end
    else % define options input for fminsearchbnd
        options = optimset('Display', 'off', 'MaxIter', opts.maxIter);
    end

%% Create objects for storing results
    thetaHat = zeros(p, H); Yout = zeros(3, m, H);
    if(opts.Ret)
        optInfo = cell(1, 4);
        optInfo{1} = theta_init;
        optInfo{2} = zeros(p, q, H); % one p-vector at convergence per starting value per bandwidth
        optInfo{3} = zeros(q, H); % one flag per staring value per bandwidth
        optInfo{4} = cell(q, H); % one set of output info per starting value per bandwidth
    else
        optInfo = [];
    end

    

%% Execute Estimation

    for l = 1:H

        % Define cost function for optimization
        hh = h(l);
        
            
        cost_fun = @(eta) WnCost(eta, Y, x, hh, opts.useAltCost);
            
        WnVec = zeros(1, q); % for storing the criterion value
        thetaMat = zeros(p, q);

        for j = 1:q

            etaCur = cart2polar(theta_init(:, j));
            if(strcmp(opts.Algorithm, 'NM'))
                
                [etaOpt, WnVal, fl, op] = fminsearchbnd(cost_fun, etaCur, -pi/2*ones(1, p - 1), pi/2*ones(1, p - 1), options);
    
            else
                
                [etaOpt, WnVal, fl, op] = fmincon(cost_fun, etaCur, [], [], [], [], -pi/2*ones(1, p - 1), pi/2*ones(1, p - 1), [], options);
            
            end
            
            % Store current optimizer and optimum
            thetaMat(:, j) = polar2cart(etaOpt);
            WnVec(j) = WnVal;

            if(opts.Ret)
        
                optInfo{2}(:, j, l) = thetaMat(:, j);
                optInfo{3}(j, l) = fl;
                optInfo{4}{j, l} = op;

            end

        end

        [~, ind] = min(WnVec);
        thetaHat(:, l) = thetaMat(:, ind);
        Yout(:, :, l) = get_sphere_fit_LF(Y, x*thetaHat(:, l), hh, xout*thetaHat(:, l));

    end

    fsi_fit = struct('thetaHat', thetaHat, 'Yout', Yout, 'optInfo', {optInfo});

end