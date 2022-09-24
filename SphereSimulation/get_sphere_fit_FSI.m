%% Fits the FSI model for spherical response data

% Inputs: Y = response data 3xn
%         x = nxp matrix of predictors corresponding to observations in Y
%         h = 1xH vector of bandwidths
%         xout = mxp matrix of predictors to be used for computing fitted
%                values (if missing, set xout = x so m = n)
%
%         opts = a structure array with the following fields:
%
%           Ret = binary flag for returning optimizer info (1 = yes [default], 0 = no)
%           maxIter = integer indicating maximum number of iterations to be
%                     executed [default = 100]
%           thetaInit = pxq matrix of q starting values for initial optimization over the
%                      index parameter OR an integer indicating the number
%                      of random starting points to draw [default = 10]
%           thetaChoose = integer indicating how many of the values returned in the
%                          initial optimization step should be used in the final round
%                          using Nelder-Mead [default = 5]; Note that iterations that
%                          did not converge in the initial optimization will be discarded
%           Algorithm = 'IP', 'SQP', or 'NM' [default]; except for 'NM', 
%                       will use fmincon with interior-point ('IP') or
%                       sequential quadratic programming ('SQP') algorithms
%                       in the final optimization; if 'NM', Nelder-Mead 
%                       will be used (fminsearchbnd)

% Output: fsi_fit = a structure array with the following fields:
%           thetaHat = pxH vector of the index parameter estimates, each
%                      corresponding to one of the bandwidths in h
%           Yout = 3xmxH array of fitted values using thetaHat estimates and local Frechet
%                  regression for the predictor values in xout for all
%                  bandwidths in h
%           optInfo = if opts.Ret = 1, then this is a cell array with five
%                     elements.
%               optInfo{1} = pxq matrix of starting values for the initial optimization
%               optInfo{2} = 1xH cell array, optInfo{2}{l} is a matrix of starting values
%                            used in the final optimization for bandwidth h(l)
%               optInfo{3} = 1xH cell array, optInfo{3}{l} is a matrix of values
%                            returned by the final optimization corresponding to
%                             starting values in optInfo{2}{l}
%               optInfo{4} = 1xH cell array, optInfo{4}{l} contains the exit flags
%                            for each starting value in the final optimization
%               optInfo{5} = 1xH cell array, optInfo{5}{l} contains additional information
%                            for each final optimization problem
%
% Details: For each bandwidth, the estimate for theta is produced in two
% steps.  First, an approximation to the estimation criterion is optimized
% using starting points specified by opts.thetaInit.  This step uses
% local linear estimates (normalized to be on the sphere) of the conditional
% Euclidean mean in lieu of local Frechet estimates of the conditional
% Frechet mean, and is very fast since these estimates are known in closed
% form, also allowing for an analytic gradient and Hessian.  Out of the
% optimizing values from the first step, the best few (determined by
% opts.thetaChoose) are identified by having the lowest approximate cost.
% In the second step, this reduced set of theta values
% are used as starting values to the optimizer using the true cost
% function.  Depending on opts.Algorithm, this second step optimization,
% which has no analytic gradient or Hessian, can be conducted by
% Nelder-Mead or by a trust region algorithm that uses the analytic
% gradient and Hessian from the approximate cost function instead of
% approximating the true gradient and Hessian numerically.

function fsi_fit = get_sphere_fit_FSI(Y, x, h, xout, opts)


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

    optsDef = struct('Ret', 1, 'maxIter', 100, 'thetaInit', 10, 'thetaChoose', 5, 'Algorithm', 'NM');

    if(nargin < 5 || isempty(opts) || ~isstruct(opts))
      opts = optsDef;
      clear optsDef;
    else
      if(~isfield(opts, 'Ret') || isempty(opts.Ret))
        opts.Ret = optsDef.Ret;
      end
      if(~isfield(opts, 'maxIter') || isempty(opts.maxIter))
        opts.maxIter = optsDef.maxIter;
      end
      if(~isfield(opts, 'thetaInit') || isempty(opts.thetaInit))
        opts.Ret = optsDef.thetaInit;
      end
      if(~isfield(opts, 'thetaChoose') || isempty(opts.thetaChoose))
        opts.maxIter = optsDef.thetaChoose;
      end
      if(~isfield(opts, 'Algorithm') || isempty(opts.Algorithm))
        opts.maxIter = optsDef.Algorithm;
      end
      clear optsDef;
    end

    if(length(opts.thetaInit) == 1)
        q = opts.thetaInit;
        eta_init = pi*rand(p - 1, q) - pi/2;
        opts.thetaInit = cell2mat(arrayfun(@(j) polar2cart(eta_init(:, j)), 1:q, 'UniformOutput', false));
    end

    q = size(opts.thetaInit, 2);

    % Define optimizer options

    optionsPrelim = optimoptions('fmincon', 'Display', 'off', 'Algorithm', 'trust-region-reflective', 'MaxIterations', opts.maxIter, 'SpecifyObjectiveGradient', true, 'HessianFcn', 'objective');
    
    switch opts.Algorithm
        
        case 'NM'
        
            optionsFinal = optimset('Display', 'off', 'MaxIter', opts.maxIter);    
            
        case 'IP'
            
            optionsFinal = optimoptions('fmincon', 'Display', 'off', 'Algorithm', 'interior-point', 'MaxIterations', opts.maxIter);
    
        case 'SQP'
            
            optionsFinal = optimoptions('fmincon', 'Display', 'off', 'Algorithm', 'sqp', 'MaxIterations', opts.maxIter);

        otherwise
            
            error('Invalid choice for Algorithm option')
            
    end

%% Create objects for storing results
    thetaHat = zeros(p, H); Yout = zeros(3, m, H);
    if(opts.Ret)
        optInfo = cell(1, 5);
        optInfo{1} = opts.thetaInit;
        optInfo{2} = cell(1, H); % each element will contain the starting values for final optimization, which may differ by bandwidth
        optInfo{3} = cell(1, H); % each element will contain the values returned by the final optimization, which may differ by bandwidth
        optInfo{4} = cell(1, H); % one flag per starting value per bandwidth
        optInfo{5} = cell(1, H); % one set of output info per starting value per bandwidth
    else
        optInfo = [];
    end



%% Execute Estimation

    for l = 1:H
        
%        disp(['Bandwidth ' num2str(l) ' of ' num2str(H)])
        % Define cost function for initial optimization
        hh = h(l);
        costFunAlt = @(eta) WnCost(eta, Y, x, hh, 1);

        WnAltVec = zeros(1, q); % for storing the criterion value
        thetaMat = zeros(p, q);

        for j = 1:q

            etaCur = cart2polar(opts.thetaInit(:, j));
            [etaOpt, WnVal, fl] = fmincon(costFunAlt, etaCur, [], [], [], [], -pi/2*ones(1, p - 1), pi/2*ones(1, p - 1), [], optionsPrelim);

            % Store current optimizer and optimum
            thetaMat(:, j) = polar2cart(etaOpt);
            if(ismember(fl, [0, -1, -2]))
              WnAltVec(j) = Inf; % if optimization fails to converge, set cost to infinity to avoid selecting this value
            else
              WnAltVec(j) = WnVal;
            end

        end

        [~, ind] = sort(WnAltVec);
        ind = ind(1:min(size(thetaMat, 2), opts.thetaChoose));

        thetaMat = round(thetaMat, 4); % round to pool together nearby convergent points
        thetaInit2 = unique(thetaMat(:, ind)', 'rows')'; % Starting values for final optimization
        qq = size(thetaInit2, 2);

        if(opts.Ret)
          optInfo{2}{l} = thetaInit2;
          optInfo{3}{l} = zeros(p, qq);
          optInfo{4}{l} = zeros(1, qq);
          optInfo{5}{l} = cell(1, qq);
        end

        % Execute final optimization on reduced set of starting values

        costFun = @(eta) WnCost(eta, Y, x, hh);
        WnVec = zeros(1, qq);
        thetaMat = zeros(p, qq);

        for j = 1:qq
%            disp(['Starting Value ' num2str(j) ' of ' num2str(qq)])

            etaCur = cart2polar(thetaInit2(:, j)/norm(thetaInit2(:, j)));
            if(strcmp(opts.Algorithm, 'NM'))
              [etaOpt, WnVal, fl, op] = fminsearchbnd(costFun, etaCur, -pi/2*ones(1, p - 1), pi/2*ones(1, p - 1), optionsFinal);
            else
              [etaOpt, WnVal, fl, op] = fmincon(costFun, etaCur, [], [], [], [], -pi/2*ones(1, p - 1), pi/2*ones(1, p - 1), [], optionsFinal);
            end

            % Store current optimizer and optimum
            thetaMat(:, j) = polar2cart(etaOpt);
            if(ismember(fl, [-2, -1]))
              WnVec(j) = Inf; % if optimization fails, set cost to infinity to avoid selecting this value
            else
              WnVec(j) = WnVal;
            end

            if(opts.Ret)

                optInfo{3}{l}(:, j) = thetaMat(:, j);
                optInfo{4}{l}(j) = fl;
                optInfo{5}{l}{j} = op;

            end

        end

        if(any(WnVec < Inf))
          [~, ind] = min(WnVec);
          thetaHat(:, l) = thetaMat(:, ind);
          Yout(:, :, l) = get_sphere_fit_LF(Y, x*thetaHat(:, l), hh, xout*thetaHat(:, l));
        else
          thetaHat(:, l) = NaN(p, 1);
          Yout(:, :, l) = NaN(3, n);
        end
    end

    fsi_fit = struct('thetaHat', thetaHat, 'Yout', Yout, 'optInfo', {optInfo});

end
