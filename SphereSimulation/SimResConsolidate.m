% Consolidate mat files by simulation setting - this code works on UNIX systems only

% Run this script after running all simulation batches to put together all
% results for the same setting

nval = [50 100 200];
pval = [2 5 10];
noise = {'Low', 'High'};

nsim = 200; % all settings have 200 runs

for j = 1:length(nval)

  for k = 1:length(pval)

    n = nval(j); p = pval(k);

    for l = 1:length(noise)

      % Initialize all variables that will be saved to final result file
      h = []; % empty - will read in bandwidths used
      theta0 = []; % empty - will read in true parameter
      tau = []; % empty - will read in noise variance
      x = zeros(n, nsim, p); % store covariates
      Y = cell(1, nsim); % store responses
      mreg = cell(1, nsim); % store true conditional Frechet means
      LFpcovFitAll = cell(1, nsim); % store local Frechet results
      fsiFitAll = cell(1, nsim); % store FSI results

      fls = dir(['Sphere_results_n' num2str(n) '*_p' num2str(p) '_noise' noise{l} '*.mat']);

      curSim = 0; % index to track the simulation numbers being read in

      for m = 1:length(fls)

        dt = load(fls(m).name);

        if(m == 1) % for the first file, store shared parameters

          h = getfield(dt, 'h');
          theta0 = getfield(dt, 'theta0');
          tau = getfield(dt, 'tau');

        end

        ind = (curSim + 1):(curSim + getfield(dt, 'nsim'));

        x(:, ind, :) = getfield(dt, 'x');
        Y(ind) = getfield(dt, 'Y');
        mreg(ind) = getfield(dt, 'mreg');
        LFpcovFitAll(ind) = getfield(dt, 'LFpcovFitAll');
        fsiFitAll(ind) = getfield(dt, 'fsiFitAll');

        curSim = curSim + getfield(dt, 'nsim');

      end % m

      fnm = strcat('FinalSimResults/Sphere_results_n', num2str(n), '_nsim', num2str(nsim), '_p', num2str(p), '_noise', noise{l}, '_Final.mat');

      save(fnm, 'n', 'p', 'theta0', 'mreg', 'tau', 'h', 'x', 'Y', 'LFpcovFitAll', 'fsiFitAll')

    end % l
  end % k
end % j
