% For each value of p and n, read in simulation results.  For each bandwidth, compute:
%
% 1) The SE in estimating theta0 for each data set using geodesic distance - also compute mean and standard deviation of these
% 2) The MSPE (observation compared to estimated means) for each data set using geodesic distance - these are the achieved minimum criterion value - then compute mean and standard deviation
% 3) The MSEE (estimated means compared to true means) for each data set using geodesic distance - this requires known the covariates for each data set
%
% Also, compute 2) and 3) for local Frechet estimates, but not 1)
%
% Then, for each n and p, find h^* that minimizes the mean SE from 1), and report results for this bandwidth only

%% Simulation parameters

nsim = 200;
nVal = [50 100 200];
pVal = [2 5 10];

% Storage Objects

FSI_resAll_LN = cell(1, length(pVal));
FSI_resAll_HN = FSI_resAll_LN;

FSI_resSum_LN = cell(1, length(pVal)); % just summaries
FSI_resSum_HN = FSI_resSum_LN;

LF_resAll_LN = cell(1, length(pVal));
LF_resAll_HN = LF_resAll_LN;

LF_resSum_LN = cell(1, length(pVal));
LF_resSum_HN = LF_resSum_LN;

hVec_LN = cell(1, length(pVal));
hVec_HN = hVec_LN;

%% Compile Results for all Simulations and their Averages and Standard Deviations

for i = 1:length(pVal)

    % Storage

    FSI_resAll_LN{i} = cell(1, length(nVal));
    FSI_resAll_HN{i} = FSI_resAll_LN{i};
    FSI_resSum_LN{i} = cell(1, length(nVal));
    FSI_resSum_HN{i} = FSI_resSum_LN{i};

    LF_resAll_LN{i} = cell(1, length(nVal));
    LF_resAll_HN{i} = LF_resAll_LN{i};
    LF_resSum_LN{i} = cell(1, length(nVal));
    LF_resSum_HN{i} = LF_resSum_LN{i};
    
    hVec_LN{i} = cell(1, length(nVal));

    for j = 1:length(nVal)

        flnmLN = ['FinalSimResults/Sphere_results_n', num2str(nVal(j)), '_nsim' num2str(nsim) '_p' num2str(pVal(i)) '_noiseLow_Final.mat'];
        flnmHN = ['FinalSimResults/Sphere_results_n', num2str(nVal(j)), '_nsim' num2str(nsim) '_p' num2str(pVal(i)) '_noiseHigh_Final.mat'];

        % Results for low noise

        load(flnmLN)
        hVec_LN{i}{j} = h;
        H = length(h);
        %[Y, mreg] = getDataVals(theta0, n, p, nsim, s_dt, tau);

        FSI_resAll_LN{i}{j} = zeros(nsim, H, 3); % For FSI, three different summaries are computed for each simulated data set, bandwidth, and sample size
        LF_resAll_LN{i}{j} = zeros(nsim, H, 2); % For LF, only two summaries are computed since no index parameter is estimated
        
        for k = 1:nsim

            FSI_resAll_LN{i}{j}(k, :, 1) = min(acos(theta0'*fsiFitAll{k}.thetaHat).^2, acos(-theta0'*fsiFitAll{k}.thetaHat).^2);
            FSI_resAll_LN{i}{j}(k, :, 2) = arrayfun(@(l) mean(arrayfun(@(g) acos(Y{k}(:, g)'*fsiFitAll{k}.Yout(:, g, l))^2, 1:n)), 1:H); % prediction error on sphere
            FSI_resAll_LN{i}{j}(k, :, 3) = arrayfun(@(l) mean(arrayfun(@(g) acos(mreg{k}(:, g)'*fsiFitAll{k}.Yout(:, g, l))^2, 1:n)), 1:H); % estimation error on sphere

            LF_resAll_LN{i}{j}(k, :, 1) = arrayfun(@(l) mean(arrayfun(@(g) acos(Y{k}(:, g)'*LFpcovFitAll{k}(:, g, l))^2, 1:n)), 1:H); % prediction error on sphere
            LF_resAll_LN{i}{j}(k, :, 2) = arrayfun(@(l) mean(arrayfun(@(g) acos(mreg{k}(:, g)'*LFpcovFitAll{k}(:, g, l))^2, 1:n)), 1:H); % estimation error on sphere

        end

        FSI_resSum_LN{i}{j}.mean = squeeze(mean(FSI_resAll_LN{i}{j}));
        FSI_resSum_LN{i}{j}.std = squeeze(std(FSI_resAll_LN{i}{j}));

        LF_resSum_LN{i}{j}.mean = squeeze(mean(LF_resAll_LN{i}{j}));
        LF_resSum_LN{i}{j}.std = squeeze(std(LF_resAll_LN{i}{j}));        

        % Results for high noise

        load(flnmHN)
        hVec_HN{i}{j} = h;
        H = length(h);
        %[Y, mreg] = getDataVals(theta0, n, p, nsim, s_dt, tau);

        FSI_resAll_HN{i}{j} = zeros(nsim, H, 3); % For FSI, three different summaries are computed for each simulated data set, bandwidth, and sample size
        LF_resAll_HN{i}{j} = zeros(nsim, H, 2); % For LF, only two summaries are computed since no index parameter is estimated
        
        for k = 1:nsim

            FSI_resAll_HN{i}{j}(k, :, 1) = min(acos(theta0'*fsiFitAll{k}.thetaHat).^2, acos(-theta0'*fsiFitAll{k}.thetaHat).^2);
            FSI_resAll_HN{i}{j}(k, :, 2) = arrayfun(@(l) mean(arrayfun(@(g) acos(Y{k}(:, g)'*fsiFitAll{k}.Yout(:, g, l))^2, 1:n)), 1:H); % prediction error on sphere
            FSI_resAll_HN{i}{j}(k, :, 3) = arrayfun(@(l) mean(arrayfun(@(g) acos(mreg{k}(:, g)'*fsiFitAll{k}.Yout(:, g, l))^2, 1:n)), 1:H); % estimation error on sphere

            LF_resAll_HN{i}{j}(k, :, 1) = arrayfun(@(l) mean(arrayfun(@(g) acos(Y{k}(:, g)'*LFpcovFitAll{k}(:, g, l))^2, 1:n)), 1:H); % prediction error on sphere
            LF_resAll_HN{i}{j}(k, :, 2) = arrayfun(@(l) mean(arrayfun(@(g) acos(mreg{k}(:, g)'*LFpcovFitAll{k}(:, g, l))^2, 1:n)), 1:H); % estimation error on sphere

        end

        FSI_resSum_HN{i}{j}.mean = squeeze(mean(FSI_resAll_HN{i}{j}));
        FSI_resSum_HN{i}{j}.std = squeeze(std(FSI_resAll_HN{i}{j}));

        LF_resSum_HN{i}{j}.mean = squeeze(mean(LF_resAll_HN{i}{j}));
        LF_resSum_HN{i}{j}.std = squeeze(std(LF_resAll_HN{i}{j}));             
        
    end

end

%% Create Plots

% Theta Estimation Plots
col = {'-or', '-ob', '-ok'};
thetaPlot = tiledlayout(2, length(pVal), 'TileSpacing', 'tight');
title(thetaPlot, 'MSE for Theta Estimation')

for i = 1:length(pVal)
        
    nexttile
    hold on
    title(['p = ' num2str(pVal(i)) ', Low Noise'])
        
    for j = 1:length(nVal)

        plot(hVec_LN{i}{j}, FSI_resSum_LN{i}{j}.mean(:, 1), col{j})

    end
    
end

for i = 1:length(pVal)
    
	nexttile
    hold on
    title(['p = ' num2str(pVal(i)) ', High Noise'])
        
    for j = 1:length(nVal)

        plot(hVec_HN{i}{j}, FSI_resSum_HN{i}{j}.mean(:, 1), col{j})

    end

end

leg = legend('n = 50', 'n = 100', 'n = 200')
leg.Layout.Tile = 'east';

% Regression Estimation Plots
col2 = {'--xr', '--xb', '--xk'};
figure;
RegPlot = tiledlayout(2, length(pVal), 'TileSpacing', 'tight');
title(RegPlot, 'MSE for Regression Estimation')

for i = 1:length(pVal)
        
    nexttile
    hold on;
    title(['p = ' num2str(pVal(i)) ', Low Noise'])
        
    for j = 1:length(nVal)

        plot(hVec_LN{i}{j}, FSI_resSum_LN{i}{j}.mean(:, 3), col{j})

    end

    for j = 1:length(nVal)

        plot(hVec_LN{i}{j}, LF_resSum_LN{i}{j}.mean(:, 2), col2{j})

    end

end

for i = 1:length(pVal)
    
    nexttile 
    hold on;
    title(['p = ' num2str(pVal(i)) ', High Noise'])
        
    for j = 1:length(nVal)

        plot(hVec_HN{i}{j}, FSI_resSum_HN{i}{j}.mean(:, 3), col{j})
        
    end
 
    for j = 1:length(nVal)

        plot(hVec_HN{i}{j}, LF_resSum_HN{i}{j}.mean(:, 2), col2{j})
        
    end

end

leg = legend('n = 50', 'n = 100', 'n = 200');
leg.Layout.Tile = 'east';


%% Boxplots