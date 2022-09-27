%% Code to create the figure of sphere plots in the document

% clear the working directory 
clear all

% create path to manifold optimization folder
addpath(genpath('manopt'))

% Define some global parameters
p = 5; pInd = 2; % we will look at data sets for covariate dimension 5 only
i = 1; % we will work with the first simulated data set in each setting
tmp = load('FinalSimResults/Sphere_results_n200_nsim200_p5_noiseHigh_Final.mat', 'theta0'); 
theta0 = tmp.theta0; clear tmp; % load true parameter value
thf0 = [0.7422, 0.1967, 0.2693, 0.0084, 0.5812]'; % arbitrary parameter value considerably far from the true parameter.

%% For each setting, read in a data set, compute estimated regression curves and plot with data and true curve.

nVal = [50, 100, 200];
noise = {'Low', 'High'};
load hopt_MSEE_indices

for nInd = 1:length(nVal)
    for nlev = 1:length(noise)
    
        flnm = ['FinalSimResults/Sphere_results_n' num2str(nVal(nInd)) '_nsim200_p5_noise' noise{nlev} '_Final.mat'];
        tmp = load(flnm, 'x', 'Y', 'h', 'fsiFitAll');
        x = squeeze(tmp.x(:, i, :)); Y = tmp.Y{i}; 
        if nlev == 1
            hopt_MSEE = tmp.h(hopt_MSEE_ln(nInd, pInd));
            thetaHat = tmp.fsiFitAll{i}.thetaHat(:, hopt_MSEE_ln(nInd, pInd));     
        else
            hopt_MSEE = tmp.h(hopt_MSEE_hn(nInd, pInd));
            thetaHat = tmp.fsiFitAll{i}.thetaHat(:, hopt_MSEE_hn(nInd, pInd)); 
        end
        clear tmp;

        z = x*theta0./sqrt(p); zout = linspace(min(z), max(z), 101);
        zHat = x*thetaHat./sqrt(p); zHatout = linspace(min(zHat), max(zHat), 101)';
        zFar = x*thf0./sqrt(p); zFarout = linspace(min(zFar), max(zFar), 101)';

        % Get true regression on zout
        ut = [sqrt(1-zout'.^2).*cos(pi*zout'), sqrt(1-zout'.^2).*sin(pi*zout'), zout']';
        m_tr = arrayfun(@(j) ut(j, :), 1:3, 'UniformOutput', false);

        % Get estimates on zHatout and zFarout
        YHat = get_sphere_fit_LF(Y, zHat, hopt_MSEE, zHatout);
        YFar = get_sphere_fit_LF(Y, zFar, hopt_MSEE, zFarout);

        [s1, s2, s3] = sphere(128);
        lightGrey = 0.9*ones(1, 3);

        set(gcf, 'PaperUnits', 'inches', 'PaperSize', [2.4, 2.6], 'PaperPosition', [0 0 2.4 2.6])
        surf(s1, s2, s3, 'FaceColor', lightGrey, 'EdgeColor', 'none', 'FaceAlpha', 0.4)
        view([-1, 1, 1])
        hold on;
        plot3(m_tr{1}, m_tr{2}, m_tr{3}, 'k', 'LineWidth', 0.5)
        plot3(Y(1, :), Y(2, :), Y(3, :), 'r.', 'MarkerSize', 2.5)

        hold on;

        plot3(YHat(1, :), YHat(2, :), YHat(3, :), 'b.', 'MarkerSize', 2.5)

        hold on;
        plot3(YFar(1, :), YFar(2, :), YFar(3, :), 'g.', 'MarkerSize', 2.5)

        % save the image in an .eps file
        saveas(gcf,['noise' noise{nlev} '_n' num2str(nVal(nInd)) '_p' num2str(p)'],'epsc')

        % finally save the outputs in an extended file to be reproduced further.
%        save([flnm(1:(end - 4)) '_extended.mat'], 'zHat', 'YHat', 'thf0', 'zFar', 'YFar', 'zout', 'm_tr', 'hopt_MSEE')
        close all
    end
end