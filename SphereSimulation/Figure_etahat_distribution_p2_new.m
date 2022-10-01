%% Codes to create figure 3 in document

clear all

% Parameter values
nsim = 200;
p=2;
nVal = [50 100 200];
noise = {'Low', 'High'};

% create the vector for MSE values
eta_sim = zeros(200, length(nVal)*length(noise));

for nlev = 1:length(noise)
    for nInd = 1:length(nVal)

        cur = (nlev - 1)*(1+length(noise)) + nInd; % indexes the column we want to populate

        flnm = ['FinalSimResults/Sphere_results_n' num2str(nVal(nInd)) '_nsim' num2str(nsim) '_p' num2str(p) '_noise' cell2mat(noise(nlev)) '_Final.mat'];
        tmp = load(flnm, 'h', 'Y', 'fsiFitAll', 'theta0');

        theta0 = tmp.theta0;
        MSE_h_sim = zeros(200, length(tmp.h));

        for i=1:200
            for l=1:length(tmp.h)
                MSE_h_sim(i, l) = acos(abs(tmp.fsiFitAll{i}.thetaHat(:,l)'*theta0))^2;
            end
        end

        MSE_h = mean(MSE_h_sim);

        xh = find(MSE_h== min(MSE_h));

        for i=1:200
            [~, eth] = cart2polar(tmp.fsiFitAll{i}.thetaHat(:,xh)); 
            eta_sim(i, cur) = eth;
        end

    end
end


fig=figure;

subplot(2,3,1)
histogram(eta_sim(:,1), 'NumBins', 20, 'Normalization','pdf');
grid on
[~, eta0]=cart2polar(theta0);
xline(eta0, 'rx', 'LineWidth', 2.0);
set(findall(gcf,'-property','FontSize'),'FontSize',15);
xlabel('$\hat{{\eta}}$','interpreter','latex', 'FontWeight','bold');
xlim([-1.25,-.1]);
ylim([0, 6.5]);
title(['n=50', ', $\sigma^2=0.4$'],'interpreter','latex', 'FontWeight','bold');

subplot(2,3,2)
histogram(eta_sim(:,2), 'NumBins', 20, 'Normalization','pdf');
grid on
[~, eta0]=cart2polar(theta0);
xline(eta0, 'rx', 'LineWidth', 2.0);
set(findall(gcf,'-property','FontSize'),'FontSize',15);
xlabel('$\hat{{\eta}}$','interpreter','latex', 'FontWeight','bold');
xlim([-1.25,-.1]);
ylim([0, 6.5]);
title(['n=100', ', $\sigma^2=0.4$'],'interpreter','latex', 'FontWeight','bold');

subplot(2,3,3)
histogram(eta_sim(:,3), 'NumBins', 20, 'Normalization','pdf');
grid on
[~, eta0]=cart2polar(theta0);
xline(eta0, 'rx', 'LineWidth', 2.0);
set(findall(gcf,'-property','FontSize'),'FontSize',15);
xlabel('$\hat{{\eta}}$','interpreter','latex', 'FontWeight','bold');
xlim([-1.25,-.1]);
ylim([0, 6.5]);
title(['n=200', ', $\sigma^2=0.4$'],'interpreter','latex', 'FontWeight','bold');

subplot(2,3,4)
histogram(eta_sim(:,4), 'NumBins', 20, 'Normalization','pdf');
grid on
[~, eta0]=cart2polar(theta0);
xline(eta0, 'rx', 'LineWidth', 2.0);
set(findall(gcf,'-property','FontSize'),'FontSize',15);
xlabel('$\hat{{\eta}}$','interpreter','latex', 'FontWeight','bold');
xlim([-1.58,1.58]);
ylim([0, 3.5]);
title(['n=50', ', $\sigma^2=0.8$'],'interpreter','latex', 'FontWeight','bold');

subplot(2,3,5)
histogram(eta_sim(:,5), 'NumBins', 20, 'Normalization','pdf');
grid on
[~, eta0]=cart2polar(theta0);
xline(eta0, 'rx', 'LineWidth', 2.0);
set(findall(gcf,'-property','FontSize'),'FontSize',15);
xlabel('$\hat{{\eta}}$','interpreter','latex', 'FontWeight','bold');
xlim([-1.58,1.58]);
ylim([0, 3.5]);
title(['n=100', ', $\sigma^2=0.8$'],'interpreter','latex', 'FontWeight','bold');

subplot(2,3,6)
histogram(eta_sim(:,6), 'NumBins', 20, 'Normalization','pdf');
grid on
[~, eta0]=cart2polar(theta0);
xline(eta0, 'rx', 'LineWidth', 2.0);
set(findall(gcf,'-property','FontSize'),'FontSize',15);
xlabel('$\hat{{\eta}}$','interpreter','latex', 'FontWeight','bold');
xlim([-1.58,1.58]);
ylim([0, 3.5]);
title(['n=200', ', $\sigma^2=0.8$'],'interpreter','latex', 'FontWeight','bold');

han=axes(fig,'visible','off'); 
han.Title.Visible='on';
han.XLabel.Visible='on';
han.YLabel.Visible='on';

set(findall(gcf,'-property','FontSize'),'FontSize',15);
ylabel('Density','interpreter','latex', 'FontWeight','bold');

saveas(gcf,'eta_hat_histogram_p2','epsc')
