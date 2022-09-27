
%% Codes to create histograms of eta hats for p=2 for various sample sizes and noise levels

fig=figure;

subplot(2,3,1)
load('Sphere_results_n50_nsim200_p2_noiseLow_Final.mat')

% choose bandwidth
MSE_h_sim = zeros(200, length(h));

for i=1:200
    for l=1:length(h)
        MSE_h_sim(i,l) = acos(abs(fsiFitAll{i}.thetaHat(:,l)'*theta0))^2;
    end
end

MSE_h = mean(MSE_h_sim);

xh = find(MSE_h== min(MSE_h));

eta_hat = NaN(200,1);

for i=1:200
    [~, eth] = cart2polar(fsiFitAll{i}.thetaHat(:,xh)); 
    eta_hat(i) = eth;
end

histogram(eta_hat, 'NumBins', 20, 'Normalization','pdf');
grid on
[~, eta0]=cart2polar(theta0);
xline(eta0, 'rx', 'LineWidth', 2.0);
set(findall(gcf,'-property','FontSize'),'FontSize',15);
xlabel('$\hat{{\eta}}$','interpreter','latex', 'FontWeight','bold');
%ylabel('Density','interpreter','latex', 'FontWeight','bold');
xlim([-1.25,-.1]);
ylim([0, 6.5]);
title('n=50, $\sigma^2=0.4$','interpreter','latex', 'FontWeight','bold');


subplot(2,3,2)
load('Sphere_results_n100_nsim200_p2_noiseLow_Final.mat')

% bandwidth
MSE_h_sim = zeros(200, length(h));

for i=1:200
    for l=1:length(h)
        MSE_h_sim(i,l) = acos(abs(fsiFitAll{i}.thetaHat(:,l)'*theta0))^2;
    end
end

MSE_h = mean(MSE_h_sim);

xh = find(MSE_h== min(MSE_h));

eta_hat = NaN(200,1);

for i=1:200
    [~, eth] = cart2polar(fsiFitAll{i}.thetaHat(:,xh)); 
    eta_hat(i) = eth;
end

histogram(eta_hat, 'NumBins', 20, 'Normalization','pdf');
grid on
[~, eta0]=cart2polar(theta0);
xline(eta0, 'rx', 'LineWidth', 2.0);
set(findall(gcf,'-property','FontSize'),'FontSize',15);
xlabel('$\hat{{\eta}}$','interpreter','latex', 'FontWeight','bold');
%ylabel('Density','interpreter','latex', 'FontWeight','bold');
xlim([-1.25,-.1]);
ylim([0, 6.5]);
title('n=100, $\sigma^2=0.4$','interpreter','latex', 'FontWeight','bold');



subplot(2,3,3)
load('Sphere_results_n200_nsim200_p2_noiseLow_Final.mat')

% bandwidth
MSE_h_sim = zeros(200, length(h));

for i=1:200
    for l=1:length(h)
        MSE_h_sim(i,l) = acos(abs(fsiFitAll{i}.thetaHat(:,l)'*theta0))^2;
    end
end

MSE_h = mean(MSE_h_sim);

xh = find(MSE_h== min(MSE_h));

eta_hat = NaN(200,1);

for i=1:200
    [~, eth] = cart2polar(fsiFitAll{i}.thetaHat(:,xh)); 
    eta_hat(i) = eth;
end

histogram(eta_hat, 'NumBins', 20, 'Normalization','pdf');
grid on
[~, eta0]=cart2polar(theta0);
xline(eta0, 'rx', 'LineWidth', 2.0);
set(findall(gcf,'-property','FontSize'),'FontSize',15);
xlabel('$\hat{{\eta}}$','interpreter','latex', 'FontWeight','bold');
%ylabel('Density','interpreter','latex', 'FontWeight','bold');
xlim([-1.25,-.1]);
ylim([0, 6.5]);
title('n=200, $\sigma^2=0.4$','interpreter','latex', 'FontWeight','bold');



subplot(2,3,4)
load('Sphere_results_n50_nsim200_p2_noiseHigh_Final.mat')

% bandwidth
MSE_h_sim = zeros(200, length(h));

for i=1:200
    for l=1:length(h)
        MSE_h_sim(i,l) = acos(abs(fsiFitAll{i}.thetaHat(:,l)'*theta0))^2;
    end
end

MSE_h = mean(MSE_h_sim);

xh = find(MSE_h== min(MSE_h));

eta_hat = NaN(200,1);

for i=1:200
    [~, eth] = cart2polar(fsiFitAll{i}.thetaHat(:,xh)); 
    eta_hat(i) = eth;
end

histogram(eta_hat, 'NumBins', 20, 'Normalization','pdf');
grid on
[~, eta0]=cart2polar(theta0);
xline(eta0, 'rx', 'LineWidth', 2.0);
set(findall(gcf,'-property','FontSize'),'FontSize',15);
xlabel('$\hat{{\eta}}$','interpreter','latex', 'FontWeight','bold');
%ylabel('Density','interpreter','latex', 'FontWeight','bold');
xlim([-1.58,1.58]);
ylim([0, 3.5]);
title('n=50, $\sigma^2=0.8$','interpreter','latex', 'FontWeight','bold');


subplot(2,3,5)
load('Sphere_results_n100_nsim200_p2_noiseHigh_Final.mat')

% bandwidth
MSE_h_sim = zeros(200, length(h));

for i=1:200
    for l=1:length(h)
        MSE_h_sim(i,l) = acos(abs(fsiFitAll{i}.thetaHat(:,l)'*theta0))^2;
    end
end

MSE_h = mean(MSE_h_sim);

xh = find(MSE_h== min(MSE_h));

eta_hat = NaN(200,1);

for i=1:200
    [~, eth] = cart2polar(fsiFitAll{i}.thetaHat(:,xh)); 
    eta_hat(i) = eth;
end

histogram(eta_hat, 'NumBins', 20, 'Normalization','pdf');
grid on
[~, eta0]=cart2polar(theta0);
xline(eta0, 'rx', 'LineWidth', 2.0);
set(findall(gcf,'-property','FontSize'),'FontSize',15);
xlabel('$\hat{{\eta}}$','interpreter','latex', 'FontWeight','bold');
%ylabel('Density','interpreter','latex', 'FontWeight','bold');
xlim([-1.58,1.58]);
ylim([0, 3.5]);
title('n=100, $\sigma^2=0.8$','interpreter','latex', 'FontWeight','bold');



subplot(2,3,6)
load('Sphere_results_n200_nsim200_p2_noiseHigh_Final.mat')

% bandwidth
MSE_h_sim = zeros(200, length(h));

for i=1:200
    for l=1:length(h)
        MSE_h_sim(i,l) = acos(abs(fsiFitAll{i}.thetaHat(:,l)'*theta0))^2;
    end
end

MSE_h = mean(MSE_h_sim);

xh = find(MSE_h== min(MSE_h));

eta_hat = NaN(200,1);

for i=1:200
    [~, eth] = cart2polar(fsiFitAll{i}.thetaHat(:,xh)); 
    eta_hat(i) = eth;
end

histogram(eta_hat, 'NumBins', 20, 'Normalization','pdf');
grid on
[~, eta0]=cart2polar(theta0);
xline(eta0, 'rx', 'LineWidth', 2.0);
set(findall(gcf,'-property','FontSize'),'FontSize',15);
xlabel('$\hat{{\eta}}$','interpreter','latex', 'FontWeight','bold');
%ylabel('Density','interpreter','latex', 'FontWeight','bold');
xlim([-1.58,1.58]);
ylim([0, 3.5]);
title('n=200, $\sigma^2=0.8$','interpreter','latex', 'FontWeight','bold');


%sgtitle('Distribution of estimates $\hat{\eta}$ for p=2','interpreter','latex', 'FontWeight','bold') 

han=axes(fig,'visible','off'); 
han.Title.Visible='on';
han.XLabel.Visible='on';
han.YLabel.Visible='on';

set(findall(gcf,'-property','FontSize'),'FontSize',15);
%xlabel('${\eta}$','interpreter','latex', 'FontWeight','bold');
ylabel('Density','interpreter','latex', 'FontWeight','bold');
%title(han,'yourTitle');

saveas(gcf,'eta_hat_histogram_p2','epsc')
