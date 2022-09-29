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
for cur = 1:size(eta_sim, 2)
    subplot(2,3,cur)
    histogram(eta_sim(:,cur), 'NumBins', 20, 'Normalization','pdf');
    grid on
    [~, eta0]=cart2polar(theta0);
    xline(eta0, 'rx', 'LineWidth', 2.0);
    set(findall(gcf,'-property','FontSize'),'FontSize',15);
    xlabel('$\hat{{\eta}}$','interpreter','latex', 'FontWeight','bold');
    if cur <= 3
        xlim([-1.25,-.1]);
        ylim([0, 6.5]);
        title(['n=', num2str(nVal(nInd)), ', $\sigma^2=0.4$'],'interpreter','latex', 'FontWeight','bold');
    else
        xlim([-1.58,1.58]);
        ylim([0, 3.5]);
        title('n=', num2str(nVal(nInd)), ', $\sigma^2=0.8$','interpreter','latex', 'FontWeight','bold');
    end
end 

han=axes(fig,'visible','off'); 
han.Title.Visible='on';
han.XLabel.Visible='on';
han.YLabel.Visible='on';

set(findall(gcf,'-property','FontSize'),'FontSize',15);
ylabel('Density','interpreter','latex', 'FontWeight','bold');

saveas(gcf,'eta_hat_histogram_p2','epsc')
