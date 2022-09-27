
%% Codes to generate figure 5 in the document

clear all

% Parameter values
nsim = 200;
nVal = [50 100 200];
pVal = [2 5 10];


%% Compute numeric summaries for high noise

hopt_MSEE_hn = zeros(length(nVal), length(pVal)); % Store best MSEE bandwidth indices for high noise
MSEE_hn = zeros(nsim, length(nVal)*length(pVal));

for pInd = 1:length(pVal)
    for nInd = 1:length(nVal)
        
        cur = (pInd - 1)*length(pVal) + nInd; % indexes the column of MSEE_hn we want to populate
        
        flnm = ['FinalSimResults/Sphere_results_n' num2str(nVal(nInd)) '_nsim' num2str(nsim) '_p' num2str(pVal(pInd)) '_noiseHigh_Final.mat'];
        load(flnm);
        
        MSEE_h_sim = zeros(200, length(h));
        for i=1:nsim
            for l=1:length(h)
                MSEE_h_sim(i,l) = mean(arrayfun(@(j) (acos(fsiFitAll{i}.Yout(:,j,l)'*mreg{i}(:,j)))^2, 1:n));
            end
        end
        
        MSEE_h = mean(MSEE_h_sim);
        xh = find(MSEE_h== min(MSEE_h));
        
        MSEE_hn(:, cur) = MSEE_h_sim(:, xh);
        hopt_MSEE_hn(nInd, pInd) = xh;
        
    end
end


% Boxplots for MSEE, high noise

    % Define global control variables for all boxplots
    group = [ones(1,200),2*ones(1,200),3*ones(1,200), 4*ones(1,200),5*ones(1,200) ,6*ones(1,200),7*ones(1,200),8*ones(1,200),9*ones(1,200)];
    positions = [1 1.25 1.50 2 2.25 2.50 3 3.25 3.50];
    color = ['c', 'y','b', 'c', 'y','b', 'c', 'y','b'];

    % MSEE boxplot
    boxplot(MSEE_hn, group, 'positions', positions);
    ylabel('${MSEE_{\oplus}}$','interpreter','latex', 'FontWeight','bold');
    set(gca,'xtick',[mean(positions(1:3)) mean(positions(4:6)) mean(positions(7:9))])
    set(gca,'xticklabel',{'p=2','p=5','p=10'})

    h = findobj(gca,'Tag','Box');
    for j=1:length(h)
       patch(get(h(j),'XData'),get(h(j),'YData'),color(j),'FaceAlpha',.5);
    end

    c = get(gca, 'Children');

    legend(c(1:3), 'n=50', 'n=100', 'n=200','Location','northwest');
    set(findall(gca,'-property','FontSize'),'FontSize',17)

    grid on
    saveas(gca,'msee_boxplots_hn','epsc')

    % Log MSEE boxplot
    boxplot(log(MSEE_hn), group, 'positions', positions);
    ylabel('${log(MSEE_{\oplus})}$','interpreter','latex', 'FontWeight','bold');
    set(gca,'xtick',[mean(positions(1:3)) mean(positions(4:6)) mean(positions(7:9))])
    set(gca,'xticklabel',{'p=2','p=5','p=10'})

    h = findobj(gca,'Tag','Box');
    for j=1:length(h)
       patch(get(h(j),'XData'),get(h(j),'YData'),color(j),'FaceAlpha',.5);
    end

    c = get(gca, 'Children');

    legend(c(1:3), 'n=50', 'n=100', 'n=200','Location','southeast');
    set(findall(gca,'-property','FontSize'),'FontSize',17)

    grid on
    saveas(gca,'log_msee_boxplots_hn','epsc')


%% Compute numeric summaries for low noise

hopt_MSEE_ln = zeros(length(nVal), length(pVal)); % Store best MSEE bandwidth indices for low noise
MSEE_ln = zeros(nsim, length(nVal)*length(pVal));

for pInd = 1:length(pVal)
    for nInd = 1:length(nVal)
        
        cur = (pInd - 1)*length(pVal) + nInd; % indexes the column of MSEE_hn we want to populate
        
        flnm = ['FinalSimResults/Sphere_results_n' num2str(nVal(nInd)) '_nsim' num2str(nsim) '_p' num2str(pVal(pInd)) '_noiseLow_Final.mat'];
        load(flnm);
        
        MSEE_h_sim = zeros(200, length(h));
        for i=1:nsim
            for l=1:length(h)
                MSEE_h_sim(i,l) = mean(arrayfun(@(j) (acos(fsiFitAll{i}.Yout(:,j,l)'*mreg{i}(:,j)))^2, 1:n));
            end
        end
        
        MSEE_h = mean(MSEE_h_sim);
        xh = find(MSEE_h== min(MSEE_h));
        
        MSEE_ln(:, cur) = MSEE_h_sim(:, xh);    
        hopt_MSEE_ln(nInd, pInd) = xh;
        
    end
end


% Boxplots for MSEE, high noise

    % Define global control variables for all boxplots
    group = [ones(1,200),2*ones(1,200),3*ones(1,200), 4*ones(1,200),5*ones(1,200) ,6*ones(1,200),7*ones(1,200),8*ones(1,200),9*ones(1,200)];
    positions = [1 1.25 1.50 2 2.25 2.50 3 3.25 3.50];
    color = ['c', 'y','b', 'c', 'y','b', 'c', 'y','b'];

    % MSEE boxplot
    boxplot(MSEE_ln, group, 'positions', positions);
    ylabel('${MSEE_{\oplus}}$','interpreter','latex', 'FontWeight','bold');
    set(gca,'xtick',[mean(positions(1:3)) mean(positions(4:6)) mean(positions(7:9))])
    set(gca,'xticklabel',{'p=2','p=5','p=10'})

    h = findobj(gca,'Tag','Box');
    for j=1:length(h)
       patch(get(h(j),'XData'),get(h(j),'YData'),color(j),'FaceAlpha',.5);
    end

    c = get(gca, 'Children');

    legend(c(1:3), 'n=50', 'n=100', 'n=200','Location','northwest');
    set(findall(gca,'-property','FontSize'),'FontSize',17)

    grid on
    saveas(gca,'msee_boxplots_ln','epsc')

    % Log MSEE boxplot
    boxplot(log(MSEE_ln), group, 'positions', positions);
    ylabel('${log(MSEE_{\oplus})}$','interpreter','latex', 'FontWeight','bold');
    set(gca,'xtick',[mean(positions(1:3)) mean(positions(4:6)) mean(positions(7:9))])
    set(gca,'xticklabel',{'p=2','p=5','p=10'})

    h = findobj(gca,'Tag','Box');
    for j=1:length(h)
       patch(get(h(j),'XData'),get(h(j),'YData'),color(j),'FaceAlpha',.5);
    end

    c = get(gca, 'Children');

    legend(c(1:3), 'n=50', 'n=100', 'n=200','Location','southeast');
    set(findall(gca,'-property','FontSize'),'FontSize',17)

    grid on
    saveas(gca,'log_msee_boxplots_ln','epsc')

    % Save optimal MSEE bandwidths
    
    save hopt_MSEE_indices.mat hopt_MSEE_ln hopt_MSEE_hn
