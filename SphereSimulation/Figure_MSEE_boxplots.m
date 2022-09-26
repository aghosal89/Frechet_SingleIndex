
%% Codes to generate figure 5 in the document

clear
% Plot of Boxplots for MSEE, high noise

load Sphere_results_n50_nsim200_p2_noiseHigh_Final.mat

MSEE_h_sim = zeros(200, length(h));

for i=1:200
    for l=1:length(h)
        MSEE_h_sim(i,l) = mean(arrayfun(@(j) (acos(fsiFitAll{i}.Yout(:,j,l)'*mreg{i}(:,j)))^2, 1:n));
    end
end

MSEE_h = mean(MSEE_h_sim);

xh = find(MSEE_h== min(MSEE_h));

MSEE_2_50_hn= zeros(200,1);

for i=1:200
    MSEE_2_50_hn(i) = mean(arrayfun(@(j) (acos(fsiFitAll{i}.Yout(:,j,xh)'* mreg{i}(:,j)))^2, 1:n));
end



load Sphere_results_n100_nsim200_p2_noiseHigh_Final.mat

MSEE_h_sim = zeros(200, length(h));

for i=1:200
    for l=1:length(h)
        MSEE_h_sim(i,l) = mean(arrayfun(@(j) (acos(fsiFitAll{i}.Yout(:,j,l)'*mreg{i}(:,j)))^2, 1:n));
    end
end

MSEE_h = mean(MSEE_h_sim);

xh = find(MSEE_h== min(MSEE_h));

MSEE_2_100_hn= zeros(200,1);

for i=1:200
    MSEE_2_100_hn(i) = mean(arrayfun(@(j) (acos(fsiFitAll{i}.Yout(:,j,xh)'* mreg{i}(:,j)))^2, 1:n));
end


load Sphere_results_n200_nsim200_p2_noiseHigh_Final.mat

MSEE_h_sim = zeros(200, length(h));

for i=1:200
    for l=1:length(h)
        MSEE_h_sim(i,l) = mean(arrayfun(@(j) (acos(fsiFitAll{i}.Yout(:,j,l)'*mreg{i}(:,j)))^2, 1:n));
    end
end

MSEE_h = mean(MSEE_h_sim);

xh = find(MSEE_h== min(MSEE_h));

MSEE_2_200_hn= zeros(200,1);

for i=1:200
    MSEE_2_200_hn(i) = mean(arrayfun(@(j) (acos(fsiFitAll{i}.Yout(:,j,xh)'* mreg{i}(:,j)))^2, 1:n));
end



load Sphere_results_n50_nsim200_p5_noiseHigh_Final.mat

MSEE_h_sim = zeros(200, length(h));

for i=1:200
    for l=1:length(h)
        MSEE_h_sim(i,l) = mean(arrayfun(@(j) (acos(fsiFitAll{i}.Yout(:,j,l)'*mreg{i}(:,j)))^2, 1:n));
    end
end

MSEE_h = mean(MSEE_h_sim);

xh = find(MSEE_h== min(MSEE_h));

MSEE_5_50_hn= zeros(200,1);

for i=1:200
    MSEE_5_50_hn(i) = mean(arrayfun(@(j) (acos(fsiFitAll{i}.Yout(:,j,xh)'* mreg{i}(:,j)))^2, 1:n));
end



load Sphere_results_n100_nsim200_p5_noiseHigh_Final.mat

MSEE_h_sim = zeros(200, length(h));

for i=1:200
    for l=1:length(h)
        MSEE_h_sim(i,l) = mean(arrayfun(@(j) (acos(fsiFitAll{i}.Yout(:,j,l)'*mreg{i}(:,j)))^2, 1:n));
    end
end

MSEE_h = mean(MSEE_h_sim);

xh = find(MSEE_h== min(MSEE_h));

MSEE_5_100_hn= zeros(200,1);

for i=1:200
    MSEE_5_100_hn(i) = mean(arrayfun(@(j) (acos(fsiFitAll{i}.Yout(:,j,xh)'* mreg{i}(:,j)))^2, 1:n));
end



load Sphere_results_n200_nsim200_p5_noiseHigh_Final.mat

MSEE_h_sim = zeros(200, length(h));

for i=1:200
    for l=1:length(h)
        MSEE_h_sim(i,l) = mean(arrayfun(@(j) (acos(fsiFitAll{i}.Yout(:,j,l)'*mreg{i}(:,j)))^2, 1:n));
    end
end

MSEE_h = mean(MSEE_h_sim);

xh = find(MSEE_h== min(MSEE_h));

MSEE_5_200_hn= zeros(200,1);

for i=1:200
    MSEE_5_200_hn(i) = mean(arrayfun(@(j) (acos(fsiFitAll{i}.Yout(:,j,xh)'* mreg{i}(:,j)))^2, 1:n));
end



load Sphere_results_n50_nsim200_p10_noiseHigh_Final.mat

MSEE_h_sim = zeros(200, length(h));

for i=1:200
    for l=1:length(h)
        MSEE_h_sim(i,l) = mean(arrayfun(@(j) (acos(fsiFitAll{i}.Yout(:,j,l)'*mreg{i}(:,j)))^2, 1:n));
    end
end

MSEE_h = mean(MSEE_h_sim);

xh = find(MSEE_h== min(MSEE_h));

MSEE_10_50_hn= zeros(200,1);

for i=1:200
    MSEE_10_50_hn(i) = mean(arrayfun(@(j) (acos(fsiFitAll{i}.Yout(:,j,xh)'* mreg{i}(:,j)))^2, 1:n));
end



load Sphere_results_n100_nsim200_p10_noiseHigh_Final.mat

MSEE_h_sim = zeros(200, length(h));

for i=1:200
    for l=1:length(h)
        MSEE_h_sim(i,l) = mean(arrayfun(@(j) (acos(fsiFitAll{i}.Yout(:,j,l)'*mreg{i}(:,j)))^2, 1:n));
    end
end

MSEE_h = mean(MSEE_h_sim);

xh = find(MSEE_h== min(MSEE_h));

MSEE_10_100_hn= zeros(200,1);

for i=1:200
    MSEE_10_100_hn(i) = mean(arrayfun(@(j) (acos(fsiFitAll{i}.Yout(:,j,xh)'* mreg{i}(:,j)))^2, 1:n));
end



load Sphere_results_n200_nsim200_p10_noiseHigh_Final.mat

MSEE_h_sim = zeros(200, length(h));

for i=1:200
    for l=1:length(h)
        MSEE_h_sim(i,l) = mean(arrayfun(@(j) (acos(fsiFitAll{i}.Yout(:,j,l)'*mreg{i}(:,j)))^2, 1:n));
    end
end

MSEE_h = mean(MSEE_h_sim);

xh = find(MSEE_h== min(MSEE_h));

MSEE_10_200_hn= zeros(200,1);

for i=1:200
    MSEE_10_200_hn(i) = mean(arrayfun(@(j) (acos(fsiFitAll{i}.Yout(:,j,xh)'* mreg{i}(:,j)))^2, 1:n));
end



% Creating multiple boxplots for comparing MSEE distributions

x = [MSEE_2_50_hn, MSEE_2_100_hn, MSEE_2_200_hn,MSEE_5_50_hn, MSEE_5_100_hn,MSEE_5_200_hn,MSEE_10_50_hn,MSEE_10_100_hn,MSEE_10_200_hn];
group = [ones(1,200),2*ones(1,200),3*ones(1,200), 4*ones(1,200),5*ones(1,200) ,6*ones(1,200),7*ones(1,200),8*ones(1,200),9*ones(1,200)];
positions = [1 1.25 1.50 2 2.25 2.50 3 3.25 3.50];

boxplot(x, group, 'positions', positions);
ylabel('${MSEE_{\oplus}}$','interpreter','latex', 'FontWeight','bold');
set(gca,'xtick',[mean(positions(1:3)) mean(positions(4:6)) mean(positions(7:9))])
set(gca,'xticklabel',{'p=2','p=5','p=10'})

color = ['c', 'y','b', 'c', 'y','b', 'c', 'y','b'];
h = findobj(gca,'Tag','Box');
for j=1:length(h)
   patch(get(h(j),'XData'),get(h(j),'YData'),color(j),'FaceAlpha',.5);
end

c = get(gca, 'Children');

hleg1 = legend(c(1:3), 'n=50', 'n=100', 'n=200','Location','northwest');
set(findall(gca,'-property','FontSize'),'FontSize',17)

grid on
saveas(gca,'msee_boxplots_hn','epsc')



% Creating multiple boxplots for comparing MSEE distributions

x = [log(MSEE_2_50_hn), log(MSEE_2_100_hn), log(MSEE_2_200_hn), log(MSEE_5_50_hn), log(MSEE_5_100_hn), log(MSEE_5_200_hn), log(MSEE_10_50_hn), log(MSEE_10_100_hn), log(MSEE_10_200_hn)];
group = [ones(1,200),2*ones(1,200),3*ones(1,200), 4*ones(1,200),5*ones(1,200) ,6*ones(1,200),7*ones(1,200),8*ones(1,200),9*ones(1,200)];
positions = [1 1.25 1.50 2 2.25 2.50 3 3.25 3.50];

boxplot(x, group, 'positions', positions);
ylabel('${log(MSEE_{\oplus})}$','interpreter','latex', 'FontWeight','bold');
set(gca,'xtick',[mean(positions(1:3)) mean(positions(4:6)) mean(positions(7:9))])
set(gca,'xticklabel',{'p=2','p=5','p=10'})

color = ['c', 'y','b', 'c', 'y','b', 'c', 'y','b'];
h = findobj(gca,'Tag','Box');
for j=1:length(h)
   patch(get(h(j),'XData'),get(h(j),'YData'),color(j),'FaceAlpha',.5);
end

c = get(gca, 'Children');

hleg1 = legend(c(1:3), 'n=50', 'n=100', 'n=200','Location','southeast');
set(findall(gca,'-property','FontSize'),'FontSize',17)

grid on
saveas(gca,'log_msee_boxplots_hn','epsc')










% For low noise

load Sphere_results_n50_nsim200_p2_noiseLow_Final.mat

MSEE_h_sim = zeros(200, length(h));

for i=1:200
    for l=1:length(h)
        MSEE_h_sim(i,l) = mean(arrayfun(@(j) (acos(fsiFitAll{i}.Yout(:,j,l)'*mreg{i}(:,j)))^2, 1:n));
    end
end

MSEE_h = mean(MSEE_h_sim);

xh = find(MSEE_h== min(MSEE_h));

MSEE_2_50_ln= zeros(200,1);

for i=1:200
    MSEE_2_50_ln(i) = mean(arrayfun(@(j) (acos(fsiFitAll{i}.Yout(:,j,xh)'* mreg{i}(:,j)))^2, 1:n));
end



load Sphere_results_n100_nsim200_p2_noiseLow_Final.mat

MSEE_h_sim = zeros(200, length(h));

for i=1:200
    for l=1:length(h)
        MSEE_h_sim(i,l) = mean(arrayfun(@(j) (acos(fsiFitAll{i}.Yout(:,j,l)'*mreg{i}(:,j)))^2, 1:n));
    end
end

MSEE_h = mean(MSEE_h_sim);

xh = find(MSEE_h== min(MSEE_h));

MSEE_2_100_ln= zeros(200,1);

for i=1:200
    MSEE_2_100_ln(i) = mean(arrayfun(@(j) (acos(fsiFitAll{i}.Yout(:,j,xh)'* mreg{i}(:,j)))^2, 1:n));
end


load Sphere_results_n200_nsim200_p2_noiseLow_Final.mat

MSEE_h_sim = zeros(200, length(h));

for i=1:200
    for l=1:length(h)
        MSEE_h_sim(i,l) = mean(arrayfun(@(j) (acos(fsiFitAll{i}.Yout(:,j,l)'*mreg{i}(:,j)))^2, 1:n));
    end
end

MSEE_h = mean(MSEE_h_sim);

xh = find(MSEE_h== min(MSEE_h));

MSEE_2_200_ln= zeros(200,1);

for i=1:200
    MSEE_2_200_ln(i) = mean(arrayfun(@(j) (acos(fsiFitAll{i}.Yout(:,j,xh)'* mreg{i}(:,j)))^2, 1:n));
end



load Sphere_results_n50_nsim200_p5_noiseLow_Final.mat

MSEE_h_sim = zeros(200, length(h));

for i=1:200
    for l=1:length(h)
        MSEE_h_sim(i,l) = mean(arrayfun(@(j) (acos(fsiFitAll{i}.Yout(:,j,l)'*mreg{i}(:,j)))^2, 1:n));
    end
end

MSEE_h = mean(MSEE_h_sim);

xh = find(MSEE_h== min(MSEE_h));

MSEE_5_50_ln= zeros(200,1);

for i=1:200
    MSEE_5_50_ln(i) = mean(arrayfun(@(j) (acos(fsiFitAll{i}.Yout(:,j,xh)'* mreg{i}(:,j)))^2, 1:n));
end



load Sphere_results_n100_nsim200_p5_noiseLow_Final.mat

MSEE_h_sim = zeros(200, length(h));

for i=1:200
    for l=1:length(h)
        MSEE_h_sim(i,l) = mean(arrayfun(@(j) (acos(fsiFitAll{i}.Yout(:,j,l)'*mreg{i}(:,j)))^2, 1:n));
    end
end

MSEE_h = mean(MSEE_h_sim);

xh = find(MSEE_h== min(MSEE_h));

MSEE_5_100_ln= zeros(200,1);

for i=1:200
    MSEE_5_100_ln(i) = mean(arrayfun(@(j) (acos(fsiFitAll{i}.Yout(:,j,xh)'* mreg{i}(:,j)))^2, 1:n));
end



load Sphere_results_n200_nsim200_p5_noiseLow_Final.mat

MSEE_h_sim = zeros(200, length(h));

for i=1:200
    for l=1:length(h)
        MSEE_h_sim(i,l) = mean(arrayfun(@(j) (acos(fsiFitAll{i}.Yout(:,j,l)'*mreg{i}(:,j)))^2, 1:n));
    end
end

MSEE_h = mean(MSEE_h_sim);

xh = find(MSEE_h== min(MSEE_h));

MSEE_5_200_ln= zeros(200,1);

for i=1:200
    MSEE_5_200_ln(i) = mean(arrayfun(@(j) (acos(fsiFitAll{i}.Yout(:,j,xh)'* mreg{i}(:,j)))^2, 1:n));
end



load Sphere_results_n50_nsim200_p10_noiseLow_Final.mat

MSEE_h_sim = zeros(200, length(h));

for i=1:200
    for l=1:length(h)
        MSEE_h_sim(i,l) = mean(arrayfun(@(j) (acos(fsiFitAll{i}.Yout(:,j,l)'*mreg{i}(:,j)))^2, 1:n));
    end
end

MSEE_h = mean(MSEE_h_sim);

xh = find(MSEE_h== min(MSEE_h));

MSEE_10_50_ln= zeros(200,1);

for i=1:200
    MSEE_10_50_ln(i) = mean(arrayfun(@(j) (acos(fsiFitAll{i}.Yout(:,j,xh)'* mreg{i}(:,j)))^2, 1:n));
end



load Sphere_results_n100_nsim200_p10_noiseLow_Final.mat

MSEE_h_sim = zeros(200, length(h));

for i=1:200
    for l=1:length(h)
        MSEE_h_sim(i,l) = mean(arrayfun(@(j) (acos(fsiFitAll{i}.Yout(:,j,l)'*mreg{i}(:,j)))^2, 1:n));
    end
end

MSEE_h = mean(MSEE_h_sim);

xh = find(MSEE_h== min(MSEE_h));

MSEE_10_100_ln= zeros(200,1);

for i=1:200
    MSEE_10_100_ln(i) = mean(arrayfun(@(j) (acos(fsiFitAll{i}.Yout(:,j,xh)'* mreg{i}(:,j)))^2, 1:n));
end



load Sphere_results_n200_nsim200_p10_noiseLow_Final.mat

MSEE_h_sim = zeros(200, length(h));

for i=1:180
    for l=1:length(h)
        MSEE_h_sim(i,l) = mean(arrayfun(@(j) (acos(fsiFitAll{i}.Yout(:,j,l)'*mreg{i}(:,j)))^2, 1:n));
    end
end

MSEE_h = mean(MSEE_h_sim);

xh = find(MSEE_h== min(MSEE_h));

MSEE_10_200_ln= zeros(200,1);

for i=1:200
    MSEE_10_200_ln(i) = mean(arrayfun(@(j) (acos(fsiFitAll{i}.Yout(:,j,xh)'* mreg{i}(:,j)))^2, 1:n));
end



% Creating multiple boxplots for comparing MSEE distributions

x = [MSEE_2_50_ln, MSEE_2_100_ln, MSEE_2_200_ln,MSEE_5_50_ln, MSEE_5_100_ln,MSEE_5_200_ln,MSEE_10_50_ln,MSEE_10_100_ln,MSEE_10_200_ln];
group = [ones(1,200),2*ones(1,200),3*ones(1,200), 4*ones(1,200),5*ones(1,200) ,6*ones(1,200),7*ones(1,200),8*ones(1,200),9*ones(1,200)];
positions = [1 1.25 1.50 2 2.25 2.50 3 3.25 3.50];

boxplot(x, group, 'positions', positions);
ylabel('${MSEE_{\oplus}}$','interpreter','latex', 'FontWeight','bold');

set(gca,'xtick',[mean(positions(1:3)) mean(positions(4:6)) mean(positions(7:9))])
set(gca,'xticklabel',{'p=2','p=5','p=10'})

color = ['c', 'y','b', 'c', 'y','b', 'c', 'y','b'];
h = findobj(gca,'Tag','Box');
for j=1:length(h)
   patch(get(h(j),'XData'),get(h(j),'YData'),color(j),'FaceAlpha',.5);
end

c = get(gca, 'Children');

hleg1 = legend(c(1:3), 'n=50', 'n=100', 'n=200','Location','northwest');
set(findall(gca,'-property','FontSize'),'FontSize',17)

grid on
saveas(gca,'msee_boxplots_ln','epsc')





% Creating multiple boxplots for comparing MSEE distributions

x = [log(MSEE_2_50_ln), log(MSEE_2_100_ln), log(MSEE_2_200_ln), log(MSEE_5_50_ln), log(MSEE_5_100_ln), log(MSEE_5_200_ln), log(MSEE_10_50_ln), log(MSEE_10_100_ln), log(MSEE_10_200_ln)];
group = [ones(1,200),2*ones(1,200),3*ones(1,200), 4*ones(1,200),5*ones(1,200) ,6*ones(1,200),7*ones(1,200),8*ones(1,200),9*ones(1,200)];
positions = [1 1.25 1.50 2 2.25 2.50 3 3.25 3.50];

boxplot(x, group, 'positions', positions);
ylabel('${log(MSEE_{\oplus})}$','interpreter','latex', 'FontWeight','bold');
set(gca,'xtick',[mean(positions(1:3)) mean(positions(4:6)) mean(positions(7:9))])
set(gca,'xticklabel',{'p=2','p=5','p=10'})

color = ['c', 'y','b', 'c', 'y','b', 'c', 'y','b'];
h = findobj(gca,'Tag','Box');
for j=1:length(h)
   patch(get(h(j),'XData'),get(h(j),'YData'),color(j),'FaceAlpha',.5);
end

c = get(gca, 'Children');

hleg1 = legend(c(1:3), 'n=50', 'n=100', 'n=200','Location','southeast');
set(findall(gca,'-property','FontSize'),'FontSize',17)

grid on
saveas(gca,'log_msee_boxplots_ln','epsc')

