
% addpath gridfitdir
% addpath mySphere

%% Code to create the figure 1 in the document
clear all
addpath(genpath('manopt'))

load NM_Sphere_results_n50_nsim200_p4_HN.mat

% since we have to work with a single dataset, we set 
i=1;

% following is an arbitrary parameter value considerably far from the true parameter.
b0f=[0.087, 0.3332, 0.487, -sqrt(1-0.087^2-0.3332^2-0.487^2)];

% following is the estimate of parameter from first dataset
bopt = beta_opt_s(i,:);

z_t2 = (bopt(1)*x(:,i,1)+bopt(2)*x(:,i,2))./sqrt(p);
z_t3 = (b0f(1)*x(:,i,1)+b0f(2)*x(:,i,2))./sqrt(p);

% read from the results
mreg2 = reg_curve(z_t2)'; 
mreg3 = reg_curve(z_t3)'; 

% read from the results
Y_t =Y{i};
Yhat= get_sphere_fit_FSI(Y_t, z_t2, h_opt);
Y0= get_sphere_fit_FSI(Y_t, z_t3, h_opt);

xout1 = linspace(min(z(:,i)), max(z(:,i)), 51);
xout2 = linspace(min(z_t2), max(z_t2), 51);
xout3 = linspace(min(z_t3), max(z_t3), 51);

yt1 = [sqrt(1-xout1'.^2).*cos(pi*xout1'), sqrt(1-xout1'.^2).*sin(pi*xout1'), xout1']';
yt2 = [sqrt(1-xout2'.^2).*cos(pi*xout2'), sqrt(1-xout2'.^2).*sin(pi*xout2'), xout2']';
yt3 = [sqrt(1-xout3'.^2).*cos(pi*xout3'), sqrt(1-xout3'.^2).*sin(pi*xout3'), xout3']';

m_tr1 = arrayfun(@(j) yt1(j, :), 1:3, 'UniformOutput', false);
m_tr2 = arrayfun(@(j) yt2(j, :), 1:3, 'UniformOutput', false);
m_tr3 = arrayfun(@(j) yt3(j, :), 1:3, 'UniformOutput', false);

[s1 s2 s3] = sphere(128);
lightGrey = 0.9*ones(1, 3);

set(gcf, 'PaperUnits', 'inches', 'PaperSize', [2.4, 2.6], 'PaperPosition', [0 0 2.4 2.6])
surf(s1, s2, s3, 'FaceColor', lightGrey, 'EdgeColor', 'none', 'FaceAlpha', 0.4)
view([1, 1, 1])
hold on;
plot3(m_tr1{1}, m_tr1{2}, m_tr1{3}, 'k', 'LineWidth', 0.5)
plot3(Y_t(1, :), Y_t(2, :), Y_t(3, :), 'r.', 'MarkerSize', 2.5)

hold on;

%plot3(m_tr2{1}, m_tr2{2}, m_tr2{3}, 'b', 'LineWidth', 0.5)
plot3(Yhat(1, :), Yhat(2, :), Yhat(3, :), 'b.', 'MarkerSize', 2.5)

hold on;
%plot3(m_tr3{1}, m_tr3{2}, m_tr3{3}, 'g', 'LineWidth', 0.5)
plot3(Y0(1, :), Y0(2, :), Y0(3, :), 'g.', 'MarkerSize', 2.5)

saveas(gcf,'highnoise_n50_p4','epsc')









clear all
addpath(genpath('manopt'))

load NM_Sphere_results_n50_nsim200_p4_LN.mat

% since we have to work with a single dataset, we set 
i=1;

% following is an arbitrary parameter value considerably far from the true parameter.
b0f=[0.087, 0.3332, 0.487, -sqrt(1-0.087^2-0.3332^2-0.487^2)];

% following is the estimate of parameter from first dataset
bopt = beta_opt_s(i,:);

z_t2 = (bopt(1)*x(:,i,1)+bopt(2)*x(:,i,2))./sqrt(p);
z_t3 = (b0f(1)*x(:,i,1)+b0f(2)*x(:,i,2))./sqrt(p);

% read from the results
Y_t =Y{i};
Yhat= get_sphere_fit_FSI(Y_t, z_t2, h_opt);
Y0= get_sphere_fit_FSI(Y_t, z_t3, h_opt);

xout1 = linspace(min(z(:,i)), max(z(:,i)), 51);
xout2 = linspace(min(z_t2), max(z_t2), 51);
xout3 = linspace(min(z_t3), max(z_t3), 51);

yt1 = [sqrt(1-xout1'.^2).*cos(pi*xout1'), sqrt(1-xout1'.^2).*sin(pi*xout1'), xout1']';
yt2 = [sqrt(1-xout2'.^2).*cos(pi*xout2'), sqrt(1-xout2'.^2).*sin(pi*xout2'), xout2']';
yt3 = [sqrt(1-xout3'.^2).*cos(pi*xout3'), sqrt(1-xout3'.^2).*sin(pi*xout3'), xout3']';

m_tr1 = arrayfun(@(j) yt1(j, :), 1:3, 'UniformOutput', false);
m_tr2 = arrayfun(@(j) yt2(j, :), 1:3, 'UniformOutput', false);
m_tr3 = arrayfun(@(j) yt3(j, :), 1:3, 'UniformOutput', false);

[s1 s2 s3] = sphere(128);
lightGrey = 0.9*ones(1, 3);

set(gcf, 'PaperUnits', 'inches', 'PaperSize', [2.4, 2.6], 'PaperPosition', [0 0 2.4 2.6])
surf(s1, s2, s3, 'FaceColor', lightGrey, 'EdgeColor', 'none', 'FaceAlpha', 0.4)
view([1, 1, 1])
hold on;
plot3(m_tr1{1}, m_tr1{2}, m_tr1{3}, 'k', 'LineWidth', 0.5)
plot3(Y_t(1, :), Y_t(2, :), Y_t(3, :), 'r.', 'MarkerSize', 2.5)

hold on;

%plot3(m_tr2{1}, m_tr2{2}, m_tr2{3}, 'b', 'LineWidth', 0.5)
plot3(Yhat(1, :), Yhat(2, :), Yhat(3, :), 'b.', 'MarkerSize', 2.5)

hold on;
%plot3(m_tr3{1}, m_tr3{2}, m_tr3{3}, 'g', 'LineWidth', 0.5)
plot3(Y0(1, :), Y0(2, :), Y0(3, :), 'g.', 'MarkerSize', 2.5)

saveas(gcf,'lownoise_n50_p4','epsc')
















clear all
addpath(genpath('manopt'))

load NM_Sphere_results_n100_nsim200_p4_HN.mat

% since we have to work with a single dataset, we set 
i=1;

% following is an arbitrary parameter value considerably far from the true parameter.
b0f=[0.087, 0.3332, 0.487, -sqrt(1-0.087^2-0.3332^2-0.487^2)];

% following is the estimate of parameter from first dataset
bopt = beta_opt_s(i,:);

z_t2 = (bopt(1)*x(:,i,1)+bopt(2)*x(:,i,2))./sqrt(p);
z_t3 = (b0f(1)*x(:,i,1)+b0f(2)*x(:,i,2))./sqrt(p);

% read from the results
Y_t =Y{i};
Yhat= get_sphere_fit_FSI(Y_t, z_t2, h_opt);
Y0= get_sphere_fit_FSI(Y_t, z_t3, h_opt);

xout1 = linspace(min(z(:,i)), max(z(:,i)), 51);
xout2 = linspace(min(z_t2), max(z_t2), 51);
xout3 = linspace(min(z_t3), max(z_t3), 51);

yt1 = [sqrt(1-xout1'.^2).*cos(pi*xout1'), sqrt(1-xout1'.^2).*sin(pi*xout1'), xout1']';
yt2 = [sqrt(1-xout2'.^2).*cos(pi*xout2'), sqrt(1-xout2'.^2).*sin(pi*xout2'), xout2']';
yt3 = [sqrt(1-xout3'.^2).*cos(pi*xout3'), sqrt(1-xout3'.^2).*sin(pi*xout3'), xout3']';

m_tr1 = arrayfun(@(j) yt1(j, :), 1:3, 'UniformOutput', false);
m_tr2 = arrayfun(@(j) yt2(j, :), 1:3, 'UniformOutput', false);
m_tr3 = arrayfun(@(j) yt3(j, :), 1:3, 'UniformOutput', false);

[s1 s2 s3] = sphere(128);
lightGrey = 0.9*ones(1, 3);

set(gcf, 'PaperUnits', 'inches', 'PaperSize', [2.4, 2.6], 'PaperPosition', [0 0 2.4 2.6])
surf(s1, s2, s3, 'FaceColor', lightGrey, 'EdgeColor', 'none', 'FaceAlpha', 0.4)
view([1, 1, 1])
hold on;
plot3(m_tr1{1}, m_tr1{2}, m_tr1{3}, 'k', 'LineWidth', 0.5)
plot3(Y_t(1, :), Y_t(2, :), Y_t(3, :), 'r.', 'MarkerSize', 2.5)

hold on;

%plot3(m_tr2{1}, m_tr2{2}, m_tr2{3}, 'b', 'LineWidth', 0.5)
plot3(Yhat(1, :), Yhat(2, :), Yhat(3, :), 'b.', 'MarkerSize', 2.5)

hold on;
%plot3(m_tr3{1}, m_tr3{2}, m_tr3{3}, 'g', 'LineWidth', 0.5)
plot3(Y0(1, :), Y0(2, :), Y0(3, :), 'g.', 'MarkerSize', 2.5)

saveas(gcf,'highnoise_n100_p4','epsc')














clear all
addpath(genpath('manopt'))

load NM_Sphere_results_n100_nsim200_p4_LN.mat

% since we have to work with a single dataset, we set 
i=1;

% following is an arbitrary parameter value considerably far from the true parameter.
b0f=[0.087, 0.3332, 0.487, -sqrt(1-0.087^2-0.3332^2-0.487^2)];

% following is the estimate of parameter from first dataset
bopt = beta_opt_s(i,:);

z_t2 = (bopt(1)*x(:,i,1)+bopt(2)*x(:,i,2))./sqrt(p);
z_t3 = (b0f(1)*x(:,i,1)+b0f(2)*x(:,i,2))./sqrt(p);

% read from the results
Y_t =Y{i};
Yhat= get_sphere_fit_FSI(Y_t, z_t2, h_opt);
Y0= get_sphere_fit_FSI(Y_t, z_t3, h_opt);

xout1 = linspace(min(z(:,i)), max(z(:,i)), 51);
xout2 = linspace(min(z_t2), max(z_t2), 51);
xout3 = linspace(min(z_t3), max(z_t3), 51);

yt1 = [sqrt(1-xout1'.^2).*cos(pi*xout1'), sqrt(1-xout1'.^2).*sin(pi*xout1'), xout1']';
yt2 = [sqrt(1-xout2'.^2).*cos(pi*xout2'), sqrt(1-xout2'.^2).*sin(pi*xout2'), xout2']';
yt3 = [sqrt(1-xout3'.^2).*cos(pi*xout3'), sqrt(1-xout3'.^2).*sin(pi*xout3'), xout3']';

m_tr1 = arrayfun(@(j) yt1(j, :), 1:3, 'UniformOutput', false);
m_tr2 = arrayfun(@(j) yt2(j, :), 1:3, 'UniformOutput', false);
m_tr3 = arrayfun(@(j) yt3(j, :), 1:3, 'UniformOutput', false);

[s1 s2 s3] = sphere(128);
lightGrey = 0.9*ones(1, 3);

set(gcf, 'PaperUnits', 'inches', 'PaperSize', [2.4, 2.6], 'PaperPosition', [0 0 2.4 2.6])
surf(s1, s2, s3, 'FaceColor', lightGrey, 'EdgeColor', 'none', 'FaceAlpha', 0.4)
view([1, 1, 1])
hold on;
plot3(m_tr1{1}, m_tr1{2}, m_tr1{3}, 'k', 'LineWidth', 0.5)
plot3(Y_t(1, :), Y_t(2, :), Y_t(3, :), 'r.', 'MarkerSize', 2.5)

hold on;

%plot3(m_tr2{1}, m_tr2{2}, m_tr2{3}, 'b', 'LineWidth', 0.5)
plot3(Yhat(1, :), Yhat(2, :), Yhat(3, :), 'b.', 'MarkerSize', 2.5)

hold on;
%plot3(m_tr3{1}, m_tr3{2}, m_tr3{3}, 'g', 'LineWidth', 0.5)
plot3(Y0(1, :), Y0(2, :), Y0(3, :), 'g.', 'MarkerSize', 2.5)

saveas(gcf,'lownoise_n100_p4','epsc')




% n200 high noise

clear all
addpath(genpath('manopt'))

load NM_Sphere_results_n200_nsim200_p4_HN.mat

% since we have to work with a single dataset, we set 
i=1;

% following is an arbitrary parameter value considerably far from the true parameter.
b0f=[0.087, 0.3332, 0.487, -sqrt(1-0.087^2-0.3332^2-0.487^2)];

% following is the estimate of parameter from first dataset
bopt = beta_opt_s(i,:);

z_t2 = (bopt(1)*x(:,i,1)+bopt(2)*x(:,i,2))./sqrt(p);
z_t3 = (b0f(1)*x(:,i,1)+b0f(2)*x(:,i,2))./sqrt(p);

% read from the results
Y_t =Y{i};
Yhat= get_sphere_fit_FSI(Y_t, z_t2, h_opt);
Y0= get_sphere_fit_FSI(Y_t, z_t3, h_opt);

xout1 = linspace(min(z(:,i)), max(z(:,i)), 51);
xout2 = linspace(min(z_t2), max(z_t2), 51);
xout3 = linspace(min(z_t3), max(z_t3), 51);

yt1 = [sqrt(1-xout1'.^2).*cos(pi*xout1'), sqrt(1-xout1'.^2).*sin(pi*xout1'), xout1']';
yt2 = [sqrt(1-xout2'.^2).*cos(pi*xout2'), sqrt(1-xout2'.^2).*sin(pi*xout2'), xout2']';
yt3 = [sqrt(1-xout3'.^2).*cos(pi*xout3'), sqrt(1-xout3'.^2).*sin(pi*xout3'), xout3']';

m_tr1 = arrayfun(@(j) yt1(j, :), 1:3, 'UniformOutput', false);
m_tr2 = arrayfun(@(j) yt2(j, :), 1:3, 'UniformOutput', false);
m_tr3 = arrayfun(@(j) yt3(j, :), 1:3, 'UniformOutput', false);

[s1 s2 s3] = sphere(128);
lightGrey = 0.9*ones(1, 3);

set(gcf, 'PaperUnits', 'inches', 'PaperSize', [2.4, 2.6], 'PaperPosition', [0 0 2.4 2.6])
surf(s1, s2, s3, 'FaceColor', lightGrey, 'EdgeColor', 'none', 'FaceAlpha', 0.4)
view([1, 1, 1])
hold on;
plot3(m_tr1{1}, m_tr1{2}, m_tr1{3}, 'k', 'LineWidth', 0.5)
plot3(Y_t(1, :), Y_t(2, :), Y_t(3, :), 'r.', 'MarkerSize', 2.5)

hold on;

%plot3(m_tr2{1}, m_tr2{2}, m_tr2{3}, 'b', 'LineWidth', 0.5)
plot3(Yhat(1, :), Yhat(2, :), Yhat(3, :), 'b.', 'MarkerSize', 2.5)

hold on;
%plot3(m_tr3{1}, m_tr3{2}, m_tr3{3}, 'g', 'LineWidth', 0.5)
plot3(Y0(1, :), Y0(2, :), Y0(3, :), 'g.', 'MarkerSize', 2.5)

saveas(gcf,'highnoise_n200_p4','epsc')










clear all
addpath(genpath('manopt'))

load NM_Sphere_results_n200_nsim200_p4_LN.mat

% since we have to work with a single dataset, we set 
i=1;

% following is an arbitrary parameter value considerably far from the true parameter.
b0f=[0.087, 0.3332, 0.487, -sqrt(1-0.087^2-0.3332^2-0.487^2)];

% following is the estimate of parameter from first dataset
bopt = beta_opt_s(i,:);

z_t2 = (bopt(1)*x(:,i,1)+bopt(2)*x(:,i,2))./sqrt(p);
z_t3 = (b0f(1)*x(:,i,1)+b0f(2)*x(:,i,2))./sqrt(p);

% read from the results
Y_t =Y{i};
Yhat= get_sphere_fit_FSI(Y_t, z_t2, h_opt);
Y0= get_sphere_fit_FSI(Y_t, z_t3, h_opt);

xout1 = linspace(min(z(:,i)), max(z(:,i)), 51);
xout2 = linspace(min(z_t2), max(z_t2), 51);
xout3 = linspace(min(z_t3), max(z_t3), 51);

yt1 = [sqrt(1-xout1'.^2).*cos(pi*xout1'), sqrt(1-xout1'.^2).*sin(pi*xout1'), xout1']';
yt2 = [sqrt(1-xout2'.^2).*cos(pi*xout2'), sqrt(1-xout2'.^2).*sin(pi*xout2'), xout2']';
yt3 = [sqrt(1-xout3'.^2).*cos(pi*xout3'), sqrt(1-xout3'.^2).*sin(pi*xout3'), xout3']';

m_tr1 = arrayfun(@(j) yt1(j, :), 1:3, 'UniformOutput', false);
m_tr2 = arrayfun(@(j) yt2(j, :), 1:3, 'UniformOutput', false);
m_tr3 = arrayfun(@(j) yt3(j, :), 1:3, 'UniformOutput', false);

[s1 s2 s3] = sphere(128);
lightGrey = 0.9*ones(1, 3);

set(gcf, 'PaperUnits', 'inches', 'PaperSize', [2.4, 2.6], 'PaperPosition', [0 0 2.4 2.6])
surf(s1, s2, s3, 'FaceColor', lightGrey, 'EdgeColor', 'none', 'FaceAlpha', 0.4)
view([1, 1, 1])
hold on;
plot3(m_tr1{1}, m_tr1{2}, m_tr1{3}, 'k', 'LineWidth', 0.5)
plot3(Y_t(1, :), Y_t(2, :), Y_t(3, :), 'r.', 'MarkerSize', 2.5)

hold on;

%plot3(m_tr2{1}, m_tr2{2}, m_tr2{3}, 'b', 'LineWidth', 0.5)
plot3(Yhat(1, :), Yhat(2, :), Yhat(3, :), 'b.', 'MarkerSize', 2.5)

hold on;
%plot3(m_tr3{1}, m_tr3{2}, m_tr3{3}, 'g', 'LineWidth', 0.5)
plot3(Y0(1, :), Y0(2, :), Y0(3, :), 'g.', 'MarkerSize', 2.5)

saveas(gcf,'lownoise_n200_p4','epsc')







%% Codes to generate figure 2 in document:
% Plotting the W_n function
% for p=2

fig= figure;
subplot(2,3,1)
load('Grid_Sphere_results_n50_nsim200_p2_LN.mat')

mspe_mat = zeros(size(beta,1), nsim, length(h));
msee_mat = zeros(size(beta,1), nsim, length(h));

% initializing sum over nsim simulations~
mspe_matsum = real(mspe{1});
msee_matsum = real(msee{1});

for i=1:nsim
    for r=1:length(beta)
        for l=1:length(h)
            mspe_mat(r, i, l) = (mspe{i}(r,l));
            msee_mat(r, i, l) = (msee{i}(r,l));
        end
    end
end

plot(f, mspe_mat(:,:,1))
grid on
xline(-.704, 'b');
xlabel('${\eta}$','interpreter','latex', 'FontWeight','bold');
%ylabel('$W_{n}^{(s)}(\theta_{\eta})$','interpreter','latex', 'FontWeight','bold');
xlim([-1.4451,1.4451]);
ylim([0, 0.84]);
title('n=50, $\sigma^2=0.2$','interpreter','latex', 'FontWeight','bold');


subplot(2,3,2)
load('Grid_Sphere_results_n100_nsim200_p2_LN.mat')

mspe_mat = zeros(size(beta,1), nsim, length(h));
msee_mat = zeros(size(beta,1), nsim, length(h));

% initializing sum over nsim simulations~
mspe_matsum = real(mspe{1});
msee_matsum = real(msee{1});

for i=1:nsim
    for r=1:length(beta)
        for l=1:length(h)
            mspe_mat(r, i, l) = (mspe{i}(r,l));
            msee_mat(r, i, l) = (msee{i}(r,l));
        end
    end
end

plot(f, mspe_mat(:,:,1))
grid on
xline(-.704,'b');
xlabel('${\eta}$','interpreter','latex', 'FontWeight','bold');
%ylabel('$W_{n}^{(s)}(\theta_{\eta})$','interpreter','latex', 'FontWeight','bold');
xlim([-1.4451,1.4451]);
ylim([0, 0.84]);
title('n=100, $\sigma^2=0.2$','interpreter','latex', 'FontWeight','bold');


subplot(2,3,3)
load('Grid_Sphere_results_n200_nsim200_p2_LN.mat')

mspe_mat = zeros(size(beta,1), nsim, length(h));
msee_mat = zeros(size(beta,1), nsim, length(h));

% initializing sum over nsim simulations~
mspe_matsum = real(mspe{1});
msee_matsum = real(msee{1});

for i=1:nsim
    for r=1:length(beta)
        for l=1:length(h)
            mspe_mat(r, i, l) = (mspe{i}(r,l));
            msee_mat(r, i, l) = (msee{i}(r,l));
        end
    end
end

plot(f, mspe_mat(:,:,1))
grid on
xline(-.704, 'b');
xlabel('${\eta}$','interpreter','latex', 'FontWeight','bold');
%ylabel('$W_{n}^{(s)}(\theta_{\eta})$','interpreter','latex', 'FontWeight','bold');
xlim([-1.4451,1.4451]);
ylim([0, 0.84]);
title('n=200, $\sigma^2=0.2$','interpreter','latex', 'FontWeight','bold');


subplot(2,3,4)

load('Grid_Sphere_results_n50_nsim200_p2_HN.mat')

mspe_mat = zeros(size(beta,1), nsim, length(h));
msee_mat = zeros(size(beta,1), nsim, length(h));

% initializing sum over nsim simulations~
mspe_matsum = real(mspe{1});
msee_matsum = real(msee{1});

for i=1:nsim
    for r=1:length(beta)
        for l=1:length(h)
            mspe_mat(r, i, l) = (mspe{i}(r,l));
            msee_mat(r, i, l) = (msee{i}(r,l));
        end
    end
end

plot(f, mspe_mat(:,:,1))
grid on
xline(-.704, 'b');
xlabel('${\eta}$','interpreter','latex', 'FontWeight','bold');
%ylabel('$W_{n}^{(s)}(\theta_{\eta})$','interpreter','latex', 'FontWeight','bold');
xlim([-1.4451,1.4451]);
ylim([0.4, 1.65]);
title('n=50, $\sigma^2=0.6$','interpreter','latex', 'FontWeight','bold');


subplot(2,3,5)
load('Grid_Sphere_results_n100_nsim200_p2_HN.mat')

mspe_mat = zeros(size(beta,1), nsim, length(h));
msee_mat = zeros(size(beta,1), nsim, length(h));

% initializing sum over nsim simulations~
mspe_matsum = real(mspe{1});
msee_matsum = real(msee{1});

for i=1:nsim
    for r=1:length(beta)
        for l=1:length(h)
            mspe_mat(r, i, l) = (mspe{i}(r,l));
            msee_mat(r, i, l) = (msee{i}(r,l));
        end
    end
end

plot(f, mspe_mat(:,:,1))
grid on
xline(-.704, 'b');
xlabel('${\eta}$','interpreter','latex', 'FontWeight','bold');
%ylabel('$W_{n}^{(s)}(\theta_{\eta})$','interpreter','latex', 'FontWeight','bold');
xlim([-1.4451,1.4451]);
ylim([0.4, 1.65]);
title('n=100, $\sigma^2=0.6$','interpreter','latex', 'FontWeight','bold');


subplot(2,3,6)
load('Grid_Sphere_results_n200_nsim200_p2_HN.mat')

mspe_mat = zeros(size(beta,1), nsim, length(h));
msee_mat = zeros(size(beta,1), nsim, length(h));

% initializing sum over nsim simulations~
mspe_matsum = real(mspe{1});
msee_matsum = real(msee{1});

for i=1:nsim
    for r=1:length(beta)
        for l=1:length(h)
            mspe_mat(r, i, l) = (mspe{i}(r,l));
            msee_mat(r, i, l) = (msee{i}(r,l));
        end
    end
end

plot(f, mspe_mat(:,:,1))
grid on
xline(-.704, 'b');
xlabel('${\eta}$','interpreter','latex', 'FontWeight','bold');
%ylabel('$W_{n}^{(s)}(\theta_{\eta})$','interpreter','latex', 'FontWeight','bold');
xlim([-1.4451,1.4451]);
ylim([0.4, 1.65]);
title('n=200, $\sigma^2=0.6$','interpreter','latex', 'FontWeight','bold');

%sgtitle('MSEE variation across grid of $\eta$ for p=2','interpreter','latex', 'FontWeight','bold') 

han=axes(fig,'visible','off'); 
han.Title.Visible='on';
han.XLabel.Visible='on';
han.YLabel.Visible='on';

set(findall(gcf,'-property','FontSize'),'FontSize',15);
%xlabel('${\eta}$','interpreter','latex', 'FontWeight','bold');
ylabel('$W_{n}^{(s)}(\theta_{\eta})$','interpreter','latex', 'FontWeight','bold');
%title(han,'yourTitle');

saveas(gcf,'Wn_variation','epsc')














%% Codes to create figure 3 in document

fig=figure;

subplot(2,3,1)
load('NM_Sphere_results_n50_nsim200_p2_LN.mat')

mspe_nm = NaN(length(h), nsim);
eta_hat = NaN(length(h), p-1, nsim);

for i=1:nsim
    for j=1:length(h)
        [~, a2] = find(mspe{i}(j,:)== min(mspe{i}(j,:))); 
        eta_hat(j,:,i) = eta_opt{i}(j, :, a2);
        mspe_nm(j,i) = mspe{i}(j, a2);
    end
end

% to find the average and standard deviation of the MSPE over simulations
mspe_nm_sum = sum(mspe_nm,2)

% obtaining the eta estimates
eta_opt_s = zeros(nsim, p-1);
for i=1:nsim
    eta_opt_s(i,:) = eta_hat(find(mspe_nm_sum == min(mspe_nm_sum(:))),:,i);
end

histogram(eta_opt_s, 'NumBins', 20, 'Normalization','pdf');
grid on
[~, eta]=cart2polar(b);
xline(eta, 'rx');
set(findall(gcf,'-property','FontSize'),'FontSize',15);
xlabel('$\hat{{\eta}}$','interpreter','latex', 'FontWeight','bold');
%ylabel('Density','interpreter','latex', 'FontWeight','bold');
xlim([-0.85,-0.58]);
ylim([0, 22]);
title('n=50, $\sigma^2=0.2$','interpreter','latex', 'FontWeight','bold');


subplot(2,3,2)
load('NM_Sphere_results_n100_nsim200_p2_LN.mat')

mspe_nm = NaN(length(h), nsim);
eta_hat = NaN(length(h), p-1, nsim);

for i=1:nsim
    for j=1:length(h)
        [~, a2] = find(mspe{i}(j,:)== min(mspe{i}(j,:))); 
        eta_hat(j,:,i) = eta_opt{i}(j, :, a2);
        mspe_nm(j,i) = mspe{i}(j, a2);
    end
end

% to find the average and standard deviation of the MSPE over simulations
mspe_nm_sum = sum(mspe_nm,2)

% obtaining the eta estimates
eta_opt_s = zeros(nsim, p-1);
for i=1:nsim
    eta_opt_s(i,:) = eta_hat(find(mspe_nm_sum == min(mspe_nm_sum(:))),:,i);
end

histogram(eta_opt_s, 'NumBins', 12, 'Normalization','pdf');
grid on
[~, eta]=cart2polar(b);
xline(eta, 'rx');
set(findall(gcf,'-property','FontSize'),'FontSize',15);
xlabel('$\hat{{\eta}}$','interpreter','latex', 'FontWeight','bold');
%ylabel('Density','interpreter','latex', 'FontWeight','bold');
xlim([-0.85,-0.58]);
ylim([0, 22]);
title('n=100, $\sigma^2=0.2$','interpreter','latex', 'FontWeight','bold');



subplot(2,3,3)
load('NM_Sphere_results_n200_nsim200_p2_LN.mat')

mspe_nm = NaN(length(h), nsim);
eta_hat = NaN(length(h), p-1, nsim);

for i=1:nsim
    for j=1:length(h)
        [~, a2] = find(mspe{i}(j,:)== min(mspe{i}(j,:))); 
        eta_hat(j,:,i) = eta_opt{i}(j, :, a2);
        mspe_nm(j,i) = mspe{i}(j, a2);
    end
end

% to find the average and standard deviation of the MSPE over simulations
mspe_nm_sum = sum(mspe_nm,2)

% obtaining the eta estimates
eta_opt_s = zeros(nsim, p-1);
for i=1:nsim
    eta_opt_s(i,:) = eta_hat(find(mspe_nm_sum == min(mspe_nm_sum(:))),:,i);
end

histogram(eta_opt_s, 'NumBins', 12, 'Normalization','pdf');
grid on
[~, eta]=cart2polar(b);
xline(eta, 'rx');
set(findall(gcf,'-property','FontSize'),'FontSize',15);
xlabel('$\hat{{\eta}}$','interpreter','latex', 'FontWeight','bold');
%ylabel('Density','interpreter','latex', 'FontWeight','bold');
xlim([-0.85,-0.58]);
ylim([0, 22]);
title('n=200, $\sigma^2=0.2$','interpreter','latex', 'FontWeight','bold');



subplot(2,3,4)
load('NM_Sphere_results_n50_nsim200_p2_HN.mat')

mspe_nm = NaN(length(h), nsim);
eta_hat = NaN(length(h), p-1, nsim);

for i=1:nsim
    for j=1:length(h)
        [~, a2] = find(mspe{i}(j,:)== min(mspe{i}(j,:))); 
        eta_hat(j,:,i) = eta_opt{i}(j, :, a2);
        mspe_nm(j,i) = mspe{i}(j, a2);
    end
end

% to find the average and standard deviation of the MSPE over simulations
mspe_nm_sum = sum(mspe_nm,2)

% obtaining the eta estimates
eta_opt_s = zeros(nsim, p-1);
for i=1:nsim
    eta_opt_s(i,:) = eta_hat(find(mspe_nm_sum == min(mspe_nm_sum(:))),:,i);
end

histogram(eta_opt_s, 'NumBins', 20, 'Normalization','pdf');
grid on
[~, eta]=cart2polar(b);
xline(eta, 'rx');
set(findall(gcf,'-property','FontSize'),'FontSize',15);
xlabel('$\hat{{\eta}}$','interpreter','latex', 'FontWeight','bold');
%ylabel('Density','interpreter','latex', 'FontWeight','bold');
xlim([-1.4,-0.20]);
ylim([0, 7.5]);
title('n=50, $\sigma^2=0.6$','interpreter','latex', 'FontWeight','bold');



subplot(2,3,5)
load('NM_Sphere_results_n100_nsim200_p2_HN.mat')

mspe_nm = NaN(length(h), nsim);
eta_hat = NaN(length(h), p-1, nsim);

for i=1:nsim
    for j=1:length(h)
        [~, a2] = find(mspe{i}(j,:)== min(mspe{i}(j,:))); 
        eta_hat(j,:,i) = eta_opt{i}(j, :, a2);
        mspe_nm(j,i) = mspe{i}(j, a2);
    end
end

% to find the average and standard deviation of the MSPE over simulations
mspe_nm_sum = sum(mspe_nm,2)

% obtaining the eta estimates
eta_opt_s = zeros(nsim, p-1);
for i=1:nsim
    eta_opt_s(i,:) = eta_hat(find(mspe_nm_sum == min(mspe_nm_sum(:))),:,i);
end

histogram(eta_opt_s, 'NumBins', 15, 'Normalization','pdf');
grid on
[~, eta]=cart2polar(b);
xline(eta, 'rx');
set(findall(gcf,'-property','FontSize'),'FontSize',15);
xlabel('$\hat{{\eta}}$','interpreter','latex', 'FontWeight','bold');
%ylabel('Density','interpreter','latex', 'FontWeight','bold');
xlim([-1.4,-0.20]);
ylim([0, 7.5]);
title('n=100, $\sigma^2=0.6$','interpreter','latex', 'FontWeight','bold');


subplot(2,3,6)
load('NM_Sphere_results_n200_nsim200_p2_HN.mat')

mspe_nm = NaN(length(h), nsim);
eta_hat = NaN(length(h), p-1, nsim);

for i=1:nsim
    for j=1:length(h)
        [~, a2] = find(mspe{i}(j,:)== min(mspe{i}(j,:))); 
        eta_hat(j,:,i) = eta_opt_mspe{i}(j, :, a2);
        mspe_nm(j,i) = mspe{i}(j, a2);
    end
end

% to find the average and standard deviation of the MSPE over simulations
mspe_nm_sum = sum(mspe_nm,2)

% obtaining the eta estimates
eta_opt_s = zeros(nsim, p-1);
for i=1:nsim
    eta_opt_s(i,:) = eta_hat(find(mspe_nm_sum == min(mspe_nm_sum(:))),:,i);
end

histogram(eta_opt_s, 'NumBins', 15, 'Normalization','pdf');
grid on
[~, eta]=cart2polar(b);
xline(eta, 'rx');
set(findall(gcf,'-property','FontSize'),'FontSize',15);
xlabel('$\hat{{\eta}}$','interpreter','latex', 'FontWeight','bold');
%ylabel('Density','interpreter','latex', 'FontWeight','bold');
xlim([-1.4,-0.20]);
ylim([0, 7.5]);
title('n=200, $\sigma^2=0.6$','interpreter','latex', 'FontWeight','bold');


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














%% Codes to generate figure 4 in document

addpath hist2d
fig = figure;

subplot(2,3,1)
load('NM_Sphere_results_n50_nsim200_p3_LN.mat')

mspe_nm = NaN(length(h), nsim);
eta_hat = NaN(length(h), p-1, nsim);

for i=1:nsim
    for j=1:length(h)
        [~, a2] = find(mspe{i}(j,:)== min(mspe{i}(j,:))); 
        eta_hat(j,:,i) = eta_opt{i}(j, :, a2);
        mspe_nm(j,i) = mspe{i}(j, a2);
    end
end

% to find the average and standard deviation of the MSPE over simulations
mspe_nm_sum = sum(mspe_nm,2)

% obtaining the eta estimates
eta_opt_s = zeros(nsim, p-1);
for i=1:nsim
    eta_opt_s(i,:) = eta_hat(find(mspe_nm_sum == min(mspe_nm_sum(:))),:,i);
end

hist2d(eta_opt_s(:,1), eta_opt_s(:,2), 'NBins', 15,'pdf');
grid on
zlim([0,150]);
hold on

[~, eta]=cart2polar(b);
xa = [eta(1), eta(1)];
ya = [eta(2), eta(2)];
za = [0, 150];

plot3(xa',ya',za', 'r-', 'linewidth', 2)
set(findall(gcf,'-property','FontSize'),'FontSize',12);
ylabel('Probability')
xlabel('${\eta_1}$','interpreter','latex', 'FontWeight','bold')
ylabel('${\eta_2}$','interpreter','latex', 'FontWeight','bold')
%zlabel('Density','interpreter','latex', 'FontWeight','bold')
xlim([-1.3,-1]);
ylim([0, 1]);
title('n=50, $\sigma^2=0.2$','interpreter','latex', 'FontWeight','bold');



subplot(2,3,2)
load('NM_Sphere_results_n100_nsim200_p3_LN.mat')

mspe_nm = NaN(length(h), nsim);
eta_hat = NaN(length(h), p-1, nsim);

for i=1:nsim
    for j=1:length(h)
        [~, a2] = find(mspe{i}(j,:)== min(mspe{i}(j,:))); 
        eta_hat(j,:,i) = eta_opt{i}(j, :, a2);
        mspe_nm(j,i) = mspe{i}(j, a2);
    end
end

% to find the average and standard deviation of the MSPE over simulations
mspe_nm_sum = sum(mspe_nm,2);

% obtaining the eta estimates
eta_opt_s = zeros(nsim, p-1);
for i=1:nsim
    eta_opt_s(i,:) = eta_hat(find(mspe_nm_sum == min(mspe_nm_sum(:))),:,i);
end

hist2d(eta_opt_s(:,1), eta_opt_s(:,2), 'NBins', 15,'pdf');
grid on
zlim([0,150]);
hold on

[~, eta]=cart2polar(b);
xa = [eta(1), eta(1)];
ya = [eta(2), eta(2)];
za = [0, 150];

plot3(xa',ya',za', 'r-', 'linewidth', 2)
set(findall(gcf,'-property','FontSize'),'FontSize',12);
ylabel('Probability')
xlabel('${\eta_1}$','interpreter','latex', 'FontWeight','bold')
ylabel('${\eta_2}$','interpreter','latex', 'FontWeight','bold')
%zlabel('Density','interpreter','latex', 'FontWeight','bold')
xlim([-1.3,-1]);
ylim([0, 1]);
title('n=100, $\sigma^2=0.2$','interpreter','latex', 'FontWeight','bold');



subplot(2,3,3)
load('NM_Sphere_results_n200_nsim200_p3_LN.mat')

mspe_nm = NaN(length(h), nsim);
eta_hat = NaN(length(h), p-1, nsim);

for i=1:nsim
    for j=1:length(h)
        [~, a2] = find(mspe{i}(j,:)== min(mspe{i}(j,:))); 
        eta_hat(j,:,i) = eta_opt{i}(j, :, a2);
        mspe_nm(j,i) = mspe{i}(j, a2);
    end
end

% to find the average and standard deviation of the MSPE over simulations
mspe_nm_sum = sum(mspe_nm,2);

% obtaining the eta estimates
eta_opt_s = zeros(nsim, p-1);
for i=1:nsim
    eta_opt_s(i,:) = eta_hat(find(mspe_nm_sum == min(mspe_nm_sum(:))),:,i);
end

hist2d(eta_opt_s(:,1), eta_opt_s(:,2), 'NBins', 15,'pdf');
grid on
zlim([0,150]);
hold on

[~, eta]=cart2polar(b);
xa = [eta(1), eta(1)];
ya = [eta(2), eta(2)];
za = [0, 150];

plot3(xa',ya',za', 'r-', 'linewidth', 2)
set(findall(gcf,'-property','FontSize'),'FontSize',12);
ylabel('Probability')
xlabel('${\eta_1}$','interpreter','latex', 'FontWeight','bold')
ylabel('${\eta_2}$','interpreter','latex', 'FontWeight','bold')
%zlabel('Density','interpreter','latex', 'FontWeight','bold')
xlim([-1.3,-1]);
ylim([0, 1]);
title('n=200, $\sigma^2=0.2$','interpreter','latex', 'FontWeight','bold');



subplot(2,3,4)
load('NM_Sphere_results_n50_nsim200_p3_HN.mat')

mspe_nm = NaN(length(h), nsim);
eta_hat = NaN(length(h), p-1, nsim);

for i=1:nsim
    for j=1:length(h)
        [~, a2] = find(mspe{i}(j,:)== min(mspe{i}(j,:))); 
        eta_hat(j,:,i) = eta_opt{i}(j, :, a2);
        mspe_nm(j,i) = mspe{i}(j, a2);
    end
end

% to find the average and standard deviation of the MSPE over simulations
mspe_nm_sum = sum(mspe_nm,2);

% obtaining the eta estimates
eta_opt_s = zeros(nsim, p-1);
for i=1:nsim
    eta_opt_s(i,:) = eta_hat(find(mspe_nm_sum == min(mspe_nm_sum(:))),:,i);
end

hist2d(eta_opt_s(:,1), eta_opt_s(:,2), 'NBins', 500,'pdf');
grid on
zlim([0, 15]);
hold on

[~, eta]=cart2polar(b);
xa = [eta(1), eta(1)];
ya = [eta(2), eta(2)];
za = [0, 15];

plot3(xa',ya',za', 'r-', 'linewidth', 2)
set(findall(gcf,'-property','FontSize'),'FontSize',12);
ylabel('Probability')
xlabel('${\eta_1}$','interpreter','latex', 'FontWeight','bold')
ylabel('${\eta_2}$','interpreter','latex', 'FontWeight','bold')
%zlabel('Density','interpreter','latex', 'FontWeight','bold')
xlim([-1.5,-0.5]);
ylim([-1, 1.5]);
title('n=50, $\sigma^2=0.6$','interpreter','latex', 'FontWeight','bold');



subplot(2,3,5)
load('NM_Sphere_results_n100_nsim200_p3_HN.mat')

mspe_nm = NaN(length(h), nsim);
eta_hat = NaN(length(h), p-1, nsim);

for i=1:nsim
    for j=1:length(h)
        [~, a2] = find(mspe{i}(j,:)== min(mspe{i}(j,:))); 
        eta_hat(j,:,i) = eta_opt{i}(j, :, a2);
        mspe_nm(j,i) = mspe{i}(j, a2);
    end
end

% to find the average and standard deviation of the MSPE over simulations
mspe_nm_sum = sum(mspe_nm,2);

% obtaining the eta estimates
eta_opt_s = zeros(nsim, p-1);
for i=1:nsim
    eta_opt_s(i,:) = eta_hat(find(mspe_nm_sum == min(mspe_nm_sum(:))),:,i);
end

hist2d(eta_opt_s(:,1), eta_opt_s(:,2), 'NBins', 50,'pdf');
grid on
zlim([0, 15]);
hold on

[~, eta]=cart2polar(b);
xa = [eta(1), eta(1)];
ya = [eta(2), eta(2)];
za = [0, 15];

plot3(xa',ya',za', 'r-', 'linewidth', 2)
set(findall(gcf,'-property','FontSize'),'FontSize',12);
ylabel('Probability')
xlabel('${\eta_1}$','interpreter','latex', 'FontWeight','bold')
ylabel('${\eta_2}$','interpreter','latex', 'FontWeight','bold')
%zlabel('Density','interpreter','latex', 'FontWeight','bold')
xlim([-1.5,-0.5]);
ylim([-1, 1.5]);
title('n=100, $\sigma^2=0.6$','interpreter','latex', 'FontWeight','bold');




subplot(2,3,6)
load('NM_Sphere_results_n200_nsim200_p3_HN.mat')

mspe_nm = NaN(length(h), nsim);
eta_hat = NaN(length(h), p-1, nsim);

for i=1:nsim
    for j=1:length(h)
        [~, a2] = find(mspe{i}(j,:)== min(mspe{i}(j,:))); 
        eta_hat(j,:,i) = eta_opt{i}(j, :, a2);
        mspe_nm(j,i) = mspe{i}(j, a2);
    end
end

% to find the average and standard deviation of the MSPE over simulations
mspe_nm_sum = sum(mspe_nm,2);

% obtaining the eta estimates
eta_opt_s = zeros(nsim, p-1);
for i=1:nsim
    eta_opt_s(i,:) = eta_hat(find(mspe_nm_sum == min(mspe_nm_sum(:))),:,i);
end

hist2d(eta_opt_s(:,1), eta_opt_s(:,2), 'NBins', 28,'pdf');
grid on
zlim([0,15]);
hold on

[~, eta]=cart2polar(b);
xa = [eta(1), eta(1)];
ya = [eta(2), eta(2)];
za = [0, 15];

plot3(xa',ya',za', 'r-', 'linewidth', 2)
set(findall(gcf,'-property','FontSize'),'FontSize',12);
ylabel('Probability')
xlabel('${\eta_1}$','interpreter','latex', 'FontWeight','bold')
ylabel('${\eta_2}$','interpreter','latex', 'FontWeight','bold')
%zlabel('Density','interpreter','latex', 'FontWeight','bold')
xlim([-1.5,-0.5]);
ylim([-1, 1.5]);
title('n=200, $\sigma^2=0.6$','interpreter','latex', 'FontWeight','bold');

%sgtitle('Distribution of estimates $\hat{\eta}$ for p=3','interpreter','latex', 'FontWeight','bold') 

han=axes(fig,'visible','off'); 
han.Title.Visible='on';
han.XLabel.Visible='on';
han.YLabel.Visible='on';
zlabel(han,'Density','interpreter','latex', 'FontWeight','bold','FontSize',15);
%xlabel(han,'${\eta}$','interpreter','latex', 'FontWeight','bold','FontSize',16);
%title(han,'yourTitle');

saveas(gcf,'eta_hat_distribution_p3','epsc')






















%% Codes to generate figure 5 in the document

% Plot of Boxplots for MSEE, high noise

load NM_Sphere_results_n50_nsim200_p2_HN.mat
dv_2_50= msee;

load NM_Sphere_results_n100_nsim200_p2_HN.mat
dv_2_100= msee;

load NM_Sphere_results_n200_nsim200_p2_HN.mat
dv_2_200= msee;

load NM_Sphere_results_n50_nsim200_p3_HN.mat
dv_3_50= msee;

load NM_Sphere_results_n100_nsim200_p3_HN.mat
dv_3_100= msee;

load NM_Sphere_results_n200_nsim200_p3_HN.mat
dv_3_200= msee;

load NM_Sphere_results_n50_nsim200_p4_HN.mat
dv_4_50= msee;

load NM_Sphere_results_n100_nsim200_p4_HN.mat
dv_4_100= msee;

load NM_Sphere_results_n200_nsim200_p4_HN.mat
dv_4_200= msee;


% Creating multiple boxplots for comparing MSEE distributions

x = [dv_2_50, dv_2_100, dv_2_200,dv_3_50, dv_3_100,dv_3_200,dv_4_50,dv_4_100,dv_4_200];
group = [ones(1,200),2*ones(1,200),3*ones(1,200), 4*ones(1,200),5*ones(1,200) ,6*ones(1,200),7*ones(1,200),8*ones(1,200),9*ones(1,200)];
positions = [1 1.25 1.50 2 2.25 2.50 3 3.25 3.50];

boxplot(x, group, 'positions', positions);

set(gca,'xtick',[mean(positions(1:3)) mean(positions(4:6)) mean(positions(7:9))])
set(gca,'xticklabel',{'p=2','p=3','p=4'})

color = ['c', 'y','b', 'c', 'y','b', 'c', 'y','b'];
h = findobj(gca,'Tag','Box');
for j=1:length(h)
   patch(get(h(j),'XData'),get(h(j),'YData'),color(j),'FaceAlpha',.5);
end

c = get(gca, 'Children');

hleg1 = legend(c(1:3), 'n=50', 'n=100', 'n=200','Location','northwest');
set(findall(gca,'-property','FontSize'),'FontSize',15)

grid on
saveas(gca,'msee_boxplots_hn','epsc')













% For low noise

load NM_Sphere_results_n50_nsim200_p2_LN.mat
dv_2_50= msee;

load NM_Sphere_results_n100_nsim200_p2_LN.mat
dv_2_100= msee;

load NM_Sphere_results_n200_nsim200_p2_LN.mat
dv_2_200= msee;

load NM_Sphere_results_n50_nsim200_p3_LN.mat
dv_3_50= msee;

load NM_Sphere_results_n100_nsim200_p3_LN.mat
dv_3_100= msee;

load NM_Sphere_results_n200_nsim200_p3_LN.mat
dv_3_200= msee;

load NM_Sphere_results_n50_nsim200_p4_LN.mat
dv_4_50= msee;

load NM_Sphere_results_n100_nsim200_p4_LN.mat
dv_4_100= msee;

load NM_Sphere_results_n200_nsim200_p4_LN.mat
dv_4_200= msee;

% Creating multiple boxplots for comparing MSEE distributions

x = [dv_2_50, dv_2_100,dv_2_200,dv_3_50,dv_3_100,dv_3_200,dv_4_50,dv_4_100,dv_4_200];
group = [ones(1,200),2*ones(1,200),3*ones(1,200), 4*ones(1,200),5*ones(1,200) ,6*ones(1,200),7*ones(1,200),8*ones(1,200),9*ones(1,200)];
positions = [1 1.25 1.50 2 2.25 2.50 3 3.25 3.50];

boxplot(x, group, 'positions', positions);

set(gca,'xtick',[mean(positions(1:3)) mean(positions(4:6)) mean(positions(7:9))])
set(gca,'xticklabel',{'p=2','p=3','p=4'})

color = ['c', 'y','b', 'c', 'y','b', 'c', 'y','b'];
h = findobj(gca,'Tag','Box');
for j=1:length(h)
   patch(get(h(j),'XData'),get(h(j),'YData'),color(j),'FaceAlpha',.5);
end

c = get(gca, 'Children');

hleg1 = legend(c(1:3), 'n=50', 'n=100', 'n=200','Location','northwest');
set(findall(gca,'-property','FontSize'),'FontSize',15)

grid on

saveas(gca,'msee_boxplots_ln','epsc')














%% Codes to generate figure 6 in the document

% Plot of Boxplots for MSEE, high noise

load NM_Sphere_results_n50_nsim200_p2_HN.mat
dv_2_50= mean(msee_lfp_sum);

load NM_Sphere_results_n100_nsim200_p2_HN.mat
dv_2_100= msee;

load NM_Sphere_results_n200_nsim200_p2_HN.mat
dv_2_200= msee;

load NM_Sphere_results_n50_nsim200_p3_HN.mat
dv_3_50= msee;

load NM_Sphere_results_n100_nsim200_p3_HN.mat
dv_3_100= msee;

load NM_Sphere_results_n200_nsim200_p3_HN.mat
dv_3_200= msee;

load NM_Sphere_results_n50_nsim200_p4_HN.mat
dv_4_50= msee;

load NM_Sphere_results_n100_nsim200_p4_HN.mat
dv_4_100= msee;

load NM_Sphere_results_n200_nsim200_p4_HN.mat
dv_4_200= msee;


% Creating multiple boxplots for comparing MSEE distributions

x = [dv_2_50, dv_2_100, dv_2_200,dv_3_50, dv_3_100,dv_3_200,dv_4_50,dv_4_100,dv_4_200];
group = [ones(1,200),2*ones(1,200),3*ones(1,200), 4*ones(1,200),5*ones(1,200) ,6*ones(1,200),7*ones(1,200),8*ones(1,200),9*ones(1,200)];
positions = [1 1.25 1.50 2 2.25 2.50 3 3.25 3.50];

boxplot(x, group, 'positions', positions);

set(gca,'xtick',[mean(positions(1:3)) mean(positions(4:6)) mean(positions(7:9))])
set(gca,'xticklabel',{'p=2','p=3','p=4'})

color = ['c', 'y','b', 'c', 'y','b', 'c', 'y','b'];
h = findobj(gca,'Tag','Box');
for j=1:length(h)
   patch(get(h(j),'XData'),get(h(j),'YData'),color(j),'FaceAlpha',.5);
end

c = get(gca, 'Children');

hleg1 = legend(c(1:3), 'n=50', 'n=100', 'n=200','Location','northwest');
set(findall(gca,'-property','FontSize'),'FontSize',15)

grid on
saveas(gca,'msee_boxplots_hn','epsc')













% For low noise

load NM_Sphere_results_n50_nsim200_p2_LN.mat
dv_2_50= msee;

load NM_Sphere_results_n100_nsim200_p2_LN.mat
dv_2_100= msee;

load NM_Sphere_results_n200_nsim200_p2_LN.mat
dv_2_200= msee;

load NM_Sphere_results_n50_nsim200_p3_LN.mat
dv_3_50= msee;

load NM_Sphere_results_n100_nsim200_p3_LN.mat
dv_3_100= msee;

load NM_Sphere_results_n200_nsim200_p3_LN.mat
dv_3_200= msee;

load NM_Sphere_results_n50_nsim200_p4_LN.mat
dv_4_50= msee;

load NM_Sphere_results_n100_nsim200_p4_LN.mat
dv_4_100= msee;

load NM_Sphere_results_n200_nsim200_p4_LN.mat
dv_4_200= msee;

% Creating multiple boxplots for comparing MSEE distributions

x = [dv_2_50, dv_2_100,dv_2_200,dv_3_50,dv_3_100,dv_3_200,dv_4_50,dv_4_100,dv_4_200];
group = [ones(1,200),2*ones(1,200),3*ones(1,200), 4*ones(1,200),5*ones(1,200) ,6*ones(1,200),7*ones(1,200),8*ones(1,200),9*ones(1,200)];
positions = [1 1.25 1.50 2 2.25 2.50 3 3.25 3.50];

boxplot(x, group, 'positions', positions);

set(gca,'xtick',[mean(positions(1:3)) mean(positions(4:6)) mean(positions(7:9))])
set(gca,'xticklabel',{'p=2','p=3','p=4'})

color = ['c', 'y','b', 'c', 'y','b', 'c', 'y','b'];
h = findobj(gca,'Tag','Box');
for j=1:length(h)
   patch(get(h(j),'XData'),get(h(j),'YData'),color(j),'FaceAlpha',.5);
end

c = get(gca, 'Children');

hleg1 = legend(c(1:3), 'n=50', 'n=100', 'n=200','Location','northwest');
set(findall(gca,'-property','FontSize'),'FontSize',15)

grid on

saveas(gca,'msee_boxplots_ln','epsc')

