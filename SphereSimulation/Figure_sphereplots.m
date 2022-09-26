%% Code to create the figure of sphere plots in the document

% clear the working directory 
clear 

% create path to manifold optimization folder
addpath(genpath('manopt'))

% load the file for a simulation output e.g. 
% here I load the estimation results for n=50, p=5, high noise setting 
load Sphere_results_n50_nsim200_p5_noiseHigh_Final.mat

% since we have to work with a single dataset, we set 
i=1;

% following is an arbitrary parameter value considerably far from the true parameter.

thf0 = [0.7422, 0.1967, 0.2693, 0.0084, 0.5812]';

% create the matrix of MSEE values corresponding to each simulation and
% bandwidth
MSEE_h_sim = zeros(200, length(h));

for i1=1:200
    for l=1:length(h)
        MSEE_h_sim(i1,l) = mean(arrayfun(@(j) (acos(fsiFitAll{i1}.Yout(:,j,l)'*mreg{i1}(:,j)))^2, 1:n));
    end
end

MSEE_h = mean(MSEE_h_sim);

xh = find(MSEE_h== min(MSEE_h));
h_opt = h(xh);

% following is the estimate of parameter from first dataset
th_opt = fsiFitAll{i}.thetaHat(:,xh);

xtmp = squeeze(x(:,i,:));

z = (xtmp*theta0)./sqrt(p);
z_t2 = (xtmp*th_opt)./sqrt(p);
z_t3 = (xtmp*thf0)./sqrt(p);

% read from the results
Y_t = Y{i};

% obtain the estimated response for estimated (previously) index parameter
Yhat= fsiFitAll{i}.Yout(:,:,xh);

% obtain the estimated response for an index parameter far from the truth
Yhatf= get_sphere_fit_LF(Y_t, z_t3, h_opt);

% create equidistant grid of points over range of single index for the true
% single index parameter
xout1 = linspace(min(z), max(z), 101);

% create the regression function 
yt1 = [sqrt(1-xout1'.^2).*cos(pi*xout1'), sqrt(1-xout1'.^2).*sin(pi*xout1'), xout1']';
m_tr1 = arrayfun(@(j) yt1(j, :), 1:3, 'UniformOutput', false);

[s1 s2 s3] = sphere(128);
lightGrey = 0.9*ones(1, 3);

set(gcf, 'PaperUnits', 'inches', 'PaperSize', [2.4, 2.6], 'PaperPosition', [0 0 2.4 2.6])
surf(s1, s2, s3, 'FaceColor', lightGrey, 'EdgeColor', 'none', 'FaceAlpha', 0.4)
view([-1, 1, 1])
hold on;
plot3(m_tr1{1}, m_tr1{2}, m_tr1{3}, 'k', 'LineWidth', 0.5)
plot3(Y_t(1, :), Y_t(2, :), Y_t(3, :), 'r.', 'MarkerSize', 2.5)

hold on;

plot3(Yhat(1, :), Yhat(2, :), Yhat(3, :), 'b.', 'MarkerSize', 2.5)

hold on;
plot3(Yhatf(1, :), Yhatf(2, :), Yhatf(3, :), 'g.', 'MarkerSize', 2.5)

% save the image in an .eps file
saveas(gcf,'highnoise_n50_p5','epsc')

% finally save the outputs in an extended file to be reproduced further.
save Sphere_results_n50_nsim200_p5_noiseHigh_Final_extended.mat thf0 Yhatf m_tr1 xh h_opt
