%% Here we fit the Local Frechet regression for more than one covariate

clear all;   % clearning all working environment
% generating path for contents of folder 'Manopt_6.0'
addpath(genpath('manopt'))

% read the data
load NM_Sphere_results_n50_nsim200_p2_HN.mat

% Creating the grid of bandwidths:
bw_min = zeros(nsim, 1);
bw_max = zeros(nsim, 1);

for i=1:nsim
    if(p==2)
        xtemp = [x(:,i,1) x(:,i,2)]./sqrt(p);
    elseif(p==3)
        xtemp = [x(:,i,1) x(:,i,2) x(:,i,3)]./sqrt(p);
    elseif(p==4)
        xtemp = [x(:,i,1) x(:,i,2) x(:,i,3) x(:,i,4)]./sqrt(p);
    end 
    ntemp = vecnorm(xtemp')';
    tempd = abs(repmat(ntemp, 1, n) - repmat(ntemp', n, 1));
    tempc = sort(tempd,2);
    bw_min(i)= max(tempc(:,4));
    bw_max(i)= max(tempd(:));
end

h_lfp=exp(linspace(log(max(bw_min)*1.25),log(max(bw_max)),10));

% Creating the array with necessary dimensions to store MSPE and MSEE respectively
mspe_lfp = arrayfun(@(i) zeros(n, length(h_lfp)), 1:nsim, 'UniformOutput', false);
msee_lfp = arrayfun(@(i) zeros(n, length(h_lfp)), 1:nsim, 'UniformOutput', false);

pool=parpool(6);   % starting computation in a parallel pool with 6 workers
parfor i = 1:nsim  % starting the parallel for loop
%for i = 1:nsim
%disp(i)

if(p==2)
    xtemp = [x(:,i,1) x(:,i,2)]./sqrt(p);
elseif(p==3)
    xtemp = [x(:,i,1) x(:,i,2) x(:,i,3)]./sqrt(p);
elseif(p==4)
    xtemp = [x(:,i,1) x(:,i,2) x(:,i,3) x(:,i,4)]./sqrt(p);
end    


for j=1:n
    
for l = 1:length(h_lfp)
    
   disp(['LF_n' num2str(n,'%01d') '_nsim' '_' num2str(i,'%03d') '_' 'Obs' '_' num2str(j,'%03d') '_' 'bandwidth' '_' num2str(l, '%02d') '_p' num2str(p, '%01d')])
   Y_fr = get_sphere_fit_pcov(Y_true{i}, xtemp, xtemp(j,:)', repelem(h_lfp(l),p)); % getting the frechet mean for given theta, covariates
   mspe_lfp{i}(j,l) = acos(Y_true{i}(:, j)'*Y_fr)^2; % compute the MSPE
   msee_lfp{i}(j,l) = acos(mreg{i}(:, j)'*Y_fr)^2;    % compute the MSEE

end
  
end

end
delete(pool) % shutting down parallel pool 


% further computations 
mspe_lfp_sum= zeros(nsim, n, length(h_lfp));
msee_lfp_sum= zeros(nsim, n, length(h_lfp));

for i=1:nsim
    for j=1:n
        for l=1:length(h_lfp)
            mspe_lfp_sum(i, j,l) = mspe_lfp{i}(j,l);
            msee_lfp_sum(i, j,l) = msee_lfp{i}(j,l);
        end
    end
end

% best bandwidth
c3 = find(sum(sum(mspe_lfp_sum, 1),2) == min(sum(sum(mspe_lfp_sum, 1),2)));

h_opt_lfp = h_lfp(c3);

mean(mean(mspe_lfp_sum(:,:,c3), 1))
std(mean(mspe_lfp_sum(:,:,c3), 1))

mean(mean(msee_lfp_sum(:,:,c3), 1))
std(mean(msee_lfp_sum(:,:,c3), 1))

save NM_Sphere_results_n50_nsim200_p3_LN.mat * 

