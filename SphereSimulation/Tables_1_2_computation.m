
% Computation of best bandwidths 

% for any simulation setting load the corresponding '.mat' file from the
% folder 'FinalSimResults' then run the following set of codes 


% Get noise level

if(tau == sqrt(0.4))
    noise = 'low';
else
    noise = 'high';
end

% create the matrix of SE values corresponding to each simulation and
% bandwidth for estimation of theta
SE_h_sim = zeros(200, length(h));

for i=1:200
    for l=1:length(h)
        SE_h_sim(i,l) = acos(abs(fsiFitAll{i}.thetaHat(:,l)'*theta0))^2;
    end
end

MSE_h = mean(SE_h_sim);

xh = find(MSE_h== min(MSE_h));
h(xh) % report best bandwidth

SE_sim = SE_h_sim(:,xh); % values only for minimizing bandwidth

% Computation of average and standard deviation of MSE
disp(['The average and standard deviation of the squared estimation error for theta when p = ' num2str(p) ', n = ' num2str(n), ' and under ' ...
    noise ' noise are:'])
mean(SE_sim)
std(SE_sim)

% create the matrix of MSEE values corresponding to each simulation and
% bandwidth
MSEE_h_sim = zeros(200, length(h));

for i=1:200
    for l=1:length(h)
        MSEE_h_sim(i,l) = mean(arrayfun(@(j) (acos(fsiFitAll{i}.Yout(:,j,l)'*mreg{i}(:,j)))^2, 1:n));
    end
end

MSEE_h = mean(MSEE_h_sim);

xh = find(MSEE_h== min(MSEE_h));

h(xh) % report best bandwidth

MSEE_sim = MSEE_h_sim(:,xh);

% mean and standard deviation of MSEE values
disp(['The average and standard deviation of the squared FSI regression estimation error when p = ' num2str(p) ', n = ' num2str(n), ' and under ' ...
    noise ' noise are:'])
mean(MSEE_sim)
std(MSEE_sim)


% create the matrix of MSPE, MSEE values corresponding to each simulation and
% bandwidth for Local Frechet with p covariates 

LFp_MSEE_h_sim = zeros(200, length(h));

for i=1:200
    for l=1:length(h)
        LFp_MSEE_h_sim(i,l) = mean(arrayfun(@(j) (acos(LFpcovFitAll{i}(:,j,l)'*mreg{i}(:,j)))^2, 1:n));
    end
end

LFp_MSEE_h = mean(LFp_MSEE_h_sim);

xh_LFp_MSEE = find(LFp_MSEE_h== min(LFp_MSEE_h));

h(xh_LFp_MSEE) % report best bandwidth

LFpMSEE_sim = LFp_MSEE_h_sim(:,xh_LFp_MSEE);

mean(LFpMSEE_sim)
std(LFpMSEE_sim)

% mean and standard deviation of MSEE values
disp(['The average and standard deviation of the squared mLF regression estimation error when p = ' num2str(p) ', n = ' num2str(n), ' and under ' ...
    noise ' noise are:'])
mean(LFpMSEE_sim)
std(LFpMSEE_sim)
