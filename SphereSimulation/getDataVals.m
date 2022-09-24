function [Y, mreg] = getDataVals(theta0, n, p, nsim, s_dt, tau)

  reg_curve = @(t) [sqrt(1-t.^2).*cos(pi*t), sqrt(1-t.^2).*sin(pi*t), t];

  rng(s_dt);
  x = rand(n, nsim, p);  % generate the explanatory variables
  z = cell2mat(arrayfun(@(j) squeeze(x(:, j, :))*theta0, 1:nsim, 'UniformOutput', false))/sqrt(p);
  mreg = arrayfun(@(j) reg_curve(z(:, j))', 1:nsim, 'UniformOutput', false);

  v1 = tau*randn(n, nsim);
  v2 = tau*randn(n, nsim);

  Y = arrayfun(@(i) cell2mat(arrayfun(@(j) add_noise(mreg{i}(:, j), [v1(j, i) v2(j, i)]), 1:n, 'UniformOutput', false)), 1:nsim, 'UniformOutput', false);

end