%% Compute gradient for optimization, using weights w, data Y and input y on manifold M

function eg = get_egrad(w, Y, y, M)

	gi = -2*cell2mat(arrayfun(@(k) w(k)*M.dist(Y(:, k), y)*(1 - (Y(:, k)'*y)^2)^(-1/2)*Y(:, k), 1:size(Y, 2), 'UniformOutput', false)); % element-wise gradient
	eg = sum(gi')';
	
end