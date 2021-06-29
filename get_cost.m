%% Compute cost function for optimization, using weights w, data Y and input y on manifold M

function c = get_cost(w, Y, y, M)

	c = dot(w, arrayfun(@(j) (M.dist(y, Y(:, j)))^2, 1:size(Y, 2)));
	
end