%% Compute gradient for optimization, using weights w, data Y and input y on manifold M

function eg = get_egrad(wm, Y, y, M)
    
    mult = abs(Y'*y);
    nz = find(wm ~= 0 & mult <1); % find non-zero weights
	
    gi = -2*cell2mat(arrayfun(@(k) wm(k)*M.dist(Y(:, k), y)*(1 - (Y(:, k)'*y)^2)^(-1/2)*Y(:, k), nz', 'UniformOutput', false)); % element-wise gradient
	
    %ng = find(abs(mean(gi)) ~= Inf);
    eg = sum(gi,2);
	
end