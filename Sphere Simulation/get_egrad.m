function eg = get_egrad(w, Y, y, M)
    
    mult = abs(Y'*y);
    nz = find(w ~= 0 & mult <0.999); % find non-zero weights
	
    gi = -2*cell2mat(arrayfun(@(k) w(k)*M.dist(Y(:, k), y)*(1 - (Y(:, k)'*y)^2)^(-1/2)*Y(:, k), nz', 'UniformOutput', false)); % element-wise gradient
	
    %ng = find(abs(mean(gi)) ~= Inf);
    eg = sum(gi,2);
	
end
