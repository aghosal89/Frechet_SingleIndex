function eg = get_egrad(w, Y, y, M, thresh)

    if(nargin < 5 || isempty(thresh))
        thresh = 0.999;
    end
    
    mult = abs(y'*Y);
    nz = find(w ~= 0 & mult < thresh); % find weights that are non-zero and sufficiently far away from +/- 1
	
    gi = -2*cell2mat(arrayfun(@(k) w(k)*M.dist(Y(:, k), y)*(1 - (Y(:, k)'*y)^2)^(-1/2)*Y(:, k), nz', 'UniformOutput', false)); % element-wise gradient
	
    %ng = find(abs(mean(gi)) ~= Inf);
    eg = sum(gi,2);
	
end
