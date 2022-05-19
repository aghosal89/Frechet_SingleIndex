%% Compute Hessian for optimization, using weights w, data Y and input y on manifold M, along vect u

function eh = get_ehess(w, Y, y, M, u, thresh)

    if(nargin < 6 || isempty(thresh))
        thresh = 0.999;
    end

	hi = zeros(3, 3, size(Y, 2));
	ip = arrayfun(@(k) Y(:, k)'*y, 1:size(Y, 2));
	
	for i = 1:size(Y, 2)
        if(w(i) ~= 0 && abs(ip(i))< 0.999) % use weights that are non-zero and sufficiently far away from +/- 1
    		hi(:, :, i) = 2*w(i)*(sqrt(1 - ip(i)^2) - ip(i)*acos(ip(i)))/((1 - ip(i)^2)^(3/2))*Y(:, i)*Y(:, i)';	
        end
    end	
	tmp = sum(hi, 3);
    eh = tmp*u;
	
end
