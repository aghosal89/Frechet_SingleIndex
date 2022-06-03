function [WnCostAlt, WnGrad, WnHess] = getAltCost(eta, Y, x, h)
%getWnGrad computes the analytic gradient of an approximation to the true
%cost function W_n.  The approximation is due to replacing the local
%Frechet fits in the expression with ordinary local linear fits, normalized
%to have unit norm. The derivative of the normalized local linear fits with
%respect to the parameter theta can be calculated in closed form

theta = polar2cart(eta);
n = size(Y, 2); p = size(x, 2);
z = x*theta;

Y_LL = zeros(size(Y));
%ttil = zeros(p, 3, n);
%Ttil = zeros(p, 3, p, n);

if(nargout >= 2)
    grad1 = zeros(p, 1); 
end
if(nargout == 3)
    hess1 = zeros(p, p);
end

for i = 1:n
    
    % Get normalized local linear fits (these hopefully approximate the LF
    % fits)
        w = getLFRweights(z, z(i), h);
        Ystar = sum(cell2mat(arrayfun(@(j) w(j)*Y(:, j), 1:n, 'UniformOutput', false)), 2);
        Ytil = Ystar/norm(Ystar);
        Y_LL(:, i) = Ytil;
     
    if(nargout >= 2)    
    
    % Get derivatives of local linear weights at x(i, :)
        
        % Kernel weights/derivatives and averages associated with local
        % linear estimation
        
        zdiff = z - z(i);    
        xdiff = x - repmat(x(i, :), n, 1);
        [Kvec, KvecDer, KvecDer2] = K(zdiff, h);
        mu0 = mean(Kvec);  mu1 = mean(Kvec.*zdiff); mu2 = mean(Kvec.*(zdiff.^2));
        sig0sq = mu0*mu2 - mu1^2;
    
        % Derivatives (w.r.t. theta) of muj (for fitting at x(i, :))
        
        h0 = mean(repmat(KvecDer, 1, p).*xdiff)';
        h1 = mean(repmat(KvecDer.*zdiff + Kvec, 1, p).*xdiff)'; 
        h2 = mean(repmat(KvecDer.*(zdiff.^2) + 2*Kvec.*zdiff, 1, p).*xdiff)';
        
        % More derivatives
        
        l = mu2*h0 + mu0*h2 - 2*mu1*h1; % deriv of sig0sq
        v = mu2 - mu1*zdiff; 
        m = repmat(h2, 1, n) - h1*zdiff' - mu1*xdiff'; % derivs of v = (mu2 - m1*zdiff)
        u = repmat((KvecDer.*v)', p, 1).*xdiff' + repmat(Kvec', p, 1).*m; % derivs of numerator of local linear weight
        q = (u/sig0sq - l*(Kvec.*v)'/(sig0sq^2))/n; % derivs of local linear weights
   
    % Get derivative of normalized local linear fit at x(i, :)    
        
        t = 0;
        for j = 1:n
            t = t + q(:, j)*Y(:, j)';
        end % t = derivative of local linear fit at x(i, :)
    
        b = (eye(3) - Ystar*Ystar'/(norm(Ystar)^2))/norm(Ystar); % derivative of normalization
        ttil = t*b; % derivative of normalized local linear fit at x(i, :)
        
        ip = Y(:, i)'*Ytil;
        grad1 = grad1 - 2*acos(ip)*(1 - ip^2)^(-1/2)*ttil*Y(:, i);
    
    end
    
    if(nargout == 3)
        
    % Get Hessians of local linear weights at x(i, :)
    
        % Hessians of muj for fitting at x(i, :)
        
        H0 = zeros(p, p); H1 = H0; H2 = H0;
        for j = 1:n
            
            H0 = H0 + KvecDer2(j)*xdiff(j, :)'*xdiff(j, :)/n;
            H1 = H1 + (KvecDer2(j)*zdiff(j) + 2*KvecDer(j))*xdiff(j, :)'*xdiff(j, :)/n;
            H2 = H2 + (KvecDer2(j)*zdiff(j)^2 + 4*KvecDer(j)*zdiff(j) + 2*Kvec(j))*xdiff(j, :)'*xdiff(j, :)/n;
            
        end
        
        % More Hessians
        
        L = h0*h2' + mu2*H0 + h2*h0' + mu0*H2 - 2*h1*h1' - 2*mu1*H1; % Hess of sig0sq
        M = arrayfun(@(j) H2 - h1*xdiff(j, :) - zdiff(j)*H1 - xdiff(j, :)'*(h1'), 1:n, 'UniformOutput', false); % Hessians of v = (mu2 - m1*zdiff)
        U = arrayfun(@(j) KvecDer2(j)*v(j)*xdiff(j, :)'*xdiff(j, :) + KvecDer(j)*(xdiff(j, :)'*m(:, j)' + m(:, j)*xdiff(j, :)) + Kvec(j)*M{j}, 1:n, 'UniformOutput', false); % Hessians of numerator of local linear weight
        Q = arrayfun(@(j) (U{j}/sig0sq - u(:, j)*l'/(sig0sq^2) - l*u(:, j)'/(sig0sq^2) - Kvec(j)*v(j)*L/(sig0sq^2) + 2*Kvec(j)*v(j)*l*l'/(sig0sq^3))/n, 1:n, 'UniformOutput', false); % Hessians of local linear weights
        
    % Get Hessian of normalized local linear fit at x(i, :)
    
        T = zeros(p, 3, p);  Ttil = T;
        for j = 1:n
            for k = 1:p
                T(:, :, k) = T(:, :, k) + Q{j}(:, k)*Y(:, j)';
            end
        end % Hessian of local linear fit at x(:, i)
        B = zeros(3, 3, 3);
        for k = 1:3
            ek = zeros(3, 1); ek(k) = 1;
            B(:, :, k) = -Ystar(k)*eye(3)/(norm(Ystar)^3) - (ek*Ystar' + Ystar*ek')/(norm(Ystar)^3) + 3*Ystar*Ystar'*Ystar(k)/(norm(Ystar)^5);
        end % Hessian of normalization
        
        for k = 1:p
            
            tmp = zeros(3, 3);
            for kk = 1:3
                for ll = 1:3
                    tmp(kk, ll) = t(k, :)*squeeze(B(kk, ll, :));
                end
            end
            
            Ttil(:, :, k) = squeeze(T(:, :, k))*b + t*tmp;
            
        end
        % Hessian of normalized local linear fit at x(:, i)
         
    tmp = cell2mat(arrayfun(@(k) squeeze(Ttil(:, :, k)*Y(:, i)), 1:p, 'UniformOutput', false));
    hess1 = hess1 - 2*(ttil*Y(:, i)*Y(:, i)'*ttil')*(1 - ip^2)^(-1)*(1 + acos(ip)*(1 - ip^2)^(-1/2)) + acos(ip)*(1 - ip^2)^(-1/2)*tmp;
    hess1 = (hess1 + hess1')/2;
    
    end

end

if(nargout >= 2)

    % Adjust for transformation to polar coordinates;

    J = zeros(p - 1, p); % Jacobian matrix
    etaA = [eta; pi/2]; % augmented eta vector

    for j = 1:(p - 1)
        for k = 1:p
        
            if j > (p - k + 1)    
                J(j, k) = 0; % eta_j not involved in computing theta_k
            elseif j == (p - k + 1)
                J(j, k) = prod(cos(etaA(1:(p - k + 1))));
            else
                tmp = cos(etaA); tmp(j) = sin(etaA(j)); tmp(p - k + 1) = sin(etaA(p - k + 1));
                J(j, k) = -prod(tmp);
            end
        end
    end

    WnGrad = J*grad1/n;

    if(nargout == 3)
        
        WnHess = J*hess1*J'/n;
        
    end
    
end

WnCostAlt = mean((arrayfun(@(g) acos(Y(:, g)'*Y_LL(:, g)), 1:n)).^2); % prediction error on sphere

end
