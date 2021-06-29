%% Takes a 2-vector u as coefficients to create a vector in the tangent space at p and compute exponential map (spherical data)

function y = add_noise(p, u)

    if(any(p == 0))
        error('Still need to add functionality for points on axes')
    end
    
    p1 = p(1); p2 = p(2); p3 = p(3);
    
    e11 = sqrt(p2^2/(p1^2 + p2^2));
    e1 = [e11, -e11*p1/p2, 0]';
    
    e21 = 1/sqrt(1 + p2^2/p1^2 + (p1^2 + p2^2)^2/(p1^2*p3^2));
    e2 = [e21, e21*p2/p1, -(p1^2 + p2^2)/(p1*p3)*e21]';
    
    v = u(1)*e1 + u(2)*e2;
    
    y = cos(norm(v))*p + sin(norm(v))*v/norm(v);

end