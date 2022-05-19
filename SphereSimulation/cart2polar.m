% The following function computes cartesian coordinates from polar coordinates
% Inputs :  x: coordinates in euclidean space R^p.
       
% Outputs:  r  : radius of polar coordinate
%           eta: Angle of polar coordinate

function [r,eta] = cart2polar(x)
    r = norm(x);        % find the norm
    p = length(x);      % dimension of the euclidean space
  eta =  zeros(p-1,1);  % dimension of the polar coordinate
  if(p==2)
      eta(1) = atan(x(2)/x(1));
  elseif(p==3) 
    eta(1) = atan(x(3)/norm([x(1),x(2)]));
    eta(2) = atan(x(2)/x(1));
  elseif(p==4) 
    eta(1) = atan(x(4)/norm([x(1),x(2),x(3)]));
    eta(2) = atan(x(3)/norm([x(1),x(2)]));
    eta(3) = atan(x(2)/x(1));
  elseif(p==5) 
    eta(1) = atan(x(5)/norm([x(1),x(2),x(3),x(4)]));
    eta(2) = atan(x(4)/norm([x(1),x(2),x(3)]));
    eta(3) = atan(x(3)/norm([x(1),x(2)]));
    eta(4) = atan(x(2)/x(1));
  end
  
end
