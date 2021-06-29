% The following function computes cartesian coordinates from polar coordinates
% Inputs:  eta: angles in radian
%            r: radius

function theta = polar2cart(eta, r)
  s = length(eta);
  theta = zeros(s+1,1);
  
  if(s==3) 
  theta(1)= r*cos(eta(1))*cos(eta(2))*cos(eta(3));
  theta(2)= r*cos(eta(1))*cos(eta(2))*sin(eta(3));
  theta(3)= r*cos(eta(1))*sin(eta(2));
  theta(4)= r*sin(eta(1));
  
  end
  
  if(s==2)
    theta(1) = r*cos(eta(1))*cos(eta(2));
    theta(2) = r*cos(eta(1))*sin(eta(2));
    theta(3) = r*sin(eta(1));
  end
  
  if(s==1)
    theta(1)= r*cos(eta(1));
    theta(2)= r*sin(eta(1));
  end
  
  if(s==4) 
  theta(1)= r*cos(eta(1))*cos(eta(2))*cos(eta(3))*cos(eta(4));
  theta(2)= r*cos(eta(1))*cos(eta(2))*cos(eta(3))*sin(eta(4));
  theta(3)= r*cos(eta(1))*cos(eta(2))*sin(eta(3));
  theta(4)= r*cos(eta(1))*sin(eta(2));
  theta(5)= r*sin(eta(1));
  end
  
end