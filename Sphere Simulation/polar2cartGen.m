%% The following function computes cartesian coordinates from polar coordinates in general dimensions
% Inputs:  eta: angles in radian
%            r: radius (default = 1)

% Output is the cartesian coordinates.

function theta = polar2cartGen(eta, r)

  if nargin < 2 % assume unit sphere if radius is not provided
      r = 1;
  end

  s = length(eta);
  theta = zeros(s+1,1);

  theta(1) = r*prod(cos(eta));
  theta(end) = r*sin(eta(1));
  for j = 2:s
      
      theta(j) = r*prod(cos(eta(1:(s-j+1))))*sin(eta(s-j+2));
      
  end
  
end
