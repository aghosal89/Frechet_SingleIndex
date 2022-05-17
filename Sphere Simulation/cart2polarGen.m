%% The following function computes cartesian coordinates from polar coordinates for general dimensions
% Inputs :  x: coordinates in euclidean space R^p.
       
% Outputs:  r  : radius of polar coordinate
%           eta: Angle of polar coordinate

function [eta,r] = cart2polar(x)
    r = norm(x);        % find the norm
    p = length(x);      % dimension of the euclidean space
    
    if p == 1
        error('Hyperspherical coordinates not meaningful for p = 1')
    end
    
  eta =  zeros(p-1,1);  % dimension of the polar coordinate
  
  for j = 1:(p-1)
  
    eta(j) = atan(x(p-j+1)/norm(x(1:p-j)));
  
  end
  
end
