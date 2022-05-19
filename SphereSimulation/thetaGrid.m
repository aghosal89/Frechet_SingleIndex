function tg = thetaGrid(nsp, s)
% Creates a grid of equispaced hyperspherical coordinates and converts to a
% Cartesian grid on the s-1 dimensional hypersphere embedded in R^s
% Inputs
%
%   nsp - number of equally spaced points in each dimension
%   s - dimension of hypersphere coordinates

spc= pi/nsp;  % increment
f = linspace(-pi/2+spc/2,pi/2-spc/2,((pi-spc)/spc)+1)'; % equidistant starting points btetween -pi/2 and pi/2

% create s - 1 output arrays [S1 ... S_{s-1}] represented as strings
outStr1 = '[S1'; outStr2 = strcat(outStr1, '(:)');
if s > 2

    for j = 1:(s - 1)
        
        outStr1 = strcat(outStr1, ' S', num2str(j));
        outStr2 = strcat(outStr2, ' S', num2str(j), '(:)');

    end
end

outStr1 = strcat(outStr1, ']');
outStr2 = strcat(outStr2, ']');

% Create arrays
eval(strcat(outStr1, '= ndgrid(f);')); 

% Reform into a 2-D array
eval(strcat('etaGrid = ', outStr2, ';'));

% Transform to Cartesian coordinates
tg = cell2mat(arrayfun(@(j) polar2cartGen(etaGrid(j, :)'), 1:size(etaGrid, 1), 'UniformOutput', false));
    

end

