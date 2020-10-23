function [R, markEdge, locationEdge] = ED_Harris(x, y, interval, dx, dy, ...
    empConst, thresholdFactor, filterName, winSize)
%ED_HARRIS processes the edge detection with the normalized Harris filter (NHF)
%
%   Notes:   
%         It is recommended to set empirical constant 'empConst' to 1.
%         
%   Input:
%         x                  matrix, coordinates of the data in x (north) direction
%         y                  matrix, coordinates of the data in y (east) direction
%         dx                 matrix, gradient of the data in x direction 
%         dy                 matrix, gradient of the data in y direction 
%         intX               scalar, interval of of the data in x direciton
%         intY               scalar, interval of of the data in y direciton
%         empConst           scalar, empirical constant [0-1]
%         thresholdFactor    scalar, threshold factor [0-1]
%         filterName         string, filter name ['gaussian' | 'sum' | 'average']
%         winSize            scalar, window size of the filter
%
%   Output:
%         R                  Harris response for input data
%         markEdge           edges info located in origanl location
%         locationEdge       edges info located in the extremum location
%
%   Reference:
%         2020, Chen and Zhang, Edge enhancement of the potential field data using the normalized Harris filter	.
%
%
%   T.Chen 14-Jan-2018.
%   Last Revision: 30-Sep-2019.
%   Copyright Tao Chen. 
%   chentaosx@hotmail.com

% arguments check
narginchk(9, 9);
validateattributes(dx,{'numeric'},{'2d'},'ed_Harris','dx');
validateattributes(dx,{'numeric'},{'2d'},'ed_Harris','dx');

% pre-process
[nr, nc] = size(dx);
winSize = abs(winSize); if winSize < 3; winSize = 3; end
winSize = 1 + 2 * floor(winSize/2);
halfWinSize = floor(winSize/2);

% 1) Compute products of derivatives
Ixx = dx.*dx;
Iyy = dy.*dy;
Ixy = dx.*dy;

% 2) Compute derivatives used in Harris filter with different method
switch lower(filterName)
    case 'gaussian'   % Gaussian lowpass filter
        sigma = 1;
        h = fspecial('gaussian', [winSize winSize], sigma);
    case 'sum'   % Harris(1988)
        h = ones(winSize);
    case 'average'
        h = fspecial('average', [winSize winSize]);
    otherwise
        error('ED_Harris: FunctionInput: Invalid_Arg: ''filterName''');
end

Ixx = filter2(h, Ixx);
Ixy = filter2(h, Ixy);
Iyy = filter2(h, Iyy);

% 3) calculate the response 'R' of the Harris detector
% the quadratic terms in the Taylor expansion of 'M' is
%                       [ Ixx   Ixy]
%                       [ Ixy   Iyy]
% det(M) = Ixx.*Iyy - Ixy.*Ixy      trace(M) = Ixx + Iyy
R = Ixx.*Iyy - Ixy.*Ixy + empConst * (Ixx + Iyy).^2;

% 4) extract the result
% for Harris edge response 'R':
% 1. A positive value corresponds to a corner region or edge region;
% 2. A small value (|R|¡Ö0) corresponds to a flat region.
threshold = thresholdFactor * max(abs(R(:)));
index = find(R <= threshold);

markEdge = zeros(nr, nc);
locationEdge = [];
for i = halfWinSize + 1 : nr - halfWinSize
    for j = halfWinSize + 1 : nc - halfWinSize
        
        if ismember(sub2ind(size(R), i, j), index); continue; end
        
        % neighbors in four direction
        M1 = R(i, j - 1 : j + 1);
        M2 = [R(i - 1, j + 1) R(i, j) R(i + 1, j - 1)];
        M3 = R(i - 1 : i + 1, j);
        M4 = [R(i - 1, j - 1) R(i, j) R(i + 1, j + 1)];
        
        % Initialization
        flag = [0 0 0 0];   % flag in four direction for extreme point
        xmax = [0 0 0 0];   % x for extreme point
        ymax = [0 0 0 0];   % y for extreme point
        rmax = [0 0 0 0];   % R(x,y) for extreme point
        
        if R(i, j) == max(M1)   % West- East direction
            flag(1) = 1;
            [xmax(1), ymax(1), rmax(1)] = locationMax(M1, interval, x(i, j), y(i, j), 90);
        end
        
        if R(i, j) == max(M2)   % SouthWest- NorthEast direction
            flag(2) = 1;
            [xmax(2), ymax(2), rmax(2)] = locationMax(M2, sqrt(2) * interval, x(i, j), y(i, j), 45);
        end
        
        if R(i, j) == max(M3)   % South - North direction
            flag(3) = 1;
            [xmax(3), ymax(3), rmax(3)] = locationMax(M3, interval, x(i, j), y(i, j), 0);
        end
        
        if R(i, j) == max(M4)   % SouthEast - NorthWest direction
            flag(4) = 1;
            [xmax(4), ymax(4), rmax(4)] = locationMax(M4, sqrt(2) * interval, x(i, j), y(i, j), -45);
        end
        
        markEdge(i, j) = sum(flag);   % sum the flag vaule
        if markEdge(i, j) == 2 || markEdge(i, j) == 3 || markEdge(i, j) == 4
            [~, indexRmax] = max(rmax);
            locationEdge = cat(1, locationEdge, [xmax(indexRmax) ymax(indexRmax) rmax(indexRmax)]);
        end

    end
end

end

function [xmax, ymax, rmax] = locationMax(M, d, x, y, theta)
% Quadratic function fitting to determine the extreme point position

a = (M(1) - 2 * M(2) + M(3))/(2 * d * d);
b = (M(3) - M(1))/(2 * d);
dmax = -b/(2 * a);
theta = theta*pi/180;

ymax = y + dmax * sin(theta);
xmax = x + dmax * cos(theta);
rmax = a * dmax * dmax + b * dmax + M(2);

end