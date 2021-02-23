function [NHFS, R, upperEnv, ind1, ind2] = NHF(x, y, dx, dy, relTol, k)
%NHF processes the edge detection with the Normalized Harris Filter.
%
%   Notes:
%         local peak is found based on Blakely's (1986) method
%
%   Input:
%         x          required, double, matrix      coordinates of the data along x (north) direction
%         y          required, double, matrix      coordinates of the data along y (east) direction
%         dx         required, double, matrix      The horizontal gradient along the x-axis
%         dy         required, double, matrix      The horizontal gradient along the y-axis
%         relTol     required, double, scalar      relative threshold, [0-1e-5]
%         k          required, double, scalar      empirical constant, range from 0 to 1, 1 is recommended
%   Output:
%         NHFS       double, matrix      normalized Harris filter response
%         R          double, matrix      Harris filter response (without normalization)
%         upperEnv   double, matrix      upper envelope used for normalization
%         ind1       double, vector      index of local maxima before thresholding
%         ind2       double, vector      index of local maxima after thresholding
%
%   Reference:
%         1986, Blakely, R. J., and R. W. Simpson, Approximating edges of source 
%               bodies from magnetic or gravity anomalies.
%
%   Update records:
%         author        date             features
%         Tao Chen      14-Jan-2018      standard implementation
%         Tao Chen      23-Feb-2021      add detailed description

% arguments check
narginchk(5, 6);
if nargin == 6
    k = 1;
end
validateattributes(x, {'numeric'}, {'2d'}, 'ED_NHF','x');
validateattributes(y, {'numeric'}, {'2d'}, 'ED_NHF','y');
validateattributes(dx, {'numeric'}, {'2d'}, 'ED_NHF','dx');
validateattributes(dy, {'numeric'}, {'2d'}, 'ED_NHF','dy');
validateattributes(relTol, {'numeric'}, {'scalar'}, 'ED_NHF','relTol');
validateattributes(k, {'numeric'}, {'scalar'}, 'ED_NHF','k');

% preliminary preparation
[nr, nc] = size(dx);
winSize = 3; 
halfWinSize = floor(winSize/2);

% 1) Calculate derivatives of local structure matrix M (eq.9)
A = dx .* dx;
B = dy .* dy;
C = dx .* dy;

h = ones(winSize);
A = filter2(h, A);   % (eq.10)
C = filter2(h, C);   % (eq.11)
B = filter2(h, B);   % (eq.12)

% 2) calculate the Harris filter response R (eq. 13)
R = A .* B - C .* C + k * (A + B).^2;

% 3) amplitude balance
% 3.1) pick out the local maxima of R using Blakely’s algorithm (1986)
indexOfLocalMaxs = zeros(nr, nc);
for i = halfWinSize + 1 : nr - halfWinSize
    for j = halfWinSize + 1 : nc - halfWinSize
        
        % neighbors of R(i, j)
        M1 = R(i, j - 1 : j + 1);
        M2 = [R(i - 1, j + 1) R(i, j) R(i + 1, j - 1)];
        M3 = R(i - 1 : i + 1, j);
        M4 = [R(i - 1, j - 1) R(i, j) R(i + 1, j + 1)];
        
        % Initialization
        flag = [0 0 0 0];   % flag of maximum
        
        if R(i, j) == max(M1)   % West- East direction
            flag(1) = 1;
        end
        
        if R(i, j) == max(M2)   % SouthWest- NorthEast direction
            flag(2) = 1;
        end
        
        if R(i, j) == max(M3)   % South - North direction
            flag(3) = 1;
        end
        
        if R(i, j) == max(M4)   % SouthEast - NorthWest direction
            flag(4) = 1;
        end
        
        indexOfLocalMaxs(i, j) = sum(flag);

    end
end

indexOfLocalMaxs = find(indexOfLocalMaxs > 1);
ind1 = indexOfLocalMaxs;   % index of local maxima without thresholding

% 3.2) remove local maximum below the threshold
absTol = relTol * max(abs(R(:)));
indexOfLocalMaxs(R(indexOfLocalMaxs) <= absTol) = [];
ind2 = indexOfLocalMaxs;   % index of local maxima after thresholding

% 3.3) calculate the upper envelope surface of the R
localMaxs = cat(2, x(indexOfLocalMaxs), y(indexOfLocalMaxs), R(indexOfLocalMaxs));
localMaxs = cat(1, localMaxs, ...
                   [x(1, :)', y(1, :)', R(1, :)' * 100], ...   % first row
                   [x(end, :)', y(end, :)', R(end, :)' * 100], ...   % last row
                   [x(:, 1), y(:, 1), R(:, 1) * 100], ...   % first column
                   [x(:, end), y(:, end), R(:, end) * 100]);   % last column
localMaxs = unique(localMaxs, 'rows');
upperEnv = griddata(localMaxs(:, 2), ...
                    localMaxs(:, 1), ...
                    localMaxs(:, 3), ...
                    y, x, 'natural');

% theoretically, the 'upperEnv' should be greater than 'R'. However, due to the 
% inevitable shorcomings of the 'natural' method, we need to correct some points where
% 'upperEnv' is less than 'R'
index = find(upperEnv < R);
index = setdiff(index, indexOfLocalMaxs);
localMaxs = cat(1, localMaxs, ...
                   [x(index), y(index), R(index) * 1.01]);
upperEnv = griddata(localMaxs(:, 2), ...
                    localMaxs(:, 1), ...
                    localMaxs(:, 3), ...
                    y, x, 'natural');

% 3.4) normalize R using upper envelope surface                     
NHFS = R ./ upperEnv;

end