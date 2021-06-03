function [NHFS, R, upperEnv, ind1, ind2] = NHF(x, y, v, relTol, k)
%NHF processes the edge detection with the Normalized Harris Filter.
%
%   Notes:
%         local peak is found based on Blakely's (1986) method
%
%   Input:
%         x                   required, double, matrix      coordinates of the data along x (north) direction
%         y                   required, double, matrix      coordinates of the data along y (east) direction
%         v                   required, double, matrix      potential field data
%         relTol              required, double, scalar      relative threshold, [0-1], initial guess 1e-3 is recommended
%         k                   optional, double, scalar      empirical constant, [0-1], 1 is recommended
%   Output:
%         NHFS                double, matrix                normalized Harris filter response
%         R                   double, matrix                Harris filter response (without normalization)
%         upperEnv            double, matrix                upper envelope used for normalization
%         ind1                double, vector                index of local maxima before thresholding
%         ind2                double, vector                index of local maxima after thresholding
%
%   Reference:
%         1986, Blakely, R. J., and R. W. Simpson, Approximating edges of source 
%               bodies from magnetic or gravity anomalies.
%
%   Update records:
%         author        date             features
%         Tao Chen      14-Jan-2018      standard implementation
%         Tao Chen      03-Jun-2021      add detailed description

% arguments check
narginchk(4, 5);
if nargin == 4
    k = 1;
end
validateattributes(x, {'numeric'}, {'2d'}, 'NHF','x');
validateattributes(y, {'numeric'}, {'2d'}, 'NHF','y');
validateattributes(v, {'numeric'}, {'2d'}, 'NHF','v');
validateattributes(relTol, {'numeric'}, {'scalar'}, 'NHF','relTol');
validateattributes(k, {'numeric'}, {'scalar'}, 'NHF','k');

% before we get the NHF response, we need calculate the derivatives of v
intX = (max(x(:)) - min(x(:)))/(size(x, 1) - 1);
intY = (max(y(:)) - min(y(:)))/(size(x, 2) - 1); 
[dx, dy, ~] = Gradients(v, 1, intX, intY);
[nr, nc] = size(dx);

% 1) Calculate derivatives of local structure matrix M
A = dx .* dx;
B = dy .* dy;
C = dx .* dy;

winSize = 3; 
halfWinSize = floor(winSize/2);
h = ones(winSize);

A = filter2(h, A);   % (eq.10), we utilize filter2 to calculate eq.10.
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
                   [x(1, :)', y(1, :)', R(1, :)' + 0.1 * max(R(:))], ...   % first row
                   [x(end, :)', y(end, :)', R(end, :)' + 0.1 * max(R(:))], ...   % last row
                   [x(:, 1), y(:, 1), R(:, 1) + 0.1 * max(R(:))], ...   % first column
                   [x(:, end), y(:, end), R(:, end) + 0.1 * max(R(:))]);   % last column
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
                   [x(index), y(index), R(index) + 0.1 * max(R(:))]);
upperEnv = griddata(localMaxs(:, 2), ...
                    localMaxs(:, 1), ...
                    localMaxs(:, 3), ...
                    y, x, 'natural');

% 3.4) normalize R using upper envelope surface                     
NHFS = R ./ upperEnv;

end

% =================================================================
% sub-functions: Gradients, GradientWavenumberDomain, taper2spline
% =================================================================

function [dx, dy, dz] = Gradients(g, order, Xint, Yint)
%GRADIENTS compute the gradients of g in 3D space.
%   GRADIENTS compute the gradients of g in three directions using FFT based method.
%
%   Notes:
%         Square area is used here for expansion
%
%   Input:
%         g          required, double, matrix      matrix of data
%         order      required, double, matrix      order of differentiation, it can be a fractional.
%         Xint       required, double, scalar      x (north) interval of the g 
%         Yint       required, double, scalar      y (east) interval of the g 
%   Output:
%         dx         double, matrix      the horizontal (x drection) gradients of the g
%         dy         double, matrix      the horizontal (y drection) gradients of the g
%         dz         double, matrix      the vertical gradients of the g
%
%   Update records:
%         author        date             features
%         Tao Chen      30-Apr-2015      standard implementation
%         Tao Chen      23-Feb-2021      add detailed description

narginchk(4, 4);

[nr, nc] = size(g);
nmax = max([nr nc]);
npts = 2^nextpow2(nmax + 1);
[dx, dy, dz] = GradientWavenumberDomain(g, npts, nc, nr, Xint, Yint, order);

if 1 == order
    [dy, dx] = gradient(g, Yint, Xint);   % For noisy data, the horizontal gradients in the spatial domain can achieve higher accuracy
end

end

function [dx, dy, dz] = GradientWavenumberDomain(data, npts, nc, nr, Xint, Yint, order)
% this function is adapted from http://software.seg.org/2008/0001 .

cdiff=floor((npts-nc)/2); rdiff=floor((npts-nr)/2);   % Determine the number of extrapolation points
data1=taper2spline(data, npts, nc, nr, cdiff, rdiff);  % Spline extrapolation
f=fft2(data1); fx=f; fy=f; fz=f;
f=fftshift(f);   % Shift zero-frequency component to center of spectrum
wnx=2.0*pi/(Xint*npts);   % Angular frequency，omega = 2*pi*f，wave-length is xint*npts
wny=2.0*pi/(Yint*npts);                                           

cx=npts/2+1; cy=cx;
for I=1:npts
    freqx=(I-cx)*wnx;   % wavenum in x axis
    for J=1:npts
        freqy=(J-cy)*wny;   % wavenum in y axis
        freq=sqrt(freqx*freqx+freqy*freqy);   % wavenum in z axis
        freq=(freq)^order;
        if freqx>0
            freq_x=(abs(freqx))^order*exp(1i*order*pi*0.5);
        else
            freq_x=(abs(freqx))^order*exp(-1i*order*pi*0.5);
        end
        if freqy>0
            freq_y=(abs(freqy))^order*exp(1i*order*pi*0.5);
        else
            freq_y=(abs(freqy))^order*exp(-1i*order*pi*0.5);
        end
        fx(I,J) = f(I,J)*freq_x;
        fy(I,J) = f(I,J)*freq_y;
        fz(I,J) = f(I,J)*freq;
    end
end

fz=fftshift(fz); fx=fftshift(fx); fy=fftshift(fy); 
fzinv=ifft2(fz); fxinv=ifft2(fx); fyinv=ifft2(fy);

dx=real(fxinv(1+rdiff:nr+rdiff,1+cdiff:nc+cdiff));
dy=real(fyinv(1+rdiff:nr+rdiff,1+cdiff:nc+cdiff));
dz=real(fzinv(1+rdiff:nr+rdiff,1+cdiff:nc+cdiff));

end

function gt=taper2spline(g,npts,nc,nr,cdiff,rdiff)
% this function is available from http://software.seg.org/2008/0001 .
% *** grcooper@iafrica.com
% *** www.wits.ac.za/science/geophysics/gc.htm

gt=zeros(npts); gt(rdiff+1:rdiff+nr,cdiff+1:cdiff+nc)=g;
gp=g(:,1:3);     [gpx1,~]=gradient(gp);  % sides
gp=g(:,nc-2:nc); [gpx2,~]=gradient(gp);
x1=0; x2=(2*cdiff)+1;
x=[1 1 0 0;x1 x2 1 1; x1^2 x2^2 2*x1 2*x2; x1^3 x2^3 3*x1^2 3*x2^2]; 
for I=1:nr
    y=[g(I,nc) g(I,1) gpx2(I,3) gpx1(I,1)];
    c=y/x;
    for J=1:cdiff
        gt(I+rdiff,J)=c(1)+(J+cdiff)*c(2)+c(3)*(J+cdiff)^2+c(4)*(J+cdiff)^3;
        gt(I+rdiff,J+nc+cdiff)=c(1)+J*c(2)+c(3)*J^2+c(4)*J^3;
    end
end

gp=g(1:3,:);     [~,gpx1]=gradient(gp);  % top and bottom
gp=g(nr-2:nr,:); [~,gpx2]=gradient(gp);
x1=0; x2=(2*rdiff)+1;
x=[1 1 0 0;x1 x2 1 1; x1^2 x2^2 2*x1 2*x2; x1^3 x2^3 3*x1^2 3*x2^2]; 
for J=1:nc
    y=[g(nr,J) g(1,J) gpx2(3,J) gpx1(1,J)];
    c=y/x;
    for I=1:rdiff
        gt(I,J+cdiff)=c(1)+(I+rdiff)*c(2)+c(3)*(I+rdiff)^2+c(4)*(I+rdiff)^3;
        gt(I+rdiff+nr,J+cdiff)=c(1)+I*c(2)+c(3)*I^2+c(4)*I^3;
    end
end

for I=rdiff+nr+1:npts % Corners
    for J=cdiff+nc+1:npts
        if (I-nr-rdiff)>(J-nc-cdiff); gt(I,J)=gt(I,nc+cdiff); else gt(I,J)=gt(nr+rdiff,J); end
    end
end

for I=1:rdiff
    for J=1:cdiff
        if I>J; gt(I,J)=gt(rdiff+1,J); else gt(I,J)=gt(I,cdiff+1); end
    end
end

for I=1:rdiff % bottom right
    for J=cdiff+nc+1:npts
        if I>(npts-J); gt(I,J)=gt(rdiff+1,J); else gt(I,J)=gt(I,cdiff+nc); end
    end
end

for I=rdiff+nr+1:npts % top left
    for J=1:cdiff
        if (npts-I)>J; gt(I,J)=gt(rdiff+nr,J); else gt(I,J)=gt(I,cdiff+1); end
    end
end

end