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

cdiff=floor((npts-nc)/2); rdiff=floor((npts-nr)/2);   % 确定行列扩边点数
data1=taper2spline(data, npts, nc, nr, cdiff, rdiff);  % 样条插值扩边
f=fft2(data1); fx=f; fy=f; fz=f;
f=fftshift(f);   % 将零频移到中间
wnx=2.0*pi/(Xint*npts);   % 空间频率换算为角频率，omega=2*pi*f，这里npts点的波长为xint*npts
wny=2.0*pi/(Yint*npts);                                           

cx=npts/2+1; cy=cx;
for I=1:npts
    freqx=(I-cx)*wnx;   % x方向波数（和频率域的频率对应）
    for J=1:npts
        freqy=(J-cy)*wny;   % y方向波数（和频率域的频率对应）
        freq=sqrt(freqx*freqx+freqy*freqy);   % z垂向波数（和频率域的频率对应）
        freq=(freq)^order;   % z垂向导数滤波算子（通过调节order计算分数阶导数）
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

dx=real(fxinv(1+rdiff:nr+rdiff,1+cdiff:nc+cdiff));   % 缩边
dy=real(fyinv(1+rdiff:nr+rdiff,1+cdiff:nc+cdiff));   % 缩边
dz=real(fzinv(1+rdiff:nr+rdiff,1+cdiff:nc+cdiff));   % 缩边

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