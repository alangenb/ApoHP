function h = hist3d_fast_wrapper(x,y,z,xmin,xmax,ymin,ymax,zmin,zmax,output_type)
% hist3d_fast_wrapper(x,y,xmin,xmax,ymin,ymax,zmin,zmax,output_type)
%
%  min/max arguments are optional
%
%  2020-07-26 added optional last argument "output_type" which can be:
%  1   double (default)
%  2   uint8 (no bin's counts go over 255)
%  3   uint16 (no bin's counts go over 65535)
%  4   uint32 (no bin's counts go over 4.2 billion)

if nargin~=3 && nargin~=4 && nargin~=9 && nargin~=10
  error('wrong number of inputs, was expecting:  x, y [, firstx, lastx, firsty, lasty, firstz, lastz] [, output_type]');
end

if ~exist('output_type','var'), output_type=1; end

x = real(double(x));
y = real(double(y));
z = real(double(z));

if any(isnan(x)) || any(isnan(y)) || any(isnan(z))
  if nargin<9
    error('Can''t handle NaNs properly without specifying xmin,xmax,ymin,ymax,zmin,zmax');
  else
    xmin = real(double(xmin));
    xmax = real(double(xmax));
    ymin = real(double(ymin));
    ymax = real(double(ymax));
    zmin = real(double(zmin));
    zmax = real(double(zmax));
    x(isnan(x)) = xmin-1;
    y(isnan(y)) = ymin-1;
    z(isnan(z)) = zmin-1;
  end
end

if nargin<9
  h = hist3d_fast(x,y,z,output_type);
else
  h = hist3d_fast(x,y,z,xmin,xmax,ymin,ymax,zmin,zmax,output_type);
end
