function h = hist2d_fast_wrapper(x,y,xmin,xmax,ymin,ymax,output_type)
% hist2d_fast_wrapper(x,y,xmin,xmax,ymin,ymax,output_type)
%  
%  min/max arguments are optional
%
%  2020-07-26 added optional last argument "output_type" which can be:
%  1   double (default)
%  2   uint8 (no bin's counts go over 255)
%  3   uint16 (no bin's counts go over 65535)
%  4   uint32 (no bin's counts go over 4.2 billion)

if nargin~=2 && nargin~=3 && nargin~=6 && nargin~=7
  error('wrong number of inputs, was expecting:  x, y [, firstx, lastx, firsty, lasty] [, output_type]');
end

if ~exist('output_type','var'), output_type=1; end
if ischar(output_type), output_type = listmap({output_type},{'double','uint8','uint16','uint32'}); end
if ~(output_type>=1 && output_type<=4), error('invalid output_type'); end

  x = real(double(x));
y = real(double(y));

if any(isnan(x)) || any(isnan(y))
  if nargin<6
    error('Can''t handle NaNs properly without specifying xmin,xmax,ymin,ymax');
  else
    xmin = real(double(xmin));
    xmax = real(double(xmax));
    ymin = real(double(ymin));
    ymax = real(double(ymax));
    x(isnan(x)) = xmin-1;
    y(isnan(y)) = ymin-1;
  end
end

if nargin<6
  h = hist2d_fast(x,y,output_type);
else
  h = hist2d_fast(x,y,xmin,xmax,ymin,ymax,output_type);
end

