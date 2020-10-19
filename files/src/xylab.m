function xylab(xlab,ylab,fs,varargin)

if ~exist('fs','var'), fs=20; end

xlabel(xlab,'fontsize',fs,'interp','none',varargin{:});
ylabel(ylab,'fontsize',fs,'interp','none',varargin{:});
