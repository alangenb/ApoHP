function b = logjit(a,width)

  if ~exist('width','var'), width=1; end
  
  b = log10(1+a+0.5*(1-width)+width*rand(size(a)));

  
