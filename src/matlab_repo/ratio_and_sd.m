function [ratio sd ci95lo ci95hi] = ratio_and_sd(n,N,method)
% [ratio sd ci] = ratio_and_sd(n,N,method)
%
% returns ratio = n/N
%         sd = estimated uncertainty ("standard dev") on ratio
%         ci95lo, ci95hi = 95% confidence interval on ratio
%
% method1 = consider integer counts n to carry sqrt(n) of noise,
%           then propagate this error
% method2 = from binomial fit (Matlab binofit)
%
% Mike Lawrence

if ~exist('method','var'), method=1; end

ratio = n./N;

if method==1   % in use ca. 2008-2019

  sn = n .^ 0.5;
  sN = N .^ 0.5;
  sd = ratio .* ((sn./n).^2+(sN./N).^2).^0.5;
  ci95lo = ratio-1.96*sd; % can produce negative numbers (weakness of method #1)
  ci95hi = ratio+1.96*sd;

elseif method==2  % introduced 2019-03-10

  [ratio ci] = binofit(n,N);
  ci95lo=ci(:,1); ci95hi=ci(:,2);

end


