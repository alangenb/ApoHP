function [r sd] = relative_ratio_and_sd(n,N,r_baseline)
% compute relative rate and sd for confidence interval
% --> normalize to rate provided as r_baseline

[r sd] = ratio_and_sd(n,N);
r = bsxfun(@rdivide, r, r_baseline);
sd = bsxfun(@rdivide, sd, r_baseline);

