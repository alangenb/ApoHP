function X = compute_pcthp_metrics(X)

besthp = ((X.site.looplen==3 & X.site.looppos==3 & X.site.ss>=10) | (X.site.looplen==4 & X.site.looppos==4 & X.site.ss>=12));
damp_hp=100;

% TpC
midx = find(X.mut.f==2 & X.mut.l==4);
h = hist2d_fast_wrapper(X.mut.pat_idx(midx),1+nansub(besthp,X.mut.site_idx(midx)),1,slength(X.pat),1,2);
X.pat.pcthp_tpc = 100*(h(:,2))./(damp_hp+sum(h,2));

% VpC
midx = find(X.mut.f==2 & X.mut.l~=4);
h = hist2d_fast_wrapper(X.mut.pat_idx(midx),1+nansub(besthp,X.mut.site_idx(midx)),1,slength(X.pat),1,2);
X.pat.pcthp_vpc = 100*(h(:,2))./(damp_hp+sum(h,2));

  
