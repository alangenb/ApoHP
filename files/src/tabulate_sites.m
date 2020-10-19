function X = tabulate_sites(X)
% X = tabulate_sites(X)
%
% --> creates table of unique site classes (by hairpin characteristics)
% --> counts number of sites and mutations for each class
% --> evaluate three nested APOBEC+ cohorts: >=10% , >=50%, >=90% APOBECC
% --> calculates rate and confidence interval, normalized to the global TpC baseline for each cohort

demand_fields(X,{'site','pat','mut'});
demand_fields(X.site,{'zone','looplen','looppos','minus2','minus1','minus0','plus1','plus2','ss'});
demand_fields(X.pat,{'frac_apobec','msupe_neg'});
demand_fields(X.mut,{'pat_idx','site_idx'});

fprintf('Tabulating sites...\n');

suse = (X.site.zone<3);    % exclude coding regions (exons+spliceflanks) from the modeling (to minimize impact of driver mutations)
puse = (X.pat.msupe_neg);  % exclude patients with >10% of MSI, smoking, UV, POLE, or ESO (to minimize impact of other hypermutation processes)

% stemstrength bins
X.ssbin.min = [0 4:22]';
X.ssbin.max = [3 4:21 inf]';
X.site.ssbin = min(max(1,X.site.ss-3),slength(X.ssbin));

% define APOBEC cohorts
X.cohort=[]; X.cohort.name={}; X.cohort.puse=logical([]);
X.cohort.name{end+1,1}='APOBEC>=90%, MSUPE-';          X.cohort.puse(end+1,:)=(puse & X.pat.frac_apobec>=0.9);
X.cohort.name{end+1,1}='APOBEC>=50% and <90%. MSUPE-'; X.cohort.puse(end+1,:)=(puse & X.pat.frac_apobec>=0.5 & X.pat.frac_apobec<0.9);
X.cohort.name{end+1,1}='APOBEC>=10% and <50%, MSUPE-'; X.cohort.puse(end+1,:)=(puse & X.pat.frac_apobec>=0.1 & X.pat.frac_apobec<0.5);
X.cohort.name{end+1,1}='APOBEC>=10% and <50%, MSUPE+'; X.cohort.puse(end+1,:)=(~puse & X.pat.frac_apobec>=0.1 & X.pat.frac_apobec<0.5);
X.cohort.min_frac_apobec = [0.9;0.5;0.1;0.1];
X.cohort.include_msupe = [0;0;0;1];
X.cohort.npat = sum(X.cohort.puse,2);
X.pat.ci = nan(slength(X.pat),1);
for ci=1:length(X.cohort.name)
  X.cohort.nmut(ci,1) = sum(nansub(suse,X.mut.site_idx) & as_column(X.cohort.puse(ci,X.mut.pat_idx)));
  X.pat.ci(X.cohort.puse(ci,:))=ci;
end

fprintf('Counting mutations\n');

% count mutations at each site
X.site.n = hist2d_fast_wrapper(X.mut.site_idx,X.pat.ci(X.mut.pat_idx),1,slength(X.site),1,slength(X.cohort),3);

fprintf('Getting site classes\n');

% get list of site classes
[u ui uj] = unique(uint32(X.site.looplen-3)+9*uint32(X.site.looppos)+9*13*uint32(X.site.minus2-1)+9*13*4*uint32(X.site.minus1-1)+9*13*4*4*uint32(X.site.minus0-1)+...
                   9*13*4*4*4*uint32(X.site.plus1-1)+9*13*4*4*4*4*uint32(X.site.plus2-1)+9*13*4*4*4*4*4*uint32(X.site.ssbin-1));
S = reorder_struct(X.site,ui);
S.N = accumarray(uj(suse),1);
S.n = nan(slength(S),slength(X.cohort)); for ci=1:slength(X.cohort), S.n(:,ci) = accumarray(uj(suse),X.site.n(suse,ci)); end

% convert cohorts to cumulative (nested)
X.cohort=rmfield(X.cohort,'puse');
X.cohort.name{1} = 'APOBEC>=90%, MSUPE-';
X.cohort.name{2} = 'APOBEC>=50%, MSUPE-';
X.cohort.name{3} = 'APOBEC>=10%, MSUPE-';
X.cohort.name{4} = 'APOBEC>=10%, all';
X.cohort.npat = cumsum(X.cohort.npat);
X.cohort.nmut = cumsum(X.cohort.nmut);
S.n = cumsum(S.n,2);

fprintf('Tabulation step #1\n');

% FIRST STEP OF TABULATION
maxlooplen=11;
N = nan(maxlooplen-2,maxlooplen+2,4,4,4,4,4,slength(X.ssbin));
n = nan(maxlooplen-2,maxlooplen+2,4,4,4,4,4,slength(X.ssbin),slength(X.cohort));
for looplen=3:maxlooplen, s1=find(S.looplen==looplen);
  for looppos=0:looplen+1, s2=s1(S.looppos(s1)==looppos);
    for minus2=1:4, s3=s2(S.minus2(s2)==minus2);
      for minus1=1:4, s4=s3(S.minus1(s3)==minus1);
        for minus0=1:4, s4a=s4(S.minus0(s4)==minus0);
          for plus1=1:4, s5=s4a(S.plus1(s4a)==plus1);
            for plus2=1:4, s6=s5(S.plus2(s5)==plus2);
              for ssbin=1:slength(X.ssbin), s7=s6(S.ssbin(s6)==ssbin);
                N(looplen-2,looppos+1,minus2,minus1,minus0,plus1,plus2,ssbin) = sum(S.N(s7));
                n(looplen-2,looppos+1,minus2,minus1,minus0,plus1,plus2,ssbin,:) = sum(S.n(s7,:),1);
end,end,end,end,end,end,end,end
              
% SECOND LEVEL OF TABULATION
% --> create table of "loop types" defined by looplen, looppos, and local sequence context (1/2/3/4 = A/C/G/T)
% --> include some "degenerate" loop types with 0 at one or more positions:
%     at minus2, minus1, plus1, plus2:  0 = N (any base)
%                           at minus0:  0 = V (A/C/G)   

fprintf('Tabulation step #2\n');

X.loop=[];
nl = round(1.1*sum([3:maxlooplen]+2)*(5.^5));  % preallocate space so arrays aren't growing during the loops
fs={'looplen','looppos','minus2','minus1','minus0','plus1','plus2'}; for fi=1:length(fs),f=fs{fi};X.loop.(f)=nan(nl,1); end;
X.loop.N = nan(nl,slength(X.ssbin));
X.loop.n = nan(nl,slength(X.ssbin),slength(X.cohort));

li=0;
for looplen=3:maxlooplen, N1=squeeze(N(looplen-2,:,:,:,:,:,:,:)); n1=squeeze(n(looplen-2,:,:,:,:,:,:,:,:));
  for looppos=0:looplen+1, N2=squeeze(N1(looppos+1,:,:,:,:,:,:)); n2=squeeze(n1(looppos+1,:,:,:,:,:,:,:));
    for minus2=0:4, if minus2==0, N3=squeeze(sum(N2,1)); n3=squeeze(sum(n2,1));
      else N3=squeeze(N2(minus2,:,:,:,:,:)); n3=squeeze(n2(minus2,:,:,:,:,:,:)); end
      for minus1=0:4, if minus1==0, N4=squeeze(sum(N3,1)); n4=squeeze(sum(n3,1));
        else N4=squeeze(N3(minus1,:,:,:,:)); n4=squeeze(n3(minus1,:,:,:,:,:)); end
        for minus0=0:4, if minus0==0, N4a=squeeze(sum(N4(1:3,:,:,:),1)); n4a=squeeze(sum(n4(1:3,:,:,:,:),1));
          else N4a=squeeze(N4(minus0,:,:,:)); n4a=squeeze(n4(minus0,:,:,:,:)); end
          for plus1=0:4, if plus1==0, N5=squeeze(sum(N4a,1)); n5=squeeze(sum(n4a,1));
            else N5=squeeze(N4a(plus1,:,:,:)); n5=squeeze(n4a(plus1,:,:,:,:)); end
            for plus2=0:4, if plus2==0, N6=squeeze(sum(N5,1)); n6=squeeze(sum(n5,1));
              else N6=squeeze(N5(plus2,:,:)); n6=squeeze(n5(plus2,:,:,:)); end
              li=li+1; for fi=1:length(fs),f=fs{fi};X.loop.(f)(li)=eval(f);end;
              X.loop.N(li,:)=N6; X.loop.n(li,:,:)=n6;
end,end,end,end,end,end,end
X.loop=reorder_struct(X.loop,1:li);   % trim list to final length
 
% compute global TCA baseline of each cohort
tca = find(X.loop.minus0==4 & X.loop.plus1==1 & X.loop.minus1+X.loop.minus2+X.loop.plus2==0);
tca_rate = bsxfun(@rdivide,nansum(nansum(X.loop.n(tca,:,:),1),2),nansum(nansum(X.loop.N(tca,:),1),2));
X.cohort.tca_rate = squeeze(tca_rate);
X.cohort.tca_baseline_muts_per_mb = 1e6*X.cohort.tca_rate./X.cohort.npat;  % overall mutations per megabase per patient

% compute relative mutation rate in each bin (and confidence intervals)
[X.loop.rate X.loop.sd] = relative_ratio_and_sd(X.loop.n,X.loop.N,tca_rate);
% --> these are the observed rates shown in Figs.2,6,S2








