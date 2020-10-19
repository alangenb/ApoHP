function X = apobec_quantitative_modeling(X,cohort)
% X = apobec_quantitative_modeling(X,cohort)
% --> call after tabulate_sites()
% --> does curve fitting to solve for the global and local parameters of the A3A hairpin optimality metric
% --> assigns X.site.relrate_exp, expected relative mutation rate,
%     expressed relative to the cohort's non-hairpin TpC baseline of 1
%
% <cohort> is which cohort to use for the modeling

demand_fields(X,{'cohort','loop','ssbin','site'});
demand_fields(X.cohort,'npat');
demand_fields(X.loop,{'looplen','looppos','minus1','minus0','plus1','N','n'});
demand_fields(X.ssbin,{'min','max'});
demand_fields(X.site,{'looplen','looppos','minus1','minus0','plus1','ssbin'});
                    
if ~exist('cohort','var'), cohort=1; end

fprintf('APOBEC quantitative modeling...\n');

Q = X.loop;
Q.n = Q.n(:,:,cohort);
Q.rate = Q.rate(:,:,cohort);
Q.sd = Q.sd(:,:,cohort);
npat = X.cohort.npat(cohort);
Q = reorder_struct(Q,Q.minus2==0 & Q.plus2==0); Q=rmfield(Q,{'minus2','plus2'});

% convert rates to mutations per megabase per patient
rescale = X.cohort.tca_baseline_muts_per_mb(cohort);
Q.rate = rescale*Q.rate;
Q.sd = rescale*Q.sd;

% CURVE FITTING

fprintf('Fitting curves.\n');

x=X.ssbin.min;                      % domain of stem strength
measurepoint=find(x==20);           % stem strength at which to quantify hairpin effect
minmut=10;                          % minimum number of mutations per bin: cap predictions where data runs out
mhps=[0.44:0.01:0.50];              % domain of possible values for global parameter m_hp = hairpin slope
xas=[0:0.5:50 51:1:100 102:2:200];  % domain of values for local parameter xa = stem strength where rate has doubled from no-hairpin

W=[]; W.mhp=as_column(mhps);
for wi=1:slength(W)
  mhp=W.mhp(wi); r0 = 1e6*sum(Q.n)/sum(Q.N)/npat;
  rhp = nan(length(x),length(xas)); for xai=1:length(xas), rhp(:,xai) = r0*2.^(mhp*(x-xas(xai))); end
  for qi=slength(Q):-1:1
    robs = as_column(Q.rate(qi,:)); sdobs = as_column(Q.sd(qi,:)); robs(Q.n(qi,:)<minmut)=nan;
    rinit = robs(1); Q.rinit(qi,1)=rinit; rexp = rinit+rhp; yobs = log10(0.0001+robs); yexp = log10(0.0001+rexp);
    yerr = nansum(bsxfun(@minus,yexp,yobs).^2,1); [err xai] = min(yerr);
    Q.err(qi,1)=err; Q.xa(qi,1)=xas(xai); Q.rexp(qi,:) = rexp(:,xai);
    lastdata = find(~isnan(robs),1,'last'); if isempty(lastdata), lastdata=1; end
    Q.rexp(qi,lastdata+1:end)=Q.rexp(qi,lastdata); Q.relrate_exp(qi,:) = Q.rexp(qi,:)/r0; Q.rr = Q.relrate_exp(:,measurepoint);
    Q.rexp(qi,lastdata+1:end)=nan; Q.lastdata(qi,1)=lastdata;
  end,W.err(wi,1)=sum(Q.err); W.Q{wi,1}=Q;
end,[tmp wi] = min(W.err); Q=W.Q{wi};

% convert back to the global TCA baseline
all_qi = find(Q.minus1==0 & Q.minus0~=0 & Q.plus1~=0); all_ssbin = 1;
all_r = sum(Q.n(all_qi,all_ssbin))/sum(Q.N(all_qi,all_ssbin));
tca_qi = find(Q.minus1==0 & Q.minus0==4 & Q.plus1==1); tca_ssbin = 1;
tca_r = nansum(Q.n(tca_qi,tca_ssbin))/nansum(Q.N(tca_qi,tca_ssbin));
ratio = tca_r/all_r;
fs = {'rate','sd','rinit','rexp','relrate_exp','rr'}; for fi=1:length(fs),f=fs{fi};Q.(f)=Q.(f)/ratio; end

% apply model to all C's in genome

fprintf('Applying model genomewide.\n');

z=nan(slength(X.site),1,'single'); X.site.q0=z; X.site.q1=z; X.site.relrate_exp0=z; X.site.relrate_exp1=z;
for looplen=3:11, i1=find(X.site.looplen==looplen);q1=find(Q.looplen==looplen);
  for looppos=0:looplen+1, i2=i1(X.site.looppos(i1)==looppos);q2=q1(Q.looppos(q1)==looppos);
    for minus0=1:4, i2a=i2(X.site.minus0(i2)==minus0);q2a=q2(Q.minus0(q2)==minus0);
      for ssbin=1:slength(X.ssbin), i3=i2a(X.site.ssbin(i2a)==ssbin);q3=q2a(Q.minus1(q2a)==0 & Q.plus1(q2a)==0); x=ssbin;
        y=q3; X.site.q0(i3)=y; X.site.relrate_exp0(i3)=Q.relrate_exp(y,x);
        for minus1=1:4, i4=i3(X.site.minus1(i3)==minus1);q4=q2a(Q.minus1(q2a)==minus1);
          for plus1=1:4, i5=i4(X.site.plus1(i4)==plus1);q5=q4(Q.plus1(q4)==plus1);
            y=q5; X.site.q1(i5)=y; X.site.relrate_exp1(i5)=Q.relrate_exp(y,x);
end,end,end,end,end,end

% use minus1/plus1 to stratify where there is suffient data, otherwise use only minus0/looplen/looppos/ssbin
X.site.relrate_exp = X.site.relrate_exp1; X.site.qi = X.site.q1; idx=find(isnan(X.site.relrate_exp1));
X.site.relrate_exp(idx)=X.site.relrate_exp0(idx); X.site.qi(idx) = X.site.q0(idx);
X.site = rmfield(X.site,{'q0','q1','relrate_exp0','relrate_exp1'});

% save model in main data object
X.model = Q;
X.site = rename_field(X.site,'qi','model_idx');








