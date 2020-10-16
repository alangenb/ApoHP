function VpC_process_data(mutfile,outdir,topdir)

srcdir = [topdir 'matlab_repo/'];
refdir = [topdir 'ref/'];
datdir = [topdir 'data/processed/'];
initdir = [topdir 'data/init/'];
figdir = [topdir 'figs/'];

addpath(srcdir);
maxNumCompThreads(8);

if ~exist('outdir','var') outdir = './'; end
ede(outdir);

load([datdir 'MAX.pat_object.v1.0.mat'],'pat');
load([datdir 'MAX.mut_object.v1.0.mat'],'mut');
load([datdir 'MAX.ttype_object.v1.0.mat'],'ttype');

X = struct(); 
X.pat = pat;
X.mut = mut;
X.ttype = ttype;

%mutfile = [datdir 'test.mutfile_input.v1.0.txt'];

m = makeapnq(load_struct_noheader(mutfile));
m = rename_fields(m,{'col1','col2','col3','col4','col5'},{'patient','chr','pos','ref','alt'});
if iscellstr(m.chr) m.chr = convert_chr(m.chr); end
if iscellstr(m.ref) m.ref = mapacross(m.ref,{'A','C','G','T'},[1:4]); end
if iscellstr(m.alt) m.alt = mapacross(m.alt,{'A','C','G','T'},[1:4]); end
ia = m.chr>=1 & m.chr<=24 & ~isnan(m.ref) & ~isnan(m.alt); m = reorder_struct(m,ia);
m.pat_idx = grp2idx(m.patient);
m.context65 = get_context(m.chr,m.pos,[refdir 'context65']);

p = struct();
p.name = unique(m.patient,'stable');
p.ttype = cell(slength(p),1); p.ttype(:) = {'NA'};
p.cohort = cell(slength(p),1); p.cohort(:) = {'NA'};
p.ttype_long = cell(slength(p),1); p.ttype_long(:) = {'NA'};
p.nmut = accumarray(m.pat_idx,repmat(1,slength(m),1));
p.ttype_idx = repmat(NaN,slength(p),1);

m = rmfield(m,{'patient'});
m = orderfields(m,[1 2 5 3 4 6]);
m.pat_idx = m.pat_idx + max(X.mut.pat_idx);

X.pat = concat_structs({X.pat,p});
X.mut = concat_structs({X.mut,m});

X.mut.patient = nansub(X.pat.name,X.mut.pat_idx);
P=[]; P.coding_only=false; X.nmf_input = make_input_for_nmf(X.mut,P);
X.mut = rmfield(X.mut,'patient');
X.nmf_input_rates = X.nmf_input;
load([initdir 'breakdowns.v1.1.mat'],'K');
for i=1:slength(X.nmf_input_rates.chan)
  X.nmf_input_rates.chan.terr(i,1)=sum(K.Ng(K.f==X.nmf_input_rates.chan.f(i)&K.l==X.nmf_input_rates.chan.l(i)&K.r==X.nmf_input_rates.chan.r(i)));
end
X.nmf_input_rates.pat.nchan_counts = X.nmf_input_rates.pat.nchan;
X.nmf_input_rates.pat.nchan = bsxfun(@rdivide,X.nmf_input_rates.pat.nchan,X.nmf_input_rates.chan.terr');

sig = makeapnq(load_struct([datdir 'APOBEC.NMF_signature.COSMIC_2_and_13.static.v1.0.txt']));
sig.nmf(95) = 0;       %remove T(C->T)G from signature

maxcos = []; maxkcor = []; maxpcor = []; seccos = []; apoidx = []; rseeds = [];
tic;
revstr = '';
randseed = 1234;
for k=3:14
X.nmf_rates = perform_nmf(X.nmf_input_rates,k,randseed);
cos = []; pcor = []; kcor = [];
for j=1:k
  tmp = X.nmf_rates.chan.nmf(:,j);
  tmp(95,:) = 0;
  cos(j) = sum(tmp.*sig.nmf)/sqrt(sum(tmp.^2)*sum(sig.nmf.^2));
  [pcor(j),~] = corr(tmp,sig.nmf,'type','pearson');
  [kcor(j),~] = corr(tmp,sig.nmf,'type','kendall');
end
[maxs,mia] = sort(cos,'descend');
maxcos(k-2) = maxs(1); apoidx(k-2) = mia(1); seccos(k-2) = maxs(2);
maxpcor(k-2) = max(pcor); maxkcor(k-2)  = max(kcor);
%[maxcos(k-2) apoidx(k-2)] = max(cos);
pct = 100*(k-2)/14;
msg = sprintf('Percent done: %3.1f',pct);
fprintf([revstr,msg]);
revstr = repmat(sprintf('\b'),1,length(msg));
end
toc     % ~15 secs for one set

i1 = find(maxcos'>=0.97,1);
k = i1+2;

X.nmf_rates = perform_nmf(X.nmf_input_rates,k,randseed);
clf; display_nmf_legos(X.nmf_rates);
ord = [apoidx(i1) setdiff([1:size(X.nmf_rates.chan.nmf,2)],apoidx(i1))];
names = [{'APOBEC'};strcat('sig',arrayfun(@num2str,[2:size(X.nmf_rates.chan.nmf,2)]','uni',0))];
X.nmf.factor.name = as_column(names);
X.nmf.chan = X.nmf_rates.chan;
X.nmf.chan.nmf = X.nmf_rates.chan.nmf(:,ord);
pord = listmap(X.pat.name,X.nmf_rates.pat.name);
X.pat.nchan_counts = X.nmf_rates.pat.nchan_counts(pord,:);
X.pat.nchan_rates = X.nmf_rates.pat.nchan(pord,:);
X.pat.nmf = X.nmf_rates.pat.nmf(pord,ord);
X.pat.nmf_norm = X.nmf_rates.pat.nmf_norm(pord,ord);
X=rmfield(X,{'nmf_input','nmf_input_rates','nmf_rates'});
X.pat=keep_fields(X.pat,{'name','cohort','ttype','ttype_long','ttype_idx','nmut','nchan_counts','nchan_rates','nmf','nmf_norm'});

X.pat.frac_apobec = X.pat.nmf_norm(:,1);

tic;X.mut = add_llftrr(X.mut);toc; % ~3min
c2gt = (X.mut.f==2 & X.mut.t~=1);
X.pat.nC = histc(X.mut.pat_idx(c2gt),1:slength(X.pat));
tca2gt = (c2gt & X.mut.l==4 & X.mut.r==1);
X.pat.nRTCA = histc(X.mut.pat_idx(tca2gt & (X.mut.ll==1 | X.mut.ll==3)),1:slength(X.pat));
X.pat.nYTCA = histc(X.mut.pat_idx(tca2gt & (X.mut.ll==2 | X.mut.ll==4)),1:slength(X.pat));
X.pat.RTCA_C = X.pat.nRTCA./X.pat.nC; X.pat.YTCA_C = X.pat.nYTCA./X.pat.nC;

%

% map mutations to all cytosines w/ model
ia = X.mut.chr>=1 & X.mut.chr<=24; X.mut = reorder_struct(X.mut,ia);
X.mut.gpos = chrpos2gpos(X.mut.chr,X.mut.pos);
[~,ia] = sort(X.mut.gpos,'ascend'); X.mut = reorder_struct(X.mut,ia);
tmp=struct(); tmp.X=X;

cfile = [datdir 'all.Cs_only.w_MAX_model_90.v1.0.mat'];
tic;load(cfile,'X');toc % ~2min
tic; X.site.gpos = chrpos2gpos(X.site.chr,X.site.pos); toc                              % ~30 sec
fs={'pat','ttype','nmf','mut'};for fi=1:length(fs),f=fs{fi};X.(f)=tmp.X.(f);end
tic; locb = ismembc2(X.mut.gpos,X.site.gpos); toc                                       % ~40 sec
X.mut.sidx = locb; X.mut.sidx(X.mut.sidx==0) = NaN;

%

rrmax=8; X.bin=[]; X.bin.num = (1:rrmax)'; X.bin.relrate = (1:rrmax)';
tic;X.site.bin = min(rrmax,max(1,round((X.site.relrate_exp))));toc % ~6 min
X.mut.bin = nansub(X.site.bin,X.mut.sidx);

% all C's
tic;X.pat.rrbin = hist2d_fast_wrapper(X.mut.pat_idx,X.mut.bin,1,slength(X.pat),1,slength(X.bin));toc % ~30 sec
X.pat.frrbin = bsxfun(@rdivide,X.pat.rrbin,sum(X.pat.rrbin,2)+100);
% TpC's
tpc = nansub(X.site.minus0,X.mut.sidx)==4;
tic;X.pat.rrbin_tpc = hist2d_fast_wrapper(X.mut.pat_idx(tpc),X.mut.bin(tpc),1,slength(X.pat),1,slength(X.bin));toc % ~30 sec
X.pat.frrbin_tpc = bsxfun(@rdivide,X.pat.rrbin_tpc,sum(X.pat.rrbin_tpc,2)+100);
% VpC's
tic;X.pat.rrbin_vpc = hist2d_fast_wrapper(X.mut.pat_idx(~tpc),X.mut.bin(~tpc),1,slength(X.pat),1,slength(X.bin));toc % ~30 sec
X.pat.frrbin_vpc = bsxfun(@rdivide,X.pat.rrbin_vpc,sum(X.pat.rrbin_vpc,2)+100);

isnew = ismember(X.pat.ttype,'NA'); X.pat.ttype_clr = repmat(NaN,slength(X.pat),3);
X.pat.ttype_clr(~isnew,:) = get_ttype_colors(X.pat.ttype(~isnew));
X.pat.ttype_clr(isnew,:) = repmat([66 245 102]./255,sum(isnew),1);

X0=X;

X = rmfield(X,{'mut','site'});
%save([datdir 'tmp.final_object.v1.2.mat'],'X');

%

%%%
%make plots

pat = X.pat;

ppp=pat;
ppp.ay = log10(0.01+0.01*rand(slength(ppp),1)+ppp.frac_apobec);
ppp.ax = 0.016+ppp.YTCA_C-2.1*ppp.RTCA_C;

tags = {'A','B','C','D'}';

figure(1);clf,fs=20;randinit(1234);
sz = repmat(30,slength(ppp),1); sz(isnew) = 75;
for i=1:4,cla,hold on,set(gca,'fontsize',fs);
  if i==1, ttl='Tumor types';
    x = ppp.ax; y = ppp.ay; xlab='A3B \leftarrow mutation character \rightarrow A3A'; ylab='% APOBEC mutations';
    show = true(slength(ppp),1);
    clr = ppp.ttype_clr;
  elseif i==2, ttl = 'APOBEC mutations at TpC in DNA hairpins';
    x = ppp.ax; y = ppp.ay; xlab='A3B \leftarrow mutation character \rightarrow A3A'; ylab='% APOBEC mutations';
    show = find(~isnan(ppp.frrbin_tpc(:,1)));
    damp=500; ppp.frrbin_tpc = bsxfun(@rdivide,ppp.rrbin_tpc,sum(ppp.rrbin_tpc,2)+damp); pcthp_tpc = 100*sum(ppp.frrbin_tpc(:,5:end),2); hmin=0.3; hmax=0.7;
    clr = max(0,min(1,(pcthp_tpc-hmin)/(hmax-hmin)));
    clr = clr*[1 0 0]; clr(isnew,:) = repmat([66 245 102]./255,sum(isnew),1);
  elseif i==3, ttl = 'APOBEC mutations at VpC in DNA hairpins';
    x = ppp.ax; y = ppp.ay; xlab='A3B \leftarrow mutation character \rightarrow A3A'; ylab='% APOBEC mutations';
    show = find(~isnan(ppp.frrbin_vpc(:,1)));
    damp=500; ppp.frrbin_vpc = bsxfun(@rdivide,ppp.rrbin_vpc,sum(ppp.rrbin_vpc,2)+damp); pcthp_vpc = 100*sum(ppp.frrbin_vpc(:,5:end),2); hmin=0.3; hmax=0.7;
    clr = max(0,min(1,(pcthp_vpc-hmin)/(hmax-hmin)));
    clr = clr*[0.3 0.3 1]; clr(isnew,:) = repmat([66 245 102]./255,sum(isnew),1);
  elseif i==4, ttl = 'TpC and VpC hairpin mutations';
    show = find(~isnan(ppp.frrbin_tpc(:,1)));
    damp=500; ppp.frrbin_tpc = bsxfun(@rdivide,ppp.rrbin_tpc,sum(ppp.rrbin_tpc,2)+damp); pcthp_tpc = 100*sum(ppp.frrbin_tpc(:,5:end),2);
    damp=500; ppp.frrbin_vpc = bsxfun(@rdivide,ppp.rrbin_vpc,sum(ppp.rrbin_vpc,2)+damp); pcthp_vpc = 100*sum(ppp.frrbin_vpc(:,5:end),2);
    x = pcthp_tpc; xlab = 'TpC % hairpin mutations';
    y = pcthp_vpc; ylab = 'VpC % hairpin mutations';
    clr = ppp.ttype_clr;
  else error('?');
  end
  clf; scatter(x(~isnew),y(~isnew),sz(~isnew),clr(~isnew,:),'filled','markerfacealpha',0.3,'markeredgealpha',0.3); hold on;
  scatter(x(isnew),y(isnew),sz(isnew),clr(isnew,:),'filled','markerfacealpha',1,'markeredgealpha',1);
  ff;set(gca,'fontsize',fs);%title(ttl,'fontsize',fs);
  xlabel(xlab,'fontsize',fs); ylabel(ylab,'fontsize',fs);
  if ismember(i,[1 2 3])
    xlim([-0.18 0.14]); set(gca,'xtick',[]);
    yt=[0 0.03 0.10 0.3 0.5 0.7 1];set(gca,'ytick',log10(0.01+yt),'yticklabel',yt*100); ylim([-2.1 0.1]);
  end
  if ismember(i,[4])
    xlim([min(x) max(x)]); ylim([min(x) max(x)]);
    set(gca,'xtick',[1 2],'ytick',[1 2]);
  end
  set(gcf,'papersize',[20 20],'paperposition',[0.2 0.2 20-0.4 20-0.4]);
  print_to_file([outdir 'FIGURE.3' tags{i} '.superimposed.png']);
  tmp = struct(); tmp.name = ppp.name; tmp.x = x; tmp.y = y;
  save_struct(tmp,[outdir tags{i} '.xydata.table.txt']);
end


