function superimpose_data(refdir, mutationfile, ttypefile, sigfile, k, allCsfile, fig3file, datoutdir, figoutdir)
% superimpose_data(refdir, mutationfile, ttypefile, sigfile, k, allCsfile, fig3file, datoutdir, figoutdir)
%
% --> re-creates "Bird Plots" from the publication
%     and superimposes an external mutation dataset onto them

if nargin~=9, error('requires: refdir, mutationfile, ttypefile, sigfile, k, allCsfile, fig3file, datoutdir, figoutdir'); end
demand_files({refdir,mutationfile,ttypefile,sigfile,allCsfile,fig3file});
if ~isnumeric(k), k = str2double(k); end
if isnan(k) || k<2 || k>1000, error('k should be a number between 2 and 1000'); end
    
% LOAD FIGURE3 FILE ("OLD" DATA)
load(fig3file,'pat');
pat.new = false(slength(pat),1);

% LOAD MUTATION DATA ("NEW DATA" from user)
% will save as compact MAT file for quicker re-loading
ede(datoutdir);
ede(figoutdir);
mutname = regexprep(mutationfile,'^.*\/([^/]+)$','$1');
mutationmat = [datoutdir '/' mutname '.mat'];
X = load_mutation_data(mutationfile,refdir,mutationmat,ttypefile);

% SIGNATURE ANALYSIS
randseed=1234;
X = analyze_mutation_signatures(X,k,randseed,sigfile);

% LOAD allCs file
tmp=load(allCsfile,'X');tmp=tmp.X;
fs = fieldnames(tmp);
for i=1:length(fs),f=fs{i};X.(f)=tmp.(f);end
clear tmp

% MAP MUTATIONS TO Cs
X.mut.site_idx = map_mutations_to_sites(X.mut,X.site);

% COMPUTE %HP METRICS
X = compute_pcthp_metrics(X);

% CREATE UNIFIED TABLE (OLD DATA + NEW DATA)
X.pat.new = true(slength(X.pat),1);
% remove NMF weightings because they were not done in the same run, and might be of different <k>
X.pat = rmfield_if_exist(X.pat,{'nmf','nmf_norm','nmf_nmut','ci'});
pat = rmfield_if_exist(pat,{'nmf','nmf_norm','nmf_nmut','ci'});
ppp = concat_structs_keep_all_fields({pat,X.pat});
ppp.ttype_clr_ref = round(255*ppp.ttype_clr(:,1));
ppp.ttype_clr_green = round(255*ppp.ttype_clr(:,2));
ppp.ttype_clr_blue = round(255*ppp.ttype_clr(:,3));
ppp.y = ppp.logjit_frac_apo;
ppp.x = ppp.apochar;
ppp.sz = repmat(30,slength(ppp),1); ppp.sz(ppp.new) = 75;
save_struct(rmfield(ppp,'ttype_clr'),[datoutdir '/figure3.superimposed.txt']);

% MAKE BIRD PLOTS ("FIGURE 3")

abcd='abcd';
for i=1:4
  figure(310+i);clf,randinit(1238); fs=35; hmin=0.3; hmax=0.7;
  hold on,set(gca,'fontsize',fs);
  if i==1, ttl='Tumor types';
    x = ppp.x; y = ppp.y; xlab='A3B \leftarrow mutation character \rightarrow A3A'; ylab='% APOBEC mutations';
    show = true(slength(ppp),1);
    clr = ppp.ttype_clr;
  elseif i==2, ttl = 'APOBEC mutations at TpC in DNA hairpins';
    x = ppp.x; y = ppp.y; xlab='A3B \leftarrow mutation character \rightarrow A3A'; ylab='% APOBEC mutations';
    show = find(~isnan(ppp.pcthp_tpc));
    clr = max(0,min(1,(ppp.pcthp_tpc-hmin)/(hmax-hmin)));
    clr = clr*[1 0 0];
  elseif i==3, ttl = 'APOBEC mutations at VpC in DNA hairpins';
    x = ppp.x; y = ppp.y; xlab='A3B \leftarrow mutation character \rightarrow A3A'; ylab='% APOBEC mutations';
    show  = find(~isnan(ppp.pcthp_vpc));
    clr = max(0,min(1,(ppp.pcthp_vpc-hmin)/(hmax-hmin)));
    clr = clr*[1 0 1];
  elseif i==4, ttl = 'TpC and VpC hairpin mutations';
    x = ppp.pcthp_tpc; xlab = 'TpC % hairpin mutations';
    y = ppp.pcthp_vpc; ylab = 'VpC % hairpin mutations';
    show = find(~isnan(x)&~isnan(y));
    clr = ppp.ttype_clr; 
  else error('?');
  end
  idx=find(~ppp.new); scatter(x(idx),y(idx),2*ppp.sz(idx),0.5*repmat([1 1 1],length(idx),1)+0.5*clr(idx,:),'filled','markerfacealpha',0.3,'markeredgealpha',0.3);
  idx=find(ppp.new);  scatter(x(idx),y(idx),ppp.sz(idx),clr(idx,:),'filled','markerfacealpha',1,'markeredgealpha',1,'markeredgecolor',[0 0 0],'linewidth',1);
  xlabel(xlab,'fontsize',fs); ylabel(ylab,'fontsize',fs);
  if i~=4, xlim([-0.18 0.14]); set(gca,'xtick',[]); yt=[0 0.03 0.10 0.3 0.5 0.7 1];set(gca,'ytick',log10(0.01+yt),'yticklabel',yt*100); ylim([-2.1 0.1]);
  else xlim([0 2.6]);ylim([0 2.6]); set(gca,'xtick',[0 1 2],'ytick',[0 1 2]); end
  if i==2, w=9; h=7.5; elseif i==4, w=7.5; h=7.5; else w=8; h=6.6; end
  set(gcf,'papersize',[w h],'paperposition',[0.2 0.2 w-0.4 h-0.4]);
  print_to_file([figoutdir '/figure3' abcd(i) '.superimposed.pdf']);
end

% tumor types colors legend
figure(319),clf,ff
[u ui uj] = unique(ppp.ttype_longname);
clrs = ppp.ttype_clr(ui,:);
x=1;y=0.1;
idx=find(strcmp(u,'')); u(idx) = repmat({'---'},length(idx),1);
x = ones(length(u),1);y=0.1+0.018*[1:length(u)];
scatter(x,y,100,clrs,'filled');
textextra(x+5*xsp,y,u,'fontsize',15,'interp','none','color',clrs);
set(gca,'visible','off','ydir','rev');
w=6;h=10;set(gcf,'papersize',[w h],'paperposition',[0.2 0.2 w-0.4 h-0.4]);
print_to_file([figoutdir '/figure3.superimposed.ttype_legend.pdf']);

% lego plot of new samples
figure(101);clf;display_nmf_legos(X.nmf)
set(gcf,'papersize',[12 9.5],'paperposition',[0.2 0.2 12-0.4 9.5-0.4]);
print_to_file([figoutdir '/new_data.lego_plots.pdf']);






