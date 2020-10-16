function superimpose_data(refdir, mutationfile, ttypefile, sigfile, k, allCsfile, fig3file, datoutdir, figoutdir)
% superimpose_data(refdir, mutationfile, ttypefile, sigfile, k, allCsfile, fig3file, datoutdir, figoutdir)
%
% --> re-creates "Bird Plots" from the publication
%     and superimposes an external mutation dataset onto them
%

if nargin~=9, error('requires: refdir, mutationfile, ttypefile, sigfile, k, allCsfile, fig3file, datoutdir, figoutdir');
demand_files({refdir,mutationfile,ttypefile,sigfile,allCsfile,fig3file});

% LOAD FIGURE3 FILE
ppp = load(fig3file,'X');
ppp = ppp.X.pat;

% LOAD MUTATION DATA
% will save as compact MAT file for quicker re-loading
ede(datoutdir);
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

% MAKE LEGO PLOT

% MAKE BIRD PLOTS ("FIGURE 3")

ppp = load(fig3file,'X'); ppp=ppp.pat;
ppp.ay = log10(0.01+0.01*rand(slength(ppp),1)+ppp.frac_apobec);
ppp.ax = ppp.apochar;
abcd='abcd';
for i=1:4
  figure(300+i);clf,randinit(1238); fs=35; sz=60; hmin=0.3; hmax=0.7;
  hold on,set(gca,'fontsize',fs);
  if i==2, ttl='Tumor types';
    x = ppp.ax; y = ppp.ay; xlab='A3B \leftarrow mutation character \rightarrow A3A'; ylab='% APOBEC mutations';
    show = true(slength(ppp),1); clr = ppp.ttype_clr;
  elseif i==1, ttl = 'APOBEC mutations at TpC in DNA hairpins';
    x = ppp.ax; y = ppp.ay; xlab='A3B \leftarrow mutation character \rightarrow A3A'; ylab='% APOBEC mutations';
    clr = max(0,min(1,(ppp.pcthp_tpc-hmin)/(hmax-hmin)));
    show = ~isnan(clr); clr = clr*[1 0 0]; % red
  elseif i==3, ttl = 'APOBEC mutations at VpC in DNA hairpins';
    x = ppp.ax; y = ppp.ay; xlab='A3B \leftarrow mutation character \rightarrow A3A'; ylab='% APOBEC mutations';
    clr = max(0,min(1,(ppp.pcthp_vpc-hmin)/(hmax-hmin)));
    show = ~isnan(clr); clr = clr*[1 0 1]; % magenta
  elseif i==4, ttl = 'TpC and VpC hairpin mutations';
    x = ppp.pcthp_tpc; xlab = 'TpC % hairpin mutations';
    y = ppp.pcthp_vpc; ylab = 'VpC % hairpin mutations';
    clr = ppp.ttype_clr; show = (~isnan(x)&~isnan(y));
  else error('?');
  end
  %scatter(x(show),y(show),sz,clr(show,:),'filled');ff;set(gca,'fontsize',fs);%title(ttl,'fontsize',fs);

  clf; scatter(x(~isnew),y(~isnew),sz(~isnew),clr(~isnew,:),'filled','markerfacealpha',0.3,'markeredgealpha',0.3); hold on;
  scatter(x(isnew),y(isnew),sz(isnew),clr(isnew,:),'filled','markerfacealpha',1,'markeredgealpha',1);
  ff;set(gca,'fontsize',fs);%title(ttl,'fontsize',fs);
  


  xlabel(xlab,'fontsize',fs); ylabel(ylab,'fontsize',fs);
  if i~=4, xlim([-0.18 0.14]); set(gca,'xtick',[]); yt=[0 0.03 0.10 0.3 0.5 0.7 1];set(gca,'ytick',log10(0.01+yt),'yticklabel',yt*100); ylim([-2.1 0.1]);
  else xlim([0 2.6]);ylim([0 2.6]); set(gca,'xtick',[0 1 2],'ytick',[0 1 2]); end
  if i==2, w=9; h=7.5; elseif i==4, w=7.5; h=7.5; else w=8; h=6.6; end
  set(gcf,'papersize',[w h],'paperposition',[0.2 0.2 w-0.4 h-0.4]);
  print_to_file([figoutdir '/figure3' abcd(i) '.superimposed.pdf']);

  tmp = struct(); tmp.name = ppp.name; tmp.x = x; tmp.y = y;
  save_struct(tmp,[figoutdir abcd(i) '.xydata.table.txt']);
end








