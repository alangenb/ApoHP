function make_figures(indir,outdir)

ede(outdir);

fprintf('Generating figures.\n');

% LEGO PLOTS

load([indir '/figure.lego_plots.mat'],'nmf');

figure(1);clf;display_nmf_legos(nmf)
set(gcf,'papersize',[12 9.5],'paperposition',[0.2 0.2 12-0.4 9.5-0.4]);
print_to_file([outdir '/figure.lego_plots.pdf']);

% FIGURE 2A

load([indir '/figure2a.mat'],'ssbin');

figure(201),clf,fs=24; yts = [0 0.5 1;0 0.3 0.6;0 0.2 0.4;0 2.5 5];
looplen = size(ssbin.rate,2);
cmap = green_orange_colormap(looplen);
for minus0=1:4, subplot(4,1,minus0),hold on,x=0;xt=[];
  rate = ssbin.rate(:,:,minus0); sd = ssbin.sd(:,:,minus0); ratehi=rate+1.96*sd; ratelo=max(0,rate-1.96*sd);
  x=0; xt=[]; ymax=0;
  for ssi=1:slength(ssbin)
    x=x+2; xt(end+1)=x+0.5+looplen/2;
    for looppos=1:looplen
      y = rate(ssi,looppos); ylo=ratelo(ssi,looppos); yhi=ratehi(ssi,looppos);
      x=x+1; bar(x,y,'facecolor',cmap(looppos,:),'linewidth',1,'clipping','off'); ymax=max(ymax,yhi);
      ebw=0.25;line([x x],[ylo yhi],'color',[0 0 0],'linewidth',0.75,'clipping','off');
      line(x+[-1 1]*ebw,[1;1]*[ylo yhi],'color',[0 0 0],'linewidth',1,'clipping','off');
    end
  end
  ff;caxis([1 looplen]);colormap(cmap);set(gca,'fontsize',fs,'xtick',[]);  
  xlim([1 x+2]); ylim([0 0.8*ymax]); set(gca,'ytick',yts(minus0,:));
  text(xt,-diff(ylim)*0.18*ones(length(xt),1),ssbin.label,'fontsize',fs,'hor','cen');
  acgt='ACGT';title([acgt(minus0) 'pC'],'fontsize',20);
  if minus0==1,
    c=colorbar('location','north','fontsize',fs-4,'linewidth',1);
    set(c,'position',[0.22 0.87 0.10 0.02]);
    text(10,0.45*ymax,'loop position','fontsize',fs-4,'hor','cen');
  elseif minus0==4
    text(mean(xlim),-0.40*diff(ylim),'stem strength','fontsize',fs,'hor','cen');
    text(-6.5,7.3,'relative mutation frequency','fontsize',fs,'rotation',90);
  end
end
w=10;h=12;set(gcf,'papersize',[w h],'paperposition',[0.2 0.2 w-0.4 h-0.4]);
print_to_file([outdir '/figure2a.pdf']);

% FIGURE 2B

load([indir '/figure2b.mat'],'loop');

minlooplen = 3;
maxlooplen = 8;
nlens=maxlooplen-minlooplen+1;
figure(202),clf,fs=24;yts=[0 1.5 3;0 0.6 1.2;0 0.4 0.8;0 10 20];
cmap = green_orange_colormap(nlens);
for minus0=1:4, subplot(4,1,minus0),hold on,x=0;xt=[];ymax=0;
  for looplen=3:maxlooplen
    x=x+2; xt(end+1)=x+0.5+looplen/2;
    rate=loop.rate(looplen-minlooplen+1,:,minus0); sd=loop.sd(looplen-minlooplen+1,:,minus0); ratelo=max(0,rate-1.96*sd); ratehi=rate+1.96*sd;
    for looppos=1:looplen,y=rate(looppos);ylo=ratelo(looppos);yhi=ratehi(looppos);
      x=x+1; bar(x,y,1,'facecolor',cmap(1+floor((size(cmap,1)-1)*(looppos-1)/(looplen-1)),:),'linewidth',1,'clipping','off');ymax=max(ymax,max(yhi));
      ebw=0.25;line([x x],[ylo yhi],'color',[0 0 0],'linewidth',0.75,'clipping','off');
      line(x+[-1 1]*ebw,[1;1]*[ylo yhi],'color',[0 0 0],'linewidth',1,'clipping','off');
  end,end
  xlim([0.52 x+1.7]);ylim([0 0.86*ymax]);set(gca,'xtick',[],'ytick',yts(minus0,:),'fontsize',fs);ff
  text(xt,-diff(ylim)*0.18*ones(length(xt),1),num2cellstr(minlooplen:maxlooplen),'fontsize',fs,'hor','cen'); 
  acgt='ACGT';title([acgt(minus0) 'pC'],'fontsize',fs);
  if minus0==1,
    colormap(cmap); c=colorbar('location','north','fontsize',fs-4,'linewidth',1);
    set(c,'position',[0.76 0.87 0.10 0.02],'ticks',[]);
    text(0.89*max(xlim),0.48*ymax,'loop position','fontsize',fs-4,'hor','cen');
  elseif minus0==4
    text(mean(xlim),-0.55*diff(ylim),'loop length','fontsize',fs,'hor','cen');
    text(-0.14*diff(xlim),1.4*diff(ylim),'relative mutation frequency','fontsize',fs,'rotation',90);
  end
end
w=10;h=12;set(gcf,'papersize',[w h],'paperposition',[0.2 0.2 w-0.4 h-0.4]);
print_to_file([outdir '/figure2b.pdf']);


% FIGURE 3

load([indir '/figure3.mat'],'pat');
pat.ay = pat.logjit_frac_apo;
pat.ax = pat.apochar;
abcd='abcd';
for i=1:4
  figure(300+i);clf, fs=35; sz=60; hmin=0.3; hmax=0.7;
  hold on,set(gca,'fontsize',fs);
  if i==1, ttl='Tumor types';
    x = pat.ax; y = pat.ay; xlab='A3B \leftarrow mutation character \rightarrow A3A'; ylab='% APOBEC mutations';
    show = true(slength(pat),1); clr = pat.ttype_clr;
  elseif i==2, ttl = 'APOBEC mutations at TpC in DNA hairpins';
    x = pat.ax; y = pat.ay; xlab='A3B \leftarrow mutation character \rightarrow A3A'; ylab='% APOBEC mutations';
    clr = max(0,min(1,(pat.pcthp_tpc-hmin)/(hmax-hmin)));
    show = ~isnan(clr); clr = clr*[1 0 0]; % red
  elseif i==3, ttl = 'APOBEC mutations at VpC in DNA hairpins';
    x = pat.ax; y = pat.ay; xlab='A3B \leftarrow mutation character \rightarrow A3A'; ylab='% APOBEC mutations';
    clr = max(0,min(1,(pat.pcthp_vpc-hmin)/(hmax-hmin)));
    show = ~isnan(clr); clr = clr*[1 0 1]; % magenta
  elseif i==4, ttl = 'TpC and VpC hairpin mutations';
    x = pat.pcthp_tpc; xlab = 'TpC % hairpin mutations';
    y = pat.pcthp_vpc; ylab = 'VpC % hairpin mutations';
    clr = pat.ttype_clr; show = (~isnan(x)&~isnan(y));
  else error('?');
  end
  scatter(x(show),y(show),sz,clr(show,:),'filled');ff;set(gca,'fontsize',fs);%title(ttl,'fontsize',fs);
  xlabel(xlab,'fontsize',fs); ylabel(ylab,'fontsize',fs);
  if i~=4, xlim([-0.18 0.14]); set(gca,'xtick',[]); yt=[0 0.03 0.10 0.3 0.5 0.7 1];set(gca,'ytick',log10(0.01+yt),'yticklabel',yt*100); ylim([-2.1 0.1]);
  else xlim([0 2.6]);ylim([0 2.6]); set(gca,'xtick',[0 1 2],'ytick',[0 1 2]); end
  if i==2, w=9; h=7.5; elseif i==4, w=7.5; h=7.5; else w=8; h=6.6; end
  set(gcf,'papersize',[w h],'paperposition',[0.2 0.2 w-0.4 h-0.4]);
  print_to_file([outdir '/figure3' abcd(i) '.pdf']);
end

% tumor types colors legend
figure(399),clf,hold on,ff
u = {'Breast TNBC';'Breast non-TNBC';'Breast';'Bladder';'Cervical';'Lung squamous cell carcinoma';'Lung adenocarcinoma';'Head and neck';'Sarcoma';'Thyroid';
     'Endometrial';'Colorectal';'Esophageal';'Stomach';'Ovarian';'Liver';'Oral';'Bone';'Glioblastoma multiforme';'Prostate';'Pancreatic';'Melanoma'};
ui = listmap(u,pat.ttype_longname);
keep=~isnan(ui);u=u(keep);ui=ui(keep);
clrs = pat.ttype_clr(ui,:);
x=1;y=0.1;
for i=1:length(u)
  scatter(x,y,100,clrs(i,:),'filled');
  text(x+0.03*diff(xlim),y,u{i},'fontsize',15);
  y=y-0.015;
  if i==11, x=x+12;y=0.1; end
end
w=6;h=10;set(gcf,'papersize',[w h],'paperposition',[0.2 0.2 w-0.4 h-0.4]);
print_to_file([outdir '/figure3.ttype_legend.pdf']);

% FIGURE 5

load([indir '/figure5.mat'],'s');

figure(500),clf,hold on,ff,set(gca,'fontsize',16);
set(gcf,'position',[400 0 926 802])
pause(3);
x = log10(s.relrate_exp);
y = logjit(s.max_apo10_ct);
clr = nansub([0 0 0;1 0.5 0],1+(s.minus0~=4));
scatter(x,y,1);ylim([0.7 max(y)*1.03]);xlim([-1.7 2.7]); xl=xlim;yl=ylim;
rectangle('position',[log10(4) yl(1) xl(end)-log10(4) diff(yl)],'linestyle','none','facecolor',[1 0.95 0.95]);
idx=find(s.known_driver); scatter(x(idx),y(idx),200,[0.4 0.6 1],'filled');
scatter(x,y,40,clr,'filled')
line(log10(4)*[1 1],ylim,'color',[1 0 0],'linestyle','--','linewidth',2);
line(xlim,yl(1)*[1 1],'color',[0 0 0]); line(xl(1)*[1 1],ylim,'color',[0 0 0]);
xlabel('substrate optimality','fontsize',24);ylabel('# mutated patients','fontsize',24);
xpts=[0.05 0.2 0.5 1  4 10  20 50 100 300]; set(gca,'xtick',log10(xpts),'xticklabel',prefix(num2cellstr(xpts),' '));
ypts=[  4 6 10 15 20 30 50 100 200]; set(gca,'ytick',log10(1.5+ypts),'yticklabel',ypts);
lab=find(y>log10(1+6));textfit(x(lab)+1.5*xsp,y(lab),s.gene(lab),'fontsize',10);
w=12; h=9.5; set(gcf,'papersize',[w h],'paperposition',[0.2 0.2 w-0.4 h-0.4]);
print_to_file([outdir '/figure5.pdf']);


% FIGURE 6

load([indir '/figure6.mat'],'ssbin','loops');

figure(601);clf; hold on;ff
cmap = green_orange_colormap(slength(ssbin));
x=0;
for i=1:2 for j=1:slength(ssbin)
    if i==2 & j==1, x=x+5; end
    rate = ssbin.rate(j,i); sd = ssbin.sd(j,i); ratehi=rate+1.96*sd; ratelo=max(0,rate-1.96*sd);
    bar(x,rate,1,'facecolor',cmap(j,:),'linewidth',1);
    ebw=0.25;line([x x],[ratelo ratehi],'color',[0 0 0],'linewidth',0.75,'clipping','off');
    line(x+[-1 1]*ebw,[1;1]*[ratelo ratehi],'color',[0 0 0],'linewidth',1,'clipping','off');    
    x=x+1;
end; end
xlim([-4 36]); ylim([0 30]); ylabel('relative mutation freq.','fontsize',36);
set(gca,'tickdir','out','xtick',[],'fontsize',28);
text(6.5,-0.5,loops{1},'horizontalalignment','center','verticalalignment','top','fontsize',42);
text(25.5,-0.5,loops{2},'horizontalalignment','center','verticalalignment','top','fontsize',42);
w=8; h=7.2; set(gcf,'papersize',[w h],'paperposition',[0.2 0.2 w-0.4 h-0.4]);
print_to_file([outdir '/figure6a.pdf']);


figure(602);clf; hold on;ff
cmap = green_orange_colormap(slength(ssbin));
x=0;
for i=3:4 for j=1:slength(ssbin)
    if i==4 & j==1, x=x+5; end
    rate = ssbin.rate(j,i); sd = ssbin.sd(j,i); ratehi=rate+1.96*sd; ratelo=max(0,rate-1.96*sd);
    bar(x,rate,1,'facecolor',cmap(j,:),'linewidth',1);
    ebw=0.25;line([x x],[ratelo ratehi],'color',[0 0 0],'linewidth',0.75,'clipping','off');
    line(x+[-1 1]*ebw,[1;1]*[ratelo ratehi],'color',[0 0 0],'linewidth',1,'clipping','off');    
    x=x+1;
end; end
xlim([-4 36]); ylim([0 30]); ylabel('relative mutation freq.','fontsize',36);
set(gca,'tickdir','out','xtick',[],'fontsize',28);
text(6.5,-0.5,loops{3},'horizontalalignment','center','verticalalignment','top','fontsize',42);
text(25.5,-0.5,loops{4},'horizontalalignment','center','verticalalignment','top','fontsize',42);
w=8; h=7.2; set(gcf,'papersize',[w h],'paperposition',[0.2 0.2 w-0.4 h-0.4]);
print_to_file([outdir '/figure6b.pdf']);


% SUPP FIG 2A

load([indir '/figureS2a.mat'],'ssbin');

figure(10201),clf
cmap = green_orange_colormap(slength(ssbin));
acgt='ACGT';
for minus0=1:4
  subplot(2,2,minus0),cla, hold on
  rates = ssbin.rate(:,minus0); sds = ssbin.sd(:,minus0);
  for i=1:slength(ssbin)
    rate = rates(i); sd = sds(i); ratehi=rate+1.96*sd; ratelo=max(0,rate-1.96*sd);
    b = bar(i,rate,'facecolor',cmap(i,:));
    ebw=0.25;line([i i],[ratelo ratehi],'color',[0 0 0],'linewidth',0.75,'clipping','off');
    line(i+[-1 1]*ebw,[1;1]*[ratelo ratehi],'color',[0 0 0],'linewidth',1,'clipping','off');
  end
  set(gca,'fontsize',26,'xtick',[]);xlim([-0.5 slength(ssbin)+0.5]);
  ylim([0 1.5*max(rates)]); if minus0==3, ylim([0 2.4]); end
  title([acgt(minus0) 'pC'],'fontsize',20);ff
  if minus0==3, text(-4,0,'relative mutation rate','rot',90,'fontsize',20); end
end
w=10; h=12; set(gcf,'papersize',[w h],'paperposition',[0.2 0.2 w-0.4 h-0.4]);
print_to_file([outdir '/figureS2a.pdf']);


% SUPP FIG 2B

load([indir '/figureS2b.mat'],'ssbin','loops');

figure(10202);clf
subplot(2,1,1); ci=3;  % APOBEC>=10%
bar(ssbin.rate(:,:,ci)',0.9);colormap(green_orange_colormap)
ff;ylabel('relative mutation freq.','fontsize',20);set(gca,'ticklength',[0.005 0.005]);
set(gca,'xtick',1:length(loops),'xticklabel',loops,'fontsize',20);xlim([0.3 0.7+length(loops)]);
ylim([0 30]);line(xlim,[4 4],'color',[0 0 0],'linewidth',0.5,'linestyle','--');
title('APOBEC>10%','fontsize',20);
subplot(2,1,2); ci=1;   % APOBEC>=90%
bar(ssbin.rate(:,:,ci)',0.9);colormap(green_orange_colormap)
ff;ylabel('relative mutation freq.','fontsize',20);set(gca,'ticklength',[0.005 0.005]);
set(gca,'xtick',1:length(loops),'xticklabel',loops,'fontsize',20);xlim([0.3 0.7+length(loops)]);
ylim([0 30]);line(xlim,[4 4],'color',[0 0 0],'linewidth',0.5,'linestyle','--');
title('APOBEC>90%','fontsize',20);

w=20; h=10; set(gcf,'papersize',[w h],'paperposition',[0.2 0.2 w-0.4 h-0.4]);
print_to_file([outdir '/figureS2b.pdf']);

fprintf('Finished generating figures.\n');










