function make_data_files_for_figures(X,outdir)

ede(outdir);

fprintf('Generating data for figures.\n');

% LEGO PLOTS

nmf = X.nmf;
save([outdir '/figure.lego_plots.mat'],'nmf');

% FIGURE 2A

ssbin=[];
ssbin.min = [0 7:2:19]';
ssbin.max = [6 8:2:19 inf]';
ssbin.label = num2cellstr(ssbin.min); ssbin.label{end}(end+1)='+'; ssbin.label{1}=['<' ssbin.label{2}];
looplen=4;
ci=2;
ssbin.n = nan(slength(ssbin),looplen,4);
ssbin.N = nan(slength(ssbin),looplen,4);
i1 = find(X.loop.looplen==looplen & X.loop.minus2+X.loop.minus1+X.loop.plus1+X.loop.plus2==0);  % (non-degenerate rows only)
for minus0=1:4
  i2=i1(X.loop.minus0(i1)==minus0);
  for looppos=1:looplen
    i3=i2(X.loop.looppos(i2)==looppos);
    for i=1:slength(ssbin)
      s=find(X.ssbin.min>=ssbin.min(i) & X.ssbin.max<=ssbin.max(i));
      ssbin.n(i,looppos,minus0) = sum(X.loop.n(i3,s,ci),2);
      ssbin.N(i,looppos,minus0) = sum(X.loop.N(i3,s),2);
    end
  end
end
[ssbin.rate ssbin.sd] = relative_ratio_and_sd(ssbin.n,ssbin.N,X.cohort.tca_rate(ci));
save([outdir '/figure2a.mat'],'ssbin');
% save as textfile version
q=keep_fields(ssbin,{'min','max','label'});
acgt='ACGT';
for minus0=1:4
  for looppos=1:looplen
    fld = [acgt(minus0) 'pC_looppos' num2str(looppos)];
    q.([fld '_rate']) = ssbin.rate(:,looppos,minus0);
    q.([fld '_sd']) = ssbin.sd(:,looppos,minus0);
  end
end
save_struct(q,[outdir '/figure2a.txt']);

% FIGURE 2B

minlooplen = 3;
maxlooplen = 8;
ssmin = 18;
ci=2;
nlens = maxlooplen-minlooplen+1;
loop.n = nan(nlens,maxlooplen,4);
loop.N = nan(nlens,maxlooplen,4);
s=find(X.ssbin.min>=ssmin);
i1 = find(X.loop.minus2+X.loop.minus1+X.loop.plus1+X.loop.plus2==0);  % (non-degenerate rows only)
for minus0=1:4
  i2=i1(X.loop.minus0(i1)==minus0);
  for looplen=minlooplen:maxlooplen
    i3=i2(X.loop.looplen(i2)==looplen);
    for looppos=1:looplen
      i4=i3(X.loop.looppos(i3)==looppos);
      loop.n(looplen-minlooplen+1,looppos,minus0) = sum(X.loop.n(i4,s,ci),2);
      loop.N(looplen-minlooplen+1,looppos,minus0) = sum(X.loop.N(i4,s),2);
    end
  end
end
[loop.rate loop.sd] = relative_ratio_and_sd(loop.n,loop.N,X.cohort.tca_rate(ci));
save([outdir '/figure2b.mat'],'loop');
% save as textfile version
q=[]; q.minus0=[]; q.looplen=[]; q.looppos=[]; q.rate=[]; q.sd=[];
acgt='ACGT';
for minus0=1:4
  for looplen=minlooplen:maxlooplen
    for looppos=1:looplen
      q.minus0(end+1,1) = minus0;
      q.looplen(end+1,1) = looplen;
      q.looppos(end+1,1) = looppos;
      q.rate(end+1,1) = loop.rate(looplen-minlooplen+1,looppos,minus0);
      q.sd(end+1,1) = loop.sd(looplen-minlooplen+1,looppos,minus0);
    end
  end
end
save_struct(q,[outdir '/figure2b.txt']);

% FIGURE 3

pat = X.pat;
pat = sort_struct(pat,'nmut');
pat = reorder_struct(pat,pat.nmut>=500);
save([outdir '/figure3.mat'],'pat');
% save as textfile version
ppp = pat;
ppp.ttype_clr_red = round(255*ppp.ttype_clr(:,1));
ppp.ttype_clr_green = round(255*ppp.ttype_clr(:,2));
ppp.ttype_clr_blue = round(255*ppp.ttype_clr(:,3));
ppp = rmfield(ppp,{'ttype_clr','ttype_idx','nmf','nmf_norm','nmf_nmut','ci'});
save_struct(ppp,[outdir '/figure3.txt']);                    

% FIGURE 5

if isfield(X.site,'wxs_apo10_ct')   % WXS provided
  X.site.max_apo10_pct = max(X.site.wgs_apo10_pct,X.site.wxs_apo10_pct);
  X.site.max_apo10_ct = max(X.site.wgs_apo10_ct,X.site.wxs_apo10_ct);
else                                % WXS not provided
  X.site.max_apo10_pct = X.site.wgs_apo10_pct;
  X.site.max_apo10_ct = X.site.wgs_apo10_ct;
end
s = reorder_struct(X.site,X.site.max_apo10_ct>=2);
s.known_driver = nansub(X.gene.known_driver,s.gene);
s.gene = nansub(X.gene.name,s.gene);
fs=fieldnames(s);for i=1:length(fs),f=fs{i};if isnumeric(s.(f)), s.(f)=double(s.(f)); end, end
s = sort_struct(s,'max_apo10_pct',-1);
s = reorder_struct(s,1:100);
save([outdir '/figure5.mat'],'s');
% save as textfile version
ss = rmfield(s,{'n'});
save_struct(ss,[outdir '/figure5.txt']);

% FIGURE 6

ssbin=[];
ssbin.min = [0 4:2:21]';
ssbin.max = [3 5:2:20 inf]';
ssbin.label = num2cellstr(ssbin.min); ssbin.label{end}(end+1)='+'; ssbin.label{1}=['<' ssbin.label{2}];
loops = {'GAAC','AAAC','GTTC','ATTC'}'; acgt={'A','C','G','T'};
ci=4;
ssbin.n = nan(slength(ssbin),length(loops));
ssbin.N = nan(slength(ssbin),length(loops));
looplen=4;
looppos=4;
i1 = find(X.loop.looplen==looplen & X.loop.looppos==looppos & X.loop.plus1+X.loop.plus2==0);
for loop=1:length(loops)
  i2=i1(X.loop.minus2(i1)==listmap(loops{loop}(1),acgt) & ...
        X.loop.minus1(i1)==listmap(loops{loop}(2),acgt) & ...
        X.loop.minus0(i1)==listmap(loops{loop}(3),acgt));
  for i=1:slength(ssbin)
    s=find(X.ssbin.min>=ssbin.min(i) & X.ssbin.max<=ssbin.max(i));
    ssbin.n(i,loop) = sum(X.loop.n(i2,s,ci),2);
    ssbin.N(i,loop) = sum(X.loop.N(i2,s),2);
  end
end
[ssbin.rate ssbin.sd] = relative_ratio_and_sd(ssbin.n,ssbin.N,X.cohort.tca_rate(ci));
save([outdir '/figure6.mat'],'ssbin','loops');
% save as textfile version
q=[]; q.stem_strength_bin = ssbin.label;
for loop=1:length(loops)
  q.([loops{loop} '_rate']) = ssbin.rate(:,loop);
  q.([loops{loop} '_sd']) = ssbin.sd(:,loop);
end
save_struct(q,[outdir '/figure6.txt']);

% SUPP FIGURE 2A

ssbin=[];
ssbin.min = [0 4:2:22]';
ssbin.max = [3 5:2:21 inf]';
ssbin.label = num2cellstr(ssbin.min); ssbin.label{end}(end+1)='+'; ssbin.label{1}=['<' ssbin.label{2}];
ssbin.n = nan(slength(ssbin),4);
ssbin.N = nan(slength(ssbin),4);
ci=1;
i1 = find(X.loop.looplen==3 & X.loop.looppos==3 & X.loop.plus1+X.loop.plus2+X.loop.minus1+X.loop.minus2==0);
for minus0=1:4
  i2 = i1(X.loop.minus0(i1)==minus0);
  for i=1:slength(ssbin)
    s=find(X.ssbin.min>=ssbin.min(i) & X.ssbin.max<=ssbin.max(i));
    ssbin.n(i,minus0) = fullsum(X.loop.n(i2,s,ci));
    ssbin.N(i,minus0) = fullsum(X.loop.N(i2,s));
  end
end
[ssbin.rate ssbin.sd] = relative_ratio_and_sd(ssbin.n,ssbin.N,X.cohort.tca_rate(ci));
save([outdir '/figureS2a.mat'],'ssbin');
% save as textfile version
q=[]; q.stem_strength_bin = ssbin.label;
acgt='ACGT';
for minus0=1:4, site_type = [acgt(minus0) 'pC'];
  q.([site_type '_rate']) = ssbin.rate(:,minus0);
  q.([site_type '_sd']) = ssbin.sd(:,minus0);
end
save_struct(q,[outdir '/figureS2A.txt']);


% SUPP FIGURE 2B

ssbin=[];
ssbin.min = [0 4:2:22]';
ssbin.max = [3 5:2:21 inf]';
ssbin.label = num2cellstr(ssbin.min); ssbin.label{end}(end+1)='+'; ssbin.label{1}=['<' ssbin.label{2}];
ssbin.n = nan(slength(ssbin),16,slength(X.cohort));
ssbin.N = nan(slength(ssbin),16);
i1 = find(X.loop.looplen==3 & X.loop.looppos==3 & X.loop.plus1+X.loop.plus2+X.loop.minus2==0);
loops = [];
acgt='ACGT'; li=1;
for minus0=1:4
  i2 = i1(X.loop.minus0(i1)==minus0);
  for minus1=1:4
    i3 = i2(X.loop.minus1(i2)==minus1);
    loops{li,1} = acgt([minus1 minus0 2]);
    for i=1:slength(ssbin)
      s=find(X.ssbin.min>=ssbin.min(i) & X.ssbin.max<=ssbin.max(i));
      ssbin.n(i,li,:) = nansum(nansum(X.loop.n(i3,s,:),2),1);
      ssbin.N(i,li) = nansum(nansum(X.loop.N(i3,s),2),1);
    end
    li=li+1;
  end
end
[ssbin.rate ssbin.sd] = relative_ratio_and_sd(ssbin.n,ssbin.N,shiftdim(X.cohort.tca_rate,-2));
save([outdir '/figureS2b.mat'],'ssbin','loops');
% save as textfile version
q=[]; q.cohort={}; q.site_type={}; q.stem_strength_bin={}; q.rate=[]; q.sd=[];
for ci=1:4, li=1;
  for minus0=1:4
    for minus1=1:4
      for i=1:slength(ssbin)
        q.cohort{end+1,1} = X.cohort.name{ci};
        q.site_type{end+1,1} = loops{li};
        q.stem_strength_bin{end+1,1} = ssbin.label{i};
        q.rate(end+1,1) = ssbin.rate(i,li,ci);
        q.sd(end+1,1) = ssbin.sd(i,li,ci);
      end
      li=li+1;
    end
  end
end
save_struct(q,[outdir '/figureS2B.txt']);

fprintf('Finished generating data for figures.\n');


