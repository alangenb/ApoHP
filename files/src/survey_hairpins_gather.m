function X = survey_hairpins_gather(refdir,datdir)
% gather C's from all blocks
% --> write full data object <datdir>/allCs.mat
% --> also saves <datdir>/best_hairpin_sites.mat
% --> also saves <datdir>/best_hairpin_sites.bed
  
X = load_genome_info(refdir);

X.block.mat = prefix(suffix(X.block.name,'.mat'),[datdir '/']);
demand_files(X.block.mat);
X.block.matsize = get_filesize(X.block.mat);
bad = (X.block.matsize<100);
if any(bad), error('%d/%d blocks did not finish',sum(bad),slength(X.block)); end

fprintf('Loading all blocks.\n');

for blockno=1:slength(X.block),fprintf('BLOCK %d\n',blockno);
  tmp = load(X.block.mat{blockno},'X');
  tmp = tmp.X.site;
  % trim interblock overlaps 
  tmp = reorder_struct(tmp,(1+X.block.st_trim(blockno)-X.block.st(blockno)):(1+X.block.en_trim(blockno)-X.block.st(blockno)));
  % keep C:G basepairs only
  X.site{blockno} = reorder_struct(tmp,tmp.ref==2|tmp.ref==3);
end
X.site=concat_structs(X.site);

% save list of best hairpin sites
besthp = reorder_struct(X.site,(X.site.looplen==3 & X.site.looppos==3 & X.site.ss>=10) | ((X.site.looplen==4 | X.site.looplen==5) & X.site.looppos==4 & X.site.ss>=12));
save([datdir '/best_hairpin_sites.mat'],'besthp','-v7.3');
% also save as BED file
B=[];
B.chr = X.chr.name(besthp.chr);
B.st = besthp.pos-1;
B.en = besthp.pos;
for i=slength(B):-1:1, B.name{i,1} = [num2str(besthp.looppos(i)) '/' num2str(besthp.looplen(i)) ' (ss=' num2str(besthp.ss(i)) ')']; end
save_struct_noheader(B,[datdir '/best_hairpin_sites.bed']);

% save full data
save([datdir '/allCs.mat'],'X','-v7.3');

