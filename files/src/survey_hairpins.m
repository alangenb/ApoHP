function X = survey_hairpins(refdir,outdir,blockno)
% SURVEY HAIRPINS
% surveys a block of the genomic reference sequence to characterize potential hairpin-forming sequences

minlooplen=3;      % consider loops of at least 3nt
maxlooplen=11;     % consider loops of at most 11nt
maxstem=100;       % consider stems of up to 100bp
minbulgepos=2;     % consider hairpins with bulges at least 2bp away from loop
maxbulgepos=8;     % consider hairpins with bulges up to 8bp away from loop
minmismatchpos=2;  % consider hairpins with mismatches in the stem at least 2 bp away from loop
maxintrastem=1;    % consider positions to be intrastem if they are within 1 bp of loop

if nargin~=3, error('should take three input arguments: refdir, outdir, blockno'); end

if ~ischar(refdir), error('refdir should be a string'); end
if ~ischar(outdir), error('outdir should be a string'); end
if ischar(blockno), blockno=str2double(blockno); end
if ~isnumeric(blockno) || isnan(blockno) || blockno<1 || blockno~=round(blockno), error('invalid blockno'); end

X = load_genome_info(refdir);
if blockno>slength(X.block), error('blockno > #blocks'); end

ede(outdir);
outfile = [outdir '/' X.block.name{blockno} '.mat'];
if exist(outfile,'file')
  fprintf('Output file already exists.\n');
  return
end
write_textfile('IN_PROGRESS',outfile);
fprintf('BLOCK %d\n',blockno);

% set up list of genomic sites in this block
X.site = [];
X.site.chr = repmat(uint8(X.block.chr(blockno)),X.block.len(blockno),1);
X.site.pos = uint32((X.block.st(blockno):X.block.en(blockno))');
ns = X.block.len(blockno);

% load reference sequence from FASTA
fprintf('Loading reference sequence.\n');
X.site.ref = load_genome_region(X.block.chr(blockno),X.block.st(blockno),X.block.en(blockno),refdir);

% annotate what the adjacent bases are (for later, when we take slices of the full list)
X.site.minus2  = [0;0;0;X.site.ref(1:end-3)];
X.site.minus1  = [0;0;X.site.ref(1:end-2)];
X.site.minus0  = [0;X.site.ref(1:end-1)];
X.site.plus1   = [X.site.ref(2:end);0];
X.site.plus2   = [X.site.ref(3:end);0;0];
X.site.plus3   = [X.site.ref(4:end);0;0;0];

% annotate gene zones: 0=IGR 1=intron 2=UTR/promoter 3=exon/spliceflank
fprintf('Annotating gene regions.\n');
X.site.gene = zeros(ns,1,'uint16');
X.site.zone = zeros(ns,1,'uint8');
offset = X.block.st(blockno)-1;
for ti=1:slength(X.tx), if X.tx.chr(ti)~=X.block.chr(blockno), continue; end
  dst = max(1,X.tx.gene_start(ti)-offset); den = min(ns,X.tx.gene_end(ti)-offset); X.site.gene(dst:den)=X.tx.gene_idx(ti); X.site.zone(dst:den)=2;
  dst = max(1,X.tx.code_start(ti)-offset); den = min(ns,X.tx.code_end(ti)-offset); X.site.gene(dst:den)=X.tx.gene_idx(ti); X.site.zone(dst:den)=1;
end
for ti=1:slength(X.tx), if X.tx.chr(ti)~=X.block.chr(blockno), continue; end 
  for ei=1:X.tx.n_coding_regions(ti)
    dst = max(1,X.tx.coding_starts{ti}(ei)-offset); den = min(ns,X.tx.coding_ends{ti}(ei)-offset); X.site.gene(dst:den)=X.tx.gene_idx(ti); X.site.zone(dst:den)=3;
  end
end

% HAIRPIN SURVEY

fprintf('Surveying hairpins.\n');

mask = maxstem+maxlooplen-1;   % mask block edges (interblock overlap of 1kb will allow trimming later when blocks are gathered)
gc = int8(X.site.ref==3 | X.site.ref==2); % (makes things quicker in loop)

fs = {'looplen','looppos','bulgepos','nbp','ngc','mmp','ss'};
for i=1:length(fs), X.site.(fs{i}) = -ones(ns,1,'int8'); end

fprintf('Loop length: ');
for looplen=minlooplen:maxlooplen, fprintf(' %d/%d',looplen, maxlooplen);
  for looppos=1-maxintrastem:looplen+maxintrastem
    for bulgepos=[0 (-minbulgepos):-1:-(maxbulgepos) minbulgepos:maxbulgepos]  % 0 = no bulge
      stem_open = false(ns,1); stem_open((1+mask):(ns-mask))=true; can_pair = false(ns,1);
      nbp = zeros(ns,1,'int8'); ngc = zeros(ns,1,'int8'); mmp = zeros(ns,1,'int8');
      lft=looppos; rgt=1+looplen-looppos; done=false;
      k_switch_to_list_mode = 3;  % after first two basepairs, switch from global mode to list mode (more efficient)
      for k=1:maxstem
        if k<k_switch_to_list_mode    % optimized for early steps (when all/most of genome is under consideration)
          can_pair((1+mask):(ns-mask)) = X.site.ref((1+mask-lft):(ns-mask-lft)) == 5-X.site.ref((1+mask+rgt):(ns-mask+rgt));
          extend = (stem_open & can_pair);
          nbp(extend)=nbp(extend)+1;
          tmp=[zeros(lft,1,'int8');gc(1:end-lft)]; ngc(extend) = ngc(extend) + tmp(extend);
          if k>minmismatchpos
            stem_open(~can_pair & mmp>0)=0;
            extend_over_mismatch = (stem_open & ~can_pair & mmp==0);
            mmp(extend_over_mismatch)=k;
          else
            stem_open(~can_pair)=0;
          end
          done = ~any(stem_open);
        else      % same procedure optimized for later steps (when smaller lists of positions remain)
          ii = find(stem_open);
          can_pair = X.site.ref(ii-lft) == 5-X.site.ref(ii+rgt);
          fe = ii(can_pair);
          if ~isempty(fe)
            nbp(fe)=nbp(fe)+1;
            ngc(fe)=ngc(fe)+gc(fe-lft);
          end
          if k>minmismatchpos
            stem_open(ii(~can_pair & mmp(ii)>0))=0;
            feomm = ii(~can_pair & mmp(ii)==0);
            mmp(feomm) = k;
          else
            stem_open(ii(~can_pair))=0;
          end
          done = ~any(stem_open(ii));
        end
        lft=lft+1; if k==-bulgepos, lft=lft+1; end            % bypass left bulge
        rgt=rgt+1; if k==+bulgepos, rgt=rgt+1; end            % bypass right bulge
        if done,break;end
      end
     
      % remove terminal mismatches
      mmp(mmp==nbp+1) = 0;
      % score hairpin strength
      ss = nbp + 2*ngc;                                       % base strength
      ss = ss + max(-100*mmp, -max(1,13-2*mmp));              % mismatch penalty
      if bulgepos~=0, ss = ss - 6; end                        % bulge penalty
      % assign hairpins, overwriting weaker with stronger
      idx = find(ss>X.site.ss);
      for i=1:length(fs), f=fs{i}; tmp=eval(f); if length(tmp)==1, X.site.(f)(idx) = tmp; else X.site.(f)(idx) = tmp(idx); end; end
    end
  end
end

% for G/A positions on reference strand, flip annotations to opposite strand so that we're always annotating C/T's
ga = (X.site.ref==3 | X.site.ref==1);
tmp=X.site.plus1(ga); X.site.plus1(ga)=5-X.site.minus0(ga); X.site.minus0(ga)=5-tmp;
tmp=X.site.plus2(ga); X.site.plus2(ga)=5-X.site.minus1(ga); X.site.minus1(ga)=5-tmp;
tmp=X.site.plus3(ga); X.site.plus3(ga)=5-X.site.minus2(ga); X.site.minus2(ga)=5-tmp;
X.site.looppos(ga) = X.site.looplen(ga)-X.site.looppos(ga)+1;
X.site.bulgepos(ga) = -X.site.bulgepos(ga);
  
% save
fprintf('\nSaving block %d... ',blockno);
X = keep_fields(X,'site');
save(outfile,'X','-v7.3');
fprintf('\n');




