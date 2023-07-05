function X = survey_hairpins_rna(refdir,outdir,blockno)
% SURVEY HAIRPINS (RNA MODE)
% surveys a block of genes to characterize potential hairpin-forming sequences
% (considers only spliced transcribed RNA sequences) 

minlooplen=3;      % consider loops of at least 3nt
maxlooplen=11;     % consider loops of at most 11nt
maxstem=20;        % consider stems of up to 20bp (reducing from 100bp to 20bp reduces dead zones at ends of transcripts)
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
outfile = [outdir '/' X.block.name{blockno} '.rna.mat'];
if exist(outfile,'file')
  fprintf('Output file already exists.\n');
  return
end
write_textfile('IN_PROGRESS',outfile);
fprintf('BLOCK %d\n',blockno);

% load chromosome reference sequence from FASTA
fprintf('Loading reference sequence.\n');
chr = X.block.chr(blockno);
ref = load_genome_region(chr,1,X.chr.len(chr),refdir);

% set up list of genomic sites in this block
X.site = [];
X.site.chr = repmat(uint8(chr),X.block.len(blockno),1);
X.site.pos = uint32((X.block.st(blockno):X.block.en(blockno))');
X.site.ref = ref(X.block.st(blockno):X.block.en(blockno));
ns = X.block.len(blockno);

% annotate what the adjacent bases are (for later, when we take slices of the full list)
X.site.minus2  = [0;0;0;X.site.ref(1:end-3)];
X.site.minus1  = [0;0;X.site.ref(1:end-2)];
X.site.minus0  = [0;X.site.ref(1:end-1)];
X.site.plus1   = [X.site.ref(2:end);0];
X.site.plus2   = [X.site.ref(3:end);0;0];
X.site.plus3   = [X.site.ref(4:end);0;0;0];

X.site.gene = zeros(ns,1,'uint16');
fs = {'looplen','looppos','bulgepos','nbp','ngc','mmp','ss'};
for i=1:length(fs), X.site.(fs{i}) = -ones(ns,1,'int8'); end

% IDENTIFY WHICH TRANSCRIPTS OVERLAP THIS BLOCK
block_st = X.block.st(blockno);
block_en = X.block.en(blockno);
tx = find(X.tx.chr==chr & X.tx.tx_start<=block_en & X.tx.tx_end>=block_st);
offset = block_st-1;
mask = maxstem+maxlooplen-1;   % mask block edges (interblock overlap of 1kb will allow trimming later when blocks are gathered)

% HAIRPIN SURVEY (RNA MODE)

fprintf('Surveying hairpins (RNA MODE).\n');

fprintf('Transcript: ');
for txi=1:length(tx), fprintf('%d/%d ',txi,length(tx));
  i=tx(txi);

  % BUILD TRANSCRIPT
  R=[];
  R.ref = []; R.pos = []; nframeshifts = 0;
  gene = X.tx.gene_idx(i);
  plusstrand = strcmp(X.tx.strand{i},'+');
  if plusstrand, forfrom=1; forstep=+1; forto=X.tx.n_exons(i);
  else forfrom=X.tx.n_exons(i); forstep=-1; forto=1;
  end
  for e=forfrom:forstep:forto
    st = X.tx.exon_starts{i}(e); en = X.tx.exon_ends{i}(e);
    if plusstrand, p = (st:en)'; d = ref(st:en);
    else p = (en:-1:st)'; d = 5-ref(en:-1:st)';
    end
    R.ref = [R.ref;d]; R.pos = [R.pos;p];
  end
  
  gc = int8(R.ref==3 | R.ref==2); % (makes things quicker in loop)

  fs = {'looplen','looppos','bulgepos','nbp','ngc','mmp','ss'};
  nr = slength(R);
  for i=1:length(fs), R.(fs{i}) = -ones(nr,1,'int8'); end

  for looplen=minlooplen:maxlooplen
    for looppos=1-maxintrastem:looplen+maxintrastem
      for bulgepos=[0 (-minbulgepos):-1:-(maxbulgepos) minbulgepos:maxbulgepos]  % 0 = no bulge
        stem_open = false(nr,1); stem_open((1+mask):(nr-mask))=true; can_pair = false(nr,1);
        nbp = zeros(nr,1,'int8'); ngc = zeros(nr,1,'int8'); mmp = zeros(nr,1,'int8');
        lft=looppos; rgt=1+looplen-looppos; done=false;
        k_switch_to_list_mode = 3;  % after first two basepairs, switch from global mode to list mode (more efficient)
        for k=1:maxstem
          if k<k_switch_to_list_mode    % optimized for early steps (when all/most of transcript is under consideration)
            can_pair((1+mask):(nr-mask)) = R.ref((1+mask-lft):(nr-mask-lft)) == 5-R.ref((1+mask+rgt):(nr-mask+rgt));
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
            can_pair = R.ref(ii-lft) == 5-R.ref(ii+rgt);
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
        idx = find(ss>R.ss);
        for i=1:length(fs), f=fs{i}; tmp=eval(f); if length(tmp)==1, R.(f)(idx) = tmp; else R.(f)(idx) = tmp(idx); end; end
      end
    end
  end

  % INTEGRATE INTO SITE LIST
  for ri=1:nr
    pos = R.pos(ri);
    xi = pos-offset;
    if xi>=1 & xi<=ns
      if R.ss(ri)>X.site.ss(xi)
        X.site.gene(xi) = gene;
        for fi=1:length(fs),f=fs{fi}; X.site.(f)(xi) = R.(f)(ri); end
      end
    end
  end

end  % next transcript

% save
fprintf('\nSaving block %d... ',blockno);
X = keep_fields(X,'site');
save(outfile,'X','-v7.3');
fprintf('\n');




