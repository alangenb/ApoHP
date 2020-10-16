function X = load_genome_info(refdir)
% X = load_genome_info(refdir)
%
% --> loads reference genome information from the files located in refdir:
%
% chromInfo.txt = list of chromosome names and their sizes
%                 --> tab-delimited, no header line, only the first two columns are used:
%                     first column = chromosome name, e.g. chr1 or chrX: must match fasta names and "chr" column in mutation data
%                     second column = chromosome length: must match fasta lengths
%                 e.g. as downloaded from http://hgdownload.cse.ucsc.edu/goldenpath/hg19/database/chromInfo.txt.gz
%                 --> can manually edit this to keep only the canonical 24 human chromsomes
%                 --> is also useful to sort it in canonical order, 1-22,X,Y
%
% chr*.fa = set of fasta files for each chromosome
%           e.g. as downloaded from http://hgdownload.cse.ucsc.edu/goldenpath/hg19/chromosomes/chr1.fa.gz
%           --> names must match chromInfo.txt
%
% refGene.txt = set of gene definitions for assigning gene names and coding/noncoding status
%             e.g. as downloaded from http://hgdownload.cse.ucsc.edu/goldenpath/hg19/database/refGene.txt.gz
%
% known_driver_genes.txt = list of known driver gene names
%            --> text file (no header line) with one gene name per line
%
% --> returns X with fields:
%
% X.chr   = list of chromosomes, names, and lengths: will use this for downstream numerical chromosome indexing
% X.block = list of 10Mb chromosome chunks (with 500bp overlap)
% X.gene  = list of genes
% X.tx    = list of transcripts
%
% --> also saves this object X as new file "ref.mat", which will be detected and loaded on future calls.

refmat = [refdir '/ref.mat'];
if exist(refmat,'file')
  load(refmat,'X');
  return
end

% verify required input files exist
chromInfo = [refdir '/chromInfo.txt'];
refGene = [refdir '/refGene.txt'];
knownDrivers = [refdir '/known_driver_genes.txt'];
demand_files({chromInfo,refGene,knownDrivers});

fprintf('Loading reference genome info from refdir = %s\n',refdir);

X=[];

blocksize = 10e6;        % 10Mb genomic blocks for splitting up the hairpin survey
splicesiteflank = 20;    % when defining exonic territory (for excluding coding mutations) expand exons 20bp into flanking territory
geneflank = 100;         % when defining gene territory include this much flanking territory

% CHROMOSOMES 

X.chr = load_struct_noheader(chromInfo);
X.chr = keep_fields(X.chr,{'col1','col2'});
X.chr = rename_fields(X.chr,{'col1','col2'},{'name','len'});
X.chr = make_numeric(X.chr,'len');
if slength(X.chr)==0, error('Failed to parse %s',chromInfo); end
if slength(X.chr)>255, error('To handle >255 chromosomes, need to modify code to use uint16 instead of uint8'); end

% verify that all fasta files exist and contain sequence of the correct length
fprintf('Checking FASTA files...\n');
X.chr.fasta = prefix(suffix(X.chr.name,'.fa'),[refdir '/']);
demand_files(X.chr.fasta);
for chr=1:slength(X.chr)
  fprintf('  %s  ',X.chr.name{chr});
  F = load_fasta(X.chr.fasta{chr});
  if slength(F)==0, error('Invalid FASTA file'); end
  if slength(F)>1, error('FASTA files should not contain multiple sequences.'); end
  if length(F.seq{1})~=X.chr.len(chr), error('Sequence length mismatch.'); end
  fprintf('OK\n');
end

% BLOCKS 

X.block = cell(slength(X.chr),1);
for i=1:slength(X.chr)
  X.block{i}.st = (1:blocksize:X.chr.len(i))';
  X.block{i}.en = min(X.chr.len(i),X.block{i}.st+blocksize-1);
  X.block{i}.chr = repmat(i,length(X.block{i}.st),1);
end
X.block = concat_structs(X.block);
X.block.name = prefix(num2cellstr(1:slength(X.block)),'block');
X.block = order_fields(X.block,{'name','chr','st','en'});
% make blocks overlap by 1000 bases (to avoid edge effects in the hairpin survey)
X.block.st = max(1,X.block.st-1000);
X.block.len = X.block.en-X.block.st+1;
X.block.st_trim = X.block.st+500; X.block.en_trim = X.block.en-500;
first_on_chr = (X.block.st==1); X.block.st_trim(first_on_chr)=X.block.st(first_on_chr);
last_on_chr = [X.block.st(2:end)==1;true]; X.block.en_trim(last_on_chr)=X.block.en(last_on_chr);
X.block.len_trim = X.block.en_trim-X.block.st_trim+1;

% GENES 

infile = [refdir '/refGene.txt'];
fprintf('Loading %s\n',infile);
T = read_table(infile,'%f%s%s%s%f%f%f%f%f%s%s%f%s%s%s%s%s',char(9),0,'whitespace','\b\r');
[X.tx.id,X.tx.transcript,X.tx.chr,X.tx.strand,X.tx.tx_start,X.tx.tx_end,X.tx.code_start,X.tx.code_end,...
   X.tx.n_exons,X.tx.exon_starts,X.tx.exon_ends,tmp,X.tx.gene,tmp,tmp,X.tx.exon_frames,X.tx.version]  = deal(T.dat{:});

fprintf('Parsing refGene... ');
for i=1:slength(X.tx)
  if ~mod(i,10000), fprintf('%d/%d ',i,slength(X.tx)); end
  X.tx.tx_start(i) = X.tx.tx_start(i)+1;
  X.tx.code_start(i) = X.tx.code_start(i)+1;
  X.tx.exon_starts{i} = str2double(split(X.tx.exon_starts{i}(1:end-1),','))+1;
  X.tx.exon_ends{i} = str2double(split(X.tx.exon_ends{i}(1:end-1),','));
  X.tx.exon_frames{i} = str2double(split(X.tx.exon_frames{i}(1:end-1),','));
end
fprintf('\n')

% add gene flanks
X.tx.gene_start = X.tx.tx_start;
X.tx.gene_end = X.tx.tx_end;
X.tx.gene_start = X.tx.gene_start - geneflank;
X.tx.gene_end = X.tx.gene_end + geneflank;

% find exons and add splice flanks
X.tx.tx_len = zeros(slength(X.tx),1);
X.tx.code_len = zeros(slength(X.tx),1);
X.tx.coding_starts = cell(slength(X.tx),1);
X.tx.coding_ends = cell(slength(X.tx),1);
X.tx.n_coding_regions = nan(slength(X.tx),1);
for i=1:slength(X.tx)
  st = X.tx.exon_starts{i};
  en = X.tx.exon_ends{i};
  X.tx.tx_len(i) = sum(en-st+1);  
  st = max(st,X.tx.code_start(i));
  en = min(en,X.tx.code_end(i));
  noncoding = (st>en); st(noncoding)=[]; en(noncoding)=[];
  st = st - splicesiteflank;
  en = en + splicesiteflank;
  X.tx.coding_starts{i} = st;
  X.tx.coding_ends{i} = en;
  X.tx.n_coding_regions(i) = length(st);
  X.tx.code_len(i) = sum(en-st+1);
end

% remove transcripts on noncanonical chromosomes and non-coding transcripts
X.tx.chr = listmap(X.tx.chr,X.chr.name);
X.tx = reorder_struct(X.tx,~isnan(X.tx.chr));
X.tx = reorder_struct(X.tx,X.tx.code_len>0);

% get list of unique genes
[X.gene.name ui X.tx.gene_idx] = unique(X.tx.gene);
if slength(X.gene)>65535, error('To handle >65535 genes, need to modify code to use uint32 instead of uint16'); end
  
% mark known driver genes
D = load_lines(knownDrivers);
X.gene.known_driver = ismember(X.gene.name,D);

fprintf('Saving completed genome info object to %s\n',refmat);
save(refmat,'X');


