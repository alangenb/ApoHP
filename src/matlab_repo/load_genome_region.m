function d = load_genome_region(chr,st,en,refdir)
% load_genome_region(chr,st,en,refdir)
%
% loads genome region from FASTA file
% --> returns as numerical representation
  
if nargin<4, error('requires 4 arguments: chr, start, end, refdir'); end

if ~isnumeric(chr), error('chromosome should be numeric'); end
if ~isnumeric(st), error('start should be numeric'); end
if ~isnumeric(en), error('end should be numeric'); end
if length(chr)>1 || length(st)>1 || length(en)>1, error('multiple ranges not supported'); end
if chr<1 || st<1 || en<1, error('chr, start, and end should be 1-based'); end
if ~ischar(refdir), error('refdir should be a string'); end

X = load_genome_info(refdir);

if chr>slength(X.chr), error('chr is larger than the number of chromosomes'); end
if st>X.chr.len(chr), error('start coordinate is larger than length of chromosome'); end
if en>X.chr.len(chr), error('end coordinate is larger than length of chromosome'); end

F = load_fasta(X.chr.fasta{chr});
if length(F.seq{1})~=X.chr.len(chr), error('FASTA sequence length mismatch'); end

d = uint8(listmap(upper(F.seq{1}(st:en))','ACGT'));   % convert to 1/2/3/4 = A/C/G/T   (0 = N)












  
