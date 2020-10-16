function X = load_mutation_data(infile,refdir,outfile,ttypefile)
% X = load_mutation_data(infile,refdir,outfile,ttypefile)
% --> loads and preprocesses mutations
%
% input file should be tab-delimited text file
% with one header line containing the names of columns
%
% required columns:
% patient = name of patient
% chr = chromosome of mutation, as name (e.g. chr1, chrX)
%       --> will map these to the chromosome information in refdir
% pos = position of mutation (numeric, 1-based)
% ref = reference base (A/C/G/T)
%       --> will check to see that this matches information in refdir
% alt = mutated base (A/C/G/T)
%
% optional columns:
% ttype = tumor type of patient
%       --> will look up tumortype longnames and colors from "ttypefile" if available
%
% --> also saves this object X as <outfile>, which will be detected and loaded on future calls.

if nargin<3, error('requires infile, refdir, outfile'); end
  
if exist(outfile,'file')
  load(outfile,'X');
  return
end
  
X = load_genome_info(refdir);

fprintf('Loading mutations.\n');

X.mut = load_struct(infile);
demand_fields(X.mut,{'patient','chr','pos','ref','alt'});

% patients + tumor types
[X.pat.name ui X.mut.pat_idx] = unique(X.mut.patient);
if slength(X.pat)>65535, error('To handle >65535 patients, need to change uint16 to uint32'); end
X.mut = rmfield(X.mut,'patient');
if isfield(X.mut,'ttype')
  X.pat.ttype = X.mut.ttype(ui);
  X.mut = rmfield(X.mut,'ttype');
else
  X.pat.ttype = repmat({'---'},slength(X.pat),1);
end
[X.ttype.name ui X.pat.ttype_idx] = unique(X.pat.ttype);
if exist('ttypefile','var') && exist(ttypefile,'file')
  tt = load_struct(ttypefile);
  demand_fields(tt,{'name','longname','red','green','blue'});
  tt = make_numeric(tt,{'red','green','blue'});
  tt.clr = [tt.red tt.green tt.blue];
  X.ttype = mapinto(X.ttype,tt,'name',{'longname','clr'});
  X.ttype.clr(isnan(X.ttype.clr))=0.6;
else
  X.ttype.longname = repmat({'---'},slength(X.ttype),1);
  X.ttype.clr = repmat([0.6 0.6 0.6],slength(X.ttype),1);
end  
X.pat.ttype_longname = nansub(X.ttype.longname,X.pat.ttype_idx);
X.pat.ttype_clr = nansub(X.ttype.clr,X.pat.ttype_idx);

% chromosomes + positions
X.mut.chr = listmap(X.mut.chr,X.chr.name);
bad = isnan(X.mut.chr);
if all(bad), error('Chromosome names in mutation file don''t match chromosome names in genomic reference info.'); end
if any(bad), fprintf('Warning: %d/%d mutations are on nonstandard chromosomes, will be ignored later.\n',sum(bad),slength(X.mut)); end
X.mut = make_numeric(X.mut,'pos');
bad = isnan(X.mut.pos);
if all(bad), error('Failed to parse genomic positions.'); end
if any(bad), fprintf('Warning: %d/%d mutations have invalid position information, will be ignored later.\n',sum(bad),slength(X.mut)); end

% ref + alt
X.mut.ref = listmap(upper(X.mut.ref),{'A','C','G','T'});
X.mut.alt = listmap(upper(X.mut.alt),{'A','C','G','T'});
bad = isnan(X.mut.ref)|isnan(X.mut.alt);
if all(bad), error('Failed to parse ref/alt bases: should be A/C/G/T'); end
if any(bad), fprintf('Warning: %d/%d mutations appear to be non-SNPs, will be ignored later.\n',sum(bad),slength(X.mut)); end

% mutation counts
X.pat.nmut = histc(X.mut.pat_idx,1:slength(X.pat));

% mutation context

fprintf('Annotating mutation context:');

X.mut.f = nan(slength(X.mut),1);          % "from" base = reference base
X.mut.t = nan(slength(X.mut),1);          % "to" base = alt (mutated) base
X.mut.l = nan(slength(X.mut),1);          % base immediately to the left (5')
X.mut.r = nan(slength(X.mut),1);          % base immediately to the right (3')
X.mut.ll = nan(slength(X.mut),1);         % base two positions to the left (5')
X.mut.rr = nan(slength(X.mut),1);         % base two positions to the right (3')

for chr=1:slength(X.chr), fprintf(' %s',X.chr.name{chr});
  len = X.chr.len(chr);
  d = load_genome_region(chr,1,len,refdir);
  mi=find(X.mut.chr==chr & X.mut.pos>=1 & X.mut.pos<=len);  X.mut.f(mi) = d(X.mut.pos(mi));
  mi2=mi(X.mut.pos(mi)>1);     X.mut.l(mi2)  = d(X.mut.pos(mi2)-1);
  mi2=mi(X.mut.pos(mi)<len);   X.mut.r(mi2)  = d(X.mut.pos(mi2)+1);
  mi2=mi(X.mut.pos(mi)>2);     X.mut.ll(mi2) = d(X.mut.pos(mi2)-2);
  mi2=mi(X.mut.pos(mi)<len-1); X.mut.rr(mi2) = d(X.mut.pos(mi2)+2);
end
fprintf('\n');

% compare "ref" to "f" to auto-detect wrong build
bad = (~isnan(X.mut.ref) & ~isnan(X.mut.f) & X.mut.ref~=X.mut.f);
if any(bad)
  fprintf('Warning: %d/%d mutations have the wrong reference base for this build.\n',sum(bad),slength(X.mut));
  if mean(bad)>0.5, error('*** PROBABLE BUILD MISMATCH ***'); end
  X.mut.ref(bad)=nan;
  X.mut.f(bad)=nan;
end

% compute "context65", which is used by some downstream functions like NMF and lego plots
X.mut.context65 = 16*(X.mut.f-1) + 4*(X.mut.l-1) + X.mut.r;

% copy alternate base to "t"
X.mut.t = X.mut.alt;

% for G/A positions, flip information in ll/l/f>t/r/rr fields to the opposite strand, so that we're always analyzing a C/T.
% (but leave "ref" field unchanged so that we can still know what the genomic reference base is)
ga = find(X.mut.f==3 | X.mut.f==1);
X.mut.f(ga) = 5-X.mut.f(ga);
X.mut.t(ga) = 5-X.mut.t(ga);
tmp = 5-X.mut.l(ga); X.mut.l(ga) = 5-X.mut.r(ga); X.mut.r(ga) = tmp;
tmp = 5-X.mut.ll(ga); X.mut.ll(ga) = 5-X.mut.rr(ga); X.mut.rr(ga) = tmp;

fprintf('Saving annotated mutation data object to %s\n',outfile);
save(outfile,'X','-v7.3');










