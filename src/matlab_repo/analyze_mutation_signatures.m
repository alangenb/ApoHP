function X = analyze_mutation_signatures(X,k,randseed,sigfile)
% X = analyze_mutation_signatures(X,k,randseed,sigfile)
%
% --> computes Gordenin A3A/A3B (YTCA/RTCA) metrics
% --> runs NMF and assigns names according to reference signature list
%
% k = number of factors
% randseed = random seed
% sigfile = mat file with list of reference signatures to compare to
  
demand_fields(X,{'mut','pat'});
demand_fields(X.mut,{'pat_idx','alt','context65','f','t','l','r','ll','rr'});
if ~isnumeric(X.mut.pat_idx), error('pat_idx should be numeric'); end
if ~isnumeric(X.mut.alt), error('alt should be numeric'); end
if ~isnumeric(X.mut.context65), error('context65 should be numeric'); end
if ~isnumeric(k), error('k should be numeric'); end
if k<2 || k>1000, error('invalid k'); end

% compute Gordenin A3A/A3B (YTCA/RTCA) metrics
c2gt = (X.mut.f==2 & X.mut.t~=1);
X.pat.nC = histc(X.mut.pat_idx(c2gt),1:slength(X.pat));
tca2gt = (c2gt & X.mut.l==4 & X.mut.r==1);
X.pat.nRTCA = histc(X.mut.pat_idx(tca2gt & (X.mut.ll==1 | X.mut.ll==3)),1:slength(X.pat));
X.pat.nYTCA = histc(X.mut.pat_idx(tca2gt & (X.mut.ll==2 | X.mut.ll==4)),1:slength(X.pat));
damp_ab=100;
X.pat.RTCA_C = X.pat.nRTCA./(damp_ab+X.pat.nC); X.pat.YTCA_C = X.pat.nYTCA./(damp_ab+X.pat.nC);
X.pat.apochar = 0.016 + X.pat.YTCA_C - 2.1*X.pat.RTCA_C;

% make input for NMF
Q=[];
Q.pat = X.pat;
Q.pat.nmut = hist3d_fast_wrapper(X.mut.pat_idx,X.mut.context65,X.mut.alt,1,slength(X.pat),1,65,1,4);
Q.pat.nchan = zeros(slength(Q.pat),192);
Q.chan.name = cell(192,1);
orig_list = generate_categ_context65_names();
base='ACGT';
i=1;
for from=1:4, for to=1:4, for left=1:4, for right=1:4
  if from==to, continue; end
  orig_name = [base(from) ' in ' base(left) '_' base(right)];
  Q.chan.name{i,1} = [base(left) ' (' base(from) '->' base(to) ') ' base(right)];
  cidx = find(strcmp(orig_name,orig_list.name)); if isempty(cidx), error('what?'); end
  Q.pat.nchan(:,i) = Q.pat.nmut(:,cidx,to);
  i=i+1;
end,end,end,end
Q = collapse_nmf_input_192_to_96(Q);
Q.chan.catnames = regexprep(Q.chan.name,'^(.) \((.)->(.)\) (.)$','$2 in $1_$4 ->$3');

% perform NMF
N = perform_nmf(Q,k,randseed);

% delete old results if they exist in the main object
X.pat = rmfield_if_exist(X.pat,{'nmf','nmf_norm','nmf_nmut'});
X = rmfield_if_exist(X,'nmf');
% save results in main object
X.pat = mapinto(X.pat,N.pat,'name',{'nmf','nmf_norm','nmf_nmut'});
X.nmf.chan = N.chan;

% ANALYZE signatures if list of reference signatures provided
if exist('sigfile','var')
  if ~exist(sigfile,'file'), error('Sig file does not exist: %s\n',sigfile); end
  
  % load reference signatures file
  S=load(sigfile,'X');S=S.X;
  demand_fields(S,{'chan','factor'});
  demand_fields(S.chan,{'catnames','nmf'});
  map = listmap(X.nmf.chan.catnames,S.chan.catnames);
  if any(isnan(map)) || length(unique(map))~=96, error('problem with sigfile'); end
  S.chan = reorder_struct(S.chan,map); 

  % compute all cosine distances
  for xi=k:-1:1
    us = X.nmf.chan.nmf_norm(:,xi);
    for si=slength(S.factor):-1:1
      them = S.chan.nmf(:,si); them=them/sum(them);
      S.factor.cos(si,xi) = sum(us.*them)/sqrt(sum(us.^2)*sum(them.^2));
    end
  end

  % find closest reference sequence (by cosine distance)
  [mx ord] = max(S.factor.cos,[],1);
  X.nmf.factor.name = S.factor.name(ord);
  X.nmf.factor.cos = as_column(mx);

  % signature(s) to define the APOBEC+ cohort
  X.nmf.factor.apobec = grepmi('APOBEC',X.nmf.factor.name);
  X.pat.frac_apobec = sum(X.pat.nmf_norm(:,X.nmf.factor.apobec),2);

  % signatures to exclude from APOBEC+ cohort
  X.nmf.factor.msupe = grepmi('MSI|Smoking|(UV|TMZ)|POLE|ESO',X.nmf.factor.name);
  X.pat.msupe_neg = all(X.pat.nmf_norm(:,X.nmf.factor.msupe)<0.1,2);
end


