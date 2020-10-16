function display_nmf_legos(X,P)

if ~exist('P','var'), P=[]; end
P=impose_default_value(P,'display','counts');

q={}; for i=1:size(X.chan.nmf,2),q{i} = X.chan.nmf(:,i);end

P.catnames=X.chan.catnames;

f=grep('^(factor|sig)$',fieldnames(X));
if ~isempty(f) && isfield(X.(f{1}),'name') && length(X.(f{1}).name)==length(q)
  titles = X.(f{1}).name;
else
  titles=num2cellstr(1:length(q));
end

if ~isfield(P,'fontsize')
  t = sum(cellfun('length',titles));  
  P.fontsize = 20 - min(14,t/60);
end

legos(q,titles,P);ff

