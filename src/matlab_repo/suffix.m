function x = suffix(x,suff)

if ~iscellstr(x), error('expects cellstr for first argument'); end
if ~ischar(suff), error('expects char for second argument'); end

x = regexprep(x,'^(.*)$',['$1' suff]);
idx=find(strcmp(x,''));x(idx)=repmat({suff},length(idx),1);
