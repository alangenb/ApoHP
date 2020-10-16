function d = difff(x,n,dim)
% d = difff(x,n,dim)
% size-preserving variant of diff()

if ~exist('n','var'), n = 1; end
if ~exist('dim','var')
  if size(x,2)>1 && size(x,1)==1
    dim=2;
  else
    dim=1;
  end
end

d = diff(x,n,dim);
  
if dim==1
  d = cat(1,nan(n,size(d,2)),d);
elseif dim==2
  d = cat(2,nan(size(d,1),n),d);
elseif dim==3
  fprintf('WARNING: difff with dim>2 is untested, please double-check!\n');
  d = cat(3,nan(size(d,1),size(d,2),n),d);
else
  error('diff is untested with dim>2');
end

if length(size(d))~=length(size(x)) || ~all(size(d)==size(x))
  error('difff failed');
end


