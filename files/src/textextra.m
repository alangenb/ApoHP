function h = textextra(x,y,txt,varargin)

vargs = varargin;
[vargs txtcolor] = extract_from_arglist(vargs,'color');
[vargs txtsize] = extract_from_arglist(vargs,'fontsize');
[vargs txtrot] = extract_from_arglist(vargs,'rotation');

x = double(x); y = double(y);  % causes weird error if they are single

h = text(x,y,txt,vargs{:});

if ~isempty(txtcolor)
  if size(txtcolor,1)==1, txtcolor = repmat(txtcolor,length(h),1); end
  if size(txtcolor,1)~=length(h), error('color has wrong number of rows: expected %d',length(h)); end
  if size(txtcolor,2)~=3, error('color should have three columns'); end
  for i=1:length(h)
    set(h(i),'color',txtcolor(i,:));
  end
end

if ~isempty(txtsize)
  if length(txtsize)==1, txtsize = repmat(txtsize,length(h),1); end
  if length(txtsize)~=length(h), error('fontsize has wrong number of entries: expected %d',length(h)); end
  for i=1:length(h)
    set(h(i),'fontsize',txtsize(i));
  end
end

if ~isempty(txtrot)
  if length(txtrot)==1, txtrot = repmat(txtrot,length(h),1); end
  if length(txtrot)~=length(h), error('fontsize has wrong number of entries: expected %d',length(h)); end
  for i=1:length(h)
    set(h(i),'rotation',txtrot(i));
  end
end





if nargout==0, clear h, end
