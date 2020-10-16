function [w h] = nmf_wrapper(varargin)
% wrapper to Brunet implementation of nmf()
% --> handles all-zero rows/cols
   
v = varargin{1};

xidx = (~all(v==0,1));
yidx = (~all(v==0,2));

[w1 h1] = nmf(v(yidx,xidx), varargin{2:end});

xmap = nan(size(v,2),1);
ymap = nan(size(v,1),1);

xmap(xidx)=1:size(h1,2);
ymap(yidx)=1:size(w1,1);

w = nansub(w1,ymap);
h = nansub(h1',xmap)';



