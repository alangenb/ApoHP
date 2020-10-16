function c = green_orange_colormap(nrows)

if ~exist('nrows','var'), nrows = 11; end

r1=46;  r2=255;
g1=196; g2=159;
b1=182; b2=28;

c = [(r1:((r2-r1)/(nrows-1)):r2)' (g1:((g2-g1)/(nrows-1)):g2)' (b1:((b2-b1)/(nrows-1)):b2)']/255;


