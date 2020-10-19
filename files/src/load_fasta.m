function S = load_fasta(filename)
%
% reads a FASTA file from the specified filename
% returns contents in a struct with the following fields:
%
% header (with ">" character removed)
% seq (with whitespace removed)
%
% Mike Lawrence 2008-06-20
%

L = load_lines(filename);

h = grep('^>',L,1);
ns = length(h);
S = [];
S.header = cell(ns,1);
S.seq = cell(ns,1);
for s=1:ns
  S.header{s} = L{h(s)}(2:end);
  firstline = h(s)+1;
  if s<ns
    lastline = h(s+1)-1;
  else
    lastline = length(L);
  end 
  S.seq{s} = cat(2,L{firstline:lastline});
  S.seq{s} = S.seq{s}(S.seq{s}~=char(10) & S.seq{s}~=char(13) & S.seq{s}~=' ' & S.seq{s}~=char(9));
end
