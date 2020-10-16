function r = look(s,i,P)
%
% look(s,i,P)
%
% shows the contents of all fields in s at record(s) i
%
% Mike Lawrence 2008-04-29
%      improved 2008-09-05
%      added functionality to deal with long fields 2010-03-01

if exist('P','var')
  if islogical(P) || isnumeric(P)
    tmp = P;
    P = [];
    P.suppress_long_fields = tmp;
  end
else
  P = [];
end

P = impose_default_value(P,'suppress_long_fields',true);
P = impose_default_value(P,'truncate_length',500);
P = impose_default_value(P,'pct_sign_prefix',true);   % add a prefix to output lines marking them as comments (for matlab log)  [added 2020-04-09]

fields = fieldnames(s);
nf = length(fields);

fnl = zeros(nf,1);
for f=1:nf
  fnl(f) = length(fields{f});
end
mfnl = max(fnl);

if ~exist('i','var'), i=':'; end

if numel(i)==1 && i==0, error('idx = [0]'); end

if ischar(i)
  if strcmpi(i,'all') || strcmp(i,':')
    i=1:slength(s);
  elseif strcmpi(i,'end')
    i=slength(s);
  else
    i = {i};  % handle below
  end
end

if iscellstr(i)
  for fi=1:length(fields), fld=fields{fi};
    if iscellstr(s.(fld))
      if length(i)==1
        tmp = find(strcmp(i,s.(fld)));
      else
        tmp = listmap(i,s.(fld));
        tmp(isnan(tmp)) = [];
      end
      if ~isempty(tmp)
        i = tmp;
        break
end,end,end,end

%  fis = find(ismember(fields,{'name','gene','key','id'}),1);
%  if isempty(fis), fis = 1; end
%  for ff=1:length(fis), fi=fis(ff);
%    if iscellstr(s.(fields{fi}))
%      tmp = find(ismember(s.(fields{fi}),i));
%      if isempty(tmp)
%        disp(i);
%        error('don''t know how to handle that as an index');
%      else
%        i=tmp;
% end,end,end,end

if islogical(i)
  i = find(i);
end
if ~isnumeric(i), error('last argument should be indices to print'); end

for n=1:length(i)
  idx = i(n);
  if length(i)>1 || idx>2, fprintf('[%d]\n',idx); end
  for f=1:nf
    fprintf(['%' num2str(mfnl+4) 's: '], fields{f});
    fld = getfield(s,fields{f});
    if ndims(fld)>4, fprintf('[too many dimensions to display]');
    elseif ndims(fld)==3
      fprintf('[rows=pages]\n');
      disp(squeeze(fld(idx,:,:))');
    else
      if iscell(fld), data = fld{idx,:};
      else data = fld(idx,:); end
      if P.suppress_long_fields && length(data)>P.truncate_length
        if size(data,1)==1
          data = data(1:P.truncate_length);
        elseif size(data,2)>P.truncate_length
          data = data(:,1:P.truncate_length);
        end
        truncflag = true;
      else
        truncflag = false;
      end
      if ischar(data)     
        if size(data,1)==1
          fprintf('%s', data);
        else
          disp(data)
        end
      elseif isnumeric(data)|islogical(data)
        mx = length(data);
        if P.suppress_long_fields && mx>15
          mx = min(10,mx);
          truncflag=true;
        end
        for x=1:mx

          val = data(x);
          if islogical(val) || val==round(val)
            fmt = '%d';
          elseif abs(val)<0.0001 || abs(val)>1e6
            fmt = '%.2d';
          elseif abs(val)>1
            fmt = '%.1f';
          elseif abs(val)>=0.1
            fmt = '%.2f';
          elseif abs(val)>=0.01
            fmt = '%.3f';
          elseif abs(val)>=0.001
            fmt = '%.4f';
          elseif abs(val)>=0.0001
            fmt = '%.5f';
          else
            fmt = '%f';
          end
          fprintf([fmt ' '],val);

%          if mod(data(x),1)==0
%            fprintf('%d ', data(x));
%          else
%            fprintf('%f ', data(x));
%          end

        end  
      elseif iscell(data)
        mx = length(data);
        if P.suppress_long_fields && mx>10
          mx = min(10,mx);
          truncflag=true;
        end
        for x=1:mx
          d2 = data{x};
          if ischar(d2)
            fprintf('%s ', d2);
          elseif isnumeric(d2)|islogical(data)
            if mod(d2,1)==0
              fprintf('%d ', d2);
            else
              fprintf('%f ', d2);
            end
          elseif iscell(d2)
            d2
    end,end,end,end
    if truncflag, fprintf(' ...'); end
    fprintf('\n');
  end
  fprintf('\n');
end
