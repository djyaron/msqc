function loc = findText(text, phrase, issueError)
% Returns location of a phrase in text cell array
% Input:
%   text:  (n,1) or (1,n) cell array of words
%   phrase: (n,1) or (1,n) cell array of phrase to find in words
%   issueError: bool  if true, error() is called if loc < 1
%               [defaults to false]
% Output:
%   loc:   text(loc,1) is first word in phrase, if found
%          loc = 0 if phrase is not found
%          loc = -1 if more than one location is found

% approach is adapted from
% http://arstechnica.com/civis/viewtopic.php?f=20&t=296197

if nargin < 3
   issueError = false;
end


ntext = size(text,1) * size(text,2);
nphrase = size(phrase,1) * size(phrase,2);
% find all location of the first word in the phrase
locs = find(ismember(text,phrase{1})==1);
nlocs = size(locs,1)*size(locs,2);
% set to 1 if locs(i) is the entire phrase
isFullPhrase = zeros(nlocs,1);

for iloc = 1:nlocs
   matches = true;
   if ((iloc+nphrase-1) > ntext)
      matches = false;
   end
   iphrase = 1;
   while (matches && (iphrase < nphrase))
      if (strcmp(text(locs(iloc) + iphrase), phrase(1+iphrase)) )
         iphrase = iphrase +1;
      else
         matches = false;
      end
   end
   isFullPhrase(iloc) = matches;
end

res = find(isFullPhrase == 1);
if (size(res,1) == 0)
   loc = 0;
   if (issueError)
      error(['Error in findText, phrase: ',phrase{1:nphrase},' not found']);
   end
%elseif (size(res,1) > 1)
%   loc = -1;
else
   loc = locs(res(1));
end