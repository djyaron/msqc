function res = checkForInput(varargin,inputName,default)
idx = find(cellfun(@(x)isequal(lower(x),lower(inputName)), varargin));
if (~isempty(idx))
    res = varargin{idx+1};
else
    res = default;
end
end
