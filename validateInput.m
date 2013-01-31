function parameters=validateInput(argsIn,validParameters,force)
%VALIDATEINPUT  Validate input for m-files
%   parameters - validateInput(varargin,validParameters,[force])
%   varargin - passed directly from the parent script.
%   validParameters - A cell of strings or cells with valid input arguments
%       validParameters = {{'print','p'},{'size','s'},'name'};
%       Will accept the following as valid input:
%              print, -print, p, -p
%              size, -size, s, -s
%              name, -name
%
%       If the input pararameter is specified as 'yes', 'on', or 'true' then the
%       parameter is set as true. If it is 'no', 'off', or 'false' then it
%       is returned as false. This is for when calling programs with out
%       parenthesis.
%
%   force - Force the output parameters struct to have all validParameters,
%   even if they are not given. All non-specified input will be set to
%   'false'.
%
%   parameters is a structure with each given input argument. In the case
%   that there are multiple options, the output is set to the first
%   'option'. 'size' and 's' will both set the 'parameters.size' field.
%
%   Example (This is intended to be called from within a function)
%    varargin={'p','s',10,'name','john doe'}
%    validParameters={{'size','s'},{'print','p'},'name'}
%    parameters=validateInput(varargin,validParameters)
%
%    varargin={'p','on','s',10,'name','john doe'}
%    validParameters={{'size','s'},{'print','p'},'name'}
%    parameters=validateInput(varargin,validParameters)
%
%    varargin={'p'}
%    validParameters={{'size','s'},{'print','p'},'name'}
%    parameters=validateInput(varargin,validParameters,true)
%
% Author: Jedediah Frey
% Created: Apr 2010
% Copyright 2010

error(nargchk(1, 3, nargin, 'struct'))
if nargin<3
    force=false;
else
    force=logical(force);
end
i=1; % Set loop variable
while i<=numel(argsIn) % Do until the end of
    % Determine if the current input is a valid parameter.
    [validParameter,parmName]=valid(argsIn{i},validParameters);
    % If it is not a valid input, die with errror.
    if ~validParameter
        error('validateInput:UnknownParameter',['Unknown Parameter: ' argsIn{i}]);
    end
    % If the parameter is the 'last' input or the next argument is a valid
    % input.
    if i+1>numel(argsIn)||valid(argsIn{i+1},validParameters)
        % Set the parameter to true (used for 'optional' calls)
        parameters.(parmName)=true;
        i=i+1; % Increment counter by 1
    else
        % Otherwise, use the next 'input' as the parameter's value
        parameters.(parmName)=argsIn{i+1};
        % If the value is logical and true, sit it to true.
        if islogical(parameters.(parmName))&&parameters.(parmName)==true
            parameters.(parmName)=true;
            % If it is 'yes' or 'on', set it to true.
        elseif strcmpi(parameters.(parmName),'yes')||strcmpi(parameters.(parmName),'on')||strcmpi(parameters.(parmName),'true')
            parameters.(parmName)=true;
        elseif strcmpi(parameters.(parmName),'no')||strcmpi(parameters.(parmName),'off')||strcmpi(parameters.(parmName),'false')
            parameters.(parmName)=false;
            % If it is a number (that may have been passed as a string,
            % then convert it to a number
        elseif ischar(parameters.(parmName))&&~isnan(str2double(parameters.(parmName)))
            parameters.(parmName)=str2double(parameters.(parmName));
        end
        i=i+2; % Increment counter by 2
    end
end
if ~force
    return;
end
for j=1:numel(validParameters)
    % Get the parameter name.
    if iscell(validParameters{j})
        name=validParameters{j}{1};
    else
        name=validParameters{j};
    end
    % If the parameter is not set, set it to false.
    if ~isfield(parameters,name)
        parameters.(name)=false;
    end
end
end

function [validParameter,name] = valid(parameter,validParameters)
% By default the current parameter isn't valid.
validParameter=false;
name=''; % Set the parameter name to something, in case nothing is returned.
% For each of the validParameters
for j=1:numel(validParameters)
    % If the parameter is a cell.
    if iscell(validParameters{j})
        % For each option in the cell
        for k=1:numel(validParameters{j})
            % If the parameter directly matches OR the first character is a
            % '-' and the rest of the string match.
            if strcmpi(validParameters{j}{k},parameter)||(strcmpi(parameter(1),'-')&&strcmpi(validParameters{j}{k},parameter(2:end)))
                validParameter=true; % Set as a valid parameter.
                name=validParameters{j}{1}; % Return the first cell as the variable name'
                return;
            end
        end
    else
        % If the validParameter isn't a cell, just directly compare it.
        if strcmpi(validParameters{j},parameter)||(strcmpi(parameter(1),'-')&&strcmpi(validParameters{j},parameter(2:end)))
            name=validParameters{j};
            validParameter=true;
            return;
        end
    end
end
end