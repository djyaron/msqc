function mexCompile(optArgs)

% Description:
%   Compile all of the MEX C functions together.
%
% Author:
%   Alex Cappiello (acappiel@andrew.cmu.edu)
%
% Date:
%   6-4-13
% Updated:
%   6-26-13
%
% Inputs:
%   optArgs: Cell array of strings representing additional compilation
%            options. See 'doc mex' for info.
%
% Outputs:
%   Compiled binaries.

if (nargin < 1)
    optArgs = {};
end

if ispc
   blaslib = fullfile(matlabroot, ...
      'extern', 'lib', 'win64', 'microsoft', 'libmwblas.lib');
   lapacklib = fullfile(matlabroot, ...
      'extern', 'lib', 'win64', 'microsoft', 'libmwlapack.lib');
else
   blaslib = '-lmwblas';
   lapacklib = '-lmwlapack';
end

compile('twoElecFock.c', optArgs);
compile('mm.c', [optArgs {blaslib}]);
compile('elementWiseCombine.c', optArgs);
end

function compile(filename, optArgs)
mex(filename, optArgs{:});
end