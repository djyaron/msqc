function dataMergeResetEnv( obj, nenv, envs )
%DATAMERGERESETENV Prepare a Fragment to hold the merged data.
% Resize all env data in the fragment, moving the desired data from the first
% dataSet to the beginning of the redefined Fragment.
%     obj:  The input Fragment, from the first dataSet.
%     nenv: The total number of envs in the merged dataSet.
%     envs: Array of the envs in frag to keep.

len = length(envs);

obj.nenv = nenv;

env = obj.env(envs);
obj.env = repmat(Environment, 1, nenv);
obj.env(1:len) = env;

H1Env = obj.H1Env(:, :, envs);
obj.H1Env = zeros(obj.nbasis, obj.nbasis, nenv);
obj.H1Env(:, :, 1:len) = H1Env;

HnucEnv = obj.HnucEnv(envs);
obj.HnucEnv = zeros(1, nenv);
obj.HnucEnv(1:len) = HnucEnv;

EhfEnv = obj.EhfEnv(envs);
obj.EhfEnv = zeros(1, nenv);
obj.EhfEnv(1:len) = EhfEnv;

MP2Env = obj.MP2Env(envs);
obj.MP2Env = zeros(1, nenv);
obj.MP2Env(1:len) = MP2Env;

EorbEnv = obj.EorbEnv(:, envs);
obj.EorbEnv = zeros(obj.nbasis, nenv);
obj.EorbEnv(:, 1:len) = EorbEnv;

orbEnv = obj.orbEnv(:, :, envs);
obj.orbEnv = zeros(obj.nbasis, obj.nbasis, nenv);
obj.orbEnv(:, :, 1:len) = orbEnv;

dipoleEnv = obj.dipoleEnv(:, envs);
obj.dipoleEnv = zeros(size(obj.dipoleEnv, 1), nenv);
obj.dipoleEnv(:, 1:len) = dipoleEnv;

end
