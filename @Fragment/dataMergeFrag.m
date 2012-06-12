function dataMergeFrag( obj, frag, envs, ienv )
%DATAMERGEFRAG Add env data to the merged Fragment from another Fragment.
%   Detailed explanation goes here
%     obj:  The input Fragment, from the first dataSet.
%     nenv: The total number of envs in the merged dataSet.
%     envs: Array of the envs in frag to keep.
%     ienv: The index to begin adding env data.

len = length(envs);
targetIndices = ienv + 1:ienv + len;

env = frag.env(envs);
obj.env(targetIndices) = env;

H1Env = frag.H1Env(:, :, envs);
obj.H1Env(:, :, targetIndices) = H1Env;

HnucEnv = frag.HnucEnv(envs);
obj.HnucEnv(targetIndices) = HnucEnv;

EhfEnv = frag.EhfEnv(envs);
obj.EhfEnv(targetIndices) = EhfEnv;

MP2Env = frag.MP2Env(envs);
obj.MP2Env(targetIndices) = MP2Env;

EorbEnv = frag.EorbEnv(:, envs);
obj.EorbEnv(:, targetIndices) = EorbEnv;

orbEnv = frag.orbEnv(:, :, envs);
obj.orbEnv(:, :, targetIndices) = orbEnv;

dipoleEnv = frag.dipoleEnv(:, envs);
obj.dipoleEnv(:, targetIndices) = dipoleEnv;

end
