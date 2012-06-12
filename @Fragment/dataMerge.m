function [ resHL, resLL ] = dataMerge( datFiles, envs, saveFilename, varargin )
%DATAMERGE Merge .dat files with HL, LL data in different environments.
%    Note: The fragments in HL, LL of the various datFiles must be identical
%          except for the environment data, including the order of the pars.
%          checks will be made where convenient to take this into account.
%          The first datFile in the list will be used for all non-env data.
%    datFiles:        A cell array of paths to the files to merge.
%    envs:            A cell array where each element is an array of envs to
%                     use.
%    saveFilename:    Where to put the combined .dat file.
%                     Pass in 0 to not save.
%
% Sample usage:
% [HL, LL] = Fragment.dataMerge({'datasets\ethaneDat.mat', ...
% 'data\ethane_field_MP2\ethane_field_MP2.mat'}, {20:25, 40:45}, ...
% 'datasets\merged.mat');

assert(length(datFiles) == length(envs), 'datFiles and envs not same size.');

ienv = length(envs{1});
nenv = 0;
for i = 1:length(envs)
    nenv = nenv + length(envs{i});
end

load(datFiles{1}, 'HL', 'LL');
assert(exist('HL', 'var') == 1, 'No HL loaded from .dat file.');
assert(exist('LL', 'var') == 1, 'No LL loaded from .dat file.');

nHLbasis = checkForInput(varargin, 'nHLbasis', size(HL, 2));
resHL = cell(size(HL, 1), nHLbasis);
for i = 1:size(HL, 1)
    for j = 1:nHLbasis
        resHL{i, j} = HL{i, j};
    end
end
resLL = LL;

for i = 1:size(resHL, 1)
    for j = 1:nHLbasis
        resHL{i, j}.dataMergeResetEnv(nenv, envs{1});
    end
end
for i = 1:size(LL, 1)
    for j = 1:size(LL, 2)
        resLL{i, j}.dataMergeResetEnv(nenv, envs{1});
    end
end

for dat = 2:length(datFiles)
    load(datFiles{dat}, 'HL', 'LL');
    assert(exist('HL', 'var') == 1, 'No HL loaded from .dat file.');
    assert(exist('LL', 'var') == 1, 'No LL loaded from .dat file.');
    
    for i = 1:size(resHL, 1)
        for j = 1:nHLbasis
            compare(resHL{i, j}, HL{i, j});
            resHL{i, j}.dataMergeFrag(HL{i, j}, envs{dat}, ienv);
        end
    end
    for i = 1:size(LL, 1)
        for j = 1:size(LL, 2)
            compare(resLL{i, j}, LL{i, j});
            resLL{i, j}.dataMergeFrag(LL{i, j}, envs{dat}, ienv);
        end
    end
    ienv = ienv + length(envs{dat});
end

if saveFilename ~= 0
    HL = resHL;
    LL = resLL;
    save(saveFilename, 'HL', 'LL');
end

end

function compare ( frag1, frag2 )
%COMPARE Check that certain features of 2 Fragments are identical.
%    Check that certain data from 2 input Fragments is identical, such that it
%    makes senes to merge the env data. This is not exhaustive, but should
%    catch any reasonable mistakes. There is no return value, since the
%    checks are handled using asserts.
%
%    frag1: A Fragment object.
%    frag2: A different Fragment object.

assert(strcmpi(frag1.config.template, frag2.config.template), ...
    'Fragments not using same template.');
assert(strcmpi(frag1.config.basisSet, frag2.config.basisSet), ...
    'Fragments not using same basis set.');
assert(strcmpi(frag1.config.method, frag2.config.method), ...
    'Fragments not using same method.');
assert(frag1.config.charge == frag2.config.charge, 'Charge does not match.');
assert(frag1.config.spin == frag2.config.spin, 'Spin does not match.');
assert(sum(size(frag1.config.par) == size(frag2.config.par)) == 2 && ...
    sum(sum(frag1.config.par == frag2.config.par)) == ...
    numel(frag1.config.par), 'par does not match.');

assert(strcmpi(frag1.templateText, frag2.templateText), ...
    'templateText does not match.');
assert(strcmpi(frag1.gaussianFile, frag2.gaussianFile), ...
    'gaussianFile does not match.');

assert(frag1.natom == frag2.natom, 'Discrepancy in natom.');
assert(frag1.nelec == frag2.nelec, 'Discrepancy in nelec.');
assert(length(frag1.Z) == length(frag2.Z) && sum(frag1.Z == frag2.Z) == ...
    length(frag1.Z), 'Discrepancy in Z.');
end
