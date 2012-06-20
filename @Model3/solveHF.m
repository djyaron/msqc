function solveHF(obj,envs,eps)
% envs: list of environments in which to solve HF (0=no environment)
% eps: tolerance for the HF convergence
if (nargin < 2)
   envs = 0:obj.nenv;
end
if (nargin < 3)
   eps = 1.0e-8;
end

nenv = size(envs,2);

if sum(size(obj.densitySave)) == 0
    obj.densitySave = cell(1,obj.nenv+1); 
end


for i = 1:nenv
    ienv = envs(i);
    if (ienv == 0)
        [obj.orb,obj.Eorb,obj.Ehf] = obj.hartreeFock(0,eps);
    else
        [obj.orbEnv(:,:,ienv), obj.EorbEnv(:,ienv), obj.EhfEnv(1,ienv)] = ...
            obj.hartreeFock(ienv,eps);
    end
end


%{
orbEnv = cell(1, nenv);
EorbEnv = cell(1, nenv);
EhfEnv = cell(1, nenv);
densitySave = cell(1, nenv);
H1 = cell(1, nenv);
Enuc = cell(1, nenv);

for i = 1:nenv
    ienv = envs(i);
    if ienv ~= 0
        H1{i} = obj.H1(ienv);
        Enuc{i} = obj.Hnuc(ienv);
    end
end

if (isempty(find(envs == 0, 1)) == false)
    [obj.orb,obj.Eorb,obj.Ehf] = obj.hartreeFock(0,eps);
end

H2 = obj.H2; % Huge. Is the whole thing needed?
S = obj.S;
Nelec = obj.frag.nelec;
densityIn = obj.densitySave; % Potentially huge.
density = obj.frag.density;
H2j = obj.H2j;
H2k = obj.H2k;
parfor i = 1:nenv
    ienv = envs(i);
    if ienv ~= 0
        [orbEnv{i},EorbEnv{i},EhfEnv{i}, ...
            densitySave{i}] = Model3.hartreeFock2(H1{i},H2,S,Enuc{i},Nelec, ...
            densityIn,density,H2j,H2k,ienv,eps);
    end
end

for i = 1:nenv
    ienv = envs(i);
    if ienv ~= 0
        obj.orbEnv(:,:,ienv) = orbEnv{i};
        obj.EorbEnv(:,ienv) = EorbEnv{i};
        obj.EhfEnv(1,ienv) = EhfEnv{i};
        obj.densitySave{ienv+1} = densitySave{i};
    end
end
%}
end


