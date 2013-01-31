function res = mcharge(obj, envs)
if (nargin < 2)
   envs = 0:obj.nenv;
end

% Calculates the mulliken charges on each atom
Q = zeros(length(obj.Z),length(envs));
ic = 0;
for ienv=envs
   ic = ic+1;
   GAP = zeros(1,obj.natom);
   P = obj.density(ienv).*obj.S;
   GOP = sum(P,1);
   arange = cell(obj.natom,1);
   for iatom = 1:obj.natom
      arange{iatom} = find(obj.basisAtom == iatom);
   end
   for i = 1:obj.natom
      GAP(i) = sum(GOP(1,arange{i}));
      Q(i,ic) = obj.Z(i)-GAP(i);
   end
end
% for compatibility with frag.mcharge
if (length(envs) == 1)
   res = Q';
else
   res = Q;
end
end