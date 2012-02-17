function res = mcharge(obj, ienv)
% Calculates the mulliken charges on each atom

Q = zeros(size(obj.Z));
GAP = zeros(1,obj.natom);
P = obj.density(ienv).*obj.S;
GOP = sum(P,1);
arange = cell(obj.natom,1);
for iatom = 1:obj.natom
   arange{iatom} = find(obj.basisAtom == iatom);
end
for i = 1:obj.natom
   GAP(i) = sum(GOP(1,arange{i}));
   Q(i) = obj.Z(i)-GAP(i);
end
res = Q;

