function res = Een(obj,atoms)
% energy of interaction of the molecule with the environment
% input:
%   atoms = list of atoms to include (default is all atoms)
% output:
%   res = vector with length =length(atoms). res(i) contains the total
%         energy of interaction between all electrons and the nucleus of
%         of atoms(i)
if (nargin < 2)
   atoms = 1:obj.natoms;
end
res = zeros(1,length(atoms));
ic = 0;
for iatom = atoms
   ic = ic+1;
   res(ic) = sum(sum( obj.density(obj.ienv).*obj.H1en(iatom,obj.ienv) ));
end
