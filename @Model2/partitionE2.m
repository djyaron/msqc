function E2p = partitionE2(obj, env, O2, arange)
% Input
%    env:        (integer) use density matrix from this environment
%                   0 for isolated fragment
%    operator:  (nbasis,nbasis,nbasis,nbasis) partition this operator
%                   defaults to obj.H2
%    arange:  {npartitions, 1} cell array of partitioned basis functions
%                defaults to partitions corresponding to atoms
% Output
%    E2p:    (natom,natom,natom,natom) partitioned operator

% We'll first make arrays that hold the basis functions for each atom
% arange{iatom} = list of all basis functions on iatom
%  (because the lists may be different length on each atom, we need to use a cell
%  array)

if (nargin < 2)
   env = 0;
end
if (nargin < 3)
   O2 = obj.H2;
end
if (nargin < 4)
   arange = cell(obj.natom,1);
   for iatom = 1:obj.natom
      arange{iatom} = find(obj.basisAtom == iatom);
   end
end

nb = obj.nbasis;
natom = size(arange,1); % so any type of arange (not just atoms) will work)

p2hf = obj.density2p(env);

E2p = zeros(natom,natom,natom,natom);
for a=1:natom
   for b=1:natom
      for c=1:natom
         for d=1:natom
            E2p(a,b,c,d) = sum(sum(sum(sum( ...
            p2hf(arange{a},arange{b},arange{c},arange{d}) ...
            .*O2(arange{a},arange{b},arange{c},arange{d}) ...
            ))));
         end
      end
   end
end
