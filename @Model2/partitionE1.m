function E1p = partitionE1(obj, env, operator, arange)
% Input
%    env:        (integer) use density matrix from this environment
%    operator:  (nbasis,nbasis) partition this operator
%                defaults to obj.H1(:,:,env)
%    arangeIn:  {npartitions, 1} cell array of partitioned basis functions
%                defaults to partitions corresponding to atoms
% Output
%    E1p:    (natom,natom) partitioned operator

% We'll first make arrays that hold the basis functions for each atom
% arange{iatom} = list of all basis functions on iatom
%  (because the lists may be different length on each atom, we need to use a cell
%  array)

if (nargin < 2)
   env = 0;
end
if (nargin < 3)
   operator = obj.H1(env);
end
if (nargin < 4)
   % this is for atom decomposition by default
   arange = cell(obj.natom,1);
   for iatom = 1:obj.natom
      arange{iatom} = find(obj.basisAtom == iatom);
   end
end
A = obj.density(env);
B = operator;
nb = obj.nbasis;
natom = size(arange,1); % so any type of arange (not just atoms) will work)

E1p = zeros(natom,natom);
% trace (A*B) = sum_{i,j} A(i,j)B(j,i);
% we'll first partition the sum on j:
%    t1(i,iatom) = sum_{j on iatom} A(i,j)B(j,i)
t1 = zeros(nb, natom);
for ib= 1:nb
   for iatom=1:natom
      t1(ib,iatom) = A(ib,arange{iatom})*B(arange{iatom},ib);
   end
end
% now we'll sum up the first indices on each atom
for i1 = 1:natom
   for i2 = 1:natom
      E1p(i1,i2) = sum(t1(arange{i1},i2));
   end
end
E1p = 2* E1p;
