function diffOverlap(obj)
% Applies differtial overlap approximation
% H2(i,j,k,l) = 0 if i and j are on different atoms
%                 or k and l are on different atoms

for i=1:obj.nbasis
   for j=1:obj.nbasis
      for k=1:obj.nbasis
         for l=1:obj.nbasis
            if (obj.basisAtom(i) ~= obj.basisAtom(j))
               obj.H2(i,j,k,l) = 0.0;
            elseif (obj.basisAtom(k) ~= obj.basisAtom(l))
               obj.H2(i,j,k,l) = 0.0;
            end
         end
      end
   end
end