function res = calcBO(obj, ienv)


P = obj.density(ienv)*obj.S;
% Initialize cell array that holds list of orbitals on each atom
arange = cell(obj.natom,1);
for iatom = 1:obj.natom
   arange{iatom} = find(obj.basisAtom == iatom);
end

%http://ccb.ou.edu/teaching/CHEM5623/lecture%2027.pdf
% This needs to be checked out. the trace is not efficiently implemented
res = zeros(obj.natom,obj.natom);
for iatom = 1:obj.natom
   for jatom = (iatom+1):obj.natom
      %if (obj.isBonded(iatom,jatom))
         P1 = P(arange{iatom},arange{jatom});
         P2 = P(arange{jatom},arange{iatom});
         res(iatom,jatom) = trace(P1*P2);
         res(jatom,iatom) = res(iatom,jatom);
      %end
   end
end

% % from: http://goldbook.iupac.org/BT07005.html 
% PS = obj.density(ienv) .* obj.S;
% % Initialize cell array that holds list of orbitals on each atom
% arange = cell(obj.natom,1);
% for iatom = 1:obj.natom
%    arange{iatom} = find(obj.basisAtom == iatom);
% end
% 
% %http://ccb.ou.edu/teaching/CHEM5623/lecture%2027.pdf
% % This needs to be checked out. the trace is not efficiently implemented
% res = zeros(obj.natom,obj.natom);
% for iatom = 1:obj.natom
%    for jatom = (iatom+1):obj.natom
%       %if (obj.isBonded(iatom,jatom))
%          res(iatom,jatom) = 2*sum(sum(PS(arange{iatom},arange{jatom})));
%          res(jatom,iatom) = res(iatom,jatom);
%       %end
%    end
% end

end