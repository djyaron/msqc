function [values names] = atomContext(obj,iatom1,ienv1)
% Context associated with a particular atom
% input:   iatom1: atom number
%          ienv1:  env number
% output:  values:     double: values of the context
%          names: string: description of the context parameter

if (isempty(obj.atomContextXSaved))
   % initialize delta charges
   for ienv = 1:obj.nenv
      obj.charges(:,ienv+1) = obj.charges(:,ienv+1) - obj.charges(:,1);
   end
   obj.charges(:,1) = zeros(size(obj.charges(:,1)));
   % create and store the contexts
   obj.atomContextNSaved = cell(obj.natom,1);
   obj.atomContextXSaved = cell(obj.natom,obj.nenv+1);
   x = zeros(3,1);
   for iatom = 1:obj.natom
      % items that do not depend on environment
      obj.atomContextNSaved{iatom} = {'rho','avg r','avg bo'};
      
      % bond lengths
      bonded = obj.isBonded(iatom,:);
      [~,bondedAtoms] = find(bonded == 1);
      bondLengths = zeros(length(bondedAtoms),1);
      ic = 0;
      for jatom = bondedAtoms
         ic = ic+1;
         bondLengths(ic) = norm(obj.rcart(:,iatom) - obj.rcart(:,jatom));
         Zs = sort([obj.Z(iatom),obj.Z(jatom)]);
         if (Zs(1) == 1 && Zs(2) == 1)
            bondLengths(ic) = bondLengths(ic) - 0.74;
         elseif (Zs(1) == 1 && Zs(2) == 6)
            bondLengths(ic) = bondLengths(ic) - 1.1;
         else
            bondLengths(ic) = bondLengths(ic) - 1.5;
         end
      end
      avgBondLength = mean(bondLengths);
      
      for ienv = 0:obj.nenv
         bondOrders  = zeros(length(bondedAtoms),1);
         ic = 0;
         for jatom = bondedAtoms
            ic = ic+1;
            bondOrders(ic) = obj.bondOrders(iatom,jatom,ienv+1) - 1;
         end
         avgBondOrder = mean(bondOrders);
         
         x(1) = obj.charges(iatom,ienv+1);
         x(2) = avgBondLength;
         x(3) = avgBondOrder;
         obj.atomContextXSaved{iatom,ienv+1} = x;
      end
   end
end
values = obj.atomContextXSaved{iatom1,ienv1+1};
if (nargout > 1)
   names = obj.atomContextNSaved{iatom1};
end
   
end

