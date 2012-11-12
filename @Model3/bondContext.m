function [values names] = bondContext(obj,iatom1,jatom1,ienv1)
% Context associated with a particular bond
% input:   iatom1: atom number of first atom
%          jatom1: atom number of second atom
%          ienv:  env number
% output:  x:     double: values of the context
%          names: string: description of the context parameter

if (isempty(obj.bondContextXSaved))
   % initialize delta charges
   for ienv = 1:obj.nenv
      obj.charges(:,ienv+1) = obj.charges(:,ienv+1) - obj.charges(:,1);
   end
   obj.charges(:,1) = zeros(size(obj.charges(:,1)));
   % create and store the contexts
   obj.bondContextNSaved = cell(obj.natom,obj.natom);
   obj.bondContextXSaved = cell(obj.natom,obj.natom,obj.nenv+1);
   x = zeros(3,1);
   for iatom = 1:obj.natom
      bonded = obj.isBonded(iatom,:);
      [~,bondedAtoms] = find(bonded == 1);
      for jatom = bondedAtoms
         obj.atomContextNSaved{iatom,jatom} = {'r','bo','drho'};
         bondLength = norm(obj.rcart(:,iatom) - obj.rcart(:,jatom));
         Zs = sort([obj.Z(iatom),obj.Z(jatom)]);
         if (Zs(1) == 1 && Zs(2) == 1)
            bondLength = bondLength - 0.74;
         elseif (Zs(1) == 1 && Zs(2) == 6)
            bondLength = bondLength - 1.1;
         else
            bondLength = bondLength - 1.5;
         end
         for ienv = 0:obj.nenv
            bondOrder = obj.bondOrders(iatom,jatom,ienv+1) -1;
            drho = obj.charges(iatom,ienv+1) - obj.charges(jatom,ienv+1);
            
            x(1) = bondLength;
            x(2) = bondOrder;
            x(3) = drho;
            obj.bondContextXSaved{iatom,jatom,ienv+1} = x;
         end
      end
   end      
end
values = obj.bondContextXSaved{iatom1,jatom1,ienv1+1};
if (nargout > 1)
   names = obj.bondContextNSaved{iatom1,jatom1};
end
   
end

