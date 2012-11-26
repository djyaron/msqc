function addModel(obj,model)

imod = obj.modelIndex(model);

for i = 1:length(model.mixers)
   mixer = model.mixers{i};
   if (~isempty(mixer.context))
      imix = obj.mixerIndex(mixer);
      if (mixer.isDiag)
         obj.context{imod,imix} ...
            = getAllAtomContexts(model,mixer.context);
      else
         obj.context{imod,imix} ...
            = getAllBondContexts(model,mixer.context);
      end
      mixer.cset = obj;
   end
end
end


function [res nc] = getAllAtomContexts(model,cstr)

contexts = parseContext(cstr);
nc = length(contexts);
res = zeros(nc,model.natom,model.nenv+1);
for i = 1:nc
   res(i,:,:) = getAtomContext(model,contexts{i});
end
end

function [res nc] = getAllBondContexts(model,cstr)
contexts = parseContext(cstr);
nc = length(contexts);
res = zeros(nc,model.natom,model.natom,model.nenv+1);
for i = 1:nc
   res(i,:,:,:) = getBondContext(model,contexts{i});
end
end

function res = getAtomContext(model,ctypeIn)
ctype = validatestring(ctypeIn,{'r','q','bo'});
res = zeros(model.natom,model.nenv+1);
for iatom = 1:model.natom
   switch ctype
      case 'r'
         bonded = model.isBonded(iatom,:);
         [~,bondedAtoms] = find(bonded == 1);
         bondLengths = zeros(length(bondedAtoms),1);
         ic = 0;
         for jatom = bondedAtoms
            ic = ic+1;
            bondLengths(ic) = ...
               norm(model.rcart(:,iatom) - model.rcart(:,jatom)) - ...
               refBondLength(model.Z(iatom),model.Z(jatom));
         end
         avgBondLength = mean(bondLengths);
         res(iatom,:) = avgBondLength * ones(1,size(res,2));
      case 'q'
         res(iatom,:) = model.charges(iatom,:)-model.charges(iatom,1);
      case 'bo'
         bonded = model.isBonded(iatom,:);
         [~,bondedAtoms] = find(bonded == 1);
         for ienv = 0:model.nenv
            bondOrders  = zeros(length(bondedAtoms),1);
            ic = 0;
            for jatom = bondedAtoms
               ic = ic+1;
               bondOrders(ic) = model.bondOrders(iatom,jatom,ienv+1) - 1;
            end
            avgBondOrder = mean(bondOrders);
            res(iatom,ienv+1) = avgBondOrder;
         end
   end
end
end

function res = getBondContext(model,ctypeIn)
ctype = validatestring(ctypeIn,{'r','q','bo'});
res = zeros(model.natom,model.natom,model.nenv+1);
for iatom = 1:model.natom
   for jatom = 1:model.natom
      switch ctype
         case 'r'
               bondLength = ...
                  norm(model.rcart(:,iatom) - model.rcart(:,jatom)) - ...
                  refBondLength(model.Z(iatom),model.Z(jatom));
            res(iatom,jatom,:) = bondLength * ones(1,size(res,3));
         case 'q'
            qi = model.charges(iatom,:) - model.charges(iatom,1);
            qj = model.charges(jatom,:) - model.charges(jatom,1);
            res(iatom,jatom,:) = abs(qi - qj);
         case 'bo'
            res(iatom,jatom,:) = model.bondOrders(iatom,jatom,:) - 1;
      end
   end
end
end


function res = parseContext(cstr)
res = cell(0,0);
remain = cstr;
while (~isempty(remain))
   [res{end+1}, remain] = strtok(remain, ' ');
end
end

function res = refBondLength(Z1,Z2)
% reference bond lengths
Zs = sort([Z1,Z2]);
if (Zs(1) == 1 && Zs(2) == 1)
   res =  0.74;
elseif (Zs(1) == 1 && Zs(2) == 6)
   res =  1.1;
else
   res =  1.5;
end
end
