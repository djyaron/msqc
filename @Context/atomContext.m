function res = atomContext(obj,mod,ienv,iatom)
bonded = find(mod.isBonded(:,iatom));
t1 = zeros(3, length(bonded));
for ibond = 1:length(bonded)
   t1(:,ibond) = ...
      Context.analyzeBond(mod,ienv,iatom,bonded(ibond));
end
% get unique order by sorting by charge on the bonded atom
[~,ic] = sort(t1(3,:));
t2 = t1(:,ic);
res = t2(:)';

end

