function res = bondContext(obj,mod,ienv,iatom,jatom)

bondedi = find(mod.isBonded(iatom,:));
bondedj = find(mod.isBonded(jatom,:));
nbondsi = length(bondedi);
nbondsj = length(bondedj);

tbond = Context.analyzeBond(mod,ienv,iatom,jatom)';

% make atom i have greater than or equal to the number of bonds
% of atom j (for unique ordering.. i.e. CH  not HC)
% if equal nbonds, sort based on atomic charge on the bonded atoms
if ((nbondsi < nbondsj) || (tbond(1) < tbond(3) ))
   tmp = iatom;       iatom = jatom;       jatom = tmp;
   tmp = bondedi;   bondedi = bondedj;   bondedj = tmp;
   tmp = nbondsi;   nbondsi = nbondsj;   nbondsj = tmp;
   tmp = tbond(1); tbond(1) = tbond(3); tbond(3) = tmp;
end

if (nbondsi > 1)
   t1 = zeros(3,nbondsi-1);
   ic = 0;
   for katom = setDiff(bondedi,jatom);
      ic = ic+1;
      t1(:,ic) = ...
         Context.analyzeBond(mod,ienv,iatom,katom);
   end
   [~,is] = sort(t1(3,:));
   ti = t1(:,is);
else
   ti = [];
end

if (nbondsj > 1)
   t1 = zeros(3,nbondsj-1);
   ic = 0;
   for katom = setDiff(bondedj,iatom);
      ic = ic + 1;
      t1(:,ic) = ...
         Context.analyzeBond(mod,ienv,jatom,katom);
   end
   [~,is] = sort(t1(3,:));
   tj = t1(:,is);
else
   tj = [];
end
t2 = [tbond,ti,tj];

res= t2(:)';

end

