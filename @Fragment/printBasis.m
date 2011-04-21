function printBasis(obj)

stype = {'s', 'p', 'd', 'f', 'g','h'};
for ib=1:obj.nbasis
   s1 = ['basis ',num2str(ib),' ',stype{obj.basisType(ib)+1},...
      num2str(obj.basisSubType(ib)),' on atom ', ...
      num2str(obj.basisAtom(ib))];
   disp(s1);
   prim = obj.basisPrims{ib};
   for ip = 1:obj.basisNprims(ib)
      disp(['   prim coef ',num2str(prim(1,ip)),' exp ',...
         num2str(prim(2,ip))]);
   end
end