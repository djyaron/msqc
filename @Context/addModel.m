function addModel(obj,mod,ienv,iatom,jatom)
% creates a unique ordering for the bond values, for a particular
% atom (if only iatom passed) or bond (if both iatom and jatom are
% passed), and adds these to a list for feature extraction
if (nargin < 5) % atom context
   if (obj.isBondContext ~= 0)
      error('Context.addModel: Context not created as atomic context');
   end
   obj.idata = obj.idata+1;
   obj.data(obj.idata,:) = obj.atomContext(mod,ienv,iatom);
else
   if (obj.isBondContext ~= 1)
      error('Context.addModel: Context not created as bond context');
   end
   obj.idata = obj.idata+1;
   obj.data(obj.idata,:) = obj.bondContext(mod,ienv,iatom,jatom);
end
end



