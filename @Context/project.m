function res = project(obj,mod,ienv,iatom,jatom)
% projects an atom or bond onto the features
if (nargin < 5) % atom context
   if (obj.isBondContext ~= 0)
      error('Context.addModel: Context not created as atomic context');
   end
   features = obj.atomContext(mod,ienv,iatom);
   f2 = (features-obj.mu)./obj.sigma;
   res = (f2*obj.coeff)';
else
   if (obj.isBondContext ~= 1)
      error('Context.addModel: Context not created as bond context');
   end
   features = obj.bondContext(mod,ienv,iatom,jatom);
   f2 = (features-obj.mu)./obj.sigma;
   res = (f2*obj.coeff)';
end

end