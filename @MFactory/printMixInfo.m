function printMixInfo(obj)
for i=1:length(obj.mixInfo)
   mi = obj.mixInfo{i};
   s1 = [mi.type,' ',num2str(mi.iatom)];
   if (isfield(mi,'jatom'))
      s1 = [s1,' ',num2str(mi.jatom)];
   end
   s1 = [s1,' mix ',mi.mixer.toString];
   disp(s1);
end

