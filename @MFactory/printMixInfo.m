function printMixInfo(obj)
for i=1:length(obj.mixInfo)
   mi = obj.mixInfo{i};
   s1 = [mi.type,' ',num2str(mi.iatom)];
   if (isfield(mi,'jatom'))
      s1 = [s1,' ',num2str(mi.jatom)];
   end
   switch mi.type
      case 'E2slater'
         disp(s1);
         disp(['    F0 --> mix ',mi.mixerF0.toString]);
         disp(['    G1 --> mix ',mi.mixerG1.toString]);
         disp(['    F2 --> mix ',mi.mixerF2.toString]);
      otherwise
         s1 = [s1,' --> mix',mi.mixer.toString];
         disp(s1);
   end
end

