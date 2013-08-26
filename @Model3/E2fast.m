function res = E2fast(obj,envs)
% Kinetic energy in environment list envs
if (nargin < 2)
   envs = 0:obj.nenv;
end
n = size(envs,2);
Ehf = zeros(1,n);
for i = 1:n
   ienv = envs(i);
   if (ienv == 0)
      Ehf(i) = obj.Ehf;
   else
      Ehf(i) = obj.EhfEnv(ienv);
   end
end
EelecNuc = obj.Een(1,envs);
for iatom = 2:obj.natom
   EelecNuc = EelecNuc + obj.Een(iatom,envs);
end

res = Ehf - obj.EKE(envs) - EelecNuc - obj.Eenv(envs) - obj.frag.Hnuc;

end