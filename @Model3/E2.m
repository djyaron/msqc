function res = E2(obj,envs)
% Kinetic energy in environment list envs
if (nargin < 2)
   envs = 0:obj.nenv;
end

n = length(envs);
res = zeros(1,n);
for i = 1:n
   ienv = envs(i);
   p2hf = obj.density2p(ienv);
   h2temp = obj.H2(ienv);
   res(i) = sum(sum(sum(sum( p2hf.*h2temp ))));
end

end