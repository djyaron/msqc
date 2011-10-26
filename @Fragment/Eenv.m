function res = Eenv(obj,ienv)
% energy of interaction of the molecule with the environment

res = sum(sum( obj.density(ienv).*obj.H1Env(:,:,ienv) ) );
res = res + frag.HnucEnv(ienv);
end

