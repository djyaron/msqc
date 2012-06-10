function res = Eenv(obj)
% energy of interaction of the molecule with the environment

res  = sum(sum( obj.density.*obj.frag.H1Env(:,:,obj.ienv) ) );
res = res + obj.frag.HnucEnv(obj.ienv);