% FIXME: no longer works...memoization.

fraghat_cfg = Fragment.defaultConfig();
fraghat_cfg.template = 'h2';
fraghat_cfg.basisSet = 'STO-3G';
fraghat_cfg.par = 1.1;
frag_cfg = fraghat_cfg;
frag_cfg.basisSet = '6-31G**';

size = [3,3,3];
nenv = 50;
env(1,nenv) = Environment;
for i = 1:nenv
   env(1,i) = Environment.newCube(size,1);
end

load('data1/env.mat');
experiments = [{env(1), env(2), [1,2], [1,2]};
               {env(3), env(4), [1,2], [1,2]}];

ell = loss(frag_cfg, fraghat_cfg, 'data/h2.tpl', experiments);