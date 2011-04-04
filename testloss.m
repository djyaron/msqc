fraghat_cfg = Fragment.defaultConfig();
fraghat_cfg.template = 'h2';
fraghat_cfg.basisSet = 'STO-3G';
fraghat_cfg.par = 1.1;
frag_cfg = fraghat_cfg;
frag_cfg.basisSet = '6-31G**';

load('data1/env.mat');
experiments = [{env(1), env(2), [1,2], [1,2]};
               {env(3), env(4), [1,2], [1,2]}];

ell = loss(frag_cfg, fraghat_cfg, 'data/h2.tpl', experiments);