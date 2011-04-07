tplpath = 'data/h2.tpl';
A = [1,2];
B = [1,2];

frag_cfg = Fragment.defaultConfig();
frag_cfg.template = 'h2';
frag_cfg.par = 1.1;
frag_cfg.basisSet = 'CEP-121G';

fraghat_cfgs{1} = frag_cfg;
fraghat_cfgs{1}.basisSet = 'STO-3G';
fraghat_cfgs{2} = frag_cfg;
fraghat_cfgs{2}.basisSet = '3-21G';
fraghat_cfgs{3} = frag_cfg;
fraghat_cfgs{3}.basisSet = '6-21G';

% Generate environments, pair up, and shuffle

n = 8;
env(1,2*n) = Environment;
for i = 1:2*n
  env(1,i) = Environment.newCube([3,3,3],1);
end
experiments = cell(n, 1);
for i = 1:n
  experiments{i} = {env(2*i - 1), env(2*i), A, B};
end
experiments = experiments(randperm(n))

% Instantiate and train

[frag, fragdir] = tmpfrag(frag_cfg, tplpath); % FIXME: memoization
aggfrag = AggregateFragment.train(frag_cfg, fraghat_cfgs, tplpath, experiments(1:(n/2)));

% Test

ell = loss(frag, aggfrag, experiments((n/2)+1:n));