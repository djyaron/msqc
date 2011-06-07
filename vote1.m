%% Generate some environments
boxSize = [6,6,6];
nenv = 10;
env(1,nenv) = Environment;
for i = 1:nenv
   env(1,i) = Environment.newCube(boxSize,3);
end
save('data2/env.mat', 'env');

%% Generate data
datapath = 'data2';
c1 = Fragment.defaultConfig();
c1.template = 'fhydeGen';
c1.basisSet = 'GEN';
c1.par = [1.0 1.1 1.1 1.2 1.2];

frag1 = Fragment(datapath,c1);
%%
load('data2/env.mat', 'env');
for i=1:nenv
   frag1.addEnv(env(i));
end

%% To recover all the environments currently calculated for a fragment
datapath = 'data2';
c1 = Fragment.defaultConfig();
c1.template = 'fhydeGen';
c1.basisSet = 'GEN';
c1.par = [1.0 1.1 1.1 0.9 0.9];

frag1 = Fragment(datapath,c1);
frag1.loadAllEnv;


