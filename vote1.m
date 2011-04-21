%% Generate some environments
size = [6,6,6];
nenv = 10;
env(1,nenv) = Environment;
for i = 1:nenv
   env(1,i) = Environment.newCube(size,10);
end
save('data2/env.mat', 'env');

%% Generate data
datapath = 'data2';
c1 = Fragment.defaultConfig();
c1.template = 'fhydeGen';
c1.basisSet = 'GEN';
c1.par = [0.9 0.9 0.9 1.1 1.1];

frag1 = Fragment(datapath,c1);
%%
for i=1:nenv
   frag1.addEnv(env(i));
end