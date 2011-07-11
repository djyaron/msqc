%% What is a reasonable value to use for the environment.
mag = [1.0 5.0 10.0 25.0];
nenv = 10;
cubSize = [6,6,6];
for imag=1:size(mag,2)
   for ienv = 1:nenv
      env{imag,ienv} = Environment.newCube(cubSize,mag(imag));
   end
end
% only run this once, or you will overwrite the env. It is commented out
% for this reason.
%save('data3/env1.mat','env');
%% Generate all of the needed quantum data
%load('data3/env1.mat');
config = Fragment.defaultConfig();
config.template = 'fhyde';
config.basisSet = 'sto-3g ';
frag0 = Fragment('c:\dave\apoly\msqc\data3', config);
%% These results suggest a magnitude of 2 to 3 cause significant
% perturbations
figure(100);
i1 = 0;
plot(0,norm(frag.dipole),'ro');
for imag = 1:size(env,1)
   for ienv = 1:size(env,2)
      frag.addEnv(env{imag,ienv});
      i1 = i1+1;
      figure(100);
      hold on;
      plot(imag,norm(frag.dipoleEnv(:,i1)),'ro');
   end
end

%% Will make 200 calculations with magnitude of 2.5
clear all;
nenv = 200;
cubSize = [6,6,6];
for ienv = 1:nenv
   env{ienv,1} = Environment.newCube(cubSize,2.5);
end
% only run this once, or you will overwrite the env. It is commented out
% for this reason.
%save('data3/env_mag25.mat','env');
%%
load('data3/env_mag25.mat');
config = Fragment.defaultConfig();
config.template = 'fhyde';
config.basisSet = '6-31G**';
frag = Fragment('c:\dave\apoly\msqc\data3', config);
%% These results suggest a magnitude of 2 to 3 cause significant
for ienv = 1:size(env,1)
   frag.addEnv(env{ienv,1});
end
%% Scaled sto-3g calcs
config = Fragment.defaultConfig();
config.template = 'fhydeGen';
config.basisSet = 'gen';
config.par = [1.0 1.0 1.0 1.0 1.0];
frag1 = Fragment('c:\dave\apoly\msqc\data3', config);
config.par = [0.8 0.8 0.8 0.8 0.8];
frag2 = Fragment('c:\dave\apoly\msqc\data3', config);
config.par = [1.2 1.2 1.2 1.2 1.2];
frag3 = Fragment('c:\dave\apoly\msqc\data3', config);
