%% What is a reasonable value to use for the environment.
% A run of Gaussian with this input geometry puts the first carbon at
% zero, and the other along the x axis. So we should displace the
% environments by 0.77 along the x axis, to keep things centered.
%       1 C        0.00000000      0.00000000      0.00000000
%       2 C        1.54000000      0.00000000      0.00000000
%       3 H       -0.35666667      1.00880567      0.00000000
%       4 H       -0.35666667     -0.50440284     -0.87365134
%       5 H       -0.35666667     -0.50440284      0.87365134
%       6 H        1.89666667      0.37651010      0.93591080
%       7 H        1.89666667     -0.99877758     -0.14188809
%       8 H        1.89666667      0.62226748     -0.79402271
clear classes;
mag = [1.0 5.0 10.0 25.0];
nenv = 5;
cubSize = [6,6,6];
cent = [0.77; 0; 0];
for imag=1:size(mag,2)
   for ienv = 1:nenv
      temp = Environment.newCube(cubSize,mag(imag));
      temp.displace(cent);
      env{imag,ienv} = temp;
   end
end
% only run this once, or you will overwrite the env. It is commented out
% for this reason.
%save('ethane1/env0.mat','env');
%% what magnitude of charge do we want
clear classes;
load('ethane1/env0.mat');
config = Fragment.defaultConfig();
config.template = 'ethane1';
config.basisSet = '6-31G**';
config.par = [1.54 1.12 60.0];
frag = Fragment('c:\dave\msqc\ethane1', config);
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
%% Generate environments for production runs 
%From the above, it looks like a magnitude of 15 will be good
clear classes;
mag = 15.0;
nenv = 100;
cubSize = [6,6,6];
cent = [0.77; 0; 0];
for ienv = 1:nenv
    temp = Environment.newCube(cubSize,mag);
    temp.displace(cent);
    env{ienv} = temp;
end
% only run this once, or you will overwrite the env. It is commented out
% for this reason.
%save('ethane1/env1.mat','env');
%% Generate data
clear classes;
load('ethane1/env1.mat');
basis{1} = '6-31G**';
basis{2} = 'STO-3G';
for itry=1:10
    try
for angle = 60:-15:0
    for rcc = 1.39:0.15:1.69 % 1.54
        for rch = 0.97:0.15:1.27 %1.12
            for ib = 1:2
                config = Fragment.defaultConfig();
                config.template = 'ethane1';
                config.basisSet = basis{ib};
                config.par = [rcc rch angle];
                frag = Fragment('c:\dave\msqc\ethane1', config);
                   disp(['rcc ',num2str(rcc), ...
                        ' rch ',num2str(rch), ...
                        ' angle ',num2str(angle), ...
                        ' basis ',basis{ib}]);                
                        % perturbations
%                 for ienv = 1:size(env,2)
%                     frag.addEnv(env{ienv});
%                     disp(['rcc ',num2str(rcc), ...
%                         ' rch ',num2str(rch), ...
%                         ' angle ',num2str(angle), ...
%                         ' basis ',basis{ib}, ...
%                         ' ienv ',num2str(ienv)]);
%                 end
            end
        end
    end
end
    catch
    end
end
disp(['done']);

   
