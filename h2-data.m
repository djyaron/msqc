%% Load data
clear classes;
root = 'c:\dave\apoly\msqc\';
% Generate environments for production runs
if (exist('ethylene1mp2/env2.mat','file'))
   disp('loading existing environments');
   load('ethylene1mp2/env2.mat');
else
   mag = 15.0;
   nenv = 100;
   cubSize = [6,6,6];
   cent = [0.77; 0; 0];
   for ienv = 1:nenv
      temp = Environment.newCube(cubSize,mag);
      temp.displace(cent);
      env{ienv} = temp;
   end
   save('ethylene1mp2/env2.mat','env');
end
nenv = size(env,2);
pars{1} = [0.5];
pars{2} = [0.65];
pars{3} = [0.7];
pars{4} = [0.75];
pars{5} = [0.8];
pars{6} = [0.85];
pars{7} = [0.9];
npar = size(pars,2);
HLbasis = {'6-31G' '6-31G**'};
HL = cell(npar,2);
LL = cell(npar,3);
%%
if (exist('h2/h2Dat.mat','file'))
   disp('loading existing data');
   load('h2/h2Dat.mat');
else
   for ipar = 1:size(pars,2)
      par = pars{ipar};
      disp(['rhh ',num2str(par(1))]);
      
      config = Fragment.defaultConfig();
      config.method = 'MP2';
      config.par = par;
      
      % HL
      for ihl = 1:size(HLbasis,2)
         config.template = 'h2';
         config.basisSet = HLbasis{ihl};
         disp(['ipar ',num2str(ipar),' loading HL ',num2str(ihl)]);
         frag1 = Fragment([root,'h2'], config);
         for ienv = 1:nenv
            display(['HL env ',num2str(ienv)]);
            frag1.addEnv(env{ienv});
         end
         HL{ipar,ihl} = frag1;
      end
      % LL 1
      config.basisSet = 'STO-3G';
      frag2 = Fragment([root,'h2'], config);
      disp(['ipar ',num2str(ipar),' loading LL 1']);
      for ienv = 1:nenv
         display(['LL env ',num2str(ienv)]);
         frag2.addEnv(env{ienv});
      end
      LL{ipar,1} = frag2;
      
      % LL 2
      config.template = 'h2-gen';
      config.basisSet = 'GEN';
      config.par = [par 0.9 0.9 0.9 0.9 0.9];
      frag3 = Fragment([root,'h2'], config);
      disp(['ipar ',num2str(ipar),' loading LL 2']);
      for ienv = 1:nenv
         display(['LL env ',num2str(ienv)]);
         frag3.addEnv(env{ienv});
      end
      LL{ipar,2} = frag3;
      % LL 3
      config.template = 'h2-gen';
      config.basisSet = 'GEN';
      config.par = [par 1.05 1.05 1.05 1.05 1.05];
      disp(['ipar ',num2str(ipar),' loading LL 3']);
      frag4 = Fragment([root,'h2'], config);
      for ienv = 1:nenv
         display(['LL env ',num2str(ienv)]);
         frag4.addEnv(env{ienv});
      end
      LL{ipar,3} = frag4;
   end
   
   % since even loading all the files will take time, we'll dave everything
   save('h2/h2Dat.mat');
end
