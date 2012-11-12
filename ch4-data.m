%% Load data
clear classes;
%root = 'c:\dave\apoly\msqc\';
root = 't:\msqc\';
% Generate environments for production runs
if (exist('ch4/env2.mat','file'))
   disp('loading existing environments');
   load('ch4/env2.mat');
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
   save('ch4/env2.mat','env');
end
nenv = size(env,2);
%nenv = 20;

pars{1} = [1.12 1.12 1.12 1.12 109.47 109.47 109.47 120.0 -120.0];
pars{2} = [0.97 0.97 0.97 0.97 109.47 109.47 109.47 120.0 -120.0];
pars{3} = [1.27 1.27 1.27 1.27 109.47 109.47 109.47 120.0 -120.0];

pars{4} = [0.98 1.12 1.12 1.12 109.47 109.47 109.47 120.0 -120.0];
pars{5} = [1.27 1.12 1.12 1.12 109.47 109.47 109.47 120.0 -120.0];
pars{6} = [0.98 0.98 1.12 1.12 109.47 109.47 109.47 120.0 -120.0];
pars{7} = [1.27 1.27 1.12 1.12 109.47 109.47 109.47 120.0 -120.0];
pars{8} = [0.98 1.27 1.12 1.12 109.47 109.47 109.47 120.0 -120.0];

pars{9} = [1.12 1.12 1.12 1.12 104.0 109.47 109.47 120.0 -120.0];
pars{10} = [1.12 1.12 1.12 1.12 115.0 109.47 109.47 120.0 -120.0];
pars{11} = [1.12 1.12 1.12 1.12 109.47 109.47 109.47 110.0 -120.0];

pars{12} = [0.98 0.98 1.12 1.12 109.47 109.47 109.47 120.0 -120.0];
pars{13} = [0.98 0.98 0.98 1.12 109.47 109.47 109.47 120.0 -120.0];
pars{14} = [0.98 0.98 0.98 0.98 109.47 109.47 109.47 120.0 -120.0];
pars{15} = [1.27 1.27 1.12 1.12 109.47 109.47 109.47 120.0 -120.0];
pars{16} = [1.27 1.27 1.27 1.12 109.47 109.47 109.47 120.0 -120.0];
pars{17} = [1.27 1.27 1.27 1.27 109.47 109.47 109.47 120.0 -120.0];
pars{18} = [1.27 0.98 1.27 1.12 109.47 109.47 109.47 120.0 -120.0];
pars{19} = [1.27 0.98 1.27 0.98 109.47 109.47 109.47 120.0 -120.0];

pars{20} = [1.02 1.02 1.02 1.02 109.47 109.47 109.47 120.0 -120.0];
pars{21} = [1.07 1.07 1.07 1.07 109.47 109.47 109.47 120.0 -120.0];
pars{22} = [1.17 1.17 1.17 1.17 109.47 109.47 109.47 120.0 -120.0];
pars{23} = [1.22 1.22 1.22 1.22 109.47 109.47 109.47 120.0 -120.0];



npar = size(pars,2);
HLbasis = {'6-31G'};% '6-31G*' '6-31G**'};
HL = cell(npar,3);
LL = cell(npar,3);
%%
%if (exist([root,'ch4/ch4Dat.mat'],'file'))
%   disp('loading existing data');
%   load([root,'ch4/ch4Dat.mat');
%else
   for ipar = 1:size(pars,2)
      par = pars{ipar};
      disp(['rch ',num2str(par(1))]);
      
      config = Fragment.defaultConfig();
      config.method = 'MP2';
      config.par = par;
      
      % HL
      for ihl = 1:size(HLbasis,2)
         config.template = 'ch4';
         config.basisSet = HLbasis{ihl};
         disp(['ipar ',num2str(ipar),' loading HL ',num2str(ihl)]);
         frag1 = Fragment([root,'ch4'], config);
         for ienv = 1:nenv
            display(['HL env ',num2str(ienv)]);
            frag1.addEnv(env{ienv});
         end
         HL{ipar,ihl} = frag1;
      end
      % LL 1
      config.basisSet = 'STO-3G';
      frag2 = Fragment([root,'ch4'], config);
      disp(['ipar ',num2str(ipar),' loading LL 1']);
      for ienv = 1:nenv
         display(['LL env ',num2str(ienv)]);
         frag2.addEnv(env{ienv});
      end
      LL{ipar,1} = frag2;
      
      % LL 2
      config.template = 'ch4-gen';
      config.basisSet = 'GEN';
      config.par = [par 0.9 0.9 0.9 0.9 0.9];
      frag3 = Fragment([root,'ch4'], config);
      disp(['ipar ',num2str(ipar),' loading LL 2']);
      for ienv = 1:nenv
         display(['LL env ',num2str(ienv)]);
         frag3.addEnv(env{ienv});
      end
      LL{ipar,2} = frag3;
      % LL 3
      config.template = 'ch4-gen';
      config.basisSet = 'GEN';
      config.par = [par 1.05 1.05 1.05 1.05 1.05];
      disp(['ipar ',num2str(ipar),' loading LL 3']);
      frag4 = Fragment([root,'ch4'], config);
      for ienv = 1:nenv
         display(['LL env ',num2str(ienv)]);
         frag4.addEnv(env{ienv});
      end
      LL{ipar,3} = frag4;
   end
   
   % since even loading all the files will take time, we'll dave everything
   save([root,'ch4/ch4Dat.mat']);
%end
