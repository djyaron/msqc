%% Load data
clear classes;
root = 'c:\dave\apoly\msqc\';
rmdir('data\h2play','s');
mkdir('data\h2play');
copyfile('templates\h2.tpl','data\h2play\h2.tpl');
copyfile('templates\h2-gen.tpl','data\h2play\h2-gen.tpl');
load('ethylene1mp2/env2.mat');
nenv = 2;
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
%
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
      frag1 = Fragment([root,'data\h2play'], config);
      for ienv = 1:nenv
         display(['HL env ',num2str(ienv)]);
         frag1.addEnv(env{ienv});
      end
      HL{ipar,ihl} = frag1;
   end
   % LL 1
   config.basisSet = 'STO-3G';
   frag2 = Fragment([root,'data\h2play'], config);
   disp(['ipar ',num2str(ipar),' loading LL 1']);
%   for ienv = 1:nenv
%      display(['LL env ',num2str(ienv)]);
%      frag2.addEnv(env{ienv});
%   end
   LL{ipar,1} = frag2;
   
   % LL 2
   config.template = 'h2-gen';
   config.basisSet = 'GEN';
   config.par = [par 0.9 0.9 0.9 0.9 0.9];
   frag3 = Fragment([root,'data\h2play'], config);
   disp(['ipar ',num2str(ipar),' loading LL 2']);
   %       for ienv = 1:nenv
   %          display(['LL env ',num2str(ienv)]);
   %          frag3.addEnv(env{ienv});
   %       end
   LL{ipar,2} = frag3;
   % LL 3
   config.template = 'h2-gen';
   config.basisSet = 'GEN';
   config.par = [par 1.05 1.05 1.05 1.05 1.05];
   disp(['ipar ',num2str(ipar),' loading LL 3']);
   frag4 = Fragment([root,'data\h2play'], config);
   %       for ienv = 1:nenv
   %          display(['LL env ',num2str(ienv)]);
   %          frag4.addEnv(env{ienv});
   %       end
   LL{ipar,3} = frag4;
end

