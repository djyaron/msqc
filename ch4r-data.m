%% Load data
clear classes;
reset(RandStream.getDefaultStream,sum(100*clock))

%root = 'c:\dave\apoly\msqc\dataz';
%dataroot = 'c:/dave/apoly/msqc/dataz/ch4r';
root = 'c:\matdl\data';
filename = 'cf4r-orig';
dataroot = [root,'\',filename];

CF4 = 1; 
loadResults = 0;

if (~exist(dataroot,'dir'))
   mkdir(dataroot,'s');
   if (~CF4)
      tplName = 'ch4';
      copyfile('templates/ch4.tpl',[dataroot,'/ch4.tpl']);
      copyfile('templates/ch4-gen.tpl',[dataroot,'/ch4-gen.tpl']);
   else
      tplName = 'cf4';
      copyfile('templates/cf4.tpl',[dataroot,'/cf4.tpl']);
   end
   copyfile('datasets/env2.mat',[dataroot,'/env2.mat']);
end

% Copying the 
load(['datasets/env2.mat']);
envOrig = env;
envs1 = [6     7     8    13    16    24];
envs2 = [5    10    14    17    20    25];
envsJ = [envs1,envs2];
env={envOrig{envsJ} };
nenv = length(env);

if (CF4)
   % r = 1.39 from http://en.wikipedia.org/wiki/Fluoromethane
   r1  = 1.39 - 0.15;
   r2 = 1.39 + 0.15;
else
   r1  = 1.12 - 0.15;
   r2 = 1.12 + 0.15;
end
t1 = 109.4 - 3;
t2 = 109.4 + 3;
p1 = 120 - 3;
p2 = 120 + 3;

pars = cell(0,0);
maxpars =20;
HLbasis = {'6-31G'};% '6-31G*' '6-31G**'};
HL = cell(0,0);
LL = cell(0,0);
% Load data into a *.mat file
if (loadResults)
   lfiles = dir([dataroot,'/*_cfg.mat']);
   parsIn = {};
   for i = 1:length(lfiles)
      % disp(lfiles(i).name);
      load([dataroot,'/',lfiles(i).name]);
      % disp([Cfile.template,' ',Cfile.basisSet]);
      if (strcmpi(Cfile.basisSet,HLbasis{1}))
         parsIn{end+1} = Cfile.par;
      end
   end
   maxpars = length(parsIn);
end
%maxpars = 35;
%%



for ipar = 1:maxpars
   %pars{1} = [1.12 1.12 1.12 1.12 109.47 109.47 109.47 120.0 -120.0];
   if (loadResults)
      par = parsIn{ipar};
   else
      par = [rr1(r1,r2) rr1(r1,r2) rr1(r1,r2) rr1(r1,r2) ...
         rr1(t1,t2) rr1(t1,t2) rr1(t1,t2) ...
         rr1(p1,p2) -rr1(p1,p2)];
   end
   pars{ipar} = par;
   disp([num2str(ipar),' par = ',num2str(par)]);
   
   config = Fragment.defaultConfig();
   config.method = 'MP2';
   config.par = par;
   
   % HL
   for ihl = 1:size(HLbasis,2)
      config.template = tplName;
      config.basisSet = HLbasis{ihl};
      disp(['ipar ',num2str(ipar),' loading HL ',num2str(ihl)]);
      frag1 = Fragment(dataroot, config);
      for ienv = 1:nenv
         display(['HL env ',num2str(ienv)]);
         frag1.addEnv(env{ienv});
      end
      HL{ipar,ihl} = frag1;
   end
   % LL 1
   config.basisSet = 'STO-3G';
   frag2 = Fragment(dataroot, config);
   disp(['ipar ',num2str(ipar),' loading LL 1']);
   for ienv = 1:nenv
      display(['LL env ',num2str(ienv)]);
      frag2.addEnv(env{ienv});
   end
   LL{ipar,1} = frag2;
   
   if (CF4)
      LL{ipar,2} = LL{ipar,1};
      LL{ipar,3} = LL{ipar,1};
   else
      % LL 2
      config.template = 'ch4-gen';
      config.basisSet = 'GEN';
      config.par = [par 0.9 0.9 0.9 0.9 0.9];
      frag3 = Fragment(dataroot, config);
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
      frag4 = Fragment(dataroot, config);
      for ienv = 1:nenv
         display(['LL env ',num2str(ienv)]);
         frag4.addEnv(env{ienv});
      end
      LL{ipar,3} = frag4;
   end
end
% since even loading all the files will take time, we'll dave everything
save([dataroot,'\',filename,'.mat'],'LL','HL');



