%% Load data
clear classes;
reset(RandStream.getDefaultStream,sum(100*clock))

%root = 'c:\matdl\ethylener';
dataroot = 'c:\matdl\data\ethylener2';
if (~exist(dataroot,'dir'))
   mkdir(dataroot);
   copyfile('templates/ethylener.tpl',[dataroot,'/ethylener.tpl']);
   copyfile('templates/ethylener-gen.tpl',[dataroot,'/ethylener-gen.tpl']);
   copyfile('ethane4mp2/env2.mat',[dataroot,'/env2.mat']);
end

load([dataroot,'/env2.mat']);
nenv = 50;

% structure will hold the characteristics of each of the parameters
randPars = {};
% first one is C-C bond
t1.type = 'CC bond';
t1.low = 1.33 - 0.15;
t1.high = 1.33 + 0.15;
randPars{end+1} = t1;
% next 4 are c-H bonds
t1.type = 'CH bond';
t1.low = 1.1 - 0.15;
t1.high = 1.1 + 0.15;
for i1 = 1:4
   randPars{end+1} = t1;
end
% next 4 are C-C-H bond angles
t1.type = 'CCH angle';
t1.low = 120 - 6;
t1.high = 120 + 6;
for i1 = 1:4
   randPars{end+1} = t1;
end
% next parameter is constrained dihedral
t1.type = 'constrained dihedral';
t1.low = 180-7;
t1.high = 180;
randPars{end+1} = t1;

% 1 rotation
% t1.type = 'rotation';
% t1.low = -90;
% t1.high = 90;
% randPars{end+1} = t1;
t1.type = 'constrained dihedral';
t1.low = -7;
t1.high = +7;
randPars{end+1} = t1;

% constrained dihedral
t1.low = 180-7;
t1.high = 180;
randPars{end+1} = t1;

pars = cell(0,0);
maxpars = 100;
HLbasis = {'6-31G'};% '6-31G*' '6-31G**'};
HL = cell(0,0);
LL = cell(0,0);
loadResults = 1;
% Find all pars for which a calculation exists
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
   maxpars = length(parsIn)-1; % contrl-C'd job, so last one not done
end
maxpars = 30;
%%
for ipar = 1:maxpars
   %pars{1} = [1.12 1.12 1.12 1.12 109.47 109.47 109.47 120.0 -120.0];
   if (loadResults)
      par = parsIn{ipar};
   else
      par = [];
      for i1 = 1:length(randPars)
         par = [par rr1(randPars{i1}.low,randPars{i1}.high)];
      end
      pars{ipar} = par;
   end
   disp([num2str(ipar),' par = ',num2str(par)]);
   
   config = Fragment.defaultConfig();
   config.method = 'MP2';
   config.par = par;
   
   % HL
   for ihl = 1:size(HLbasis,2)
      config.template = 'ethylener';
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
   
   % LL 2
   config.template = 'ethylener-gen';
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
   %config.template = 'ethane1r-gen';
   %config.basisSet = 'GEN';
   config.par = [par 1.05 1.05 1.05 1.05 1.05];
   disp(['ipar ',num2str(ipar),' loading LL 3']);
   frag4 = Fragment(dataroot, config);
   for ienv = 1:nenv
      display(['LL env ',num2str(ienv)]);
      frag4.addEnv(env{ienv});
   end
   LL{ipar,3} = frag4;
end

% since even loading all the files will take time, we'll dave everything
save([dataroot,'/ethylener2Dat.mat'],'LL','HL');




