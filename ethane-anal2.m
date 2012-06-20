%% Load data
clear classes;
dataroot = 'C:/dave/apoly/msqc/ethane4mp2';
datasetroot = 'c:/dave/apoly/msqc/datasets2';
% Generate environments for production runs
if (exist([dataroot,'/env2.mat'],'file'))
   disp('loading existing environments');
   load([dataroot,'/env2.mat']);
else
   error('no env found');
   mag = 15.0;
   nenv = 100;
   cubSize = [6,6,6];
   cent = [0.77; 0; 0];
   for ienv = 1:nenv
      temp = Environment.newCube(cubSize,mag);
      temp.displace(cent);
      env{ienv} = temp;
   end
   save([dataroot,'/env2.mat'],'env');
end
nenv = size(env,2);
pars{1} = [1.54 1.12 60];
pars{2} = [1.54 1.12 30];
pars{3} = [1.54 1.12 0];
pars{4} = [1.39 1.12 60];
pars{5} = [1.69 1.12 60];
pars{6} = [1.54 0.97 60];
pars{7} = [1.54 1.27 60];
npar = size(pars,2);
HLbasis = {'6-31G'};% '6-31G*' '6-31G**'};
HL = cell(npar,3);
LL = cell(npar,3);
%%
if (exist([datasetroot,'/ethaneDat.mat'],'file'))
   disp('loading existing data');
   load([datasetroot,'/ethaneDat.mat']);
else
   for ipar = 1:size(pars,2)
      par = pars{ipar};
      disp(['rcc ',num2str(par(1)), ...
         ' rch ',num2str(par(2)), ...
         ' angle ',num2str(par(3))]);
      
      config = Fragment.defaultConfig();
      config.method = 'MP2';
      config.par = par;
      
      % HL
      for ihl = 1:size(HLbasis,2)
         config.template = 'ethane1';
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
      config.template = 'ethane1-gen';
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
      config.template = 'ethane1-gen';
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
   
   % since even loading all the files will take time, we'll dave everything
   save([datasetroot,'/ethaneDat.mat']);
end
%% Single geometry fitmod = Model2(reg,nar,dif);
f1 = Fitme;
for ipar = 1:7
   mod{ipar} = Model2(LL{ipar,1},LL{ipar,2},LL{ipar,3});
   mod{ipar}.sepKE = 0;
   mod{ipar}.sepSP = 0;
   f1.addFrag(mod{ipar},HL{ipar,3});
end
% start with 4 parameter fit
f1.exactDensity = 1;
pt{1} = [0.98928     0.73149     0.60662     0.11493];
f1.updateDensity(pt{1});
err{1} = f1.err(pt{1});
nerr = size(err{1},1)*size(err{1},2);
disp([' starting rms err = ',num2str(sqrt(err{1}*err{1}')/nerr)]);
pt{2} = lsqnonlin(@f1.err, pt{1});
err{2} = f1.err(pt{2});
disp(['KE=0 SP=0 rms err = ',num2str(sqrt(err{2}*err{2}')/nerr)]);
%% 16 parameter fit
f1 = Fitme;
for ipar = 1:7
   mod{ipar} = Model2(LL{ipar,1},LL{ipar,2},LL{ipar,3});
   mod{ipar}.sepKE = 1;
   mod{ipar}.sepSP = 1;
   f1.addFrag(mod{ipar},HL{ipar,3});
end
f1.exactDensity = 1;
%pt{1} = [0 0.75 0.75 0 0.75 0.75 0.75 0.75  0  0  0 0.75 0.75 0 0 0];
%pt{1} = [0.59959     -1.4638     0.70091    0.041315     0.49848     0.50381      1.1166      0.0      0.0     0.97819     0.63104     0.79428     0.99788    -0.13183     0.12354     0.18347];
pt{1} = zeros(1,16);
f1.updateDensity(pt{1});
err{1} = f1.err(pt{1});
nerr = size(err{1},1)*size(err{1},2);
disp([' starting rms err = ',num2str(sqrt(err{1}*err{1}')/nerr)]);
limits = 3 * ones(1,16);
pt{2} = lsqnonlin(@f1.errDiffs, pt{1},-limits,limits);
err{2} = f1.err(pt{2});
disp(['KE=0 SP=0 rms err = ',num2str(sqrt(err{2}*err{2}')/nerr)]);

%%
f1.corrPlot(0);
%% Fit each geometry separately
clear classes;
load('ethane4/ethanedat.mat');
for ipar = 7:7
   disp(['starting fit on geometry ',num2str(ipar)]);
   f1 = Fitme;
   mod{ipar} = Model2(LL{ipar,1},LL{ipar,2},LL{ipar,3});
   mod{ipar}.sepKE = 1;
   mod{ipar}.sepSP = 1;
   mod{ipar}.rhodep = 1;
   mod{ipar}.mixType = 1;
   f1.addFrag(mod{ipar},HL{ipar,1});
   f1.exactDensity = 1;
   nfitpar = mod{ipar}.npar;
   %start = zeros(1,nfitpar);
   %results from fitting ipar=7 to HL{ipar,3}
   %start = [-0.1056   -0.0251   -0.1511    0.2214    0.1249   -2.1994   -0.0109   -0.0758    0.0000    0.0009    0.0153    0.2796   -0.7136   -0.0092 -0.1874    0.9148];
   start = [0.62574     0.47368      0.7247     -1.1403      0.1813     0.40036    -0.45548    -0.49231     -0.4414    -0.44616    -0.44858     0.17348      0.1549    0.045701    -0.07031    0.076608 0.0 0.0 0.0];
   limits = 3 * ones(1,nfitpar);
   pt{ipar} = lsqnonlin(@f1.err, start,-limits,limits);
   err{ipar} = f1.errDiffs(pt{ipar});
   corrPlot(f1,pt{ipar}, 0, 800+ipar);
   figure(810);
   hold on;
   [LL1{ipar}, HL1{ipar}] = corrPlot(f1,pt{ipar}, 0, 810);
end
%%
save('ethane4/fit12.mat','pt','err','LL1','HL1');
%%
load('ethane4/fit12.mat');

%% Plots versus mulliken charges
ipar = 1;
mod = Model2(LL{ipar,1},LL{ipar,2},LL{ipar,3});
hl = HL{ipar,1};
%%
clear classes;
load('ethane4/ethanedat.mat');
ipar = 1;
ll = LL{ipar,1};
hl = HL{ipar,1};
%%
m1 = zeros(ll.nenv,ll.natom);
m2 = m1;
for ienv = 1:ll.nenv
  m1(ienv,:) = ll.mcharge(ienv);
  mdiff(ienv) = m1(ienv,2) - m1(ienv,1);
  llke(ienv) = ll.EKE(ienv);
  hlke(ienv) = hl.EKE(ienv);
  d1(ienv) = norm(ll.dipoleEnv(:,ienv));
  m2(ienv,:) = hl.mcharge(ienv);
  m2diff(ienv) = m2(ienv,2) - m2(ienv,1);
  d2(ienv) = norm(hl.dipoleEnv(:,ienv));
end
