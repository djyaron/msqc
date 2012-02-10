%% Load data
clear classes;
root = 'c:\Users\Alex\Programming\msqc\';
% Generate environments for production runs
if (exist('ethane5mp2/env2.mat','file'))
   disp('loading existing environments');
   load('ethane5mp2/env2.mat');
else
    nenv = 100;
    fieldType = [ 1 0 0; 0 1 0; 0 0 1; 2 0 0; 0 2 0; 0 0 2; 1 1 0; 1 0 1; ...
                  0 1 1; 3 0 0; 0 3 0; 0 0 3; 2 1 0; 1 2 0; 2 0 1; 1 0 2; ...
                  0 2 1; 0 1 2; 1 1 1; 4 0 0; 0 4 0; 0 0 4; 3 1 0; 1 3 0; ...
                  3 0 1; 1 0 3; 0 3 1; 0 1 3; 2 2 0; 2 0 2; 0 2 2 ];
    i = 1;
    for ifield = 1:nenv
        new = Environment;
        new.nfield = 1;
        tmp = 0;
        while tmp == 0
            tmp = int16(rand * size(fieldType,1));
        end
        new.fieldType = fieldType(tmp, :);
        new.fieldMag = int16(rand * 650) + 50;
        %disp(new.fieldMag)
        new.ncharge = 0; new.rho = 0; new.r = 0;
        env{ifield} = new;
    end
    save('ethane5mp2/env2.mat','env');
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
HLbasis = {'6-31G' '6-31G*' '6-31G**'};
HL = cell(npar,3);
LL = cell(npar,3);
%%
if (exist('ethane5mp2/ethaneDat.mat','file'))
   disp('loading existing data');
   load('ethane5mp2/ethaneDat.mat');
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
         frag1 = Fragment([root,'ethane5mp2'], config);
         for ienv = 1:nenv
            display(['HL env ',num2str(ienv)]);
            frag1.addEnv(env{ienv});
         end
         HL{ipar,ihl} = frag1;
      end
      % LL 1
      config.basisSet = 'STO-3G';
      frag2 = Fragment([root,'ethane5mp2'], config);
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
      frag3 = Fragment([root,'ethane5mp2'], config);
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
      frag4 = Fragment([root,'ethane5mp2'], config);
      for ienv = 1:nenv
         display(['LL env ',num2str(ienv)]);
         frag4.addEnv(env{ienv});
      end
      LL{ipar,3} = frag4;
   end
   
   % since even loading all the files will take time, we'll dave everything
   save('ethane5mp2/ethaneDat.mat');
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
pt{1} = [0.59959     -1.4638     0.70091    0.041315     0.49848     0.50381      1.1166      0.0      0.0     0.97819     0.63104     0.79428     0.99788    -0.13183     0.12354     0.18347];
f1.updateDensity(pt{1});
err{1} = f1.err(pt{1});
nerr = size(err{1},1)*size(err{1},2);
disp([' starting rms err = ',num2str(sqrt(err{1}*err{1}')/nerr)]);
pt{2} = lsqnonlin(@f1.err, pt{1});
err{2} = f1.err(pt{2});
disp(['KE=0 SP=0 rms err = ',num2str(sqrt(err{2}*err{2}')/nerr)]);
