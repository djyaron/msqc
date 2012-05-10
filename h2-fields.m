%% Load data
clear classes;
%root = 'c:\Users\Alex\Programming\msqc\';
% Generate environments for production runs
root = 'c:\dave\apoly\msqc\';
if (exist('h2fld/env1.mat','file'))
   disp('loading existing environments');
   load('h2fld/env1.mat');
else
    nenv = 10;
    fieldType = [ 1 0 0; 0 1 0; 0 0 1; 2 0 0; 0 2 0; 0 0 2; 1 1 0; 1 0 1; ...
                  0 1 1; 3 0 0; 0 3 0; 0 0 3; 2 1 0; 1 2 0; 2 0 1; 1 0 2; ...
                  0 2 1; 0 1 2; 1 1 1; 4 0 0; 0 4 0; 0 0 4; 3 1 0; 1 3 0; ...
                  3 0 1; 1 0 3; 0 3 1; 0 1 3; 2 2 0; 2 0 2; 0 2 2 ];
    i = 1;
    for ifield = 1:nenv
        new = Environment;
        new.nfield = 1;
%         tmp = 0;
%         while tmp == 0
%             tmp = int16(rand * size(fieldType,1));
%         end
%         new.fieldType = fieldType(tmp, :);
        new.fieldType = [0 0 1];
        %new.fieldMag = int16(rand * 80) + 20;
        new.fieldMag = 100 * ifield;
        %disp(new.fieldMag)
        new.ncharge = 0; new.rho = 0; new.r = 0;
        env{ifield} = new;
    end
    save('h2fld/env1.mat','env');
end
nenv = size(env,2);
%pars{1} = [0.5];
pars{1} = [0.65];
pars{2} = [0.7];
pars{3} = [0.75];
pars{4} = [0.8];
%pars{6} = [0.85];
%pars{7} = [0.9];
npar = size(pars,2);
HLbasis = {'6-31G' '6-31G**'};
HL = cell(npar,2);
LL = cell(npar,3);
%%
if (exist('h2fld/h2fldDat1.mat','file'))
   disp('loading existing data');
   load('h2fld/h2fldDat1.mat');
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
         frag1 = Fragment([root,'h2fld'], config);
         for ienv = 1:nenv
            display(['HL env ',num2str(ienv)]);
            frag1.addEnv(env{ienv});
         end
         HL{ipar,ihl} = frag1;
      end
      % LL 1
      config.basisSet = 'STO-3G';
      frag2 = Fragment([root,'h2fld'], config);
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
      frag3 = Fragment([root,'h2fld'], config);
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
      frag4 = Fragment([root,'h2fld'], config);
      for ienv = 1:nenv
         display(['LL env ',num2str(ienv)]);
         frag4.addEnv(env{ienv});
      end
      LL{ipar,3} = frag4;
   end
   
   % since even loading all the files will take time, we'll dave everything
   save('h2fld/h2fldDat1.mat');
end
