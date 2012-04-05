%% Create the environment
% While none of this is "randomly" generated, it represents very wide
% range of fields. Next step is narrowing in more on what range of fields
% are useful.
clear classes;
fieldType = [ 1 0 0; 0 1 0; 0 0 1; 2 0 0; 0 2 0; 0 0 2; 1 1 0; 1 0 1; ...
              0 1 1; 3 0 0; 0 3 0; 0 0 3; 2 1 0; 1 2 0; 2 0 1; 1 0 2; ...
              0 2 1; 0 1 2; 1 1 1; 4 0 0; 0 4 0; 0 0 4; 3 1 0; 1 3 0; ...
              3 0 1; 1 0 3; 0 3 1; 0 1 3; 2 2 0; 2 0 2; 0 2 2 ];
i = 1;
for ifield = 1:size( fieldType, 1 ) 
   for imag=50:50:800
      new = Environment;
      new.nfield = 1;
      new.fieldType = fieldType( ifield, : );
      new.fieldMag = imag;
      new.ncharge = 0; new.rho = 0; new.r = 0;
      env{i} = new;
      i = i + 1;
   end
end
% only run this once, or you will overwrite the env. It is commented out
% for this reason.
save('data/ethane/ethane3/env0.mat','env');

%% Same parameters as Charge run.
clear classes;
root = 'C:\Users\Alex\Programming\msqc';
load('data\ethane\ethane3\env0.mat');
nenv = size(env,2);
pars{1} = [1.54 1.12 60];
pars{2} = [1.54 1.12 30];
pars{3} = [1.54 1.12 0];
pars{4} = [1.39 1.12 60];
pars{5} = [1.69 1.12 60];
pars{6} = [1.54 0.97 60];
pars{7} = [1.54 1.27 60];
npar = size(pars,2);
HLbasis = '6-31G';
HL = cell(npar,1);
LL = cell(npar,3);

%%
dataFolder = '\data\ethane\ethane3';
for itry = 1:1
    %try
        for ipar = 1:size(pars,2)
            par = pars{ipar};
            disp(['rcc ',num2str(par(1)), ...
                ' rch ',num2str(par(2)), ...
                ' angle ',num2str(par(3))]);
            
            config = Fragment.defaultConfig();
            config.par = par;
            
            % HL
            config.template = 'ethane1';
            config.basisSet = HLbasis;
            disp('loading HL');
            frag1 = Fragment([root,dataFolder], config);
            for ienv = 1:nenv
                display(['HL env ',num2str(ienv)]);
                frag1.addEnv(env{ienv});
            end
            HL{ipar,1} = frag1;
            % LL 1
            config.basisSet = 'STO-3G';
            frag2 = Fragment([root,dataFolder], config);
            disp('loading HL 1');
            for ienv = 1:nenv
                display(['LL env ',num2str(ienv)]);
                frag2.addEnv(env{ienv});
            end
            LL{ipar,1} = frag2;
            
            % LL 2
            config.template = 'ethane1-gen';
            config.basisSet = 'GEN';
            config.par = [par 0.9 0.9 0.9 0.9 0.9];
            frag3 = Fragment([root,dataFolder], config);
            for ienv = 1:nenv
                display(['LL env ',num2str(ienv)]);
                frag3.addEnv(env{ienv});
            end
            LL{ipar,2} = frag3;
            % LL 3
            config.template = 'ethane1-gen';
            config.basisSet = 'GEN';
            config.par = [par 1.05 1.05 1.05 1.05 1.05];
            frag4 = Fragment([root,dataFolder], config);
            for ienv = 1:nenv
                display(['LL env ',num2str(ienv)]);
                frag4.addEnv(env{ienv});
            end
            LL{ipar,3} = frag4;
        end
    %catch
    %end
end