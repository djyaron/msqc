%% Load data
clear classes;
load('ethane1/env1.mat');
nenv = size(env,2);
pars{1} = [1.54 1.12 60];
pars{2} = [1.54 1.12 30];
pars{3} = [1.54 1.12 0];
pars{4} = [1.39 1.12 60];
pars{5} = [1.69 1.12 60];
pars{6} = [1.54 0.97 60];
pars{7} = [1.54 1.27 60];
npar = size(pars,2);
HLbasis = '6-31G**';
HL = cell(npar,1);
LL = cell(npar,3);
for itry = 1:1
    try
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
            frag1 = Fragment('c:\dave\msqc\ethane1', config);
            for ienv = 1:nenv
                display(['HL env ',num2str(ienv)]);
                frag1.addEnv(env{ienv});
            end
            HL{ipar,1} = frag1;
            % LL 1
            config.basisSet = 'STO-3G';
            frag2 = Fragment('c:\dave\msqc\ethane1', config);
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
            frag3 = Fragment('c:\dave\msqc\ethane1-gen', config);
            for ienv = 1:nenv
                display(['LL env ',num2str(ienv)]);
                frag3.addEnv(env{ienv});
            end
            LL{ipar,2} = frag3;
            % LL 3
            config.template = 'ethane1-gen';
            config.basisSet = 'GEN';
            config.par = [par 1.05 1.05 1.05 1.05 1.05];
            frag4 = Fragment('c:\dave\msqc\ethane1-gen', config);
            for ienv = 1:nenv
                display(['LL env ',num2str(ienv)]);
                frag4.addEnv(env{ienv});
            end
            LL{ipar,3} = frag4;
        end
    catch
    end
end
%%
save('c:\dave\msqc\ethane1\data1.mat','HL','LL');
%%
clear classes;
load('c:\dave\msqc\ethane1-gen\data1.mat');
%% Select non-crazy environments
figure(101);
sym = {'ro','bo','go','ko','rx','bx','gx'};
edifft = [];
for ipar = 1:size(HL,1)
    ediff =  HL{ipar}.EhfEnv - LL{ipar,1}.EhfEnv;
    edifft = [edifft,ediff];
    hold on;
    plot(ediff,sym{ipar})
end
mean = edifft;

%%


npar = size(HL,1);
figure(100);
for ipar = 1:npar
    for ienv = 1:HL{ipar}.nenv
        hold on;
        plot(HL{ipar}.EhfEnv(ienv), LL{ipar,1}.EhfEnv(ienv),'b.');
    end
end
