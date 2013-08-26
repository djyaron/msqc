%% H2
clear;
jobname = 'h2';
pars{1} = [0.5];
pars{2} = [0.65];
pars{3} = [0.7];
pars{4} = [0.75];
pars{5} = [0.8];
pars{6} = [0.85];
pars{7} = [0.9];

h2 = Datagen(jobname, pars, 'method', 'MP2', 'envType', 'mixed', ...
    'HLbasis', {'6-31G' '6-31G**'}, 'nenv', 20);

h2.initDir();
h2.makeEnv();

diary ([h2.dataFolder, '_out.log']);
diary on;
[HL, LL] = h2.runData();
diary off;

%% Methane
clear;
jobname = 'ch4';
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

pars{12} = [0.98 0.98 0.98 1.12 109.47 109.47 109.47 120.0 -120.0];
pars{13} = [0.98 0.98 0.98 0.98 109.47 109.47 109.47 120.0 -120.0];
pars{14} = [1.27 1.27 1.27 1.12 109.47 109.47 109.47 120.0 -120.0];
pars{15} = [1.27 0.98 1.27 1.12 109.47 109.47 109.47 120.0 -120.0];
pars{16} = [1.27 0.98 1.27 0.98 109.47 109.47 109.47 120.0 -120.0];

methane = Datagen(jobname, pars, 'method', 'MP2', 'envType', 'mixed', ...
    'HLbasis', {'6-31G' '6-31G*' '6-31G**'}, 'nenv', 20);
disp(methane.nenv);
methane.initDir();
methane.makeEnv();

diary ([methane.dataFolder, '_out.log']);
diary on;
[HL, LL] = methane.runData();
diary off;

%% Ethane
clear;
jobname = 'ethane';
pars{1} = [1.54 1.12 60];
pars{2} = [1.54 1.12 30];
pars{3} = [1.54 1.12 0];
pars{4} = [1.39 1.12 60];
pars{5} = [1.69 1.12 60];
pars{6} = [1.54 0.97 60];
pars{7} = [1.54 1.27 60];

ethane = Datagen(jobname, pars, 'method', 'MP2', 'envType', 'mixed', ...
    'tplName', 'ethane1', 'tplNameGen', 'ethane1-gen', 'nenv', 20);

ethane.initDir();
ethane.makeEnv();

diary ([ethane.dataFolder, '_out.log']);
diary on;
[HL, LL] = ethane.runData();
diary off;

%% Ethylene
clear;
jobname = 'ethylene';
pars{1} = [1.32 1.08 0];
pars{2} = [1.32 1.08 45];
pars{3} = [1.32 1.08 85];
pars{4} = [1.17 1.08 0];
pars{5} = [1.47 1.08 0];
pars{6} = [1.32 0.93 0];
pars{7} = [1.32 1.23 0];

ethylene = Datagen(jobname, pars, 'method', 'MP2', 'envType', 'mixed', ...
    'nenv', 20);

ethylene.initDir();
ethylene.makeEnv();

diary ([ethylene.dataFolder, '_out.log']);
diary on;
[HL, LL] = ethylene.runData();
diary off;
