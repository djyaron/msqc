%% Parallel version
% 409.0881
iP = 1;
ftype = 2;
ke.H = Mixer(iP,1,'ke.H',ftype);
ke.Cs = Mixer(iP,1,'ke.C',ftype);
ke.Cp = ke.Cs;
ke.HH = Mixer(iP,1,'ke.HH',ftype);
ke.CsH = Mixer(iP,1,'ke.CH',ftype);
ke.CpH = ke.CsH;
ke.CsCs = Mixer(iP,1,'ke.CC',ftype);
ke.CsCp = ke.CsCs;
ke.CpCp = ke.CsCs;

en.H = Mixer(iP,1,'en.H',ftype);
en.Cs = Mixer(iP,1,'en.C',ftype);
en.Cp = en.Cs;
en.HH = Mixer(iP,1,'en.HH',ftype);
en.CsH = Mixer(iP,1,'en.CH',ftype);
en.CpH = en.CsH;
en.HCs = Mixer(iP,1,'en.HC',ftype);
en.HCp = en.HCs;
en.CsCs = Mixer(iP,1,'en.CC',ftype);
en.CsCp = en.CsCs;
en.CpCp = en.CsCs;

e2.H = Mixer(iP,1,'e2.H',ftype);

e2.C = Mixer(iP,1,'e2.C',ftype);
e2.HH = Mixer(iP,1,'e2.HH',ftype);
e2.CC = Mixer(iP,1,'e2.CC',ftype);
e2.CH = Mixer(iP,1,'e2.CH',ftype);

f2 = makeFitme('ch4',1:16,'ethane',1:7,'ethylene',1:7,'envs',0:25,...
    'enstruct1',en,'kestruct',ke,'e2struct',e2,'verify',0,...
    'datadir','../msqc/datasets');
f2.scratchDir = 'I:/';
p = f2.getPars();
f2.plot = 0;

f2.parHF = zeros(1,14);
tic
err2 = f2.err(p);
toc
