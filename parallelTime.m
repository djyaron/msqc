clear classes
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
f1 = makeFitme('ch4',1:17,'envs',1:10, 'enstruct1',en,'kestruct',ke, ...
   'e2struct',e2);
%%
p =ones(1,12);
f1.plot = 0;
%%
f1.parallel = 0;
f1.parHF = zeros(1,12);
tic
e = f1.err(p);
toc

%%
f1.parallel = 1;
f1.parHF = zeros(1,14);
tic
e2 = f1.err(p);
toc
max(abs(e2-e))