%% Fitting multiple molecules, using makeFitme
clear classes;
% Hydrogen
% >>>> START WITH CONSTANTS <<<<
dataDir = 'tmp/h2-1/';
junk = Mixer(0,1,'junk');
ke.H = Mixer(0,1,'ke.H');
ke.Cs = junk;
ke.Cp = junk;
ke.HH = Mixer(0,1,'ke.HH');
ke.CsH = junk;
ke.CpH = junk;
ke.CsCs = junk;
ke.CsCp = junk;
ke.CpCp = junk;
en = ke;
en.H = Mixer(0,1,'en.H');
en.HH = Mixer(0,1,'en.HH');

ftest = makeFitme('h2',[1 3 5],'enstruct',en,'kestruct',ke,'plot',2);
f1 = makeFitme('h2',[2 4 6 7],'enstruct',en,'kestruct',ke, ...
   'testFitme',ftest);
disp('Starting fit for H2');
limits = [];
options = optimset('DiffMinChange',1.0e-5);
if (exist(dataDir,'dir') ~= 7)
   status = mkdir(dataDir);
end
diary([dataDir,'out.diary']);
diary on;
tic
start = f1.getPars;
[pt,resnorm,residual,exitflag,output,lambda,jacobian] = ...
   lsqnonlin(@f1.err, start,-limits,limits,options);
clockTime = toc
pt
resnorm
f1.printMixers;
diary off;
figure(799); saveas(gcf,[dataDir,'error.fig']);
figure(800); saveas(gcf,[dataDir,'h2-train.fig']);
figure(810); saveas(gcf,[dataDir,'h2-test.fig']);

% >>>>> ADD CH/BO dependence KE <<<<
dataDir = 'tmp/h2-2/';
ke.H.mixType = 2;
ke.H.par(2) = 0;
ke.H.fixed(2) = 0;
ke.HH.mixType = 3;
ke.HH.par(2) = 0;
ke.HH.fixed(2) = 0;

ftest = makeFitme('h2',[1 3 5],'enstruct',en,'kestruct',ke,'plot',2);
f1 = makeFitme('h2',[2 4 6 7],'enstruct',en,'kestruct',ke, ...
   'testFitme',ftest);
if (exist(dataDir,'dir') ~= 7)
   status = mkdir(dataDir);
end
diary([dataDir,'out.diary']);
tic
start = f1.getPars;
[pt,resnorm,residual,exitflag,output,lambda,jacobian] = ...
   lsqnonlin(@f1.err, start,-limits,limits,options);
clockTime = toc
pt
resnorm
f1.printMixers;
diary off;
figure(799); saveas(gcf,[dataDir,'error.fig']);
figure(800); saveas(gcf,[dataDir,'h2-train.fig']);
figure(810); saveas(gcf,[dataDir,'h2-test.fig']);

% >>>> ADD CH/BO DEPENDENCE TO EN <<<<
dataDir = 'tmp/h2-3/';
en.H.mixType = 2;
en.H.par(2) = 0;
en.H.fixed(2) = 0;
en.HH.mixType = 3;
en.HH.par(2) = 0;
en.HH.fixed(2) = 0;
% 
ftest = makeFitme('h2',[1 3 5],'enstruct',en,'kestruct',ke,'plot',2);
f1 = makeFitme('h2',[2 4 6 7],'enstruct',en,'kestruct',ke, ...
   'testFitme',ftest);
if (exist(dataDir,'dir') ~= 7)
   status = mkdir(dataDir);
end
diary([dataDir,'out.diary']);
diary on;
tic
start = f1.getPars;
[pt,resnorm,residual,exitflag,output,lambda,jacobian] = ...
   lsqnonlin(@f1.err, start,-limits,limits,options);
clockTime = toc
pt
resnorm
f1.printMixers;
diary off;
figure(799); saveas(gcf,[dataDir,'error.fig']);
figure(800); saveas(gcf,[dataDir,'h2-train.fig']);
figure(810); saveas(gcf,[dataDir,'h2-test.fig']);

%% start with just constants
ke.H = Mixer(0,1,'ke.H');
ke.C = Mixer(0,1,'ke.C');
ke.Cp = ke.C;
ke.HH = Mixer(0,1,'ke.HH');
ke.CsH = Mixer(0,1,'ke.CH');
ke.CpH = ke.CsH;
ke.CsCs = Mixer(0,1,'ke.CC');
ke.CsCp = ke.CsCs;
ke.CpCp = ke.CsCs;

en.H = Mixer(0,1,'en.H');
en.C = Mixer(0,1,'en.C');
en.Cp = en.C;
en.HH = Mixer(0,1,'en.HH');
en.CsH = Mixer(0,1,'en.CH');
en.CpH = en.CsH;
en.CsCs = Mixer(0,1,'en.CC');
en.CsCp = en.CsCs;
en.CpCp = en.CsCs;

f1 = makeFitme('h2',[],'ch4',1:7,'enstruct',en,'kestruct',ke);
disp('Starting fit for ch4');
start = f1.getPars;
diary('tmp/ch4-1.diary');
diary on;
limits = [];
options = optimset('DiffMinChange',1.0e-5);
tic
pt = lsqnonlin(@f1.err, start,-limits,limits,options);
clockTime = toc
pt
f1.printMixers;
diary off;



