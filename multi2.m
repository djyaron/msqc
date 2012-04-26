%% Fitting multiple molecules, using makeFitme
clear classes;
%trainMols  = {'h2',[2 4 6 7]};
%testMols = {'h2',[1 3 5]};

% >>>> START WITH CONSTANTS <<<<
ke.H = Mixer(0,1,'ke.H');
ke.Cs = Mixer(0,1,'ke.C');
ke.Cp = ke.Cs;
ke.HH = Mixer(0,1,'ke.HH');
ke.CsH = Mixer(0,1,'ke.CH');
ke.CpH = ke.CsH;
ke.CsCs = Mixer(0,1,'ke.CC');
ke.CsCp = ke.CsCs;
ke.CpCp = ke.CsCs;

en.H = Mixer(0,1,'en.H');
en.Cs = Mixer(0,1,'en.C');
en.Cp = en.Cs;
en.HH = Mixer(0,1,'en.HH');
en.CsH = Mixer(0,1,'en.CH');
en.CpH = en.CsH;
en.CsCs = Mixer(0,1,'en.CC');
en.CsCp = en.CsCs;
en.CpCp = en.CsCs;

dataDir = 'tmp/c2h4-1/';
ftest = makeFitme('h2',[],'ethylene',[1 3 5],'enstruct',en,'kestruct',ke,'plot',2);
f1 = makeFitme('h2',[],'ethylene',[2 4 6 7],'enstruct',en,'kestruct',ke, ...
   'testFitme',ftest);
disp('Starting fit for H2');
limits = [];
options = optimset('DiffMinChange',1.0e-5,'TolFun',1.0e-4,'TolX',1.0e-3);
if (exist(dataDir,'dir') ~= 7)
   status = mkdir(dataDir);
end
diary([dataDir,'out.diary']);
diary on;
tic
start = f1.getPars;
%%
start = ptsave + 0.02*rand(size(ptsave));
options = optimset('DiffMinChange',1.0e-5,'TolFun',1.0e-3,'TolX',1.0e-2);
[pt,resnorm,residual,exitflag,output,lambda,jacobian] = ...
   lsqnonlin(@f1.err, start,-limits,limits,options);
clockTime = toc
pt
resnorm
f1.printMixers;
diary off;
figure(799); saveas(gcf,[dataDir,'error.fig']);
figure(800); saveas(gcf,[dataDir,'train.fig']);
figure(810); saveas(gcf,[dataDir,'test.fig']);

%% >>>>> ADD CH/BO dependence KE <<<<
dataDir = 'tmp/ch4-2/';
ke.H.mixType = 2;   ke.H.par(2) = 0;  ke.H.fixed(2) = 0;
ke.C.mixType = 2;   ke.C.par(2) = 0;  ke.C.fixed(2) = 0;
ke.HH.mixType = 3;  ke.HH.par(2) = 0; ke.HH.fixed(2) = 0;
ke.CH.mixType = 3;  ke.CH.par(2) = 0; ke.CH.fixed(2) = 0;
ke.CC.mixType = 3;  ke.CC.par(2) = 0; ke.CC.fixed(2) = 0;

ftest = makeFitme('ch4',[1 3 5],'enstruct',en,'kestruct',ke,'plot',2);
f1 = makeFitme('ch4',[2 4 6 7],'enstruct',en,'kestruct',ke, ...
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
dataDir = 'tmp/ch4-3/';
en.H.mixType = 2;   en.H.par(2) = 0;    en.H.fixed(2) = 0;
en.C.mixType = 2;   en.C.par(2) = 0;    en.C.fixed(2) = 0;
en.HH.mixType = 3;  en.HH.par(2) = 0;   en.HH.fixed(2) = 0;
en.CH.mixType = 3;  en.CH.par(2) = 0;   en.CH.fixed(2) = 0;
en.CC.mixType = 3;  en.CC.par(2) = 0;   en.CC.fixed(2) = 0;
% 
ftest = makeFitme('ch4',[1 3 5],'enstruct',en,'kestruct',ke,'plot',2);
f1 = makeFitme('ch4',[2 4 6 7],'enstruct',en,'kestruct',ke, ...
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

