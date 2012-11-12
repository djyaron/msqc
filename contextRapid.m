%function context1(iprocess)
% Playing around with ways to ramp up context
clear classes;
close all;
% choose representative environments
% load('datasets/ch4rDat.mat');
% m1 = LL{1,1};
% Ehf = m1.EhfEnv;
% [a,i] = sort(Ehf);
% plot(a,'b.');
% %ikeep = i(5:10:90); % used for train of methane
% %ikeep = i(8:10:88); % used for test of methane
% ikeep1 = i(4:3:20);
% ikeep1 = sort(ikeep1);
% ikeep2 = i(5:3:21);
% ikeep2 = sort(ikeep2);
% hold on;
% plot(ikeep1,Ehf(ikeep1),'ro');
% plot(ikeep2,Ehf(ikeep2),'go');
% %save('ch4keep.mat','ikeep');
%
topDir = 'C:/matdl/yaron/10-31-12/context-rapid/';
%topDir = '/brashear/yaron/matdl/9-2-12/context-psc-batchqueue/';
ftype = 3;
fitmeParallel = 1;
showPlots = 0;
separateSP = 0;
psc = 0; % does not use optimization toolbox

for iprocess = 11; % [3 8 6];

if (iprocess == 1)
   trainC{1} = {'h2',[2 3 4],'envs',1:5};
   testC{1} = {'h2',[6 7],'envs',10:15};
   filePrefix{1} = 'h2-cross2';
elseif (iprocess == 2)
   %   load('ch4keep.mat');
   ikeep = [ 12    26    37    41    47    48    62    92    94];
   trainC{1} = {'ch4',1:3,'envs',ikeep};
   ikeep2 = [7    20    22    23    30    54    75    81    88];
   testC{1} = {'ch4',4:19,'envs',ikeep2};
   filePrefix{1} = 'ch4-cross1';
elseif (iprocess == 3)
   ikeep = [6     7     8    13    16    24];
   trainC{1} = {'ch4r',1:10,'envs',ikeep};
   ikeep2 = [5    10    14    17    20    25];
   testC{1} = {'ch4r',11:25,'envs',ikeep2};
   filePrefix{1} = 'ch4r-cross1';
elseif (iprocess == 4)
   %   load('ch4keep.mat');
   ikeep = [6     7     8    13    16    24];
   trainC{1} = {'ch4r',1:10,'ethane',[1 3 5 7],'envs',ikeep};
   ikeep2 = [5    10    14    17    20    25];
   testC{1} = {'ch4r',11:20,'ethane',[2 4 6],'envs',ikeep2};
   filePrefix{1} = 'hcfit1';
elseif (iprocess == 5)
   %   load('ch4keep.mat');
   ikeep = [6     7     8    13    16    24];
   trainC{1} = {'ethane',[1 3 5 7],'envs',ikeep};
   ikeep2 = [5    10    14    17    20    25];
   testC{1} = {'ethane',[2 4 6],'envs',ikeep2};
   filePrefix{1} = 'ethane-cross';
elseif (iprocess == 6)
   %   load('ch4keep.mat');
   ikeep = [6     7     8    13    16    24];
   trainC{1} = {'ch4r',1:10,'ethaner',1:10,'envs',ikeep};
   ikeep2 = [5    10    14    17    20    25];
   testC{1} = {'ch4r',11:20,'ethaner',11:20,'envs',ikeep2};
   filePrefix{1} = 'ch4r-ethaner';
elseif (iprocess == 7)
   %   load('ch4keep.mat');
   ikeep = [6     7     8    13    16    24];
   trainC{1} = {'ethaner',1:10,'envs',ikeep};
   ikeep2 = [5    10    14    17    20    25];
   testC{1} = {'ethaner',11:20,'envs',ikeep2};
   filePrefix{1} = 'ethaner';
elseif (iprocess == 8)
   %   load('ch4keep.mat');
   ikeep = [6     7     8    13    16    24];
   trainC{1} = {'ethaner',1:10,'envs',ikeep};
   ikeep2 = [5    10    14    17    20    25];
   testC{1} = {'ethaner',11:20,'envs',ikeep2};
   filePrefix{1} = 'ethaner';
elseif (iprocess == 11)
   %   load('ch4keep.mat');
   trainC{1} = {'ch4r',1:10,'envs',1:10};
   ikeep2 = [5    10    14    17    20    25];
   testC{1} = {'ch4r',11:20,'envs',11:20};
   filePrefix{1} = 'ch4-env1';
end
filePre = filePrefix{1};
dataDir = [topDir,filePre];
if (exist(dataDir,'dir') ~= 7)
   status = mkdir(dataDir);
end
summaryName = [topDir,filePre,'/summary.txt'];
% if (exist(summaryName,'file'))
%    delete(summaryName);
% end
summaryFile = fopen(summaryName,'a');
diaryName = [topDir,filePre,'/cfit.diary'];
% if (exist(diaryName,'file'))
%    delete(diaryName);
% end
diary(diaryName);
diary on;
commonIn = {'silent',0,'plot',0};
%
iC = 1;
trainIn = trainC{iC};
testIn = testC{iC};
filePre = filePrefix{iC};
ke = [];
en = [];
e2 = [];
f1 = [];

% Assumes atom and bond context have 3 elements
iP2 = [1 0 0 0];
iP3 = [1 0 0 0 0];
diagNames = {'val','rho','avg r','avg bo','shift'};
bondNames = {'val','r','bo','drho'};

ke.H = Mixer(iP3,11,'ke.H',3);
if (separateSP)
   ke.Cs = Mixer(iP3,11,'ke.Cs',3);
   ke.Cp = Mixer(iP3,11,'ke.Cp',3);
else
   ke.Cs = Mixer(iP3,11,'ke.C',3);
   ke.Cp = ke.Cs;
end
ke.HH = Mixer(iP2,12,'ke.HH',2);
ke.CH = Mixer(iP2,12,'ke.CH',2);
ke.CH.hybrid = 1;
ke.CCs = Mixer(iP2,12,'ke.CCs',2);
ke.CCs.hybrid = 1;
ke.CCp = Mixer(iP2,12,'ke.CCp',2);
ke.CCp.hybrid = 2;

en.H = Mixer(iP3,11,'en.H',3);
if (separateSP)
   en.Cs = Mixer(iP3,11,'en.Cs',3);
   en.Cp = Mixer(iP3,11,'en.Cp',3);
else
   en.Cs = Mixer(iP3,11,'en.C',3);
   en.Cp = en.Cs;
end
en.HH = Mixer(iP2,12,'en.HH',2);
en.CH = Mixer(iP2,12,'en.CH',2);
en.CH.hybrid = 1;
en.HC = en.CH;
en.CCs = Mixer(iP2,12,'en.CCs',2);
en.CCs.hybrid = 1;
en.CCp = Mixer(iP2,12,'en.CCp',2);
en.CCp.hybrid = 2;

e2.H = Mixer(iP2,11,'e2.H',2);
e2.C = Mixer(iP2,11,'e2.C',2);
if (iprocess == 1) % h2
   e2.HH = Mixer(iP2,12,'e2.HH',2);
else % want only bond length dependence in H-H for methane
   e2.HH = Mixer([1 0],4,'e2.HH',2);
end
e2.CC = Mixer(iP2,12,'e2.CC',2);
e2.CH = Mixer(iP2,12,'e2.CH',2);

% Create fitme object
f1 = makeFitme(trainIn{:},commonIn{:},'enstructh',en, ...
   'kestructh',ke,'e2struct',e2);
f1.parallel = fitmeParallel;
ftest = makeFitme(testIn{:},commonIn{:},'enstructh',en, ...
   'kestructh',ke,'e2struct',e2);
ftest.parallel = fitmeParallel;
%f1 = makeFitme(trainIn{:},commonIn{:},'enmods',0, ...
%   'kestructh',ke);



% Fix all context sensitive parameters
for imix = 1:length(f1.mixers)
   mix = f1.mixers{imix};
   if ((mix.mixType == 11) && mix.funcType == 3)
      mix.fixed = [0 1 1 1 0];
   elseif ((mix.mixType == 11) && mix.funcType == 2)
      mix.fixed = [0 1 1 1];
   elseif (mix.mixType == 12)
      mix.fixed = [0 1 1 1];
   elseif (mix.mixType == 4)
      mix.fixed = [0 1];
   end
end

startName = [topDir,filePre,'/start.mat'];
toSave = {'f1','ftest','currentTrainErr','currentPar','currentErr'};
if (exist(startName,'file'))
   fprintf(1,'LOADING START \n');
   fprintf(summaryFile,'LOADING START \n');
   load(startName,toSave{:});
else
   [currentTrainErr,currentPar,currentErr] = contextFit2(f1,ftest,0,0,0,500,psc);
   save(startName);
end

str1 = 'initial error %12.5f test %12.5f \n';
fprintf(1,str1,currentTrainErr,currentErr);
fprintf(summaryFile,str1,currentTrainErr,currentErr);
ticID = tic;
for iter = 1:5
   allName = [topDir,filePre,'/all-',num2str(iter),'.mat'];
   if (exist(allName,'file'))
      fprintf(1,'LOADING ITERATION %i \n',iter);
      fprintf(summaryFile,'LOADING ITERATION %i \n',iter);
      load(allName,toSave{:});
   else
      fprintf(1,'STARTING ITERATION %i \n',iter);
      fprintf(summaryFile,'STARTING ITERATION %i \n',iter);
      % set up loop over imix and ipar, so we can do one big parfor loop
      ic = 0;
      mixes = {};
      % unfix all parameters
      for imix = 1:length(f1.mixers)
         mix = f1.mixers{imix};
         for ipar = 1:length(mix.par)
            if (mix.fixed(ipar) == 1)
               mix.fixed(ipar) = 0;
            end
         end
      end
      [currentTrainErr,currentPar,currentErr] = contextFit2(f1,ftest,0,0,0,500,psc);
      save(allName);
      
      str2 = 'context error %12.5f test %12.5f \n';
      fprintf(1,str2,currentTrainErr,currentErr);
      fprintf(summaryFile,str2,currentTrainErr,currentErr);

   end
end
runTime = toc(ticID)
diary off;
fclose(summaryFile);
%    %%
%    for i=1:length(errors)
%       mix = f1.mixers{mixes{i}.imix};
%       ipar = mixes{i}.ipar;
%       etemp = errors(i);
%       disp([mix.desc,' context ',num2str(ipar),' err ', ...
%          num2str(etemp)]);
%       disp(['pars ',num2str(pars{i})]);
%    end
end