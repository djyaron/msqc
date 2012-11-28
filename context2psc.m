function context2psc %(topDir)
%if (nargin < 1)
%   topDir = 'C:/matdl/yaron/11-26-12/psc6/';
   topDir = '/brashear/yaron/matdl/11-26-12/psc18/';
%end   
disp('got to 1');
iprocess = 1;
disp('got to 2');

runParallel = 1;
showPlots = 0;
psc = 0; % true if you do not want to use optimization toolbox

ics = 1;

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
end
disp('got to 3');

filePre = filePrefix{1};
dataDir = [topDir,filePre];
if (exist(dataDir,'dir') ~= 7)
   [status, message] = mkdir(dataDir);
   if (~status)
       error(['failed to create ',dataDir, ' with error ', message]);
   end
end
disp('got to 4');

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
ke.Cs = Mixer(iP3,11,'ke.C',3);
ke.Cp = ke.Cs;
ke.HH = Mixer(iP2,12,'ke.HH',2);
ke.CH = Mixer(iP2,12,'ke.CH',2);
ke.CH.hybrid = 1;
ke.CCs = Mixer(iP2,12,'ke.CCs',2);
ke.CCs.hybrid = 1;
ke.CCp = Mixer(iP2,12,'ke.CCp',2);
ke.CCp.hybrid = 2;

en.H = Mixer(iP3,11,'en.H',3);
en.Cs = Mixer(iP3,11,'en.C',3);
en.Cp = en.Cs;
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
f1.silent = 1;
f1.parallel = 0;
ftest = makeFitme(testIn{:},commonIn{:},'enstructh',en, ...
   'kestructh',ke,'e2struct',e2);
ftest.silent = 1;
ftest.parallel = 0;
%f1 = makeFitme(trainIn{:},commonIn{:},'enmods',0, ...
%   'kestructh',ke);
disp('got to 5');

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
disp('got to 6');

startName = [topDir,filePre,'/start.mat'];
toSave = {'f1','ftest','currentTrainErr','currentPar','currentErr'};
if (exist(startName,'file'))
   fprintf(1,'LOADING START \n');
   fprintf(summaryFile,'LOADING START \n');
   load(startName,toSave{:});
else
   [currentTrainErr,currentPar,currentErr] = contextFit2(f1,ftest,0,0,0,500,1);
   save(startName);
end

str1 = 'initial error %12.5f test %12.5f \n';
fprintf(1,str1,currentTrainErr,currentErr);
fprintf(summaryFile,str1,currentTrainErr,currentErr);
%ticID = tic;
for iter = 1:1
   allName = [topDir,filePre,'/all-',num2str(iter),'.mat'];
   if (exist(allName,'file'))
      fprintf(1,'LOADING ITERATION %i \n',iter);
      fprintf(summaryFile,'LOADING ITERATION %i \n',iter);
      load(allName,toSave{:});
   else
      fprintf(1,'STARTING ITERATION %i \n',iter);
      fprintf(summaryFile,'STARTING ITERATION %i \n',iter);
      if (runParallel)
         save('f1temp.mat','f1','ftest');
      end
      % set up loop over imix and ipar, so we can do one big parfor loop
      ic = 0;
      mixes = {};
      for imix = 1:length(f1.mixers)
         mix = f1.mixers{imix};
         for ipar = 1:length(mix.par)
            if (mix.fixed(ipar) == 1)
               temp2.imix = imix;
               temp2.ipar = ipar;
               if (mix.mixType == 11)
                  names = diagNames;
               elseif (mix.mixType == 12)
                  names = bondNames;
               elseif (mix.mixType == 4)
                  names = {'val','r'};
               end
               temp2.name = [mix.desc,' ',names{ipar}];
               ic = ic+1;
               mixes{ic} = temp2;
            end
         end
      end
      nSave = length(mixes);
      etrain = zeros(nSave,1);
      errors = zeros(nSave,1);
      pars = cell(nSave,1);
      ticInput = tic
      disp('got to 7');

      parfor ic = 1:nSave
         imix = mixes{ic}.imix;
         ipar = mixes{ic}.ipar;
         %name = mixes{ic}.name;
         %      if (runParallel)
         [etemp,ptemp,etest] = contextFit2([],[],imix,ipar,0,100,psc);
         %      else
         %         [etemp, ptemp] = contextFit(f1,imix,ipar);
         %      end
         etrain(ic) = etemp;
         errors(ic) = etest;
         pars{ic} = ptemp;
         %       temp.imix = imix;
         %       temp.ipar = ipar;
         %       temp.name = [desc,' ',names{ipar}];
         %       mixes{ic} = temp;
         %      disp([desc,' context ',names{ipar},' err ', ...
         %         num2str(etemp - currentError)]);
         %       if (~runParallel)
         %          f1.printMixers;
         %          %            input('junk');
         %          f1.mixers{imix}.fixed(ipar) = 1;
         %          f1.mixers{imix}.par(ipar) = 0;
         %          f1.setPars(currentPar);
         %       end
         %      fprintf(summaryFile, ...
         %         '%s %12.5f \n',name,etemp);
         %       fprintf(summaryFile, ' %i %12.5f \n',ic,etemp);
         
      end
      toc(ticInput)
      disp('Done with loop over all mixers');
      for ic=1:nSave
         fprintf(summaryFile, ...
            '%s train %12.5f test %12.5f \n',mixes{ic}.name,etrain(ic),errors(ic));
      end
      [lowErr,lowi] = min(errors);
      [lowErrTrain,lowiTrain] = min(etrain);
      currentError = lowErr;
      currentPar = pars{lowi};
      imix = mixes{lowi}.imix;
      ipar = mixes{lowi}.ipar;
      name = mixes{lowi}.name;
      str2 = '%s Lowest error of %12.5f and %12.5f is from %s \n';
      fprintf(1,str2, 'train', etrain(lowiTrain), errors(lowiTrain),...
         mixes{lowiTrain}.name);
      fprintf(summaryFile, str2, 'train', etrain(lowiTrain), errors(lowiTrain),...
         mixes{lowiTrain}.name);
      fprintf(1,str2, 'test', etrain(lowi), currentError,name);
      fprintf(summaryFile, str2, 'test', etrain(lowi), currentError,name);
      f1.mixers{imix}.fixed(ipar) = 0;
      f1.setPars(pars{lowi});
      f1.printMixers;
      save(allName);
   end
end
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
%runTime = toc(ticID)
end