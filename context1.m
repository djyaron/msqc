%% Playing around with ways to ramp up context
clear classes;
close all;
topDir = 'C:\matdl\yaron\8-10-12\context\';
ftype = 3;
runParallel = 0;
showPlots = 1;

ics = 1;

%trainC{1} = {'h2',2:4,'envs',1:5};
%testC{1} = []; %{'h2',4,'envs',1:100};
%filePrefix{1} = 'h2-27';

trainC{1} = {'ch4',1:17,'envs',1:20};
testC{1} = []; %{'h2',4,'envs',1:100};
filePrefix{1} = 'ch4-17';

filePre = filePrefix{1};
dataDir = [topDir,filePre];
if (exist(dataDir,'dir') ~= 7)
   status = mkdir(dataDir);
end
summaryName = [topDir,filePre,'/summary.txt'];
if (exist(summaryName,'file'))
   delete(summaryName);
end
summaryFile = fopen(summaryName,'w');
diaryName = [topDir,filePre,'/cfit.diary'];
if (exist(diaryName,'file'))
   delete(diaryName);
end
diary(diaryName);
diary on;
commonIn = {};
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

% Create fitme object
%       f1 = makeFitme(trainIn{:},commonIn{:},'enstructh',en, ...
%          'kestructh',ke,'e2struct',e2);
f1 = makeFitme(trainIn{:},commonIn{:},'enmods',0, ...
   'kestructh',ke);

% Fix all context sensitive parameters
for imix = 1:length(f1.mixers)
   mix = f1.mixers{imix};
   if (mix.funcType == 3)
      mix.fixed = [0 1 1 1 0];
   else
      mix.fixed = [0 1 1 1];
   end
end

[currentError,currentPar] = contextFit(f1,0,0);
disp(['starting error ', num2str(currentError)]);
fprintf(summaryFile,'initial error %12.5f \n',currentError);
errors = [];
mixes = {};
pars = {};
for iter = 1:8
   fprintf(1,'STARTING INTERATION %i \n',iter);
   fprintf(summaryFile,'STARTING INTERATION %i \n',iter);
   for imix = 1:length(f1.mixers)
      mix = f1.mixers{imix};
      if (mix.mixType == 11)
         names = diagNames;
      else
         names = bondNames;
      end
      for ipar = 1:length(mix.par)
         if (mix.fixed(ipar) == 1)
            [etemp, ptemp] = contextFit(f1,imix,ipar);
            errors(end+1) = etemp;
            pars{end+1} = ptemp;
            temp.imix = imix;
            temp.ipar = ipar;
            temp.name = [mix.desc,' ',names{ipar}];
            mixes{end+1} = temp;
            disp([mix.desc,' context ',names{ipar},' err ', ...
               num2str(etemp - currentError)]);
            f1.printMixers;
%            input('junk');
            f1.mixers{imix}.fixed(ipar) = 1;
            f1.mixers{imix}.par(ipar) = 0;
            f1.setPars(currentPar);
            fprintf(summaryFile, ...
               '  %s context %s %12.5f \n',mix.desc,names{ipar},etemp);
         end
      end
   end
   disp('Done with loop over all mixers');
   [lowErr,lowi] = min(errors);
   currentError = lowErr;
   currentPar = pars{lowi};
   imix = mixes{lowi}.imix;
   ipar = mixes{lowi}.ipar;
   name = mixes{lowi}.name;
   fprintf(1,'Lowest error of %12.5f is from %s \n', ...
      currentError,name);
   fprintf(summaryFile,'Lowest error of %12.5f is from %s \n', ...
      currentError,name);
   f1.mixers{imix}.fixed(ipar) = 0;
   f1.setPars(pars{lowi});
   f1.printMixers;
end
diary off;
%    %%
%    for i=1:length(errors)
%       mix = f1.mixers{mixes{i}.imix};
%       ipar = mixes{i}.ipar;
%       etemp = errors(i);
%       disp([mix.desc,' context ',num2str(ipar),' err ', ...
%          num2str(etemp)]);
%       disp(['pars ',num2str(pars{i})]);
%    end
