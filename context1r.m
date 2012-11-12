%function context1(iprocess)
% Playing around with ways to ramp up context
clear classes;
close all;
% choose representative environments
% load('datasets/ch4Dat.mat');
%
% m1 = LL{1,1};
% Ehf = m1.EhfEnv;
% [a,i] = sort(Ehf);
% plot(a,'b.');
% ikeep = i(5:10:90);
% ikeep = sort(ikeep);
% hold on;
% plot(ikeep,Ehf(ikeep),'ro');
% save('ch4keep.mat','ikeep');
%
topDir = 'C:/matdl/yaron/8-12-12/context-par/';
%topDir = '/brashear/yaron/matdl/8-12-12/context-psc/';
filePre = 'ch4-bl-r';

load([topDir,'ch4-bl/all-28.mat'],'f1','currentError','currentPar');

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
commonIn = {'silent',0};
%

diagNames = {'val','rho','avg r','avg bo','shift'};
bondNames = {'val','r','bo','drho'};

% Fix all context sensitive parameters
removeable = cell(length(f1.mixers),1);
for imix = 1:length(f1.mixers)
   mix = f1.mixers{imix};
   if ((mix.mixType == 11) && mix.funcType == 3)
      removeable{imix} = [0 1 1 1 0];
   elseif ((mix.mixType == 11) && mix.funcType == 2)
      removeable{imix} = [0 1 1 1];
   elseif (mix.mixType == 12)
      removeable{imix} = [0 1 1 1];
   elseif (mix.mixType == 4)
      removeable{imix} = [0 1];
   else
      error('unidentified mixtype');
   end
end

disp(['starting error ', num2str(currentError)]);
fprintf(summaryFile,'initial error %12.5f \n',currentError);
ticID = tic;
for iter = 1:25
   allName = [topDir,filePre,'/all-',num2str(iter),'.mat'];
   if (exist(allName,'file'))
      fprintf(1,'LOADING ITERATION %i \n',iter);
      fprintf(summaryFile,'LOADING ITERATION %i \n',iter);
      load(allName,'f1','currentError','currentPar');
   else
      fprintf(1,'STARTING ITERATION %i \n',iter);
      fprintf(summaryFile,'STARTING ITERATION %i \n',iter);
      save('f1temp.mat','f1');
      % set up loop over imix and ipar, so we can do one big parfor loop
      ic = 0;
      mixes = {};
      for imix = 1:length(f1.mixers)
         mix = f1.mixers{imix};
         remove = removeable{imix};
         for ipar = 1:length(mix.par)
            % check if removeable and not yet removed
            if ((remove(ipar) == 1) && (mix.fixed(ipar) == 0))
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
      errors = zeros(nSave,1);
      pars = cell(nSave,1);
      parfor ic = 1:nSave
         imix = mixes{ic}.imix;
         ipar = mixes{ic}.ipar;
         [etemp,ptemp] = contextFit([],imix,ipar,1,50);
         errors(ic) = etemp;
         pars{ic} = ptemp;
      end
      disp('Done with loop over all mixers');
      for ic=1:nSave
         fprintf(summaryFile, ...
            '%s %12.5f \n',mixes{ic}.name,errors(ic));
      end
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
      f1.mixers{imix}.fixed(ipar) = 1;
      f1.mixers{imix}.par(ipar) = 0.0;
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
runTime = toc(ticID)