clear classes;
close all;
topDir = 'C:/matdl/yaron/11-29-12/factory2/';
maxIter = 500;

% CREATE MODEL SETS
dataf = {'ch4rDat','ch4rDat-1c','ch4rDat-diponly','ch4rDat-linrho','ethanerDat'};
dsets = cell(1,2);
dname = cell(1,1);
for idata = 1:length(dataf);
   dname{idata} = dataf{idata};
   dfile = ['datasets/',dataf{idata},'.mat'];
   % train data
   ms = MSet;
   ms.addData(dfile, 1:10, 1:2:20 ,1,791);
   dsets{idata,1} = ms;
   % test data
   ms = MSet;
   ms.addData(dfile, 11:20, 2:2:20 ,1,791);
   dsets{idata,2} = ms;
end

% CREATE POLICIES
pname{1} = 'hybrid1';
m1 = MFactory;
m1.addPolicy('o','KE', 'f','scale', 'sp','sonly', 'i',1, 'c','q r bo');
m1.addPolicy('o','EN', 'f','scale', 'sp','sonly', 'i',1, 'c','q r bo');
m1.addPolicy('o','E2', 'f','scale', 'sp','sonly', 'i',1, 'c','q r bo');

m1.addPolicy('o','KE', 'f','scale', 'sp','hybrid', 'i',6, 'j',1, ...
   'c','r bo q');
m1.addPolicy('o','EN', 'f','scale', 'sp','hybrid', 'i',6, 'j',1, ...
   'c','r bo q');
m1.addPolicy('o','E2', 'f','scale', 'sp','hybrid', 'i',6, 'j',1, ...
   'c','r bo q');

m1.addPolicy('o','KE', 'f','scale', 'sp','combine', 'i',6, 'c','q r bo');
m1.addPolicy('o','KE', 'f','const', 'i',6, 'sp','combine');
m1.addPolicy('o','EN', 'f','scale', 'sp','combine', 'i',6, 'c','q r bo');
m1.addPolicy('o','EN', 'f','const', 'i',6, 'sp','combine');
m1.addPolicy('o','E2', 'f','scale', 'sp','combine', 'i',6, 'c','q r bo');

m1.addPolicy('o','E2', 'f','scale', 'sp','sonly', 'i',1, 'j',1, ...
   'c','r','nb',1);
policies{1} = m1.policy;

for ipol = 1:length(policies)
   for idata = 1:size(dsets,1)
      filePre=[pname{ipol},'/',dname{idata}];
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
      
      % Create fitme object
      fact  = MFactory;
      fact.policy = policies{ipol};
      fact.makeMixInfo(dsets{idata,1}.atomTypes);
      f1    = fact.makeFitme(dsets{idata,1});
      ftest = fact.makeFitme(dsets{idata,2});
      
      fprintf(summaryFile,'train and test starting error \n');
      f1.printEDetails(summaryFile);
      ftest.printEDetails(summaryFile);
      
      %
      startName = [topDir,filePre,'/start.mat'];
      toSave = {'fact','f1','ftest','currentTrainErr', ...
         'currentPar','currentErr'};
      if (exist(startName,'file'))
         fprintf(1,'LOADING START \n');
         fprintf(summaryFile,'LOADING START \n');
         load(startName,toSave{:});
      else
         [currentTrainErr,currentPar,currentErr] = ...
            contextFit3(f1,ftest,maxIter);
         save(startName,toSave{:});
      end
      
      str1 = 'initial error %12.5f test %12.5f \n';
      fprintf(1,str1,currentTrainErr,currentErr);
      fprintf(summaryFile,str1,currentTrainErr,currentErr);
      f1.printEDetails(summaryFile);
      ftest.printEDetails(summaryFile);
      
      ticID = tic;
      for iter = 1:3
         allName = [topDir,filePre,'/all-',num2str(iter),'.mat'];
         if (exist(allName,'file'))
            fprintf(1,'LOADING ITERATION %i \n',iter);
            fprintf(summaryFile,'LOADING ITERATION %i \n',iter);
            load(allName,toSave{:});
         else
            fprintf(1,'STARTING ITERATION %i \n',iter);
            fprintf(summaryFile,'STARTING ITERATION %i \n',iter);
            % unfix 1 level of context
            for imix = 1:length(f1.mixers)
               mix = f1.mixers{imix};
               for ipar = 1:length(mix.fixed)
                  if (mix.fixed(ipar) == 1)
                     mix.fixed(ipar) = 0;
                     break;
                  end
               end
            end
            [currentTrainErr,currentPar,currentErr] = ...
               contextFit3(f1,ftest,maxIter);
            save(allName,toSave{:});
         end
         str2 = 'context error %12.5f test %12.5f \n';
         fprintf(1,str2,currentTrainErr,currentErr);
         fprintf(summaryFile,str2,currentTrainErr,currentErr);
         f1.printEDetails(summaryFile);
         ftest.printEDetails(summaryFile);
         
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
end