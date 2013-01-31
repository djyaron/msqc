clear classes;
close all;
rootDir = 'C:/matdl/yaron/2010a/';
maxIter = 500;

for propWeights = 0
for EtotWeight = 1 %$[1e7 1 5 10 20 0.1 0.5 30 0.25 0.75]
topDir = [rootDir,'w',num2str(EtotWeight)];
if (propWeights)
   topDir = [topDir,'p/'];
else
   topDir = [topDir,'/'];
end

printDetailsOnLoad = 0;

dname1=cell(0,0); fname=cell(0,0);
gTrain=cell(0,0); gTest=cell(0,0);
eTrain=cell(0,0); eTest=cell(0,0); pnn=[];

h2 = 1;
fname{end+1} = 'h2Dat'; dname1{end+1}=fname{end};
gTrain{end+1}=[1 3 5 7];  eTrain{end+1}=1:2:20;
gTest{end+1} =[2 4 6];    eTest{end+1} =2:2:20; pnn(end+1) = 796;
ch4r = 2;
fname{end+1} = 'ch4rDat'; dname1{end+1}=fname{end};
gTrain{end+1}=1:10;  eTrain{end+1}=1:2:20;
gTest{end+1} =11:20; eTest{end+1} =2:2:20; pnn(end+1) = 791;
ethaner = 3;
fname{end+1} = 'ethanerDat'; dname1{end+1}=fname{end};
gTrain{end+1}=1:10;  eTrain{end+1}=1:2:20;
gTest{end+1} =11:20; eTest{end+1} =2:2:20; pnn(end+1) = 792;
ethylener =4;
fname{end+1} = 'ethylenerDat'; dname1{end+1}=fname{end};
gTrain{end+1}=1:10;  eTrain{end+1}=1:2:20;
gTest{end+1} =11:20; eTest{end+1} =2:2:20; pnn(end+1) = 793;

toFit = {[ethaner]};%{[ch4r], [ethaner], [ch4r,ethaner]};

dsets = cell(1,2);
dname = cell(1,1);
ic = 0;
for ifit = 1:length(toFit)
   ic = ic+1;
   name = '';
   for idata = toFit{ifit}
      if (isempty(name))
         name = dname1{idata};
      else
         name = [name,'_',dname1{idata}];
      end
   end
   mtrain = MSet;
   mtest  = MSet;
   for idata = toFit{ifit}
      dfile = ['datasets/',fname{idata},'.mat'];
      mtrain.addData(dfile, gTrain{idata}, eTrain{idata} ,1,pnn(idata));
      mtest.addData(dfile,  gTest{idata},  eTest{idata}  ,1,pnn(idata));
   end
   dname{ic} = name;
   dsets{ic,1} = mtrain;
   dsets{ic,2} = mtest;
end

% CREATE POLICIES
pname = {'hybridslater1'};
%
for ipol = 1:length(pname)
   for idata = 1:size(dsets,1)
      filePre=[pname{ipol},'/',dname{idata}];
      dataDir = [topDir,filePre];
      if (exist(dataDir,'dir') ~= 7)
         status = mkdir(dataDir);
      end
      copyfile('c:/dave/apoly/msqc/contextFactory2.m',...
          [dataDir,'/contextFactory2.m']);
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
      fact.setPolicies(pname{ipol});
      fact.makeMixInfo(dsets{idata,1}.atomTypes);
      f1    = fact.makeFitme(dsets{idata,1});
      ftest = fact.makeFitme(dsets{idata,2});
      input junk;
      % Add weighting
      %f1.setWeights(EtotWeight,propWeights);
      
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
         loaded = 1;
      else
         [currentTrainErr,currentPar,currentErr] = ...
            contextFit3(f1,ftest,maxIter);
         save(startName,toSave{:});
         loaded = 0;
      end
      
      str1 = 'initial error %12.5f test %12.5f \n';
      fprintf(1,str1,currentTrainErr,currentErr);
      fprintf(summaryFile,str1,currentTrainErr,currentErr);
      if (~loaded || printDetailsOnLoad)
         f1.printEDetails(summaryFile);
         ftest.printEDetails(summaryFile);
      end
      ticID = tic;
      for iter = 1:3
         allName = [topDir,filePre,'/all-',num2str(iter),'.mat'];
         if (exist(allName,'file'))
            fprintf(1,'LOADING ITERATION %i \n',iter);
            fprintf(summaryFile,'LOADING ITERATION %i \n',iter);
            load(allName,toSave{:});
            loaded = 1;
         else
            loaded = 0;
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
         if (~loaded || printDetailsOnLoad)
            f1.printEDetails(summaryFile);
            ftest.printEDetails(summaryFile);
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
end

end
end