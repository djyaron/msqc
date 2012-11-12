%% Fitting multiple molecules, using makeFitme
clear classes;
topDir = 'T:\matdl\yaron\8-3-12\quadratic-sp\';
%topDir = 'scalehybridparallel/';
ics = 1;

ftype = 2;
runParallel = 0;
showPlots = 1;

%trainC{1}  = {'h2',3,'envs',1:100};
%testC{1} = []; %{'h2',4,'envs',1:100};
%filePrefix{1} = 'h2-geom3';

for igeom =[1:3 20:23]
trainC{1}  = {'h2',[],'ch4',igeom,'envs',1:100};
testC{1} = []; %{'h2',[],'ch4',1:17,'envs',20:30};
filePrefix{1} = ['ch4-geom',num2str(igeom)];

trainC{2}  = {'h2',[],'ethane',1:7,'envs',1:10};
testC{2} = {'h2',[],'ethane',1:7,'envs',20:30};
filePrefix{2} = 'c2h6';

trainC{3}  = {'h2',[],'ethylene',1:7,'envs',1:10};
testC{3} = {'h2',[],'ethylene',1:7,'envs',20:30};
filePrefix{3} = 'c2h4z2';

trainC{4}  = {'h2',[],'ch4',1:7,'ethane',1:7,'envs',1:10};
testC{4} = {'h2',[],'ch4',1:7,'ethane',1:7,'envs',20:30};
filePrefix{4} = 'ch4-c2h6';

trainC{5}  = {'h2',[],'ch4',1:7,'ethane',1:7,'ethylene',1:7,'envs',1:10};
testC{5} = {'h2',[],'ch4',1:7,'ethane',1:7,'ethylene',1:7,'envs',20:30};
filePrefix{5} = 'ch4-c2h6-c2h4';

trainC{6}  = {'h2',[],'ch4',1:19,'ethane',1:7,'envs',1:10};
testC{6} = {'h2',[],'ch4',1:19,'ethane',1:7,'envs',20:30};
filePrefix{6} = 'ch4f-c2h6';

trainC{7}  = {'h2',[],'ch4',1:19,'ethane',1:7,'ethylene',1:7,'envs',1:10};
testC{7} = {'h2',[],'ch4',1:19,'ethane',1:7,'ethylene',1:7,'envs',20:30};
filePrefix{7} = 'ch4f-c2h6-c2h4';

trainC{8}  = {'h2',[],'propane',1:7,'envs',1:10};
testC{8} = {'h2',[],'propane',1:7,'envs',20:30};
filePrefix{8} = 'c3h8';

trainC{9}  = {'h2',[],'ch4',1:19,'ethane',1:7,'propane',1:7,'envs',1:10};
testC{9} = {'h2',[],'ch4',1:19,'ethane',1:7,'propane',1:7,'envs',20:30};
filePrefix{9} = 'ch4f-c2h6-c3h8';

commonIn = {};
%
for iC = ics
   trainIn = trainC{iC};
   testIn = testC{iC};
   filePre = filePrefix{iC};
   ke = [];
   en = [];
   e2 = [];
   f1 = [];
   if (exist([topDir,filePre],'dir') ~= 7)
      status = mkdir([topDir,filePre]);
   end
   summaryName = [topDir,filePre,'\summary.txt'];
   summaryFile = fopen(summaryName,'a');
   for iPar = 0:9
      if (iPar == 0)
         fprintf(summaryFile,' %s \n','no shift, no context');
         if (ftype == 2)
            iP = 1;
         else
            iP = 0;
         end
         
         ke.H = Mixer(iP,1,'ke.H',ftype);
         ke.Cs = Mixer(iP,1,'ke.C',ftype);
         ke.Cp = ke.Cs;
         ke.HH = Mixer(iP,1,'ke.HH',ftype);
         ke.CH = Mixer(iP,1,'ke.CH',ftype);
         ke.CH.hybrid = 1;
         ke.CCs = Mixer(iP,1,'ke.CCs',ftype);
         ke.CCs.hybrid = 1;
         ke.CCp = Mixer(iP,1,'ke.CCp',ftype);
         ke.CCp.hybrid = 2;
         
         en.H = Mixer(iP,1,'en.H',ftype);
         en.Cs = Mixer(iP,1,'en.C',ftype);
         en.Cp = en.Cs;
         en.HH = Mixer(iP,1,'en.HH',ftype);
         en.CH = Mixer(iP,1,'en.CH',ftype);
         en.CH.hybrid = 1;
         en.HC = en.CH;
         en.CCs = Mixer(iP,1,'en.CCs',ftype);
         en.CCs.hybrid = 1;
         en.CCp = Mixer(iP,1,'en.CCp',ftype);
         en.CCp.hybrid = 2;
         
         e2.H = Mixer(iP,1,'e2.H',ftype);
         e2.C = Mixer(iP,1,'e2.C',ftype);
         e2.HH = Mixer(iP,1,'e2.HH',ftype);
         e2.CC = Mixer(iP,1,'e2.CC',ftype);
         e2.CH = Mixer(iP,1,'e2.CH',ftype);
         
         if (isempty(testIn))
            f1 = makeFitme(trainIn{:},commonIn{:},'enstructh',en, ...
               'kestructh',ke,'e2struct',e2);
         else
            ftest = makeFitme(testIn{:},commonIn{:},'enstructh',en, ...
               'kestructh',ke,'e2struct',e2,'plot',2);
            ftest.parallel = runParallel;
            ftest.plot = showPlots;
            f1 = makeFitme(trainIn{:},commonIn{:},'enstructh',en, ...
               'kestructh',ke,'e2struct',e2,'testFitme',ftest);
         end
         f1.plot = showPlots;
         f1.parallel = runParallel;
      elseif (iPar == 1)
         fprintf(summaryFile,' %s \n','with shift, no context');
         ke.Cp = ke.Cs.deepCopy;
         en.Cp = en.Cs.deepCopy;         
      elseif (iPar == 2)
         fprintf(summaryFile,' %s \n','with shift, no context');
         for m1 = [ke.H ke.Cs ke.Cp en.H en.Cs en.Cp]
            if (ftype == 2)
               m1.funcType = 3;
            else
               m1.funcType = 4;
            end
            m1.par(2) = 0;
            m1.fixed(2) = 0;
         end
      elseif (iPar == 3)
         fprintf(summaryFile,' %s \n','ke diag linear');
         for m1 = [ke.H ke.Cs ke.Cp]
            m1.mixType = 2;
            m1.par(3) = m1.par(2);
            m1.par(2) = 0;
            m1.fixed(3) = 0;
         end
      elseif (iPar == 3)
         fprintf(summaryFile,' %s \n','ke diag linear');
         for m1 = [ke.H ke.Cs ke.Cp]
            m1.mixType = 2;
            m1.par(3) = m1.par(2);
            m1.par(2) = 0;
            m1.fixed(3) = 0;
         end
      elseif (iPar == 4)
         fprintf(summaryFile,' %s \n','en diag linear');
         for m1 = [en.H en.Cs en.Cp]
            m1.mixType = 2;
            m1.par(3) = m1.par(2);
            m1.par(2) = 0;
            m1.fixed(3) = 0;
         end
      elseif (iPar == 5)
         fprintf(summaryFile,' %s \n','ke diag quad');
         for m1 = [ke.H ke.Cs ke.Cp]
            m1.mixType = 22;
            m1.par(4) = m1.par(3);
            m1.par(3) = 0.0;
            m1.fixed(4) = 0;
         end
      elseif (iPar == 6)
         fprintf(summaryFile,' %s \n','en diag quad');
         for m1 = [en.H en.Cs en.Cp]
            m1.mixType = 22;
            m1.par(4) = m1.par(3);
            m1.par(3) = 0.0;
            m1.fixed(4) = 0;
         end
      elseif (iPar == 7)
         fprintf(summaryFile,' %s \n','ke off diag BO');
         for m1 = ke.CH
            m1.mixType = 3;
            m1.par(2) = 0.0;
            m1.fixed(2) = 0;
         end
      elseif (iPar == 8)
         fprintf(summaryFile,' %s \n','en off diag BO');
         for m1 = en.CH
            m1.mixType = 3;
            m1.par(2) = 0.0;
            m1.fixed(2) = 0;
         end
      elseif (iPar == 9)
         fprintf(summaryFile,' %s \n','e2 off diag BO');
         for m1 = e2.CH
            m1.mixType = 3;
            m1.par(2) = 0.0;
            m1.fixed(2) = 0;
         end
      end
      
      dataDir = [topDir,filePre,'/fit-',num2str(iPar),'/'];
      allFile = [dataDir,'all.mat'];
      if (exist(allFile,'file'))
         load(allFile,'pt');
         f1.setPars(pt);
         disp([filePre,' iC ',num2str(iC),'fit# ',num2str(iPar),...
            'loaded from file']);
      else
         options = optimset('DiffMinChange',1.0e-5,'TolFun',1.0e-3, ...
            'TolX',3.0e-3,'MaxFunEvals',500);
         if (exist(dataDir,'dir') ~= 7)
            status = mkdir(dataDir);
         end
         diary([dataDir,'out.diary']);
         diary on;
         tic
         lowLimits = zeros(f1.npar,1);
         highLimits = lowLimits;
         i1 = 1;
         for imix = 1:length(f1.mixers)
            mix = f1.mixers{imix};
            if ((mix.funcType == 2)||(mix.funcType == 3))
               lowLimits(i1) = 0.0;
               highLimits(i1) = 10.0;
               i1 = i1+1;
               for i2 = 2:mix.npar
                  lowLimits(i1) = -inf;
                  highLimits(i1) = inf;
                  i1 = i1+1;
               end
            else
               for i2 = 1:mix.npar
                  lowLimits(i1) = -inf;
                  highLimits(i1) = inf;
                  i1 = i1+1;
               end
            end
         end
         restartFile = [dataDir,'restart.mat'];
         f1.restartFile = restartFile;
         if (exist(restartFile,'file'))
            disp('loading restart file');
            load(restartFile);
            start = ptSave;
            f1.itcount = itSave;
            f1.errTrain = errTrainSave;
            f1.errTest = errTestSave;
         else
            start = f1.getPars;
         end
         [pt,resnorm,residual,exitflag,output,lambda,jacobian] = ...
            lsqnonlin(@f1.err, start,lowLimits,highLimits,options);
         %  options = LMFnlsq;
         %  options.Display =1;
         %  options.FunTol = 1.0e-6;
         %  options.XTol = 1.0e-5;
         %  [pt,resnorm, CNT, Res, XY] = LMFnlsq(@f1.err,start',options);
         clockTime = toc
         pt
         resnorm
         f1.printMixers;
         f1.printEDetails(summaryFile);
         save([dataDir,'all.mat']);
         
         diary off;
         if (f1.plot)
            figure(799); saveas(gcf,[dataDir,'error.fig']);
            if (~isempty(find(cellfun(@(x)isequal(lower(x),'ch4'),trainIn)) ))
               figure(801); saveas(gcf,[dataDir,'ch4-train.fig']);
               figure(811); saveas(gcf,[dataDir,'ch4-test.fig']);
            end
            if (~isempty(find(cellfun(@(x)isequal(lower(x),'ethane'),trainIn)) ))
               figure(802); saveas(gcf,[dataDir,'c2h6-train.fig']);
               figure(812); saveas(gcf,[dataDir,'c2h6-test.fig']);
            end
            if (~isempty(find(cellfun(@(x)isequal(lower(x),'ethylene'),trainIn)) ))
               figure(803); saveas(gcf,[dataDir,'c2h4-train.fig']);
               figure(813); saveas(gcf,[dataDir,'c2h4-test.fig']);
            end
         end
      end
   end
   fclose(summaryFile);
end
end