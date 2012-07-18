%% Fitting multiple molecules, using makeFitme
%clear classes;
topDir = 'C:\matdl\yaron\7-13-12\scalehybrid\parallel2\fixedE2\';
%topDir = 'scalehybridparallel/';
runParallel = 1;
% if (Aprocess == 1)
%    ics = [1 6];
% elseif (Aprocess == 2)
%    ics = [2 7];
% else
%    ics = [3 9];
% end
ics = [1 2 3 6 7 9];
%Aprocess = 1;
%ics = [1 2 3 6 7];
% if (Aprocess == 1)
%    runParallel = 1;
%    topDir = 'scalehybridparallel/fixedE2/';
% else
%    runParallel = 0;
%    topDir = 'scalehybridparics/tight10/';
% end
%trainC{1}  = {'h2',2:7,'envs',1:10};
%testC{1} = {'h2',2:7,'envs',20:30};
ftype = 2;
trainC{1}  = {'h2',[],'ch4',1:17,'envs',1:10};
testC{1} = {'h2',[],'ch4',1:17,'envs',20:30};
filePrefix{1} = 'ch4';

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
for iC = ics% [1 2 3 4 6 7]
   trainIn = trainC{iC};
   testIn = testC{iC};
   filePre = filePrefix{iC};
   ke = [];
   en = [];
   e2 = [];
   f1 = [];
   for iPar = 1:5
      if (iPar == 1)
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
         ke.CC.hybrid = 1;
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
         %           ftest = makeFitme(testIn{:},commonIn{:},'enstruct1',en, ...
         %              'kestruct',ke,'e2struct',e2,'plot',2);
         %           ftest.parallel = 0;
         %           ftest.plot = 0;
         f1 = makeFitme(trainIn{:},commonIn{:},'enstructh',en,'kestructh',ke, ...
            'e2struct',e2);%,'testFitme',ftest);
         f1.plot = 0;
         f1.parallel = runParallel;
         %pst = [ -1.2840    1.4139   -0.9773   -0.1648    2.9684   -1.7791    5.7310   -9.6449    8.0355  12.5867   -0.1876   -0.1118    2.0048   -0.3105];
         %f1.setPars(pst);
         %          f1.parHF = zeros(size(f1.getPars));
         %          etest1 = f1.err(f1.getPars);
         %          f1.parHF = zeros(size(f1.getPars));
         %          f1.parallel = 1;
         %          etest2 = f1.err(f1.getPars);
         %          input('hi');
      elseif (iPar == 2) % add constants
         for m1 = [ke.H ke.Cs en.H en.Cs]
            if (ftype == 2)
               m1.funcType = 3;
            else
               m1.funcType = 4;
            end
            m1.par(2) = 0;
            m1.fixed(2) = 0;
         end
      elseif (iPar == 3) % add context sensitive
         for m1 = [ke.H ke.Cs en.H en.Cs]
            m1.mixType = 2;
            m1.par(3) = m1.par(2);
            m1.fixed(3) = 0;
            m1.par(2) = 0;
         end
         for m1 = [e2.H e2.C]
            m1.mixType = 2;
            m1.par(2) = 0;
            m1.fixed(2) = 0;
         end
         for m1 = [ke.HH ke.CH ke.CCs ke.CCp en.HH en.CH en.HC en.CCs ...
               en.CCp] % e2.HH e2.CH e2.CC]
            m1.mixType = 3;
            m1.par(2) = 0;
            m1.fixed(2) = 0;
         end
      elseif (iPar == 4) % add context sensitive (bond length)
         for m1 = [ke.HH ke.CH ke.CCs ke.CCp en.HH en.CH en.HC en.CCs ...
               en.CCp] % e2.HH e2.CH e2.CC]
            m1.mixType = 4;
            m1.par(2) = 0;
         end
      elseif (iPar == 5) % add context sensitive (both)
         for m1 = [ke.HH ke.CH ke.CCs ke.CCp en.HH en.CH en.HC en.CCs ...
               en.CCp] % e2.HH e2.CH e2.CC]
            m1.mixType = 5;
            m1.par(3) = m1.par(2);
            m1.par(2) = 0.0;
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
         %f1.parallel = 0;
         %etest3 = f1.err(start);
         %f1.parallel = 1;
         [pt,resnorm,residual,exitflag,output,lambda,jacobian] = ...
            lsqnonlin(@f1.err, start,lowLimits,highLimits,options);
         %       options = LMFnlsq;
         %      options.Display =1;
         %    options.FunTol = 1.0e-6;
         %    options.XTol = 1.0e-5;
         %    [pt,resnorm, CNT, Res, XY] = LMFnlsq(@f1.err,start',options);
         clockTime = toc
         pt
         resnorm
         f1.printMixers;
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
end


