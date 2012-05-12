%% Fitting multiple molecules, using makeFitme
clear classes;
trainC{1}  = {'h2',[],'ch4',1:19,'envs',1:10};
testC{1} = {'h2',[],'ch4',1:19,'envs',20:30};
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

commonIn = {};

for iC = 3:3
   for iPar = 1:4
      trainIn = trainC{iC};
      testIn = testC{iC};
      filePre = filePrefix{iC};
      if (iPar == 1)
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
      elseif (iPar == 2)
         ke.H.mixType = 2;    ke.H.par(2) = 0;    ke.H.fixed(2) = 0;
         ke.Cs.mixType = 2;   ke.Cs.par(2) = 0;   ke.Cs.fixed(2) = 0;
         ke.HH.mixType = 3;   ke.HH.par(2) = 0;   ke.HH.fixed(2) = 0;
         ke.CsH.mixType = 3;  ke.CsH.par(2) = 0;  ke.CsH.fixed(2) = 0;
         ke.CsCs.mixType = 3; ke.CsCs.par(2) = 0; ke.CsCs.fixed(2) = 0;
      elseif (iPar == 3)
         en.H.mixType = 2;    en.H.par(2) = 0;      en.H.fixed(2) = 0;
         en.Cs.mixType = 2;   en.Cs.par(2) = 0;     en.Cs.fixed(2) = 0;
         en.HH.mixType = 3;   en.HH.par(2) = 0;     en.HH.fixed(2) = 0;
         en.CsH.mixType = 3;  en.CsH.par(2) = 0;    en.CsH.fixed(2) = 0;
         en.CsCs.mixType = 3; en.CsCs.par(2) = 0;   en.CsCs.fixed(2) = 0;
      elseif (iPar == 4)
         ke.Cp = ke.Cs.deepCopy();     ke.Cs.desc ='ke.Cs';     ke.Cp.desc = 'ke.Cp';
         ke.CpH = ke.CsH.deepCopy();   ke.CsH.desc ='ke.CsH';   ke.CpH.desc = 'ke.CpH';
         ke.CsCp = ke.CsCs.deepCopy(); ke.CsCs.desc ='ke.CsCs'; ke.CsCp.desc = 'ke.CsCp';
         ke.CpCp = ke.CsCs.deepCopy();                          ke.CpCp.desc = 'ke.CpCp';
         
         en.Cp = en.Cs.deepCopy();     en.Cs.desc ='en.Cs';     en.Cp.desc = 'en.Cp';
         en.CpH = en.CsH.deepCopy();   en.CsH.desc ='en.CsH';   en.CpH.desc = 'en.CpH';
         en.CsCp = en.CsCs.deepCopy(); en.CsCs.desc ='en.CsCs'; en.CsCp.desc = 'en.CsCp';
         en.CpCp = en.CsCs.deepCopy();                          en.CpCp.desc = 'en.CpCp';
         
      end
      
      dataDir = ['tmp/',filePre,'/fit-',num2str(iPar),'/'];
      ftest = makeFitme(trainIn{:},commonIn{:},'enstruct',en,'kestruct',ke);
      f1 = makeFitme(testIn{:},commonIn{:},'enstruct',en,'kestruct',ke, ...
         'testFitme',ftest);
      limits = [];
      options = optimset('DiffMinChange',1.0e-5,'TolFun',1.0e-4,'TolX',1.0e-3);
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
      save([dataDir,'all.mat']);
      diary off;
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