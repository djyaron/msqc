%% factoryUse
function factoryUse()
%diffEnv;
%testEDetails();
bigPlotsGenData();
%bigPlot1(); % each molecule gets own window
%bigPlot2(); % ethane train and test in one window
%bigPlot3(); % all molecules in one window
%bigTable(); % Print table of all results
%weightPlotGenData();
%weightPlotGenData2();
%weightPlot1();
%weightPlot2();
%weightPlot3();
end
%% bigPlot2
function bigPlot2()
load('D:\matdl\yaron\dec12e\bigplot.mat');
ethTrain = 1; ethTest = 2; meth = 3; prop = 4; nbut = 5; tbut = 6;
toplot = {'ke','H','C','e2','etot','all'};
pcol = {'g','c','b','m','r','k'};
lineType = {'-',':'};
close all;
er1 = cell(size(errs));
niter = size(errs,2);
for idata = [ethTrain ethTest]
   for iter = 1:niter
      disp(['idata iter ',num2str(idata),' ',num2str(iter)]);
      er1{idata,iter} = EDetails(errs{idata,iter},0,0);
   end
end
leg1 = cell(0,0);
leg2 = cell(0,0);
for idata = [ethTrain ethTest]
   for itype = 1:length(toplot) % data type
      for il = 1 % err or standard deviation
         x = zeros(niter,1);
         y = zeros(niter,1);
         for iter = 1:niter
            x(iter) = iter;
            y(iter) = getfield(er1{idata,iter},toplot{itype});
         end
         switch itype
            case {5,6}
               figure(71)
               hold on
               plot(x,y,[pcol{itype},lineType{idata}]...
                  ,'LineWidth',2);
               leg1{end+1} = toplot{itype};
            otherwise
               figure(72)
               hold on
               plot(x,y,[pcol{itype},lineType{idata}] ...
                  ,'LineWidth',2);
               leg2{end+1} = toplot{itype};
         end
      end
   end
end
%set(gca,'YGRID','on');
figure(71);
legend(leg1);
ylabel('RMS error (kcal/mol)');
xlabel('iteration');
set(gca,'YSCALE','log');
set(gca,'YTick',[2 3 4 5 6 8 10 20 40 60 80 100 200 400 600 800]);
figure(72);
legend(leg2);
ylabel('RMS error (kcal/mol)');
xlabel('iteration');
set(gca,'YSCALE','log');
set(gca,'YTick',[2 3 4 5 6 8 10 20 40 60 80 100 200 400 600 800]);
for isig = 2:length(iterSig)
   for ifig = [71 72]
      figure(ifig);
      hold on;
      xL = iterSig(isig);
      yL = get(gca,'YLim');
      line([xL xL],yL,'color','k');
   end
end

printFig(71,'ethaneTot');
printFig(72,'ethaneOpers');

end
%% printFig
function printFig(num,name)
figure(num)
box on;
set(findall(gcf,'type','text'),'fontSize',16,'fontName','Times New Roman')
input(['adjust figure ',num2str(num)]);
print('-depsc',['figs2/',name]);

end
%% bigPlot3
function bigPlot3()
load('C:\matdl\yaron\dec12e\bigplot.mat');
ethTrain = 1; ethTest = 2; meth = 3; prop = 4; nbut = 5; tbut = 6;
toplot = {'etot','oper'};
dnames = {'Ethane (train)','Ethane (test)','Methane',...
   'Propane','nButane','tButane'};
pcol = {'k','b','g','m','r','c'};
lineType = {'-',':'};
% errs{dataset, iter, err/sd}
close all;
erTemp = cell(size(errs));
er2 = cell(0,0);
niter = size(errs,2);
for iavg = 1:2
   for idata = 1:tbut
      for iter = 1:niter
         disp(['idata iter ',num2str(idata),' ',num2str(iter)]);
         erTemp{idata,iter} = EDetails(errs{idata,iter},iavg-1,0);
      end
   end
   er2{end+1} = erTemp;
end
leg1 = cell(0,0);
leg2 = cell(0,0);
leg3 = cell(0,0);
for itype = 1:2 % etot oper
   for idata = 1:tbut
      for iavg = 1:2 % is average subtracted
         x = zeros(niter,1);
         y = zeros(niter,1);
         for iter = 1:niter
            x(iter) = iter;
            er1 = er2{iavg};
            y(iter) = getfield(er1{idata,iter},toplot{itype});
         end
         if ((iavg == 1) && (itype == 1))
            figure(81)
            hold on
            plot(x,y,[pcol{idata},lineType{itype}]...
               ,'LineWidth',2);
            leg1{end+1} = dnames{idata};
         elseif ((iavg == 2) && (itype == 1))
            figure(82)
            hold on
            plot(x,y,[pcol{idata},lineType{itype}] ...
               ,'LineWidth',2);
            leg2{end+1} = dnames{idata};
         elseif ((iavg == 2) && (itype == 2))
            figure(83)
            hold on
            plot(x,y,[pcol{idata},lineType{1}] ...
               ,'LineWidth',2);
            leg3{end+1} = dnames{idata};            
         end
      end
   end
end
%set(gca,'YGRID','on');
figure(81);
title('Etot absolute')
legend(leg1);
ylabel('RMS error (kcal/mol)');
xlabel('iteration');
set(gca,'YSCALE','log');
set(gca,'YTick',[2 3 4 5 6 8 10 20 40 60 80 100 200 400 600 800]);
figure(82);
title('Etot mean substracted');
legend(leg2);
ylabel('RMS error (kcal/mol)');
xlabel('iteration');
set(gca,'YSCALE','log');
set(gca,'YTick',[2 3 4 5 6 8 10 20 40 60 80 100 200 400 600 800]);
figure(83);
title('Oper mean substracted');
legend(leg2);
ylabel('RMS error (kcal/mol)');
xlabel('iteration');
set(gca,'YSCALE','log');
set(gca,'YTick',[2 3 4 5 6 8 10 20 40 60 80 100 200 400 600 800]);
for isig = 2:length(iterSig)
   for ifig = [81 82 83]
      figure(ifig);
      hold on;
      xL = iterSig(isig);
      yL = get(gca,'YLim');
      line([xL xL],yL,'color','k');
   end
end

%printFig(81,'allMolsAbsError');
%printFig(82,'allMolsSubMean');
printFig(83,'allMolsOpers');

end
%% errOut
function errOut = getTotalError(errIn,nC,nH)
errOut.ke = errIn.ke;
errOut.H = nH * errIn.H;
errOut.C = nC * errIn.C;
errOut.e2 = errIn.e2;
errOut.etot = errIn.etot;
errOut.all = errOut.ke + errOut.H + errOut.C + errOut.e2 + errOut.etot;
end
%% diffEnv
function diffEnv()
%clear classes;
%close all;

fct = cell(0,0);
ms = cell(0,0);
fitme = cell(0,0);
fname = cell(0,0);
% % dataExt = {'','-1c','-linrho','-diponly'};
% % for iext = 1:4
% %    load(['C:\matdl\yaron\11-29-12\factory\hybrid1\ch4r',dataExt{iext},'\all-3.mat'])
% %    fct{end+1} = fact;  fitme{end+1} = f1;  ms{end+1} = MSet.fromFitme(f1);
% % end
% % fname = {'orig','1c','lin','dip'};
% 
% dataExt = {'ethanerDat','propanerDat'};
% for iext = 1:length(dataExt)
%    load(['C:\matdl\yaron\dec12a\hybridslater\',dataExt{iext},'\all-3.mat'])
%    fct{end+1} = fact;  fitme{end+1} = f1;  ms{end+1} = MSet.fromFitme(f1);
% end
% fname = dataExt;
% 

dataExt = {''};%,'-1c','-linrho','-diponly'};
for iext = 1:length(dataExt)
   load(['C:\matdl\yaron\11-29-12\factory\hybrid1\ch4r',dataExt{iext},'\all-3.mat'])
   for imix = 1:length(f1.mixers)
      f1.mixers{imix}.bonded = 1;
   end
   f1.mixers{end}.bonded = 0;
   fct{end+1} = fact;  fitme{end+1} = f1;  ms{end+1} = MSet.fromFitme(f1);
end
fname = {'orig','1c','lin','dip'};


nfact = length(fct);
ndata = length(fct);
esum = zeros(nfact,ndata);
for ifact = 1:nfact
   for idata = 1:ndata
      disp(['factory: ',fname{ifact},'   data: ',fname{idata}]);
      f2 = fct{ifact}.makeFitme(ms{idata},fitme{idata});
      f2.silent = 1;
%       res1 = f1.printEDetails;
%       disp(' ');
%       esum(ifact,idata) = res1{1}.etot;
      
         [e0 plotnum etype modelnum envnum] = f2.err(f2.getPars);
   f2.printEDetails;
   e0 =e0*627.509;
   eke = e0(etype==1);
   eH  = e0(etype==11);
   eC  = e0(etype==16);
   e2  = e0(etype==2);
   etot = e0(etype==3);
   disp(['mean ke ',num2str(mean(eke)),' H ',num2str(mean(eH)),' C ',...
      num2str(mean(eC)),' E2 ',num2str(mean(e2)),' Etot ', num2str(mean(etot))]);
   disp(['std  ke ',num2str(std(eke)),' H ',num2str(std(eH)),' C ',...
      num2str(std(eC)),' E2 ',num2str(std(e2)),' Etot ', num2str(std(etot))]);

      
   end
end
for idata = 1:nfact
   fprintf(ofile,'\t%s ',fname{idata});
end
fprintf(ofile,'\n');

for ifact = 1:nfact
   fprintf(ofile,'%s ',fname{ifact});
   for idata = 1:ndata
      fprintf(ofile,'\t %5.2f ',esum(ifact,idata));
   end
   fprintf(ofile,'\n');
end

end
%% bigPlotsGenData
function bigPlotsGenData()
dataroot = ...
   'd:\matdl\yaron\dec12e\iter100\w1\hybridslater1\ethanerDat';
%dataroot = 'C:\matdl\yaron\dec12e\c1\w1\h2fits\h2Dat';
lfiles = {'start.mat', 'all-1.mat', 'all-2.mat', 'all-3.mat'};

load([dataroot,'\start.mat']);
ms = cell(0,0);
mtemp = MSet;
mtemp.addData('datasets/ethanerDat.mat',1:10,1:2:20,1,791);
ms{end+1} = mtemp;
mtemp = MSet;
mtemp.addData('datasets/ethanerDat.mat',11:20,2:2:20,1,791);
ms{end+1} = mtemp;
mtemp = MSet;
mtemp.addData('datasets/ch4rDat.mat',1:10,1:2:20,1,791);
ms{end+1} = mtemp;
mtemp = MSet;
mtemp.addData('datasets/propaner2-orig.mat',1:10,1:2:20,1,791);
ms{end+1} = mtemp;
mtemp = MSet;
mtemp.addData('datasets/butaner2-orig.mat',1:10,1:2:20,1,791);
ms{end+1} = mtemp;
mtemp = MSet;
mtemp.addData('datasets/tbutaner2-orig.mat',1:10,1:2:20,1,791);
ms{end+1} = mtemp;

% errs{dataset, iter}
errs = cell(0,0);
iterSig = [];
iter = 0;
for i = 1:length(lfiles)
  load([dataroot,'\',lfiles{i}]);
  ngood = length(monitor.accepted);
  iterSig = [iterSig iter+1];
  for igood = 1:ngood
     iter = iter + 1;
     ipar = monitor.accepted(igood);
     pars = monitor.param{ipar};
     disp([lfiles{i},' infile ',num2str(ipar),' total ',...
        num2str(iter)]);
     % Get the parameters into the factory
     f1.setPars(pars);
     for iset = 1:length(ms)
        fm2=fact.makeFitme(ms{iset});
        saveWeights = fm2.operWeights;
        obj.operWeights = [];
        [eout plotnum etype modelnum envnum] = fm2.err(fm2.getPars);
        fm2.operWeights = saveWeights;
        tsave.err = eout;
        tsave.plotnum = plotnum;
        tsave.etype = etype;
        tsave.modelnum = modelnum;
        tsave.envnum = envnum;
        errs{iset,iter} = tsave;
     end
  end
end
save([dataroot, '\bigplot.mat'],'errs','iterSig'); %,'emeth','smeth');
disp('done');
end
%% bigPlot1
function bigPlot1()
% each molecule gets its own window
dataroot = ...
   'C:\matdl\yaron\dec12e\iter100\w1\hybridslater1\ethanerDat';
load([dataroot, '\bigplotdeb.mat']); %,'emeth','smeth');
toplot = {'ke','H','C','e2','etot'};
psym = {'ko','c^','b^','ks','ro'};

ltype = {'-','--'};

% errs{dataset, iter, err/sd}
for idata = 1:size(errs,1) % data set
   for itype = 1:length(toplot) % data type
      for il = 1:2 % err or standard deviation
         niter = size(errs,2);
         x = zeros(niter,1);
         y = zeros(niter,1);
         for iter = 1:niter
            [er ersd] = EDetails(errs{idata,iter});
            %input junk
            x(iter) = iter;
            if (il == 1)
               y(iter) = getfield(er,toplot{itype});
            else
               y(iter) = getfield(ersd,toplot{itype});
            end
         end
         figure(70+idata);
         subplot(2,1,il)
         hold on;
         plot(x,y,[psym{itype},ltype{il}]);
      end
   end
   figure(70+idata);
   subplot(2,1,1);
   set(gca,'YSCALE','log');
   set(gca,'YGRID','on');
   subplot(2,1,2);
   set(gca,'YSCALE','log');
   set(gca,'YGRID','on');
   ylabel('average error');
   xlabel('iteration');
   legend(toplot);
end
end
%% weightPlotGenData
function weightPlotGenData()
ws = [0 0.1 0.25 0.5 0.75 1 5 10 20 30 1e8];
ms = cell(0,0);
mtemp = MSet;
mtemp.addData('datasets/ethanerDat.mat',1:10,1:2:20,1,791);
ms{end+1} = mtemp;
mtemp = MSet;
mtemp.addData('datasets/ethanerDat.mat',11:20,2:2:20,1,791);
ms{end+1} = mtemp;
mtemp = MSet;
mtemp.addData('datasets/ch4rDat.mat',1:10,1:2:20,1,791);
ms{end+1} = mtemp;
mtemp = MSet;
mtemp.addData('datasets/propaner2-orig.mat',1:10,1:2:20,1,791);
ms{end+1} = mtemp;
mtemp = MSet;
mtemp.addData('datasets/butaner2-orig.mat',1:10,1:2:20,1,791);
ms{end+1} = mtemp;
mtemp = MSet;
mtemp.addData('datasets/tbutaner2-orig.mat',1:10,1:2:20,1,791);
ms{end+1} = mtemp;


% errs{dataset, w, err/sd}
errs = cell(0,0);
iterSig = [];
iter = 0;
iw = 0;
for wt = ws
   iw = iw+1;
   load(['C:\matdl\yaron\dec12e\w',num2str(wt), ...
      '\hybridslater1\ethanerDat\all-3.mat']);
   ntest = length(monitor.etest);
   errTest = zeros(ntest,1);
   for i = 1:ntest
      errTest(i) = norm(monitor.etest{i});
   end
   [emin,imin] = min(errTest);
   iterMin = monitor.accepted(imin);
   pars = monitor.param{iterMin};
   disp(['wt = ',num2str(wt),' min is ',num2str(imin),' of ', ...
      num2str(length(errTest))]);
   % Get the parameters into the factory
   fm1=f1;
   fm1.setPars(pars);
   for iset = 1:length(ms)
        fm2=fact.makeFitme(ms{iset});
        saveWeights = fm2.operWeights;
        fm2.operWeights = [];
        [eout plotnum etype modelnum envnum] = fm2.err(fm2.getPars);
        fm2.operWeights = saveWeights;
        tsave.err = eout;
        tsave.plotnum = plotnum;
        tsave.etype = etype;
        tsave.modelnum = modelnum;
        tsave.envnum = envnum;
        errs{iset,iw} = tsave;
   end
end
save(['c:\matdl\yaron\dec12e\wplot.mat'],'ws','errs');
end
%% weightPlotGenData2
function weightPlotGenData2()
wt = 5;
ms = cell(0,0);
mtemp = MSet;
mtemp.addData('datasets/ethanerDat.mat',1:10,1:2:20,1,791);
ms{end+1} = mtemp;
mtemp = MSet;
mtemp.addData('datasets/ethanerDat.mat',11:20,2:2:20,1,791);
ms{end+1} = mtemp;
mtemp = MSet;
mtemp.addData('datasets/ch4rDat.mat',1:10,1:2:20,1,791);
ms{end+1} = mtemp;
mtemp = MSet;
mtemp.addData('datasets/propaner2-orig.mat',1:10,1:2:20,1,791);
ms{end+1} = mtemp;
mtemp = MSet;
mtemp.addData('datasets/butaner2-orig.mat',1:10,1:2:20,1,791);
ms{end+1} = mtemp;
mtemp = MSet;
mtemp.addData('datasets/tbutaner2-orig.mat',1:10,1:2:20,1,791);
ms{end+1} = mtemp;


% errs{dataset, ifile}
errs = cell(0,0);
iterSig = [];
iter = 0;
iw = 0;
files = {'start.mat' 'start.mat' 'all-3.mat'}; 
for ifile = 1:3
   load(['C:\Users\yaron\Dropbox\MSQCdata\dec12e\w',num2str(wt), ...
      '\hybridslater1\ethanerDat\',files{ifile}]);
   ntest = length(monitor.etest);
   errTest = zeros(ntest,1);
   for i = 1:ntest
      errTest(i) = norm(monitor.etest{i});
   end
   [emin,imin] = min(errTest);
   iterMin = monitor.accepted(imin);
   pars = monitor.param{iterMin};
   disp(['file = ',files{ifile},' min is ',num2str(imin),' of ', ...
      num2str(length(errTest))]);
   % Get the parameters into the factory
   fm1=f1;
   if (ifile == 1)
      fm1.setPars(zeros(size(pars)));
   else
      fm1.setPars(pars);
   end
   for iset = 1:length(ms)
      fm2=fact.makeFitme(ms{iset});
      saveWeights = fm2.operWeights;
      fm2.operWeights = [];
      [eout plotnum etype modelnum envnum] = fm2.err(fm2.getPars);
      fm2.operWeights = saveWeights;
      tsave.err = eout;
      tsave.plotnum = plotnum;
      tsave.etype = etype;
      tsave.modelnum = modelnum;
      tsave.envnum = envnum;
      errs{iset,ifile} = tsave;
   end
end
save(['c:\matdl\yaron\dec12e\wplot2.mat'],'errs');
end
%% weightPlot1
function weightPlot1()
load('c:\matdl\yaron\dec12e\wplot.mat');
close all;
toplot = {'ke','H','C','e2','etot'};
psym = {'ko','c^','b^','ks','ro'};
ltype = {'-','--'};

ws(end) = 70;
% errs{dataset, iter, err/sd}
for idata = 1:size(errs,1) % data set
   for itype = 1:length(toplot) % data type
      for il = 1:2 % err or standard deviation
         nw = size(errs,2);
         x = zeros(nw,1);
         y = zeros(nw,1);
         for iw = 1:nw
            [er ersd] = EDetails(errs{idata,iw});
            x(iw) = ws(iw);
            if (il == 1)
               y(iw) = getfield(er,toplot{itype});
            else
               y(iw) = getfield(ersd,toplot{itype});
            end
         end
         figure(70+idata);
         subplot(2,1,il)
         hold on;
         plot(x,y,[psym{itype},ltype{il}]);
      end
   end
   figure(70+idata);
   subplot(2,1,1);
   set(gca,'YSCALE','log');
   set(gca,'YGRID','on');
   set(gca,'XSCALE','log');
   subplot(2,1,2);
   set(gca,'YSCALE','log');
   set(gca,'YGRID','on');
   set(gca,'XSCALE','log');
   ylabel('average error');
   xlabel('iteration');
   legend(toplot);
end

end
%% weightPlot2
function weightPlot2()
load('c:\matdl\yaron\dec12e\wplot.mat');
close all;
toplot = {'ke','H','C','e2','etot'};
psym = {'ko','c^','b^','ms','ro'};
ltype = {'-',':'};

ws(1) = 0.01;
ws(end) = 1000;
% errs{dataset, iter, err/sd}
for idata = 1:2 % data set
   for itype = 1:length(toplot) % data type
      il = 1;
      nw = size(errs,2);
      x = zeros(nw,1);
      y = zeros(nw,1);
      for iw = 1:nw
         [er ersd] = EDetails2(errs{idata,iw});
         x(iw) = ws(iw);
         y(iw) = getfield(er,toplot{itype});
      end
      figure(55);
      hold on;
      plot(x,y,[psym{itype},ltype{idata}]);
   end
end
set(gca,'YSCALE','log');
set(gca,'YTick',[0.5 1 2 5 10 20 50 100 200 500 1000 2000 5000]);
set(gca,'XTick',ws);
%set(gca,'YGRID','on');
set(gca,'XSCALE','log');
set(gca,'YLim',[0.5 12000]);
ylabel('average error (kcal/mol)');
xlabel('weight of Etot');
legend(toplot);

%printFig(55,'weightOper');

end
%% weightPlot3
function weightPlot3()
load('c:\matdl\yaron\dec12e\wplot.mat');
ws(1) = 0.05;
ws(end) = 100;
close all;
ethTrain = 1; ethTest = 2; meth = 3; prop = 4; nbut = 5; tbut = 6;
toplot = {'etot','oper'};
dnames = {'Ethane (train)','Ethane (test)','Methane',...
   'Propane','nButane','tButane'};
pcol = {'k','b','g','m','r','c'};
lineType = {'o-','o-'};
% errs{dataset, iw}
close all;
erTemp = cell(size(errs));
er2 = cell(0,0);
niter = size(errs,2);
for iw = 1:length(ws)
   for idata = 1:tbut
      disp(['idata iter ',num2str(idata),' ',num2str(ws(iw))]);
      er2{idata,iw} = EDetails(errs{idata,iw},1,0);
   end
end
leg1 = cell(0,0);
leg2 = cell(0,0);
leg3 = cell(0,0);
for itype = 1:2 % etot oper
   for idata = 1:tbut
      x = zeros(niter,1);
      y = zeros(niter,1);
      for iw = 1:length(ws);
         x(iw) = ws(iw);
         y(iw) = getfield(er2{idata,iw},toplot{itype});
      end
      figure(50+itype);
      hold on;
      plot(x,y,[pcol{idata},lineType{itype}]...
         ,'LineWidth',2);
      leg1{end+1} = dnames{idata};
   end
end
%set(gca,'YGRID','on');
figure(51);
title('Etot')
legend(leg1);
ylabel('RMS error (kcal/mol)');
xlabel('Etot weight');
set(gca,'YSCALE','log');
set(gca,'YTick',[1 2 3 4 5 6 8 10 15]);
set(gca,'XSCALE','log');
set(gca,'XTick',[0.1 0.25 0.5 1 5 10 20 30]);
figure(52);
title('oper')
legend(leg1);
ylabel('RMS error (kcal/mol)');
xlabel('Etot weight');
set(gca,'YSCALE','log');
set(gca,'YTick',[3 5 10 30 50 100 300 500 1000 3000 5000]);
set(gca,'XSCALE','log');
set(gca,'XTick',[0.1 0.25 0.5 1 5 10 20 30]);

printFig(51,'weightEtot');
printFig(52,'weightoper');

end

%% printBigTable()
function printBigTable()

load('c:\matdl\yaron\dec12e\wplot.mat');



end

%% testEDetails
function testEDetails()

load('C:\matdl\yaron\dec12e\bigplot.mat');
er = errs{6,300};
[res0 resSD0] = EDetails(er,0,1);
[res resSD] = EDetails(er,1,1);
a=1;
end

%% EDetails
function [res resSD] = EDetails(errIn, subtractAvg, ofile)
if (nargin < 2)
   ofile = 1;
end
err = errIn.err;
etype = errIn.etype;

err =err*627.509;
eke = err(etype==1);
eH  = err(etype==11);
eC  = err(etype==16);
e2  = err(etype==2);
etot = err(etype==3);
if (subtractAvg)
   eke = eke - mean(eke);
   eH  = eH  - mean(eH);
   eC  = eC  - mean(eC);
   e2  = e2  - mean(e2);
   etot = etot -mean(etot);
end
eall = [eke eH eC e2 etot];
eoper = [eke eH eC e2];
x.ke = norm(eke)/sqrt(length(eke));
x.H = norm(eH)/sqrt(length(eH));
x.C = norm(eC)/sqrt(length(eC));
x.e2 = norm(e2)/sqrt(length(e2));
x.etot = norm(etot)/sqrt(length(etot));
x.oper = norm(eoper)/sqrt(length(eoper));
x.all = norm(eall)/sqrt(length(eall));
x.err = norm(err)/sqrt(length(err));
res = x;
x2.ke = std(eke); x2.H = std(eH); x2.C = std(eC); 
x2.e2=std(e2); x2.etot = std(etot); x2.oper = std(eoper);
x2.all = std(eall); x2.err = std(err);
resSD = x2;
if (ofile ~= 0)
fprintf(ofile,...
   'x  all %5.3f ke %5.3f H %5.3f C %5.3f E2 %5.3f tot %5.3f oper %5.3f err %5.3f\n',...
   x.all, x.ke, x.H, x.C, x.e2, x.etot, x.oper, x.err);
fprintf(ofile,...
   'x2 all %5.3f ke %5.3f H %5.3f C %5.3f E2 %5.3f tot %5.3f oper %5.3f err %5.3f\n',...
   x2.all, x2.ke, x2.H, x2.C, x2.e2, x2.etot, x2.oper, x2.err);
end
end

%% EDetails2: contribution by operator type (H summed)
function [res resSD] = EDetails2(errIn, ofile)
if (nargin < 2)
   ofile = 1;
end
%[eout plotnum etype modelnum envnum]
err = errIn.err;
etype = errIn.etype;
modelnum = errIn.modelnum;
envnum = errIn.envnum;

err =err*627.509;
mods = unique(modelnum);
envs = unique(envnum);
eke = [];
eH = [];
eC = [];
e2 = [];
etot = [];
for imod = mods
   for ienv = envs
      emol = err((modelnum == imod) & (envnum == ienv) & (etype~=3));
      etot1 = sum(emol);
      eka = sum(err((modelnum == imod) & (envnum == ienv) & (etype==1)));
      eHa = sum(err((modelnum == imod) & (envnum == ienv) & (etype==11)));
      eCa = sum(err((modelnum == imod) & (envnum == ienv) & (etype==16)));
      e2a = sum(err((modelnum == imod) & (envnum == ienv) & (etype==2)));
      eta = sum(err((modelnum == imod) & (envnum == ienv) & (etype==3)));
      if ((eta-etot1) > 1e-10)
         disp(['etot error ',num2str(etot1),' ',num2str(eta),' ', ...
            num2str(eta-etot1)]);
      end
      eke(end+1) = eka;
      eH(end+1) = eHa;
      eC(end+1) = eCa;
      e2(end+1) = e2a;
      etot(end+1) = eta;
   end
end
x.ke = norm(eke)/sqrt(length(eke));
x.H = norm(eH)/sqrt(length(eH));
x.C = norm(eC)/sqrt(length(eC));
x.e2 = norm(e2)/sqrt(length(e2));
x.etot = norm(etot)/sqrt(length(etot));
x.all = norm(err)/sqrt(length(err));
x2.ke = std(eke); x2.H = std(eH); x2.C = std(eC); 
x2.e2=std(e2); x2.etot = std(etot);
if (ofile ~= 0)
   fprintf(ofile,'ke %5.3f H %5.3f C %5.3f E2 %5.3f tot %5.3f\n',...
      x.ke, x.H, x.C, x.e2, x.etot);
   fprintf(ofile,'ke %5.3f H %5.3f C %5.3f E2 %5.3f tot %5.3f\n',...
      x2.ke, x2.H, x2.C, x2.e2, x2.etot);
end
res = x;
resSD = x2;
end

