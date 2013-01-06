function factoryUse()
%diffEnv;
%bigPlot1(); % each molecule gets own window
bigPlot2(); % ethane train and test in one window

end
%% bigPlot2
function bigPlot2()
dataroot = ...
   'C:\Users\yaron\Documents\My Dropbox\MSQCdata\dec12e\iter100\w1\hybridslater1\ethanerDat';
load([dataroot, '\bigplot.mat']);
ethTrain = 1; ethTest = 2; meth = 3; prop = 4; nbut = 5; tbut = 6;
nC = [2,2,1,3,4,4];
nH = [6,6,4,8,10];
toplot = {'ke','H','C','e2','etot','all'};
pcol = {'g','c','b','m','r','k'};
lineType = {'-',':'};
% errs{dataset, iter, err/sd}
figure(70)
clf
for idata = [ethTrain ethTest]
   for itype = 1:length(toplot) % data type
      for il = 1 % err or standard deviation
         niter = size(errs,2);
         x = zeros(niter,1);
         y = zeros(niter,1);
         for iter = 1:niter
            er1 = errs{idata,iter,il};
            er = getTotalError(er1{1},nC(idata),nH(idata));
            x(iter) = iter;
            y(iter) = getfield(er,toplot{itype});
         end
         hold on;
         plot(x,y,[pcol{itype},lineType{idata}]);
      end
   end
end
set(gca,'YSCALE','log');
set(gca,'YTick',[2 4 6 8 10 20 40 60 80 100 200 400 600 800]);
for isig = 2:length(iterSig)
   xL = iterSig(isig);
   yL = get(gca,'YLim');
   line([xL xL],yL,'color','k');
end
%set(gca,'YGRID','on');
ylabel('RMS error (kcal/mol)');
xlabel('iteration');
legend(toplot);
end

function errOut = getTotalError(errIn,nC,nH)
errOut.ke = errIn.ke;
errOut.H = nH * errIn.H;
errOut.C = nC * errIn.C;
errOut.e2 = errIn.e2;
errOut.etot = errIn.etot;
errOut.all = errOut.ke + errOut.H + errOut.C + errOut.e2 + errOut.etot;
end

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

function BigPlotsGenData()
dataroot = ...
   'C:\matdl\yaron\dec12e\iter100\w1\hybridslater1\ethanerDat';
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
mtemp.addData('datasets/butaner-orig.mat',1:10,1:2:20,1,791);
ms{end+1} = mtemp;
mtemp = MSet;
mtemp.addData('datasets/tbutaner-orig.mat',1:10,1:2:20,1,791);
ms{end+1} = mtemp;

% errs{dataset, iter, err/sd}
errs = cell(0,0,0);
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
        [a b] = fm2.printEDetails;
        errs{iset,iter,1} = a;
        errs{iset,iter,2} = b;
     end
  end
end
save([dataroot, '\bigplot.mat'],'errs','iterSig'); %,'emeth','smeth');
end
function bigPlot1()
% each molecule gets its own window
dataroot = ...
   'C:\Users\yaron\Documents\My Dropbox\MSQCdata\dec12e\iter100\w1\hybridslater1\ethanerDat';
load([dataroot, '\bigplot.mat']); %,'emeth','smeth');
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
            er = errs{idata,iter,il};
            %input junk
            x(iter) = iter;
            y(iter) = getfield(er{1},toplot{itype});
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
mtemp.addData('datasets/butaner-orig.mat',1:10,1:2:20,1,791);
ms{end+1} = mtemp;
mtemp = MSet;
mtemp.addData('datasets/tbutaner-orig.mat',1:10,1:2:20,1,791);
ms{end+1} = mtemp;

fitmes = cell(0,0);
fitmes{end+1} = f1;
fitmes{end+1} = ftest;
fitmes{end+1} = [];
fitmes{end+1} = [];
fitmes{end+1} = [];
% errs{dataset, w, err/sd}
errs = cell(0,0,0);
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
      if (iset == 2)
         fm2.setWeights(wt,0);
         errCheck = fm2.err(pars);
         disp(['test error should agree: ',...
            num2str(max(abs(errCheck-monitor.etest{imin})))]);
         fm2.operWeights = [];
      end
      [a b] = fm2.printEDetails;
      errs{iset,iw,1} = a;
      errs{iset,iw,2} = b;
   end
end
save(['c:\matdl\yaron\dec12e\wplot.mat'],'ws','errs');
end
function weightPlot1()
load('c:\matdl\yaron\dec12e\wplot.mat');
close all;
toplot = {'ke','H','e2','etot'};%{'ke','H','C','e2','etot'};
psym = {'co','bo','k^','ro'};%{'co','bo','b^','k^','ro'};
ltype = {'-','--'};

ws(end) = 1000
% errs{dataset, iter, err/sd}
for idata = 1:size(errs,1) % data set
   for itype = 1:length(toplot) % data type
      for il = 1:2 % err or standard deviation
         nw = size(errs,2);
         x = zeros(nw,1);
         y = zeros(nw,1);
         ii = 0;
         for iw = [1 3 2 4:10]
            ii = ii+1;
            er = errs{idata,iw,il};
            %input junk
            x(ii) = ws(iw);
            y(ii) = getfield(er{1},toplot{itype});
         end
         figure(70+idata);
         subplot(2,1,il)
         hold on;
         plot(x,y,[psym{itype},ltype{il}]);
      end
   end
   figure(70+idata);
   subplot(2,1,1);
   %set(gca,'YSCALE','log');
   set(gca,'YGRID','on');
   set(gca,'XSCALE','log');
   subplot(2,1,2);
   %set(gca,'YSCALE','log');
   set(gca,'YGRID','on');
   set(gca,'XSCALE','log');
   ylabel('average error');
   xlabel('iteration');
   legend(toplot);
end

end