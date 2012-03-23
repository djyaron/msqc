%% 
clear classes;
reload = 1;
modelType = 'noCharge';
handFit = 0;
doFit = 0;
plotResults = 1;

useStart = 1;
% fit for modelType = 'nocharge'
pstart = [1.43099 -5.37587 1.07249 23.5346 -3.07946 13.122 -10.6423 -3.30887];
%
if (reload)
   disp('reloading data');
   load('ethane4/ethaneDat.mat');
end
if strcmp(modelType,'noCharge')
   disp('setting up model with charge independent parameters');
   mixKEh = Mixer(0,1);
   mixKEcs = Mixer(0,1);
   mixKEcp = Mixer(0,1);
   mixKEchs = Mixer(0,1);
   mixKEchp = Mixer(0,1);
   mixKEccss = Mixer(0,1);
   mixKEccpp = Mixer(0,1);
   mixKEccsp = Mixer(0,1);
   m = cell(1,7);
   ic = 0;
   for ipar = 1:7
      t1 = LL{ipar,1}.EKE;
      lke(ic+1:ic+size(t1,2)) = t1;
      hke(ic+1:ic+size(t1,2)) = HL{ipar,1}.EKE;
      
      m{ipar} = Model3(LL{ipar,1},LL{ipar,2},LL{ipar,3});
      m{ipar}.addKEmodDiag(1,1,mixKEh);
      m{ipar}.addKEmodDiag(6,1,mixKEcs);
      m{ipar}.addKEmodDiag(6,2,mixKEcp);
      m{ipar}.addKEmodBonded(1,6,1,1, mixKEchs);
      m{ipar}.addKEmodBonded(1,6,1,2, mixKEchp);
      m{ipar}.addKEmodBonded(6,6,1,1, mixKEccss);
      m{ipar}.addKEmodBonded(6,6,2,2, mixKEccpp);
      m{ipar}.addKEmodBonded(6,6,1,2, mixKEccsp);
      ic = ic+size(t1,2);
   end
elseif (strcmp(modelType,'chargeDep'))
   disp('Building charge dependent model');
   error('ChargeDep NOT YET IMPLIMENTED');
   for ipar = 1:7
      for ienv = 0:m{ipar}.nenv
         m{ipar}.charges(:,ienv+1) = m{ipar}.charges(:,ienv+1) - m{1}.charges(:,1);
      end
   end
   
end
if (handFit)
   while 0
      p(1) = input('H constant ');
      p(2) = input('H slope    ');
      p(3) = input('C constant ');
      p(4) = input('C slope    ');
      p = [1.2353   -4.7202    0.4195   -7.0583];
      ic = 0;
      for ipar = 1:7
         disp(['starting calc on ipar ',num2str(ipar)]);
         m{ipar}.setPar(p);
         m{ipar}.solveHF;
         t1 = m{ipar}.EKE;
         mke(ic+1:ic+size(t1,2)) = t1;
         ic = ic + size(t1,2);
      end
      figure(100);
      hold off;
      plot(lke,hke,'r.');
      hold on;
      plot(lke,lke,'k.');
      plot(lke,mke,'b.');
      figure(200);
      plot(mke,hke,'g.');
   end
end
if (doFit)
   disp('Starting to do parameter fitting');
   f1 = Fitme;
   for ipar = 1:7
      f1.addFrag(m{ipar},HL{ipar,1});
   end
   f1.exactDensity = 1;
   f1.includeEN = zeros(1,m{1}.natom);
   
   nfitpar = m{1}.npar;
   start = zeros(1,nfitpar);
   if (useStart)
      start = pstart;
   end
   % code from using matlabs optimization toolbox
   %limits = 3 * ones(1,nfitpar);
   %options = optimset('DiffMinChange',1.0e-5);
   %pt = lsqnonlin(@f1.err, start,-limits,limits,options);
   options = LMFnlsq;
   options.Display =1;
   options.FunTol = 1.0e-5;
   options.XTol = 1.0e-4;
   [pfit, Ssq, CNT, Res, XY] = LMFnlsq(@f1.err,start',options);
end

if (plotResults)
   disp('Starting to do plots');
   if (doFit)
      pt = pfit;
   else
      pt = pstart;
   end
   ic = 0;
   for ipar = 1:7
      disp(['starting calc on ipar ',num2str(ipar)]);
      m{ipar}.setPar(pt);
      m{ipar}.solveHF;
      t1 = m{ipar}.EKE;
      mke(ic+1:ic+size(t1,2)) = t1;
      ic = ic + size(t1,2);
   end
   figure(100);
   hold off;
   plot(lke,hke,'r.');
   hold on;
   plot(lke,lke,'k.');
   plot(lke,mke,'b.');
   figure(200);
   plot(mke,hke,'g.');
end
%%
ar = 1:700;
figure(100);
hold off;
plot(lke(ar),hke(ar),'r.');
hold on;
plot(lke(ar),lke(ar),'k.');
plot(lke(ar),mke(ar),'b.');
figure(200);
plot(mke(ar),hke(ar),'g.');
%% atomic charges
ic = 0;
% HL and LL charges on the atoms
chl=zeros(8,700);
cll=zeros(8,700);
% charges on teh carbons, with hydrogens summed into carbon
cshl=zeros(2,700);
csll=zeros(2,700);
for ipar = 1:7
   hl = HL{ipar,1};
   ll = LL{ipar,1};
   t1 = zeros(8,nenv);
   t2 = zeros(8,nenv);
   t1s = zeros(2,nenv);
   t2s = zeros(2,nenv);
   for ienv = 1:hl.nenv
      mh = hl.mcharge(ienv);
      ml = ll.mcharge(ienv);
      t1(:,ienv) = mh(1,:);
      t2(:,ienv) = ml(1,:);
      t1s(1,ienv) = sum(mh(1,[1,3,4,5]));
      t1s(2,ienv) = sum(mh(1,[2,6,7,8]));
      t2s(1,ienv) = sum(ml(1,[1,3,4,5]));
      t2s(2,ienv) = sum(ml(1,[2,6,7,8]));
   end
   rr = ic+1:ic+size(t1,2);
   chl(:,rr) = t1;
   cll(:,rr) = t2;
   cshl(:,rr) = t1s;
   csll(:,rr) = t2s;
   ic = ic + size(t1,2);
end
%%
ar = 1:700;
figure(500);
hold off;
plot(csll(1,ar),csll(2,ar),'b.');
hold on;
plot(cshl(1,ar),cshl(2,ar),'r.');
%%
ar = 1:700;
figure(501);
hold off;
plot(csll(1,ar),lke(ar),'b.');
figure(502);
hold off;
plot(csll(1,ar),hke(ar),'r.');

%%
hold on;
plot(lke(ar),lke(ar),'k.');
plot(lke(ar),mke(ar),'b.');
figure(200);
plot(mke(ar),hke(ar),'g.');

