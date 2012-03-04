%% ethaneDat.mat is from top two sections of ethane-anal2.m
clear classes;
load('ethane4/ethaneDat.mat');
%% Correlation plots
for ipar = 1
   %set(gca,'NextPlot','replacechildren');
   %ic = 0;
   x = LL{ipar,1}.EKE;
   figure(100+ipar)
   subplot(3,3,1);
   plot(x,LL{ipar,1}.EKE,'k.');
   hold on;
   plot(x,LL{ipar,2}.EKE,'b.');
   plot(x,LL{ipar,3}.EKE,'y.');
   plot(x,HL{ipar,1}.EKE,'r.');
   plot(x,HL{ipar,2}.EKE,'r.');
   plot(x,HL{ipar,3}.EKE,'r.');
   title('KE')
   %legend('KE','LL2','LL3','HL1','HL2','HL3');
   %ic = ic+1;
   %F(ic) = getframe();
   
   for iatom = 1:LL{ipar,1}.natom
      x = LL{ipar,1}.Een(iatom);
      subplot(3,3,iatom+1);
      hold off;
      plot(x,LL{ipar,1}.Een(iatom),'k.');
      hold on;
      plot(x,LL{ipar,2}.Een(iatom),'b.');
      plot(x,LL{ipar,3}.Een(iatom),'y.');
      plot(x,HL{ipar,1}.Een(iatom),'r.');
      plot(x,HL{ipar,2}.Een(iatom),'r.');
      plot(x,HL{ipar,3}.Een(iatom),'r.');
      title(['EN ',num2str(iatom)])
      %legend(['EN ',num2str(iatom)],'LL2','LL3','HL1','HL2','HL3');
      %   ic = ic+1;
      %   F(ic) = getframe();
   end
end
%movie(F,2,0.3)
%% Replacing operators
% What happens if we just replace the KE operator with the narrow and
% diffuse
ipar = 1;
m = Model3(LL{ipar,1},LL{ipar,2},LL{ipar,3});

mixKE = Mixer(0,1); % linear interpolation of entire KE operator
%m.addKEmodDiag(1,1,mixKE);
%m.addKEmodDiag(6,[1,2],mixKE);
m.addKEmodBonded(1,6,[1],[1 2],mixKE);
m.addKEmodBonded(6,6,[1 2],[1 2],mixKE);

figure(200)
x = LL{ipar,1}.EKE;
plot(x,x,'k.');
hold on;
plot(x,LL{ipar,2}.EKE,'b.');
plot(x,LL{ipar,3}.EKE,'y.');
sym = {'r.','g.','c.'};
ic = 0;

for p = -1 %:1:1
   m.setPar(p);
   disp(['solving HF for ',num2str(p)]);
   m.solveHF;
   ic = ic+1;
   plot(x,m.EKE,sym{ic});
end

hold on;
plot(x,HL{ipar,3}.EKE,'rx')

%% Checking out charge dependence
ipar = 1;
m = Model3(LL{ipar,1},LL{ipar,2},LL{ipar,3});
% replace charges with difference from having no environment
for ienv = 0:m.nenv
   m.charges(:,ienv+1) = m.charges(:,ienv+1) - m.charges(:,1);
end

mixKE = Mixer([0 0],2); % linear interpolation of entire KE operator
mixKE.fixed = [1 0];
m.addKEmodDiag(1,1,mixKE);
m.addKEmodDiag(6,[1,2],mixKE);

figure(200)
x = LL{ipar,1}.EKE;
plot(x,x,'k.');
hold on;
plot(x,LL{ipar,2}.EKE,'b.');
plot(x,LL{ipar,3}.EKE,'y.');
sym = {'r.','g.','c.'};
ic = 0;

for p = -10:10:10
   m.setPar(p);
   disp(['solving HF for ',num2str(p)]);
   m.solveHF;
   ic = ic+1;
   hold on;
   plot(x,m.EKE,sym{ic});
end

hold on;
plot(x,HL{ipar,3}.EKE,'rx')


%% hand fit charge dep mixers
clear classes;
load('ethane4/ethaneDat.mat');

%%
%mixKEh = Mixer([0,0],2); 
%mixKEc = Mixer([0,0],2);
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
   %m{ipar}.addKEmodDiag(1,1,mixKEh);
   %m{ipar}.addKEmodDiag(6,[1,2],mixKEc);
   %m{ipar}.addKEmodBonded(1,6,[1],[1 2], mixKEch);
   %m{ipar}.addKEmodBonded(6,6,[1 2],[1 2], mixKEcc);
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

% adjust charges
for ipar = 1:7
   for ienv = 0:m{ipar}.nenv
      m{ipar}.charges(:,ienv+1) = m{ipar}.charges(:,ienv+1) - m{1}.charges(:,1);
   end
end
% hand fit
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

% general fit

f1 = Fitme;
for ipar = 1:7
   f1.addFrag(m{ipar},HL{ipar,1});
end
f1.exactDensity = 1;
f1.includeEN = zeros(1,m{1}.natom);

nfitpar = m{1}.npar;
start = zeros(1,nfitpar);
%start(1) = 0.0;
%start = [1.4962 -5 0 0]; % Fit of all four gives [1.1734  -3 1.1326 2.6529]
limits = 3 * ones(1,nfitpar);
options = optimset('DiffMinChange',1.0e-5);
%pt = lsqnonlin(@f1.err, start,-limits,limits,options);
options = LMFnlsq;
options.Display =1;
options.FunTol = 1.0e-5;
options.XTol = 1.0e-4;
[pt, Ssq, CNT, Res, XY] = LMFnlsq(@f1.err,start',options);
%%
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

%%
% pt = [1.1734  -3 1.1326 2.6529];
m.setPar(pt);
m.solveHF;
mke = m.EKE;
figure(100);
hold off;
plot(lke,hke,'r.');
hold on;
plot(lke,lke,'k.');
plot(lke,mke,'b.');
figure(200);
plot(mke,hke,'g.');
%% Plots versus various parameters
clear classes;
load('ethane4/ethaneDat10env.mat');

lke = LL{ipar,1}.EKE;
hke = HL{ipar,1}.EKE;

ipar = 1;
m = Model3(LL{ipar,1},LL{ipar,2},LL{ipar,3});

%mix1 = cell(0,0);

m.addKEmodDiag(6,1, Mixer(0,1));
m.addKEmodDiag(6,2, Mixer(0,1));
m.addKEmodDiag(1,1, Mixer(0,1));

f1 = Fitme;
f1.addFrag(m,LL{ipar,1});
f1.exactDensity = 1;
f1.includeKE = 1;
f1.includeEN = zeros(1,m.natom);

% ic=0;
% for p = -1:0.25:1
%    ic = ic+1;
%    e1 = f1.errDiffs(p);
%    x1(ic) = p;
%    y1(ic) = norm(e1);
% end
%%
nfitpar = m.npar;
start = zeros(1,nfitpar);
limits = 3 * ones(1,nfitpar);
options = optimset('DiffMinChange',1.0e-5);
pt = lsqnonlin(@f1.errDiffs, start,-limits,limits,options);

%% loop through parameters, adjusting each separately
clear classes;
load('ethane4/ethaneDat10env.mat');
% KE diagonal    1: H  2: Cs  3: Cp
% Hen diagonal   4: H  5: Cs  6: Cp
% KE bonding     7: H-Cs  8: H-Cp  9: Cs-Cs 10: Cs-Cp 11: Cp-Cp
% Hen bonding   12: H-Cs 13: H-Cp 14: Cs-Cs 15: Cs-Cp 16: Cp-Cp
% KE charge     17: H    18:Cs   19:Cp

ipar = 1;
m = Model3(LL{ipar,1},LL{ipar,2},LL{ipar,3});

mix1 = cell(0,0);
mix1{1,1} = Mixer(0,1);
%mix1{1,2} = Mixer(0,1);

%mixKEC = Mixer(0,1); 
%m.addKEmodDiag(6,1,mixKEC);
%m.addKEmodDiag(6,2,mixKEC);
%mixKEH = Mixer(0,1); 
%m.addKEmodDiag(1,1,mixKEH);
%m.addKEmodBonded(1,6,1,1,mixKE);
%m.addKEmodBonded(1,6,1,2,mixKE);
%m.addKEmodBonded(6,6,1,1,mixKE);
%m.addKEmodBonded(6,6,1,2,mixKE);
%m.addKEmodBonded(6,6,2,2,mixKE);

f1 = Fitme;
f1.addFrag(m,HL{ipar,1});
f1.exactDensity = 1;
f1.includeKE = 1;
f1.includeEN = zeros(1,m.natom);

% nmix = size(m.mixers,2);
% disp('1');
% for im = 1:nmix
%    m.mixers{1,im}.fixed = 1;
% end
% disp('2');

nmix = size(mix1,2);
for im = 1:npar1
   switch im 
      case 1
         m.addKEmodDiag(6,2,mix1{1,im});
         
      case 2
         m.addKEmodDiag(6,1,mix1{1,im});
   end
   ic = 0;
   for p = -2:0.5:2
      ic = ic+1;
      e1 = f1.err(p);
      x1(ic) = p;
      y1(ic) = norm(e1);
   end
  figure(102);
  hold off;
  plot(x1,y1,'bx-');
  title(['im = ',num2str(im)]);

  fminbnd(@f1.err,-2,2,optimset('TolX',0.05));
  mix1{1,im}.fixed = 1;
end



