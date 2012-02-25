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
m.addKEmodDiag(1,1,mixKE);
m.addKEmodDiag(6,[1,2],mixKE);
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

for p = -1:1:1
   m.setPar(p);
   disp(['solving HF for ',num2str(p)]);
   m.solveHF;
   ic = ic+1;
   plot(x,m.EKE,sym{ic});
end


%%

ipar = 7;
disp(['starting fit on geometry ',num2str(ipar)]);
m = Model3(LL{ipar,1},LL{ipar,2},LL{ipar,3});

mixKE = Mixer(0,1); 
m.addKEmodDiag(1,1,mixKE);
m.addKEmodDiag(6,[1,2],mixKE);
m.addKEmodBonded(1,6,[1],[1 2],mixKE);
m.addKEmodBonded(6,6,[1 2],[1 2],mixKE);

f1 = Fitme;
f1.addFrag(m,HL{ipar,3});
f1.exactDensity = 1;
f1.includeEN = zeros(1,m.natom);

% ic = 0;
% for p = -1:0.25:1
%    ic = ic+1;
%    e1{ic} = f1.errDiffs(p);
%    x1(ic) = p;
%    y1(ic) = norm(e1{ic});
% end

nfitpar = m.npar;
start = zeros(1,nfitpar);
start(1) = 0.0;
limits = 3 * ones(1,nfitpar);
options = optimset('DiffMinChange',1.0e-5);
pt = lsqnonlin(@f1.errDiffs, start,-limits,limits,options);
err = f1.errDiffs(pt);
corrPlot(f1,pt, 0, 800+ipar);

%% Plots versus various parameters
clear classes;
load('ethane4/ethaneDat10env.mat');

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



