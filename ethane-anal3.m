%% ethaneDat.mat is from top two sections of ethane-anal2.m
clear classes;
load('ethane4/ethaneDat.mat');
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


