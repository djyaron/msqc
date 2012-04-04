function fitme = makeFitme()
% location of the data
dataDir = './datasets';
% The high level theory to fit to (second index on the HL's read from data
% directory
nhl = 1;
% plot results as fitting occurs
doPlot = 1;
% environments to include in the fit
envs = 0:20; % environments to include in fit
geomsH2 = []; %2:7;
geomsCH4 = 1:7; % allowed range is 1:19
geomsEthane = 1:7; % allowed range is 1:7

includeKEmods = 1;
includeENmods = 1;
useDeltaCharges = 1;

% Load data
LL1 = cell(0,0);
HL1 = cell(0,0);
ic = 0;
plotNumber = [];
for i = geomsCH4
   load([dataDir,'/ch4Dat.mat']);
   ic = ic+1;
   plotNumber(1,ic) = 801;
   for j = 1:size(LL,2)
      LL1{ic,j} = LL{i,j};
   end
   HL1{ic,1} = HL{i,nhl};
end
for i = geomsH2
   load([dataDir,'/h2Dat.mat']);
   ic = ic+1;
   plotNumber(1,ic) = 800;
   for j = 1:size(LL,2)
      LL1{ic,j} = LL{i,j};
   end
   HL1{ic,1} = HL{i,nhl};
end
for i = geomsEthane
   load([dataDir,'/ethaneDat.mat']);
   ic = ic+1;
   plotNumber(1,ic) = 802;
   for j = 1:size(LL,2)
      LL1{ic,j} = LL{i,j};
   end
   HL1{ic,1} = HL{i,nhl};
end
params = 1:ic;
LL = LL1;
HL = HL1;

% build models
%disp('building models');
m = cell(1,size(params,2));
for ipar = params
   m{ipar} = Model3(LL{ipar,1},LL{ipar,2},LL{ipar,3});
end
if (includeKEmods)
   mixKEdiagH = Mixer([0 0],2,'KEdiagH');
   mixKEdiagC = Mixer([0 0],2,'KEdiagC');
   mixKEdiagCp = Mixer([0 0],2,'KEdiagCp');
   mixKEbondHH = Mixer(0,1,'KEbondHH');
   mixKEbondCH  = Mixer(0,1,'KEbondCH');
   mixKEbondCHp  = Mixer(0,1,'KEbondCHp');
   mixKEbondCC  = Mixer(0,1,'KEbondCC');
   for ipar = params
      m{ipar}.addKEmodDiag(1,1,mixKEdiagH);
      m{ipar}.addKEmodDiag(6,1,mixKEdiagC);
      m{ipar}.addKEmodDiag(6,2,mixKEdiagCp);
      m{ipar}.addKEmodBonded(1,1,1,1,mixKEbondHH);
      m{ipar}.addKEmodBonded(1,6,1,1,mixKEbondCH);
      m{ipar}.addKEmodBonded(1,6,1,2,mixKEbondCHp);
      m{ipar}.addKEmodBonded(6,6,[1 2],[1 2],mixKEbondCC);
   end
end
if (includeENmods)
   mixENdiagH = Mixer([0 0],2,'ENdiagH');
   mixENdiagC = Mixer([0 0],2,'ENdiagC');
   mixENdiagCp = Mixer([0 0],2,'ENdiagCp');
   mixENbondHH = Mixer(0,1,'ENbondHH');
   mixENbondCH  = Mixer(0,1,'ENbondCH');
   mixENbondCHp  = Mixer(0,1,'ENbondCHp');
   mixENbondCC  = Mixer(0,1,'ENbondCC');
   for ipar = params
      m{ipar}.addENmodDiag(1,1,mixENdiagH);
      m{ipar}.addENmodDiag(6,1,mixENdiagC);
      m{ipar}.addENmodDiag(6,2,mixENdiagCp);
      m{ipar}.addENmodBonded(1,1,1,1,mixENbondHH);
      m{ipar}.addENmodBonded(1,6,1,1,mixENbondCH);
      m{ipar}.addENmodBonded(1,6,1,2,mixENbondCHp);
      m{ipar}.addENmodBonded(6,6,[1 2],[1 2],mixENbondCC);
   end
end
if (useDeltaCharges)
   for ipar = params
      for ienv = 1:m{ipar}.nenv
         m{ipar}.charges(:,ienv+1) = m{ipar}.charges(:,ienv+1) - m{ipar}.charges(:,1);
      end
      m{ipar}.charges(:,1) = m{ipar}.charges(:,1) - m{ipar}.charges(:,1);
   end
end

fitme = Fitme;
for ipar = params
   fitme.addFrag(m{ipar},HL{ipar,1},plotNumber(ipar));
end
fitme.includeKE = includeKEmods;
fitme.includeEN = includeENmods * ones(1,6);
fitme.setEnvs(envs);

end