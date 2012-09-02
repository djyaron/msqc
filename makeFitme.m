function fitme = makeFitme(varargin)
% Create a fitme object, using data on H2, CH4 and ethane
% Change default parameters by passing list of key, value pairs, e.g.
%   fitme('datadir','./mydir','plot',0);
% keys and default values are:
%  datadir  "./datasets"     directory with *.mat data files
%  nhl      1                high level theory (1..3)
%  plot     1                if >0 fitme will plot results when parameters change
%                            (figure numbers are H2:800 ch4:801 ethane:802
%                             if plot=1, 810/811/812 if plot=2 etc.
%  envs     0:20             range of environments to include in the fits
%  h2       2:7              range of h2 geometries to include (1..7)
%  ch4      []               range of ch4 geometries to include (1..19)
%  ethane   []               range of ethane geometries to include (1..7)
%  ethaner   []              range of random ethane geometries to include (1..50)
%  ethylene []               range of ethylene geometries to include (1..7)
%  kemods   1                include mixers to modify kinetic energy
%  enmods   1                include mixers to modify elec-nuc interaction
%  deltarho 1                base charge dependence on charges induced
%                            by field
%  enstruct  []              structure with enmods
%  enstruct1 []              structure with enmods (1 oper only)
%  enstructh []              structure with enmods hybrid
%  kestruct  []              structure with kemods
%  kestructh  []             structure with kemods hybrid
%  e2struct  []              structure with e2mods
%  testFitme []              Fitme object with test data
%  silent    0               suppress fitme displayed output

% To test parsing of input, use:
%makeFitme('datadir','./mydata','nhl',2,'plot',0,'envs',1:3,'h2',1:3, ...
%   'ch4',1:4,'ethane',1:2,'kemods',0,'enmods',0,'deltarho',1)

% location of the data sets
dataDir = checkForInput(varargin,'datadir','./datasets');
% The high level theory to fit to (second index on the HL's read from data
% directory
nhl = checkForInput(varargin,'nhl',1);
% plot results as fitting occurs
doPlot = checkForInput(varargin,'plot',1);
% environments to include in the fit
envs = checkForInput(varargin,'envs',0:20); 
geomsH2 = checkForInput(varargin,'h2',[]); % allowed range is 1:7
geomsCH4 = checkForInput(varargin,'ch4',[]); % allowed range is 1:19
geomsCH4r = checkForInput(varargin,'ch4r',[]); % allowed range is 1:25
geomsEthane = checkForInput(varargin,'ethane',[]); % allowed range is 1:7
geomsEthaner = checkForInput(varargin,'ethaner',[]); % allowed range is 1:7
geomsEthylene = checkForInput(varargin,'ethylene',[]); % allowed range is 1:7
geomsPropane = checkForInput(varargin,'propane',[]); % allowed range is 1:7
geomsPropene = checkForInput(varargin,'propene',[]); % allowed range is 1:9
geomsCH3F = checkForInput(varargin,'ch3f',[]); % allowed range is 1:19
geomsC2H5F = checkForInput(varargin,'c2h5f',[]); % allowed range is 1:9
includeKEmods = checkForInput(varargin,'kemods',1);
includeENmods = checkForInput(varargin,'enmods',1);
useDeltaCharges = checkForInput(varargin,'deltarho',1);
enstruct = checkForInput(varargin,'enstruct',[]);
enstruct1 = checkForInput(varargin,'enstruct1',[]);
enstructh = checkForInput(varargin,'enstructh',[]);
kestruct = checkForInput(varargin,'kestruct',[]);
kestructh = checkForInput(varargin,'kestructh',[]);
e2struct = checkForInput(varargin,'e2struct',[]);
testFitme = checkForInput(varargin,'testFitme',[]);
silent = checkForInput(varargin,'silent',0);

if (~isempty(enstruct) && ~isempty(enstruct1))
   error('Do not set both enstruct and enstruct1');
end

% Load data
LL1 = cell(0,0);
HL1 = cell(0,0);
ic = 0;
plotNumber = [];
if (~isempty(geomsCH4))
   load([dataDir,'/ch4Dat.mat']);
   for i = geomsCH4
      ic = ic+1;
      plotNumber(1,ic) = 801 + 10 * (doPlot-1);
      for j = 1:size(LL,2)
         LL1{ic,j} = LL{i,j};
      end
      HL1{ic,1} = HL{i,nhl};
   end
end
if (~isempty(geomsCH4r))
   load([dataDir,'/ch4rDat.mat']);
   for i = geomsCH4r
      ic = ic+1;
      plotNumber(1,ic) = 801 + 10 * (doPlot-1);
      for j = 1:size(LL,2)
         LL1{ic,j} = LL{i,j};
      end
      HL1{ic,1} = HL{i,nhl};
   end
end
if (~isempty(geomsH2))
   load([dataDir,'/h2Dat.mat']);
   for i = geomsH2
      ic = ic+1;
      plotNumber(1,ic) = 800 + 10 * (doPlot-1);
      for j = 1:size(LL,2)
         LL1{ic,j} = LL{i,j};
      end
      HL1{ic,1} = HL{i,nhl};
   end
end
if (~isempty(geomsEthane))
   for i = geomsEthane
      load([dataDir,'/ethaneDat.mat']);
      ic = ic+1;
      plotNumber(1,ic) = 802 + 10 * (doPlot-1);
      for j = 1:size(LL,2)
         LL1{ic,j} = LL{i,j};
      end
      HL1{ic,1} = HL{i,nhl};
   end
end
if (~isempty(geomsEthaner))
   for i = geomsEthaner
      load([dataDir,'/ethanerDat.mat']);
      ic = ic+1;
      plotNumber(1,ic) = 802 + 10 * (doPlot-1);
      for j = 1:size(LL,2)
         LL1{ic,j} = LL{i,j};
      end
      HL1{ic,1} = HL{i,nhl};
   end
end
if (~isempty(geomsEthylene))
   for i = geomsEthylene
      load([dataDir,'/ethyleneDat.mat']);
      ic = ic+1;
      plotNumber(1,ic) = 803 + 10 * (doPlot-1);
      for j = 1:size(LL,2)
         LL1{ic,j} = LL{i,j};
      end
      HL1{ic,1} = HL{i,nhl};
   end
end
if (~isempty(geomsPropane))
   for i = geomsPropane
      load([dataDir,'/propaneDat.mat'], 'LL', 'HL');
      ic = ic+1;
      plotNumber(1,ic) = 804 + 10 * (doPlot-1);
      for j = 1:size(LL,2)
         LL1{ic,j} = LL{i,j};
      end
      HL1{ic,1} = HL{i,nhl};
   end
end
if (~isempty(geomsPropene))
   for i = geomsPropene
      load([dataDir,'/propeneDat.mat'], 'LL', 'HL');
      ic = ic+1;
      plotNumber(1,ic) = 805 + 10 * (doPlot-1);
      for j = 1:size(LL,2)
         LL1{ic,j} = LL{i,j};
      end
      HL1{ic,1} = HL{i,nhl};
   end
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
   if (isempty(kestruct) && isempty(kestructh))
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
   elseif (size(kestruct,1) == 1)
      for ipar = params
         m{ipar}.addKEmodDiag(1,1,kestruct.H);
         m{ipar}.addKEmodDiag(6,1,kestruct.Cs);
         m{ipar}.addKEmodDiag(6,2,kestruct.Cp);
         m{ipar}.addKEmodBonded(1,1,1,1,kestruct.HH);
         m{ipar}.addKEmodBonded(1,6,1,1,kestruct.CsH);
         m{ipar}.addKEmodBonded(1,6,1,2,kestruct.CpH);
         m{ipar}.addKEmodBonded(6,6,1,1,kestruct.CsCs);
         m{ipar}.addKEmodBonded(6,6,1,2,kestruct.CsCp);
         m{ipar}.addKEmodBonded(6,6,2,2,kestruct.CpCp);
      end
   elseif (size(kestructh,1) == 1)
      for ipar = params
         m{ipar}.addKEmodDiag(1,1,kestructh.H);
         m{ipar}.addKEmodDiag(6,1,kestructh.Cs);
         m{ipar}.addKEmodDiag(6,2,kestructh.Cp);
         m{ipar}.addKEmodBondedh(1,1,kestructh.HH);
         m{ipar}.addKEmodBondedh(1,6,kestructh.CH);
         m{ipar}.addKEmodBondedh(6,6,kestructh.CCs);
         m{ipar}.addKEmodBondedh(6,6,kestructh.CCp);
      end
   end
end
if (includeENmods)
   if (isempty(enstruct) && isempty(enstruct1) && isempty(enstructh))
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
   elseif (size(enstruct,1) == 1)
      for ipar = params
         m{ipar}.addENmodDiag(1,1,enstruct.H);
         m{ipar}.addENmodDiag(6,1,enstruct.Cs);
         m{ipar}.addENmodDiag(6,2,enstruct.Cp);
         m{ipar}.addENmodBonded(1,1,1,1,enstruct.HH);
         m{ipar}.addENmodBonded(1,6,1,1,enstruct.CsH);
         m{ipar}.addENmodBonded(1,6,1,2,enstruct.CpH);
         m{ipar}.addENmodBonded(6,6,1,1,enstruct.CsCs);
         m{ipar}.addENmodBonded(6,6,1,2,enstruct.CsCp);
         m{ipar}.addENmodBonded(6,6,2,2,enstruct.CpCp);
      end
   elseif (size(enstruct1,1) == 1)
      for ipar = params
         m{ipar}.addENmodDiag(1,1,enstruct1.H);
         m{ipar}.addENmodDiag(6,1,enstruct1.Cs);
         m{ipar}.addENmodDiag(6,2,enstruct1.Cp);
         m{ipar}.addENmodBonded(1,1,1,1,enstruct1.HH);
         m{ipar}.addENmodBonded1(6,1,1,1,enstruct1.CsH);
         m{ipar}.addENmodBonded1(6,1,2,1,enstruct1.CpH);
         m{ipar}.addENmodBonded1(1,6,1,1,enstruct1.HCs);
         m{ipar}.addENmodBonded1(1,6,1,2,enstruct1.HCp);
         m{ipar}.addENmodBonded(6,6,1,1,enstruct1.CsCs);
         m{ipar}.addENmodBonded(6,6,1,2,enstruct1.CsCp);
         m{ipar}.addENmodBonded(6,6,2,2,enstruct1.CpCp);
      end
   elseif (size(enstructh,1) == 1)
      for ipar = params
         m{ipar}.addENmodDiag(1,1,enstructh.H);
         m{ipar}.addENmodDiag(6,1,enstructh.Cs);
         m{ipar}.addENmodDiag(6,2,enstructh.Cp);
         m{ipar}.addENmodBonded1h(1,1,enstructh.HH);
         m{ipar}.addENmodBonded1h(6,1,enstructh.CH);
         m{ipar}.addENmodBonded1h(1,6,enstructh.HC);
         m{ipar}.addENmodBonded1h(6,6,enstructh.CCs);
         m{ipar}.addENmodBonded1h(6,6,enstructh.CCp);
      end
   end
end
if (~isempty(e2struct) > 0)
   for ipar = params
      m{ipar}.addH2modDiag(1,e2struct.H);
      m{ipar}.addH2modDiag(6,e2struct.C);
      m{ipar}.addH2modOffDiag(1,1,e2struct.HH);
      m{ipar}.addH2modOffDiag(6,6,e2struct.CC);
      m{ipar}.addH2modOffDiag(1,6,e2struct.CH);
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
if (~isempty(e2struct))
   fitme.includeE2 = 1;
end
fitme.silent = silent;
fitme.setEnvs(envs);
% setEnvs calculates the HL values of everything we are fitting to, so the
% HLs are no longer needed. By removing these from fitme, we make the fitme
% object quite a bit smaller.
fitme.HLs = [];
if (doPlot > 0)
   fitme.plot = 1;
else
   fitme.plot = 0;
end
fitme.testFitme = testFitme;

end

function res = checkForInput(varargin,inputName,default)
idx = find(cellfun(@(x)isequal(lower(x),lower(inputName)), varargin));
if (~isempty(idx))
   res = varargin{idx+1};
else
   res = default;
end
end

