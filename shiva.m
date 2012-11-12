clear classes;
close all;
% All output will be placed in this directory:
topDir = 'C:/matdl/yaron/10-23-12/shiva/';

% Configuration variables
fitmeParallel = 1; % true = take advantage of matlabpool
showPlots = 0;  % true = plots comparisons of HL and LL Model on each
                %        call to fitme.err
separateSP = 0; % true = treat s and p orbitals differently
optRoutine = 1; % 0 = nlsqmin   1 = alternative (see below)
iprocess = 2;   % 1 = h2 (fast) 2 = ch4 (medium) 3 = ethane (slow)
                % 4 = ch4 and ethane (very slow)
maxIter = 500;  % for optRoutine = 0

% End of configuration variables


if (iprocess == 1)
   trainC{1} = {'h2',[2 3 4],'envs',1:5};
   testC{1} = {'h2',[6 7],'envs',10:15};
   filePrefix{1} = 'h2-cross2';
elseif (iprocess == 2)
   ikeep = [6     7     8    13    16    24];
   trainC{1} = {'ch4r',1:10,'envs',ikeep};
   ikeep2 = [5    10    14    17    20    25];
   testC{1} = {'ch4r',11:25,'envs',ikeep2};
   filePrefix{1} = 'ch4r-cross1';
elseif (iprocess == 3)
   ikeep = [6     7     8    13    16    24];
   trainC{1} = {'ethane',[1 3 5 7],'envs',ikeep};
   ikeep2 = [5    10    14    17    20    25];
   testC{1} = {'ethane',[2 4 6],'envs',ikeep2};
   filePrefix{1} = 'ethane-cross';
elseif (iprocess == 4)
   ikeep = [6     7     8    13    16    24];
   trainC{1} = {'ch4r',1:10,'ethane',[1 3 5 7],'envs',ikeep};
   ikeep2 = [5    10    14    17    20    25];
   testC{1} = {'ch4r',11:20,'ethane',[2 4 6],'envs',ikeep2};
   filePrefix{1} = 'hcfit1';
end

filePre = filePrefix{1};
dataDir = [topDir,filePre];
if (exist(dataDir,'dir') ~= 7)
   status = mkdir(dataDir);
end
summaryName = [topDir,filePre,'/summary.txt'];
% if (exist(summaryName,'file'))
%    delete(summaryName);
% end
summaryFile = fopen(summaryName,'a');
diaryName = [topDir,filePre,'/cfit.diary'];
% if (exist(diaryName,'file'))
%    delete(diaryName);
% end
diary(diaryName);
diary on;
commonIn = {'silent',0};

% Create the mixers
iC = 1;
trainIn = trainC{iC};
testIn = testC{iC};
filePre = filePrefix{iC};
ke = [];
en = [];
e2 = [];
iP2 = [1 0 0 0];
iP3 = [1 0 0 0 0];

ke.H = Mixer(iP3,11,'ke.H',3);
if (separateSP)
   ke.Cs = Mixer(iP3,11,'ke.Cs',3);
   ke.Cp = Mixer(iP3,11,'ke.Cp',3);
else
   ke.Cs = Mixer(iP3,11,'ke.C',3);
   ke.Cp = ke.Cs;
end
ke.HH = Mixer(iP2,12,'ke.HH',2);
ke.CH = Mixer(iP2,12,'ke.CH',2);
ke.CH.hybrid = 1;
ke.CCs = Mixer(iP2,12,'ke.CCs',2);
ke.CCs.hybrid = 1;
ke.CCp = Mixer(iP2,12,'ke.CCp',2);
ke.CCp.hybrid = 2;

en.H = Mixer(iP3,11,'en.H',3);
if (separateSP)
   en.Cs = Mixer(iP3,11,'en.Cs',3);
   en.Cp = Mixer(iP3,11,'en.Cp',3);
else
   en.Cs = Mixer(iP3,11,'en.C',3);
   en.Cp = en.Cs;
end
en.HH = Mixer(iP2,12,'en.HH',2);
en.CH = Mixer(iP2,12,'en.CH',2);
en.CH.hybrid = 1;
en.HC = en.CH;
en.CCs = Mixer(iP2,12,'en.CCs',2);
en.CCs.hybrid = 1;
en.CCp = Mixer(iP2,12,'en.CCp',2);
en.CCp.hybrid = 2;

e2.H = Mixer(iP2,11,'e2.H',2);
e2.C = Mixer(iP2,11,'e2.C',2);
if (iprocess == 1) % h2
   e2.HH = Mixer(iP2,12,'e2.HH',2);
else % want only bond length dependence in H-H for methane
   e2.HH = Mixer([1 0],4,'e2.HH',2);
end
e2.CC = Mixer(iP2,12,'e2.CC',2);
e2.CH = Mixer(iP2,12,'e2.CH',2);

% Create fitme object
f1 = makeFitme(trainIn{:},commonIn{:},'enstructh',en, ...
   'kestructh',ke,'e2struct',e2);
f1.parallel = fitmeParallel;
ftest = makeFitme(testIn{:},commonIn{:},'enstructh',en, ...
   'kestructh',ke,'e2struct',e2);
ftest.parallel = fitmeParallel;
%f1 = makeFitme(trainIn{:},commonIn{:},'enmods',0, ...
%   'kestructh',ke);

% Fix all context sensitive parameters, i.e. no context fitting
for imix = 1:length(f1.mixers)
   mix = f1.mixers{imix};
   if ((mix.mixType == 11) && mix.funcType == 3)
      mix.fixed = [0 1 1 1 0];
   elseif ((mix.mixType == 11) && mix.funcType == 2)
      mix.fixed = [0 1 1 1];
   elseif (mix.mixType == 12)
      mix.fixed = [0 1 1 1];
   elseif (mix.mixType == 4)
      mix.fixed = [0 1];
   end
end

f1.plot = showPlots; 
ftest.plot = showPlots;

% Set limits on the parameters
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
start = f1.getPars;

if (optRoutine == 0)
   options = optimset('DiffMinChange',1.0e-4,'TolFun',1.0e-3, ...
      'TolX',3.0e-3,'MaxFunEvals',maxIter);
   [pt,resnorm,residual,exitflag,output,lambda,jacobian] = ...
      lsqnonlin(@f1.err, start,lowLimits,highLimits,options);
   err = sqrt(residual*residual'/length(residual))*627.509;
elseif (optRoutine == 1)
end
testres = ftest.err(pt);
testErr = sqrt(testres*testres'/length(testres))*627.509;

str2 = 'error train %12.5f test %12.5f \n';
fprintf(1,str2,err,testErr);

% With iprocess = 1 and optRoutine = 0, should give:
% Fitme.err called with par = 1.0007    -0.29093      1.9715      1.0032      0.3377      1.2321      0.6567      1.3157
% density matrices will be recalculated Densities reset
% RMS err/ndata = 0.0034796 kcal/mol err = 15.1275
% error train      5.64975 test     15.12751 

