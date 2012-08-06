function demo4
% Function that demonstrates use of getFitme

% EXAMPLE 1 Hydrogen
f1 = getFitme(1,1,1,1,'h2',2:7);

[lowLimits,highLimits] = getLimits(f1);
start = f1.getPars;
options = optimset('DiffMinChange',1.0e-5,'TolFun',1.0e-3,'TolX',3.0e-3);

[pt,resnorm,residual,exitflag,output,lambda,jacobian] = ...
   lsqnonlin(@f1.err, start,lowLimits,highLimits,options);
display('Hydrogen results');
display(['pt = ',num2str(pt)]);
display(['resnorm = ',num2str(resnorm)]);

%pt = 0.62343  1.9482  0.67907  1.2712  0.48982  1.4852
%resnorm = 0.076006  (7.706 kcal/mol)

%
% EXAMPLE 2 methane fit: KE only for 7 geometries (geoms can go up to 1:17)
f1 = getFitme(1,0,0,1,'ch4',1:7);

[lowLimits,highLimits] = getLimits(f1);
start = f1.getPars;
options = optimset('DiffMinChange',1.0e-5,'TolFun',1.0e-3,'TolX',3.0e-3);

[pt,resnorm,residual,exitflag,output,lambda,jacobian] = ...
   lsqnonlin(@f1.err, start,lowLimits,highLimits,options);

display('Methane results');
display(['pt = ',num2str(pt)]);
display(['resnorm = ',num2str(resnorm)]);

%Methane results
%pt = 1.0211     0.85864      1.3767
%resnorm = 0.2444 (25.5857 kcal/mol)

% EXAMPLE 3: methane fit: all energies, for 7 geometries (geoms can go up to 1:17)
f1 = getFitme(1,1,1,1,'ch4',1:7);

[lowLimits,highLimits] = getLimits(f1);
start = f1.getPars;
options = optimset('DiffMinChange',1.0e-5,'TolFun',1.0e-3,'TolX',3.0e-3);

[pt,resnorm,residual,exitflag,output,lambda,jacobian] = ...
   lsqnonlin(@f1.err, start,lowLimits,highLimits,options);

display('Methane results');
display(['pt = ',num2str(pt)]);
display(['resnorm = ',num2str(resnorm)]);

%Methane results
%pt = 1.1305 0.77577 1.9704 0.79484 0.92042 1.1347 0.20348 0.95222 0.95053 0.93869
%resnorm = 0.32642 (11.1763 kcal/mol)
end

function f1 = getFitme(KE,EN,E2,scale,varargin)
% input
%   KE,EN,E2: int that specifies how to handle these energy terms
%             0 = ignore   1 = regular   2 = hybrid
%             E2 does not yet support hybrid
%   scale   : bool using scaling (with constant) if true
%   varargin : other argument to pass to makeFitme, especially 
%              for h2 use 'h2',2:7   for methane: 'ch4',1:17

args = varargin;
if (scale)
   ftype = 2;
else
   ftype = 4;
end
if (KE==0)
   args = {args{:},'kemods',0};
elseif (KE==1)
   ke = createStruct('ke',ftype,0);
   args = {args{:},'kestruct',ke};   
elseif (KE==2)
   ke = createStruct('ke',ftype,1);
   args = {args{:},'kestructh',ke};   
end
if (EN==0)
   args = {args{:},'enmods',0};
elseif (EN==1)
   en = createStruct('en',ftype,0);
   args = {args{:},'enstruct',en};   
elseif (EN==2)
   en = createStruct('en',ftype,1);
   args = {args{:},'enstructh',en};   
end

if (E2 == 1)
   if (scale == 1)
      iP = 1;
   else
      iP = 0;
   end
   e2.H = Mixer(iP,1,'e2.H',ftype);
   e2.C = Mixer(iP,1,'e2.C',ftype);
   e2.HH = Mixer(iP,1,'e2.HH',ftype);
   e2.CC = Mixer(iP,1,'e2.CC',ftype);
   e2.CH = Mixer(iP,1,'e2.CH',ftype);
   args = {args{:},'e2struct',e2};   
elseif (E2 == 2)
   error('getFitme: E2 == 2, but hybrid not supported for E2');
end

f1 = makeFitme(args{:});
end

function res = createStruct(label,ftype,hybrid)
% creates a structure with all relevant mixers for makeFitme
% Input:
%   label: string, desc field of Mixer will start with this label
%   ftype:  ftype for Mixer (see Mixer.m)
%   hybrid: bool, creates Mixers for hybrid orbitals
% Output:
%   res: structure with fields shown below, each of which holds a handle to
%        a Mixer. fields in parens have handles to same Mixer, e.g.
%        (Cs Cp) means res.Cs  and res.Cp point to same Mixer
%      hybrid = 0 : H (Cs Cp) HH (CsH CpH) (CsCs CsCp CpCp)
%      hybrid = 1 : H (Cs Cp) HH CH CCs CCp

% use 1 for initial parameter if ftype is scaling, 0 for interpolation
if ((ftype == 2) || (ftype == 3))
   iP = 1;
else
   iP = 0;
end

if (~hybrid)
   res.H    = Mixer(iP,1,[label,'.H'],ftype);
   res.Cs   = Mixer(iP,1,[label,'.C'],ftype);
   res.Cp   = res.Cs;
   res.HH   = Mixer(iP,1,[label,'.HH'],ftype);
   res.CsH  = Mixer(iP,1,[label,'.CH'],ftype);
   res.CpH  = res.CsH;
   res.CsCs = Mixer(iP,1,[label,'.CC'],ftype);
   res.CsCp = res.CsCs;
   res.CpCp = res.CsCs;
else
   res.H   = Mixer(iP,1,[label,'.H'],ftype);
   res.Cs  = Mixer(iP,1,[label,'.C'],ftype);
   res.Cp  = res.Cs;
   res.HH  = Mixer(iP,1,[label,'.HH'],ftype);
   res.CH  = Mixer(iP,1,[label,'.CH'],ftype);
   res.CH.hybrid = 1;
   res.CCs = Mixer(iP,1,[label,'.CCs'],ftype);
   res.CC.hybrid = 1;
   res.CCp = Mixer(iP,1,[label,'.CCp'],ftype);
   res.CCp.hybrid = 2;
end
end

function [lowLimits, highLimits] = getLimits(f1)
% For scale parameters, sets limits from 0 to 5
% all other parameters have limits -infinity to infinity
maxScale = 5;

lowLimits = zeros(f1.npar,1);
highLimits = lowLimits;
i1 = 1;
for imix = 1:length(f1.mixers)
   mix = f1.mixers{imix};
   if ((mix.funcType == 2)||(mix.funcType == 3))
      lowLimits(i1) = 0.0;
      highLimits(i1) = maxScale;
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
end
