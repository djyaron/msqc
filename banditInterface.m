clear classes;
ke.H = Mixer(0,1,'ke.H');
ke.Cs = Mixer(0,1,'ke.C');
ke.Cp = ke.Cs;
ke.HH = Mixer(0,1,'ke.HH');
ke.CsH = Mixer(0,1,'ke.CH');
ke.CpH = ke.CsH;
ke.CsCs = Mixer(0,1,'ke.CC');
ke.CsCp = ke.CsCs;
ke.CpCp = ke.CsCs;
en.H = Mixer(0,1,'en.H');
en.Cs = Mixer(0,1,'en.C');
en.Cp = en.Cs;
en.HH = Mixer(0,1,'en.HH');
en.CsH = Mixer(0,1,'en.CH');
en.CpH = en.CsH;
en.CsCs = Mixer(0,1,'en.CC');
en.CsCp = en.CsCs;
en.CpCp = en.CsCs;
 
f1 = makeFitme('h2',2:7,'enstruct',en,'kestruct',ke);
p = f1.getPars; % just getting current parameters to have correct
                % length vector
a = f1.randMolError(p);

%% Traditional minimization approach
start = f1.getPars;
options = optimset('DiffMinChange',1.0e-5,'TolFun',1.0e-4,'TolX',1.0e-3);
limits = [];
[pt2,resnorm,residual,exitflag,output,lambda,jacobian] = ...
   lsqnonlin(@f1.err, start,-limits,limits,options);
% pt2 is the set of parameters that minimize the energy
