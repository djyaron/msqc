clear classes;
K = 25; % number of arms
T = 1e6; % number of rounds
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
%%
f1.generateArms(K,-10,10);

% sort arms from best to worst
errArm = zeros(K,1);
for iarm = 1:K
   errArm(iarm) = f1.armError(iarm);
end
[errs,ierrs] = sort(errArm);
f1.arms = f1.arms(:,ierrs);
%for iarm = 1:K
%   errArm(iarm) = band.fullError(iarm);
%end
figure(100)
plot(errs,'bo-');
%%

gamma = 0.5;
w = ones(K,1);
%
ic = 0;
esum = 0;
for t=1:T
   p = (1-gamma) * w/sum(w) + gamma/K;
   armToPull = find(mnrnd(1,p));
   error = f1.pullArm(armToPull);
   ic = ic + 1;
   esum = esum + error;
   reward = exp(-error); % will have required range of 0 to 1
   w(armToPull) = w(armToPull) * exp(gamma*reward/(K*p(armToPull)));
   if (rem(t,100) == 0)
      figure(101);
      plot(p,'rx-');
      title(['after ',num2str(t),' iterations ','avg error ',num2str(esum/ic)]);
      input('ok');
   end
end

%%
limits = [];
options = optimset('DiffMinChange',1.0e-5,'TolFun',1.0e-4,'TolX',1.0e-3);
start = zeros(1,f1.npar);
[pt,resnorm,residual,exitflag,output,lambda,jacobian] = ...
   lsqnonlin(@f1.err, start,-limits,limits,options);
%%
start = f1.arms(:,1)';
[pt2,resnorm,residual,exitflag,output,lambda,jacobian] = ...
   lsqnonlin(@f1.err, start,-limits,limits,options);

