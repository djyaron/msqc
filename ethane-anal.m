%% Load data
clear classes;
root = 'c:\dave\apoly\msqc\';
% Generate environments for production runs
if (exist('ethane4/env2.mat','file'))
   disp('loading existing environments');
   load('ethane4/env2.mat');
else
   mag = 15.0;
   nenv = 50;
   cubSize = [6,6,6];
   cent = [0.77; 0; 0];
   for ienv = 1:nenv
      temp = Environment.newCube(cubSize,mag);
      temp.displace(cent);
      env{ienv} = temp;
   end
   save('ethane4/env2.mat','env');
end
nenv = size(env,2);
pars{1} = [1.54 1.12 60];
pars{2} = [1.54 1.12 30];
pars{3} = [1.54 1.12 0];
pars{4} = [1.39 1.12 60];
pars{5} = [1.69 1.12 60];
pars{6} = [1.54 0.97 60];
pars{7} = [1.54 1.27 60];
npar = size(pars,2);
HLbasis = {'6-31G**'}; %{'6-31G' '6-31G*' '6-31G**'};
HL = cell(npar,3);
LL = cell(npar,3);
%%
for ipar = 1:size(pars,2)
   par = pars{ipar};
   disp(['rcc ',num2str(par(1)), ...
      ' rch ',num2str(par(2)), ...
      ' angle ',num2str(par(3))]);
   
   config = Fragment.defaultConfig();
   config.method = 'MP2';
   config.par = par;
   
   % HL
   for ihl = 1:size(HLbasis,2)
      config.template = 'ethane1';
      config.basisSet = HLbasis{ihl};
      disp(['ipar ',num2str(ipar),' loading HL ',num2str(ihl)]);
      frag1 = Fragment([root,'ethane4mp2'], config);
      for ienv = 1:nenv
         display(['HL env ',num2str(ienv)]);
         frag1.addEnv(env{ienv});
      end
      HL{ipar,ihl} = frag1;
   end
   % LL 1
   config.basisSet = 'STO-3G';
   frag2 = Fragment([root,'ethane4mp2'], config);
   disp(['ipar ',num2str(ipar),' loading LL 1']);
   for ienv = 1:nenv
      display(['LL env ',num2str(ienv)]);
      frag2.addEnv(env{ienv});
   end
   LL{ipar,1} = frag2;
   
   % LL 2
   config.template = 'ethane1-gen';
   config.basisSet = 'GEN';
   config.par = [par 0.9 0.9 0.9 0.9 0.9];
   frag3 = Fragment([root,'ethane4mp2'], config);
   disp(['ipar ',num2str(ipar),' loading LL 2']);
   for ienv = 1:nenv
      display(['LL env ',num2str(ienv)]);
      frag3.addEnv(env{ienv});
   end
   LL{ipar,2} = frag3;
   % LL 3
   config.template = 'ethane1-gen';
   config.basisSet = 'GEN';
   config.par = [par 1.1 1.1 1.1 1.1 1.1];
   disp(['ipar ',num2str(ipar),' loading LL 3']);
   frag4 = Fragment([root,'ethane4mp2'], config);
   for ienv = 1:nenv
      display(['LL env ',num2str(ienv)]);
      frag4.addEnv(env{ienv});
   end
   LL{ipar,3} = frag4;
end


%% since even loading all the files will take time, we'll dave everything
save('ethane4mp2/ethaneDat.mat');


%% Aggregator fits
%clear classes;
load('ethaneDat.mat');
agg = Aggregator;
for ipar = 1:1 %size(HL,1)
   agg.addFrags(HL{ipar},LL{ipar,1},LL{ipar,2});
end
[diff, ehigh, epred]  = agg.err1([0 1]);
fit = lsqnonlin(@agg.err1, [0.5 0.5]);
%%
%trying err 2
clear classes;
load('temp1.mat');
agg = Aggregator;
for ipar = [1 4 5]%size(HL,1) %%%need this to work for multiple parameters
   %add second low level into calculation
   agg.addFrags(HL{ipar},LL{ipar,1},LL{ipar,2});
end
%[diff, ehigh, epred, elow]  = agg.err2([0 1 1 1 1 1 1 1]);%%add in list of environmet, env
fit = lsqnonlin(@agg.err2, [1 1 1 1 1 1 1 1]);
[diff, ehigh, epred elow] = agg.err2(fit);
%%
[diff, ehigh, epred]  = agg.err1(fit);
%%
save([root,'ethane2\data3.mat'],'HL','LL');
%%
clear classes;
%root = 'c:\dave\apoly\msqc\';
load(['ethane2\data3.mat']);
%% determine the energy of interaction with the environment
EHL = cell(7,1);
for ipar = 1:7
   frag = HL{ipar,1};
   Eenv = zeros(1,frag.nenv);
   for ienv = 1:frag.nenv
      Eenv(ienv) = sum(sum( frag.density(ienv).*frag.H1Env(:,:,ienv) ) );
      Eenv(ienv) = Eenv(ienv) + frag.HnucEnv(ienv);
   end
   EHL{ipar,1} = Eenv;
end
ELL = cell(7,3);
for ipar = 1:7
   for i2 = 1:3
      frag = LL{ipar,i2};
      Eenv = zeros(1,frag.nenv);
      for ienv = 1:frag.nenv
         Eenv(ienv) = sum(sum( frag.density(ienv).*frag.H1Env(:,:,ienv) ) );
         Eenv(ienv) = Eenv(ienv) + frag.HnucEnv(ienv);
      end
      ELL{ipar,i2} = Eenv;
   end
end


%% plot versus environment
figure(101);
sym = {'ro','bo','go','ko','rx','bx','gx'};
for ipar = 1:1 %size(HL,1)
   etotal =  HL{ipar,1}.EhfEnv;
   emol = etotal - EHL{ipar,1};
   hold off;
   plot(etotal,'bo');
   hold on;
   plot(emol,'ro');
   %plot(LL{ipar}.EhfEnv,HL{ipar}.EhfEnv,sym{ipar});
end
%mean = edifft;

%%
npar = size(HL,1);
figure(90);
sym={'ro','bo','go'};
for ipar = 1:3
   for ienv = 1:LL{ipar,3}.nenv
      hold on;
      %plot(HL{ipar}.EhfEnv(ienv), LL{ipar,1}.EhfEnv(ienv),'b.');
      plot(ienv,norm(LL{ipar}.dipoleEnv(:,ienv)),sym{ipar});
   end
end
%%  Look at energy versus some key variables
econv = 627.509; % kcal/mole per hartree
x = zeros(3,1);
y = zeros(3,1);
figure(100);
goodEnv = [];
for ienv = 0:HL{1,1}.nenv
   ic = 0;
   for ipar = [1 2 3] % this goes over dihedral angle
      frag = HL{ipar,1};
      Eenv = EHL{ipar,1};
      ic = ic + 1;
      x(ic) = frag.config.par(3);
      if (ienv == 0)
         y(ic) = frag.Ehf;
      else
         y(ic) = frag.EhfEnv(1,ienv) -Eenv(1,ienv);
      end
   end
   y = (y - y(1)) * econv;
   ytest = max(y)-min(y);
   if (ytest > 10)
      disp(['environment ',num2str(ienv),' has y spread of ', num2str(ytest)]);
   else
      goodEnv = [goodEnv, ienv];
      hold on;
      if (ienv == 0)
         plot(x,y,'bo-');
      else
         plot(x,y,'r-');
      end
   end
end
%% used above to find intersection of good environments for
% HL, LL{:,1} and LL{:,3}. This gave me 87 environments.
save('ethane2\realGood','realGood');
%% correlation between LL and HL data for barrier heights
for ienv =1:100
   barrierHL(ienv) = ( HL{3,1}.EhfEnv(ienv)-EHL{3,1}(ienv) ) - ...
      ( HL{1,1}.EhfEnv(ienv) - EHL{1,1}(ienv) );
   barrierLL(ienv) = ( LL{3,1}.EhfEnv(ienv)-ELL{3,1}(ienv) ) - ...
      ( LL{1,1}.EhfEnv(ienv) - ELL{1,1}(ienv) );
end
barrierLL = 627.509 * barrierLL;
barrierHL = 627.509 * barrierHL;
figure(100);
plot(-barrierLL,-barrierHL,'b.');
hold on;
plot(-barrierHL,-barrierHL,'k-');

%% Can we build a lame model to do the above
econv = 627.509; % kcal/mole per hartree
nenv = 100;
%x = zeros(2,nenv);
%y = zeros(1,nenv*3);
clear x;
clear y;
idata = 0;
for ienv = 1:100
   if ((barrierLL(ienv) < -2) || (barrierLL(ienv)> -0.8) )
      
   else
      for ipar = 1:3  % dihedral angles
         idata = idata+1;
         x(idata,1) = (LL{ipar,1}.EhfEnv(ienv)-ELL{ipar,1}(ienv)) * econv;
         x(idata,2) = (LL{ipar,2}.EhfEnv(ienv)-ELL{ipar,2}(ienv)) * econv;
         x(idata,3) = 1.0;
         y(idata,1) = (HL{ipar,1}.EhfEnv(ienv)-EHL{ipar,1}(ienv)) * econv;
      end
   end
end
%%
figure(120);
hold off;
plot(x(:,1),y(:,1),'b.');
hold on
plot(x(:,2),y(:,1),'r.');
%%

[b,bint,r,rint,stats] = regress(y,x);
yp = x*b;
%%
nenv = idata/3;
ienv = 0;
for i=1:3:nenv*3
   ienv = ienv+1;
   BHLL(ienv) = x(i+2,1)-x(i,1);
   BHHL(ienv) = y(i+2,1)-y(i,1);
   BHpred(ienv) = yp(i+2,1) - yp(i,1);
end
figure(350)
hold off;
%plot(-BHHL,-BHLL,'bo');
%hold on;
plot(-BHHL,-BHpred,'rx');
%hold on;
%plot(-BHHL,-BHHL,'k-');

%% PLaying with model2
clear classes;
load('ethanefixed.mat');

f1 = LL{1,1};
f2 = LL{1,2};
f3 = LL{1,3};
fhl = LL{1,1};

m2 = Model2(f1, f2, f3, fhl);
m2.par = [0 0 0 0];
m2.nenv = 20;
m2.solveHF;
disp(['HF energy ',num2str(m2.Ehf),' ',num2str(f1.Ehf), ...
   ' ', num2str(m2.Ehf - f1.Ehf)]);
disp(['max Eorb diff ',num2str( max(m2.Eorb-f1.Eorb))]);
%disp(['max EhfEnv diff ',num2str( max(max(m2.EhfEnv-f1.EhfEnv)))]);
%disp(['max EorbEnv diff ',num2str( max(max(m2.EorbEnv-f1.EorbEnv)))]);
par = [ 0 0 0 0];
%%
par = [ 0 0 0 0];
for i=1:5
   m2.updateDensity(par);
   par = lsqnonlin(@m2.errApprox, par);
   errFit = m2.errApprox(par);
   errReal = m2.err(par);
   disp(['iter ',num2str(i),' par ',num2str(par)]);
   disp([' errorApprox ',...
      num2str(max(errFit)),' errorReal ',num2str(max(errReal))]);
end

