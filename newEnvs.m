%% CH4 charge plots, comparison of high-level and low-level
clear all;
load('datasets/ethanerDat.mat');
qs = cell(6,2);  % {Z, LL or HL}
ch3sum = cell(2,1);
close all;
col = {'b','r'};
egs = cell(2,1);
for ilevel = 1:2
   for imod = 1:size(LL,1);
      switch ilevel
         case 1
            m1 = LL{imod,1};
         case 2
            m1 = HL{imod,1};
      end
      for ienv = 1:m1.nenv
         t1 = m1.mcharge(ienv)-m1.mcharge(0);
         egs{ilevel}(end+1) = m1.EhfEnv(ienv) - m1.Ehf;
         ch3sum{ilevel}(end+1) = sum(t1([1 3 4 5]));
         ch3sum{ilevel}(end+1) = sum(t1([2 6 7 8]));
         for iatom = 1:length(t1)
            qs{m1.Z(iatom),ilevel}(end+1) = t1(iatom);
         end
      end
   end
   figure(101)
   hold on;
   subplot(2,2,2*(ilevel-1) +1);
   hist(qs{1,ilevel},10);
   hold on;
   subplot(2,2,2*(ilevel-1) + 2);
   hist(qs{6,ilevel},10);
end
figure(300)
hist(ch3sum{1},20)
figure(301)
hist(ch3sum{2},20)
%%
disp(['mean of C charges from LL and HL ', num2str(mean(qs{6,1})), ' ', num2str(mean(qs{6,2}))]);
disp(['mean of H charges from LL and HL ', num2str(mean(qs{1,1})), ' ', num2str(mean(qs{1,2}))]);
disp(['std of C charges from LL and HL  ', num2str(std(qs{6,1})), ' ', num2str(std(qs{6,2}))]);
disp(['std of H charges from LL and HL  ', num2str(std(qs{1,1})), ' ', num2str(std(qs{1,2}))]);
disp(['std of CH3 charges from LL and HL  ', num2str(std(ch3sum{1})), ' ', num2str(std(ch3sum{2}))]);

%%  Effect of a point charge on methane
% make environments
envs = cell(1,1);
ic = 0;
% for r=2:0.25:5
%    q = 1;
%    ic = ic+1;
%    envs{ic} =Environment;
%    envs{ic}.addCharge([0 0 r], q);
% end
r = 5;
for q=1:0.5:4
   ic = ic+1;
   envs{ic} =Environment;
   envs{ic}.addCharge([0 0 r], q);
end
dataroot = 'c:\matdl\yaron\10-30-12\methane-q1\';
if (exist(dataroot,'dir') ~= 7)
   status = mkdir(dataroot);
end
copyfile('templates/ch4.tpl',[dataroot,'/ch4.tpl']);

config = Fragment.defaultConfig();
config.template = 'ch4';
config.par = [1.12 1.12 1.12 1.12 109.47 109.47 109.47 120.0 -120.0]


config.basisSet = 'STO-3G';
LL = Fragment([dataroot], config);
for ienv = 1:length(envs)
   display(['LL env ',num2str(ienv)]);
   LL.addEnv(envs{ienv});
end
config.basisSet = '6-31G';
HL = Fragment([dataroot], config);
for ienv = 1:length(envs)
   display(['HL env ',num2str(ienv)]);
   HL.addEnv(envs{ienv});
end
%%
close all
for i= 1:length(envs)
   %x = envs{i}.r(3,1);
   x = envs{i}.rho(1);
   figure(500);
   hold on;
   plot(x,(LL.EhfEnv(i) - LL.Ehf) * 627.509,'bo');
   plot(x,(HL.EhfEnv(i) - HL.Ehf) * 627.509,'ro');
   figure(600);
   hold on;
   qLL = LL.mcharge(i) - LL.mcharge(0);
   qHL = HL.mcharge(i) - HL.mcharge(0);
   plot(x,qLL(1),'bo');
   plot(x,qHL(1),'ro');
   plot(x,qLL(2),'bx');
   plot(x,qHL(2),'rx');
   plot(x,qLL(3),'b+');
   plot(x,qHL(3),'r+');
end

%% Effects of dipoles from cube with edges of 5
mags = [];
edge = 5;
ic = 0;
envs = cell(1,1);
for dip = 5:5:20
   for iavg = 1:10
      ic = ic + 1;
      envs{ic} = Environment.dipCube(edge,0.1,dip);
      mags(ic) = dip;
   end
   dataroot = 'c:\matdl\yaron\10-30-12\methane-dip';
   if (exist(dataroot,'dir') ~= 7)
      status = mkdir(dataroot);
   end
end
save([dataroot,'\envs.mat'], 'envs','mags')
copyfile('templates/ch4.tpl',[dataroot,'\ch4.tpl']);
%%
clear all
dataroot = 'c:\matdl\yaron\10-30-12\methane-dip';
load([dataroot,'\envs.mat'], 'envs','mags');
config = Fragment.defaultConfig();
config.template = 'ch4';
config.par = [1.12 1.12 1.12 1.12 109.47 109.47 109.47 120.0 -120.0];

config.basisSet = 'STO-3G';
LL = Fragment([dataroot], config);
for ienv = 1:length(envs)
   display(['LL env ',num2str(ienv)]);
   LL.addEnv(envs{ienv});
end
config.basisSet = '6-31G';
HL = Fragment([dataroot], config);
for ienv = 1:length(envs)
   display(['HL env ',num2str(ienv)]);
   HL.addEnv(envs{ienv});
end
%%
close all
xold = mags(1);
yLL = [];
yHL = [];
for i = 1:length(envs);
   x = mags(i);
   figure(500);
   hold on;
   plot(x,(LL.EhfEnv(i) - LL.Ehf) * 627.509,'bo');
   plot(x,(HL.EhfEnv(i) - HL.Ehf) * 627.509,'ro');
   figure(600);
   hold on;
   qLL = LL.mcharge(i) - LL.mcharge(0);
   qHL = HL.mcharge(i) - HL.mcharge(0);
   if ((x == xold) && i<length(envs))
      yLL(end+1) = qLL(1);
      yHL(end+1) = qHL(1);
   else
      plot(xold,std(yLL),'bo');
      plot(xold,std(yHL),'ro');
      xold = x; yLL = []; yHL = [];
   end
%    plot(x,qLL(2),'bx');
%    plot(x,qHL(2),'rx');
%    plot(x,qLL(3),'b+');
%    plot(x,qHL(3),'r+');
end
%%
clear all;
nenv = 25;
edge = 5;
env = cell(25,1);
for ienv = 1:25
   t1 = Environment.dipCube(edge,0.1,25);
   r = -1 + (2 *rand(3,1));  % each element will go from -1 to +1 
   orien = r/norm(r);
   orien = 5 * orien;
   rho = 5 * (-1 + 2 * rand(1,1));
   t1.addCharge(orien,rho);
   env{ienv} = t1;
end
save('datasets/env3.mat','env');
