% Exploratory data analysis

%% create environments
clear classes;
size = [3,3,3];
nenv = 50;
env(1,nenv) = Environment;
for i = 1:nenv
   env(1,i) = Environment.newCube(size,1);
end
save('data1/env.mat', 'env');

%% Generate the data
clear classes;
load('data1/env.mat');
nenv = size(env,2);

c1 = Fragment.defaultConfig();
c1.template = 'h2';
c1.basisSet = 'STO-3G';
c1.par = 1.1;
c2 = c1;
c2.basisSet = '6-31G**';

%%
frag1 = Fragment('data1',c1);
for i=1:25
   frag1.addEnv(env(i));
end
%%
frag2 = Fragment('data1',c2);
for i=1:25
   frag2.addEnv(env(i));
end

%%
frag1 = Fragment('data1',c1);
frag2 = Fragment('data1',c2);
for i=1:25
   frag1.addEnv(env(i));
   frag2.addEnv(env(i));
end

%% Compare Hnuc to make sure things are ok
plot(frag1.HnucEnv-frag1.Hnuc, frag2.HnucEnv-frag2.Hnuc, 'r.');

%% The difference in the energy on an atom is related to the density on 
%  that atom
nenv = 25;
x = zeros(nenv,1);
y = zeros(nenv,1);
z = zeros(nenv,1);
sym = {'r.','bo'};
for iatom = 1:2
   for i=1:nenv
      %   x(i) = sum(sum(frag1.partitionE1(i)));
      %   y(i) = sum(sum(frag2.partitionE1(i)));
      t1 = frag1.partitionE1(i);
      t2 = frag2.partitionE1(i);
      x(i) = t1(iatom,iatom);
      y(i) = t2(iatom,iatom);
      %x(i) = sum(sum(t1));
      %y(i) = sum(sum(t2));
      t3 = frag1.density(i);
      z(i) = t3(iatom,iatom);
      if (iatom == 1)
         hold off;
      else
         hold on;
      end
      plot(z,x-y,sym{iatom});
   end
end

%% The difference in the energy on an is related to the density on that
%  atom
nenv = 25;
x = zeros(nenv,1);
y = zeros(nenv,1);
z = zeros(nenv,1);
sym = {'r.','bo'};
for i=1:nenv
   %   x(i) = sum(sum(frag1.partitionE1(i)));
   %   y(i) = sum(sum(frag2.partitionE1(i)));
   t1 = frag1.partitionE1(i);
   t2 = frag2.partitionE1(i);
   x(i) = t1(1,2);
   y(i) = t2(1,2);
   t3 = frag1.density(i);
   t4 = frag2.density(i);
   z(i) = abs(t4(1,2)-t3(1,2));
end
plot(z,x-y,sym{iatom});

%% Create model classes
clear classes;
load('data1/env.mat');
nenv = size(env,2);

c1 = Fragment.defaultConfig();
c1.template = 'h2';
c1.basisSet = 'STO-3G';
c1.par = 1.1;
c2 = c1;
c2.basisSet = '6-31G**';

frag1 = Fragment('data1',c1);
frag1.loadAllEnv();
frag2 = Fragment('data1',c2);
for i=1:frag1.nenv
   frag2.addEnv(frag1.env(i));
end
%%
m1 = Model1(frag1);
m2 = Model1(frag2);
%m1.diffOverlap();
m1.solveHF;
%m2.diffOverlap();
m2.solveHF;
disp(['max HF diff c1 ', num2str(max(abs(m1.EhfEnv - frag1.EhfEnv)))]);
disp(['max HF diff c2 ', num2str(max(abs(m2.EhfEnv - frag2.EhfEnv)))]);

%%
frag = frag2;
mod  = m2;
for i=0:frag.nenv
   t1 = frag.partitionE1(i);
   t2 = mod.partitionE1(i);
   disp(['env ',num2str(i),' partitioned E1 agree to ', ...
      num2str(max(max(abs(frag.partitionE1(i)-mod.partitionE1(i))))), ...
      ' and E2 to ', num2str(max(max(max(max(abs(frag.partitionE2(i) - ...
      mod.partitionE2(i)))))))]);
end

