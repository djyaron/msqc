%% Read data into frag1 and frag2
clear classes;

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
m1.diffOverlap();
m1.solveHF;
m2.diffOverlap();
m2.solveHF;
disp(['max HF diff c1 ', num2str(max(abs(m1.EhfEnv - frag1.EhfEnv)))]);
disp(['max HF diff c2 ', num2str(max(abs(m2.EhfEnv - frag2.EhfEnv)))]);

%%
plot(m1.EhfEnv - frag1.EhfEnv, ...
   m2.EhfEnv - frag2.EhfEnv, 'r.');

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

%%
px = 0.5:0.1:1.5;
nx = size(px,2);
err = zeros(1,nx);
frag = frag1;
for ix=1:nx
   m1 = Model1(frag);
   m1.scaleH1(1.0,px(ix));
   m1.solveHF;
   t1 = m1.EhfEnv - frag2.EhfEnv;
   %t1 = t1 - mean(t1);
   for j = 1:size(t1,2)
      hold on
      plot(px(ix), t1(j),'r.');
   end
end
