%%
clear classes;
dataDir = 'D:\tmp';
LLFile = 'ethane1r-gen4H';
HLFile = 'ethane1r-gen31mult';
copyfile(['D:\dave\apoly\msqc\templates\',LLFile,'.tpl'],...
   [dataDir,'\',LLFile,'.tpl']);
copyfile(['D:\dave\apoly\msqc\templates\',HLFile,'.tpl'],...
   [dataDir,'\',HLFile,'.tpl']);

obj = FitBasis(dataDir,HLFile,LLFile);
r1 = 0.001;
t1 = 0.001;
obj.setRandomGeom('CC bond', 1.54 - r1, 1.54 + r1);
obj.setRandomGeom('CH bond', 1.12 - t1, 1.12 + t1, 6);
obj.setRandomGeom('CCH angle', 110.5 - t1, 110.5 + t1, 6);
obj.setRandomGeom('free dihedral', 0, 0+t1, 1);
obj.setRandomGeom('constrained dihedral', 120-t1, 120+t1, 2);
obj.setRandomGeom('constrained dihedral', -120-t1, -120+t1, 2);
obj.addRandomGeom(1);
% load('D:\dave\apoly\msqc\datasets\env2.mat');
% for i= 1:2
%    obj.addEnv(env{i});
% end
obj.cacheHL = false;
obj.cacheLL = false;

obj.generateHLData(1);
save('d:\tmp\verifyFitBasis1.mat');
%%
clear classes;
load('d:\tmp\verifyFitBasis1.mat');
obj.setSpin(1);
aone = ones(obj.LLnpar - length(obj.geomPars{1}),1);
ic = 0;
for sc = 0.25:0.25:2
   ic = ic+1;
  [e1{ic}, hl1{ic}, ll1{ic}] = obj.err(sc*aone);
end
%%

for i=1:length(e1)
   disp(['err = ',num2str(norm(e1{i}))]);
end
for i=1:length(e1)
   disp([num2str(e1{i}(:)')]);
end

%%
start = [0.75 0.75 0.75];
lowLimits = [0.1 0.1 0.1];
highLimits = [10 10 10];
options = optimset('DiffMinChange',1.0e-4,'TolFun',1.0e-10, ...
   'TolX',3.0e-10,'MaxFunEvals',1e5);
[pt,resnorm,residual,exitflag,output,lambda,jacobian] = ...
   lsqnonlin(@obj.err, start,lowLimits,highLimits,options);
