clear classes;
load('C:\Users\yaron\Documents\My Dropbox\MSQCdata\dec12e\w1\hybridslater1\ch4rDat\start.mat')
load('datasets/ch4rDat.mat');

%f1.setPars(zeros(size(f1.getPars)));
%f1.silent = 0;
%etest = f1.err(f1.getPars);
%%
q1 = cell(0,0);
q2 = cell(0,0);
q3 = cell(0,0);
qd = [];
for imod = 1:length(f1.models)
   f1.models{imod}.solveHF(f1.envs{imod});
   q1{imod} = f1.models{imod}.charges(:,f1.envs{imod}+1);
   q2{imod} = f1.models{imod}.mcharge(f1.envs{imod});
   qs = zeros(size(q2{imod}));
   ic = 0;
   for ienv = f1.envs{imod}
      ic = ic+1;
      qs(:,ic) = HL{imod,1}.mcharge(ienv)';
   end
   q3{imod} = qs;
   qd1 = q1{imod}-q2{imod};
   qd = [qd qd1(:)'];
end
%%
clc
for imod = 1:length(f1.models)
   rcart = f1.models{imod}.rcart;
   c1 = q1{imod};
   c2 = q2{imod};
   c3 = q3{imod};
   cd = c1-c2;
   evs = f1.envs{imod};
   for ic = 1:length(evs)
      disp(['model ',num2str(imod),' env ',num2str(evs(ic))]);
      qq = c1(:,ic);
      dip = zeros(3,1);
      for iatom =1:size(rcart,2);
         dip = dip + qq(iatom) * rcart(:,iatom);
      end
      disp(['q1 ',num2str(c1(:,ic)'),' ',num2str(norm(dip))]);
      qq = c2(:,ic);
      dip = zeros(3,1);
      for iatom =1:size(rcart,2);
         dip = dip + qq(iatom) * rcart(:,iatom);
      end
      disp(['q2 ',num2str(c2(:,ic)'),' ',num2str(norm(dip))]);
      qq = c3(:,ic);
      dip = zeros(3,1);
      for iatom =1:size(rcart,2);
         dip = dip + qq(iatom) * rcart(:,iatom);
      end
      disp(['q3 ',num2str(c3(:,ic)'),' ',num2str(norm(dip))]);
      disp(['df ',num2str(cd(:,ic)')]);
   end
end

%% Change in charges with geometry
q1 = cell(length(f1.models),f1.models{1}.nenv+1);
q2 = cell(length(f1.models),f1.models{1}.nenv+1);
for imod = 1:length(f1.models)
   for ienv = 0:f1.models{imod}.nenv
      f1.models{imod}.solveHF(ienv);
      q1{imod,ienv+1} = f1.models{imod}.frag.mcharge(ienv);
      q2{imod,ienv+1} = f1.models{imod}.mcharge(ienv);
   end
end
%% absolute charges
qstruct = [];
qenv = [];
qs = q2;
atoms = 2:5;
for imod = 1:length(f1.models)
   qstruct = [qstruct qs{imod,1}(atoms)'];
   for ienv = 1:f1.models{imod}.nenv
      qd = qs{imod,ienv+1}-qs{imod,1};
      qenv = [qenv qd(atoms)'];
   end
end
figure(1)
hist(qstruct);
figure(2)
hist(qenv);
%% how do induced charges differ between LL and model
qindH1 = [];
qindH2 = [];
qindC1 = [];
qindC2 = [];
for imod = 1:length(f1.models)
   for ienv = 1:f1.models{imod}.nenv
      qd1 = q1{imod,ienv+1}-q1{imod,1};
      qindH1 = [qindH1 qd1(2:5)];
      qindC1 = [qindC1 qd1(1)];
      qd2 = q2{imod,ienv+1}-q2{imod,1};
      qindH2 = [qindH2 qd2(2:5)];
      qindC2 = [qindC2 qd2(1)];
   end
end
figure(10)
plot(qindH1,qindH2,'r.');
figure(11)
plot(qindC1,qindC2,'b.');
