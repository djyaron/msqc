%% Testing analyze bond
% m1 is the train set;  m2 is the test set
%clear classes;
%load('datasets\ch4rDat.mat');
load('datasets\ethanerDat.mat');
m1 = Model3(LL{1,1},LL{1,1},LL{1,1});
m1.solveHF;
c1 = Context(nmod);
%%
atom1 = 1;
ic = 0;
res = zeros(1,3);
for ienv = 0:m1.nenv
   for atom2 = 2
      ic = ic + 1;
      res(ic,:) = c1.analyzeBond(m1,ienv,atom1,atom2)';
   end
end

%% Test of fillInContexts
% m1 is the train set;  m2 is the test set
clear classes;
load('datasets\ch4rDat.mat');
ftrain = LL(1:10,1);
ftest = LL(11:19,1);
envs1 = 1:10;
envs2 = 11:15;
mtrain = cell(size(ftrain));
mtest = cell(size(ftest));
envsTrain = cell(size(ftrain));
envsTest = cell(size(ftest));
for i = 1:length(ftrain)
   mtrain{i} = Model3(ftrain{i},ftrain{i},ftrain{i});
   mtrain{i}.solveHF;
   envsTrain{i} = envs1;
end
for i = 1:length(ftest)
   mtest{i} = Model3(ftest{i},ftest{i},ftest{i});
   mtest{i}.solveHF;
   envsTest{i} = envs2;
end
load('datasets\ethanerDat.mat');
ftrain = LL(1:10,1);
ftest = LL(11:19,1);
for i = 1:length(ftrain)
   mtrain{end+1} = Model3(ftrain{i},ftrain{i},ftrain{i});
   mtrain{end}.solveHF;
   envsTrain{end+1} = envs1;
end
for i = 1:length(ftest)
   mtest{end+1} = Model3(ftest{i},ftest{i},ftest{i});
   mtest{end}.solveHF;
   envsTest{end+1} = envs2;
end
save('tmp1.mat','mtrain','envsTrain','mtest','envsTest');
%%
for i = 1:length(ftest)
   mtest{i} = Model3(ftest{i},ftest{i},ftest{i});
   mtrain{i}.solveHF;
   envsTest{i} = envs2;
end
save('tmp1.mat');
%%
data = c1.data;
% data: <obs|orig>  coeff: <orig|pca>  score: <obs|pca>
%    <obs|pca> = <obs|orig><orig|pca>  or score = data*coeff
data = zscore(data);
[coeff,score,latent] = princomp(data);
test = max(max(abs(score- data*coeff)));
figure(100);
plot(latent,'rx');
figure(101);
plot(cumsum(latent)./sum(latent));

%%
Context.fillInContexts(mtrain,envsTrain,mtest,envsTest);

%% Verification of combined contexts
clear all;
load('C:\matdl\yaron\10-16-12\context-rapid\ch4r-cross1\all-1.mat');
x1 = f1; % old rapid context
load('C:\matdl\yaron\10-17-12\contextPCA\ch4r\all-3.mat');
x2 = f1; % shoudl be identical to x1
% includeAdhoc = 1;
% mtrain = cell(0,0);
% HLtrain = cell(0,0);
% envsTrain = cell(0,0);
% mtest = cell(0,0);
% HLtest = cell(0,0);
% envsTest = cell(0,0);
% files = {'datasets/ch4rDat.mat'};
% for i1 = 1:length(files)
%    load(files{i1});
%    train = 1:10;
%    test = 11:20;
%    envs1 = [6     7     8    13    16    24];
%    envs2 = [5    10    14    17    20    25];
%    for i = train
%       mtrain{end+1} = Model3(LL{i,1},LL{i,1},LL{i,1});
%       mtrain{end}.solveHF;
%       HLtrain{end+1} = HL{i,1};
%       envsTrain{1,end+1} = envs1;
%    end
%    for i = test
%       mtest{end+1} = Model3(LL{i,1},LL{i,1},LL{i,2});
%       mtest{end}.solveHF;
%       HLtest{end+1} = HL{i,1};
%       envsTest{end+1} = envs2;
%    end
% end
% [x2 ftest] = Context.makeFitme(mtrain,envsTrain,HLtrain, ...
%    mtest,envsTest,HLtest,includeAdhoc);
% x2.silent = 1;
% x2.parallel = 1;
%%
diff1 = [];
for imod = 1:x1.nmodels
   mod1 = x1.models{imod};
   mod2 = x2.models{imod};
   for iatom = 1:mod1.natom
      for ienv = x1.envs{imod};
         ac1 = mod1.atomContext(iatom,ienv);
         ac2 = mod2.atomContext(iatom,ienv);
         diff1(end+1) = max(abs(ac1-ac2(1:3)));
      end
   end
end
disp(['diff in atom contexts ',num2str(max(diff1))]);
%%
diff1 = [];
for imod = 1:x1.nmodels
   for iatom = 1:mod1.natom
      for jatom = 1:mod1.natom
         if (mod1.isBonded(iatom,jatom))
            for ienv = x1.envs{imod};
               ac1 = mod1.bondContext(iatom,jatom,ienv);
               ac2 = mod2.bondContext(iatom,jatom,ienv);
               diff1(end+1) = max(abs(ac1-ac2(1:3)));
%                if (diff1(end) < 0.1)
%                   disp([num2str(diff1(end)), ' for ', num2str(iatom), ...
%                      ' and ',num2str(jatom)]);
%                end
            end
         end
      end
   end
end
disp(['diff in bond contexts ',num2str(max(diff1))]);

%%
% reset mixers to a null state
for i=1:length(x1.mixers)
   mix1 = x1.mixers{i};
   mix2 = x2.mixers{i};
   mix1.par = [1, zeros(1,length(mix1.par)-1)];
   mix2.par = [1, zeros(1,length(mix2.par)-1)];
end
err1 = x1.err(x1.getPars);
err2 = x2.err(x2.getPars);
d1 = err1-err2;
disp(['diff in unmixed energies ',num2str(max(d1))]);

%%
% map the mixers of mix2 to those of mix1
%                1    2    3     4    5    6     7   8     9     10
% order in x2: ke.H en.H e2.H ke.CH en.CH e2.CH ke.C en.C e2.C e2.HH
% order in x1: ke.H ke.C ke.CH en.H en.C en.CH e2.H e2.C e2.HH e2.CH
%                1    7    4     2    8    5     3    9   10     6
m1 = x1.mixers;
imap = [1    7    4     2    8    5     3    9   10     6];
for i = 1:length(imap)
   m2{i} = x2.mixers{imap(i)};
   disp([m1{i}.desc,' --> ',m2{i}.desc]);
end
%%
ic = 0;
res = {};
for imix = 1:length(m1)
   mix1 = m1{imix};
   mix2 = m2{imix};
   for ipar = 1:min(4,length(mix1.par))
      mix1.par(ipar) = mix1.par(ipar) + 0.1;
      err1 = x1.err(x1.getPars);
      mix1.par(ipar) = mix1.par(ipar) - 0.1;
      
      mix2.par(ipar) = mix2.par(ipar) + 0.1;
      err2 = x2.err(x2.getPars);
      mix2.par(ipar) = mix2.par(ipar) - 0.1;
      
      d2 = max(abs(err1-err2));
      ic = ic+1;
      res{ic} = [mix1.desc,' par = ',num2str(ipar),' diff = ',num2str(d2)];
      disp(res{ic});
   end
end
%%
for ic = 1: length(res)
   disp(res{ic});
end

%%
modp{1} = x1.models{1};
modp{2} = x2.models{1};
for im = 1:2
   disp(['for im = ',num2str(im)]);
   kems = modp{im}.KEmods;
   for i= 1:length(kems)
      kem = kems{i};
      disp([num2str(kem.ilist),' to ',num2str(kem.jlist),' with ',kem.mixer.desc]);
   end
end

%%
modp{1} = x1.models{1};
modp{2} = x2.models{1};
for im = 1:2
   disp(['for im = ',num2str(im)]);
   kems = modp{im}.H2mods;
   for i= 1:length(kems)
      kem = kems{i};
      disp([num2str(kem.ilist),' , ',num2str(kem.jlist), ' | '...
         num2str(kem.klist),' , ', num2str(kem.llist),' with ',kem.mixer.desc]);
   end
end
