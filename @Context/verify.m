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



