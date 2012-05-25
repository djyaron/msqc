%%
clear classes;
load('5-3-12/ch4f-c2h6/fit-1/all.mat');
%%
dataDir = ['tmp/',filePre,'/fit-',num2str(iPar),'/'];
ftest = makeFitme('propene',[1 2 3 4 5 6 7],'enstruct',en,'kestruct',ke);
start = ftest.getPars;
%%
ftest.err(start);
%%
limits = [];
options = optimset('DiffMinChange',1.0e-5,'TolFun',1.0e-4,'TolX',1.0e-3);
if (exist(dataDir,'dir') ~= 7)
    status = mkdir(dataDir);
end
diary([dataDir,'out.diary']);
diary on;
tic
start = f1.getPars;
[pt,resnorm,residual,exitflag,output,lambda,jacobian] = ...
    lsqnonlin(@f1.err, start,-limits,limits,options);
clockTime = toc
pt
resnorm