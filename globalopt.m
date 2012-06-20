%% Fit H2 data with some environments as test and some as train
% if fastFit = 1, it does just hydrogen (h2), so it will run fast
% if fastFit = 0, it does h2 and ch4, so it will run slow, but be 
%                 more interesting data
clear classes
for calcType = 1:3
close all
fastFit = 0;  

%  'h2',[2 4 6]   means use the second fourth and sixth geometries of h2
%  'ch4',[1 3]    means the 1 and 3 geometry of methane
%  'env', 1:10    means the 1 through 10 environments around these
%                 molecules
if (fastFit)
   ftest = makeFitme('h2',[2 4 6],'env',1:5,'plot',0);
else
   ftest = makeFitme('h2',[2 4 6],'ch4',[1 3],'env',1:10,'plot',0);
end
%  ftest.npar      is the number of parameters in the model
%  ftest.err(par)  is a function that returns the disagreement
%                  between the model and the data, for the parameters par
%                  returns a long list of numbers that we want to make
%                  as small as possible. 

% We want to be able to fit the parameters to one set of data and test
% the fit on another set of set of data. 
% What this does is sets the TRAIN data to h2, geometries 3 5 7
%   ch4, geometries 2 4, and environments 20 through 30
% The test data is what ever we put into ftest above. 
% (the arguments 'testFitme', ftest, say us ftest as the training data
if (fastFit)
   f1 = makeFitme('h2',[3 5 7],'env',20:25);%, 'testFitme',ftest);
else
   f1 = makeFitme('h2',[3 5 7],'ch4',[2 4],'env',20:30, 'testFitme',ftest);
end
% The "least squares" problem is to try to find the parameters that make
% the sum of the squares of the error as small as possible. 
% To solve this problem, we are using the lsqnonlin function of matlab.

%% Levenberg Marquadt
if (calcType == 1)
% if wanted to limit the range of the parameters (like say, every value in 
% par to stay between -100 and +100, we could do it by setting limits here
% By setting limits to [], we are saying that we aren't setting any limits
% so the paramters can take on any value they want.
limits = [];
% This is how we set parameters for the optimization. We pass them as
% string, value pairs, (just like above for makefitme)
%  DiffMinChange = value it will add to each paramters to see what the
%                  effect of that parameter is on the fit
%  TolFun        = stop the fit if the change in error is smaller than this
%                  value (this is one of the convergence criteria)
%  TolX          = stop if the guess at the paramters changes by less than
%                  this value
options = optimset('DiffMinChange',1.0e-5,'TolFun',1.0e-5,'TolX',1.0e-5, ...
   'MaxFunEvals',200);
% we have to start at some guess for the parameters. f1.getPars gets the
% parameters from the model (which will start at zero)
start = f1.getPars;
% [pt resnorm etc] are the return values from lsqnonlin.
%    pt = the best guess at the parameters. These are the parameters that
%         give the lowest error for the TRAIN set. The TEST set is
%         is calculated and displayed, but not used to determine the
%         parameters.
dataDir = ['tmp/global/LM/'];
if (exist(dataDir,'dir') ~= 7)
   status = mkdir(dataDir);
end
diary([dataDir,'out.diary']);
diary on;
tic
[pt,resnorm,residual,exitflag,output,lambda,jacobian] = ...
   lsqnonlin(@f1.err, start,-limits,limits,options);
clockTime = toc
pt
resnorm
f1.printMixers;
save([dataDir,'all.mat']);
diary off;
figure(799); saveas(gcf,[dataDir,'error.fig']);

% figure 799 shows the error in the test (red plusses) and train (blue
% circles) as a function of the iteration.
end
%% To get the error of the train and test data, we do:
% trainError = f1.errTrain;
% testError = f1.errTest;
% 
% figure(100)
% plot(trainError,'bo');
% hold on;
% plot(testError,'r+');
%% Genetic algorithm
if (calcType == 2)
dataDir = ['tmp/global/ga/'];
if (exist(dataDir,'dir') ~= 7)
   status = mkdir(dataDir);
end
diary([dataDir,'out.diary']);
diary on;
tic
options = gaoptimset(@ga);
options.PopInitRange = [-7; 7];
options.Generations = 200;
[x fval] = ga(@f1.normErr, f1.npar, [],[],[],[],[],[],[],options);
clockTime = toc
x
fval
f1.printMixers;
save([dataDir,'all.mat']);
diary off;
figure(799); saveas(gcf,[dataDir,'error.fig']);
end
%% Simulated annealing
if (calcType == 3)
dataDir = ['tmp/global/sa/'];
if (exist(dataDir,'dir') ~= 7)
   status = mkdir(dataDir);
end
diary([dataDir,'out.diary']);
diary on;
tic

x0 = zeros(f1.npar,1);
lb = -10 * ones(f1.npar,1);
ub = 10 * ones(f1.npar,1);
[x fval] = simulannealbnd(@f1.normErr,x0,lb,ub);
clockTime = toc
x
fval
f1.printMixers;
save([dataDir,'all.mat']);
diary off;
figure(799); saveas(gcf,[dataDir,'error.fig']);

end
%%
%x0 = zeros(f1.npar,1);
%x = fminsearch(@f1.normErr, x0);
end