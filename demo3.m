% Only really needed if changes were made to class files
clear classes;

%% H2 fit
f1 = makeFitme;
disp('Starting fit for h2');
start = zeros(1,f1.npar);
%limits = [4 20 4 20 4 20 4 4 4 20 4 20 4 20 4 4 4 4];
limits = [];
options = optimset('DiffMinChange',1.0e-5);
diary('tmp/h2.diary');
diary on;
pt = lsqnonlin(@f1.err, start,-limits,limits,options);
pt
diary off;
% Should give pt = -0.5127    3.7062    7.6738   -0.2253    4.1360    4.3110

%% CH4 fits
clear classes;
f1 = makeFitme('h2',[],'ch4',1:7);
disp('Starting fit for ch4');
start = zeros(1,f1.npar);
diary('tmp/ch4.diary');
diary on;
limits = [];
options = optimset('DiffMinChange',1.0e-5);
tic
pt = lsqnonlin(@f1.err, start,-limits,limits,options);
clockTime = toc
pt
diary off;
% with ch4 1:7 this gives
% pt =
%    12.9908   17.3272  -14.3800  -23.6141   -0.3528   14.7917   10.9971   16.6052    1.8722
%     9.2113   -0.3798   -4.0887   -2.3191    7.6619    0.8430    2.3617
% RMS err/ndata = 0.00036669


%% H2, CH4 and ethane fit
f1 = makeFitme('ch4',1:7,'ethane',1:7);
disp('Starting fit ch4 and ethane');
start = zeros(1,f1.npar);
%limits = [4 20 4 20 4 20 4 4 4 20 4 20 4 20 4 4 4 4];
limits = [];
options = optimset('DiffMinChange',1.0e-5);
pt = lsqnonlin(@f1.err, start,-limits,limits,options);

%  Leads to:
% Fitme.err called with par = -0.397823      7.11072     -14.2068     -40.7873      3.98352       17.195      6.65837      1.50171    -0.143203      6.99015      -2.3632     -5.03394      2.61736      8.90068     0.925762     0.953655      7.29745       4.1667    -0.806963     0.310333
% solving for density matrices
% RMS err/ndata = 0.00058134
% 
% Solver stopped prematurely.
% 
% lsqnonlin stopped because it exceeded the function evaluation limit,
%options.MaxFunEvals = 2000 (the default value).


%% Code for using the LMFnlsq optimizer
%    options = LMFnlsq;
%    options.Display =1;
%    options.FunTol = 1.0e-6;
%    options.XTol = 1.0e-5;
%    [pfit, Ssq, CNT, Res, XY] = LMFnlsq(@f1.err,start',options);
