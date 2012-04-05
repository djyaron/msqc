clear classes;

%% H2 fit
f1 = makeFitme;
disp('Starting fit for h2');
start = zeros(1,f1.npar);
%limits = [4 20 4 20 4 20 4 4 4 20 4 20 4 20 4 4 4 4];
limits = [];
options = optimset('DiffMinChange',1.0e-5);
pt = lsqnonlin(@f1.err, start,-limits,limits,options);

%% CH4 and ethane fit
f1 = makeFitme('ch4',1:7,'ethane',1:7);
disp('Starting fit ch4 and ethane');
start = zeros(1,f1.npar);
%limits = [4 20 4 20 4 20 4 4 4 20 4 20 4 20 4 4 4 4];
limits = [];
options = optimset('DiffMinChange',1.0e-5);
pt = lsqnonlin(@f1.err, start,-limits,limits,options);

   
%    options = LMFnlsq;
%    options.Display =1;
%    options.FunTol = 1.0e-6;
%    options.XTol = 1.0e-5;
%    [pfit, Ssq, CNT, Res, XY] = LMFnlsq(@f1.err,start',options);
