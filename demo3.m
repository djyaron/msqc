f1 = makeFitme;
disp('Starting fit');
start = zeros(1,f1.npar);
   limits = [4 20 4 20 4 20 4 4 4 20 4 20 4 20 4 4 4 4];
   options = optimset('DiffMinChange',1.0e-5);
   pt = lsqnonlin(@f1.err, start,-limits,limits,options);

   
%    options = LMFnlsq;
%    options.Display =1;
%    options.FunTol = 1.0e-6;
%    options.XTol = 1.0e-5;
%    [pfit, Ssq, CNT, Res, XY] = LMFnlsq(@f1.err,start',options);
