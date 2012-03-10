clear classes;
reload = 1;
plotCorrelations = 0;
modelType = 'noCharge';
handFit = 0;
doFit = 1;
plotResults = 1;
useStart = 1;
pstart =  [-0.5 3 7];
params = 2:7;

if (reload)
   load('h2/h2Dat.mat');
end

if (plotCorrelations)
   chl=zeros(2,100*size(params,2));
   cll=zeros(2,100*size(params,2));
   ic = 0;
   for ipar = params
      ll = HL{ipar,1};
      hl = HL{ipar,2};
      t1 = ll.EKE;
      rr = (ic+1):(ic+size(t1,2));
      lke(ic+1:ic+size(t1,2)) = t1;
      hke(ic+1:ic+size(t1,2)) = hl.EKE;
      le1(ic+1:ic+size(t1,2)) = ll.Een(1);
      le2(ic+1:ic+size(t1,2)) = ll.Een(2);
      he1(ic+1:ic+size(t1,2)) = hl.Een(1);
      he2(ic+1:ic+size(t1,2)) = hl.Een(2);
      c1 = zeros(2,100);
      c2 = zeros(2,100);
      for ienv = 1:hl.nenv
         mh = hl.mcharge(ienv);
         ml = ll.mcharge(ienv);
         ch(:,ienv) = mh(1,:);
         cl(:,ienv) = ml(1,:);
      end
      chl(:,rr) = ch;
      cll(:,rr) = cl;
      ic = ic + size(t1,2);
   end
   CL = cll(1,:);
   CH = chl(1,:);
   figure(1)
   hold off
   plot(CL,lke,'b.');
   hold on
   plot(CL,hke,'r.');
   title('ke');
   figure(2)
   hold off
   plot(CL,le1,'b.');
   hold on
   plot(CL,he1,'r.');
   title('en1');
   figure(3)
   hold on
   plot(CL,hke-lke,'g.');
   title('ke diff');
   figure(4)
   hold on
   plot(CL,he2-le1,'g.');
end

if (doFit || handFit)
   disp('building models');
   mixKEdiag = Mixer([0 0],2);
   mixKEbond = Mixer(0,1);
   m = cell(1,size(params,2));
   for ipar = params
      m{ipar} = Model3(LL{ipar,1},LL{ipar,2},LL{ipar,3});
      m{ipar}.addKEmodDiag(1,1,mixKEdiag);
      m{ipar}.addKEmodBonded(1,1,1,1,mixKEbond);
   end
end

if (handFit)
   ic = 0;
   for ipar = params
      ll = LL{ipar,1};
      hl = HL{ipar,1};
      t1 = ll.EKE;
      rr = (ic+1):(ic+size(t1,2));
      lke(rr) = t1;
      hke(rr) = hl.EKE;
      ic = ic + size(t1,2);
   end
   while 1
      p(1) = input('diag const ');
      p(2) = input('diag charge ');
      p(3) = input('bond const ');
      ic = 0;
      for ipar = params
         disp(['starting calc on ipar ',num2str(ipar)]);
         m{ipar}.setPar(p);
         m{ipar}.solveHF;
         t1 = m{ipar}.EKE;
         mke(ic+1:ic+size(t1,2)) = t1;
         ic = ic + size(t1,2);
      end
      figure(100);
      hold off;
      plot(lke,hke,'r.');
      hold on;
      plot(lke,lke,'k.');
      plot(lke,mke,'b.');
      %figure(200);
      %plot(mke,hke,'g.');
   end
end

if (doFit)
   disp('Starting to do parameter fitting');
   f1 = Fitme;
   for ipar = params
      f1.addFrag(m{ipar},HL{ipar,1});
   end
   f1.exactDensity = 1;
   nn = params(1);
   f1.includeEN = zeros(1,m{nn}.natom);
   
   nfitpar = m{nn}.npar;
   start = zeros(1,nfitpar);
   if (useStart)
      start = pstart;
   end
   % code from using matlabs optimization toolbox
   %limits = 3 * ones(1,nfitpar);
   %options = optimset('DiffMinChange',1.0e-5);
   %pt = lsqnonlin(@f1.err, start,-limits,limits,options);
   options = LMFnlsq;
   options.Display =1;
   options.FunTol = 1.0e-6;
   options.XTol = 1.0e-5;
   [pfit, Ssq, CNT, Res, XY] = LMFnlsq(@f1.err,start',options);
end

if (plotResults)
   disp('Starting to do plots');
   if (doFit)
      pt = pfit;
   else
      pt = pstart;
   end
   ic = 0;
   for ipar = params
      ll = LL{ipar,1};
      hl = HL{ipar,1};
      t1 = ll.EKE;
      rr = (ic+1):(ic+size(t1,2));
      lke(rr) = t1;
      hke(rr) = hl.EKE;
      disp(['starting calc on ipar ',num2str(ipar)]);
      m{ipar}.setPar(pt);
      m{ipar}.solveHF;
      t1 = m{ipar}.EKE;
      mke(rr) = t1;
      ic = ic + size(t1,2);
   end
   figure(100);
   hold off;
   plot(lke,hke,'r.');
   hold on;
   plot(lke,lke,'k.');
   plot(lke,mke,'b.');
   figure(200);
   plot(mke,hke,'g.');
end
