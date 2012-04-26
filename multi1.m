clear classes;
reload = 1;
nhl = 1;
plotCorrelations = 0;
includeKEmods = 1;
includeENmods = 1;
useDeltaCharges = 1;
debugModel = 0;
handFit = 0;
doFit = 1;
plotResults = 0;
useStart = 0;
%pstart =  [7.2274 9.0707 -22.1916 -14.9135 4.2480  6.6412 22.2263 14.5396  1.2623  5.9730  -2.8560 -2.5580 1.2478 3.4480 2.0505 2.4942 0 0];
%methane and ethane fits
pstart = [0.517818  8.1441 2.13829 6.56335 0.464458 9.73133 -10.4231 ...
   7.37763  -0.0349262 7.59848 1.07518 1.11226 0.209319 4.52014  ...
   -0.594433  1.94629  0.683698  0.890883];
envs = 0:20; % environments to include in fit
geomsH2 = []; %2:7;
geomsCH4 = 1:7;% 1:3;
geomsEthane =1:7;% [];
plotNumber = [];

if (reload)
   load('h2/h2Dat.mat');
   LLh2 = LL;
   HLh2 = HL;
   load('ch4/ch4Dat.mat');
   LLch4 = LL;
   HLch4 = HL;
   load('ethane4/ethaneDat.mat');
   LLeth = LL;
   HLeth = HL;
   LL = cell(0,0);
   HL = cell(0,0);
   ic = 0;
   for i = geomsCH4
      ic = ic+1;
      plotNumber(1,ic) = 801;
      for j = 1:size(LLch4,2)
         LL{ic,j} = LLch4{i,j};
      end
      for j = 1:size(HLch4,2)
         HL{ic,j} = HLch4{i,j};
      end
   end
   for i = geomsH2
      ic = ic+1;
      plotNumber(1,ic) = 800;
      for j = 1:size(LLh2,2)
         LL{ic,j} = LLh2{i,j};
      end
      for j = 1:size(HLh2,2)
         HL{ic,j} = HLh2{i,j};
      end
   end
   for i = geomsEthane
      ic = ic+1;
      plotNumber(1,ic) = 802;
      for j = 1:size(LLeth,2)
         LL{ic,j} = LLeth{i,j};
      end
      for j = 1:size(HLeth,2)
         HL{ic,j} = HLeth{i,j};
      end
   end
   params = 1:ic;
end

% build models
if (doFit || handFit || debugModel)
   disp('building models');
   m = cell(1,size(params,2));
   for ipar = params
      m{ipar} = Model3(LL{ipar,1},LL{ipar,2},LL{ipar,3});
   end
   if (includeKEmods)
      mixKEdiagH = Mixer([0 0],2,'KEdiagH');
      mixKEdiagC = Mixer([0 0],2,'KEdiagC');
      mixKEdiagCp = Mixer([0 0],2,'KEdiagCp');
      mixKEbondHH = Mixer(0,1,'KEbondHH');
      mixKEbondCH  = Mixer([0 0],3,'KEbondCH');
      mixKEbondCHp  = Mixer([0 0],3,'KEbondCHp');
      mixKEbondCC  = Mixer([0 0],3,'KEbondCC');
      for ipar = params
         m{ipar}.addKEmodDiag(1,1,mixKEdiagH);
         m{ipar}.addKEmodDiag(6,1,mixKEdiagC);
         m{ipar}.addKEmodDiag(6,2,mixKEdiagCp);
         m{ipar}.addKEmodBonded(1,1,1,1,mixKEbondHH);
         m{ipar}.addKEmodBonded(1,6,1,1,mixKEbondCH);
         m{ipar}.addKEmodBonded(1,6,1,2,mixKEbondCHp);
         m{ipar}.addKEmodBonded(6,6,[1 2],[1 2],mixKEbondCC);         
      end
   end
   if (includeENmods)
      mixENdiagH = Mixer([0 0],2,'ENdiagH');
      mixENdiagC = Mixer([0 0],2,'ENdiagC');
      mixENdiagCp = Mixer([0 0],2,'ENdiagCp');
      mixENbondHH = Mixer(0,1,'ENbondHH');
      mixENbondCH  = Mixer([0 0],3,'ENbondCH');
      mixENbondCHp  = Mixer([0 0],3,'ENbondCHp');
      mixENbondCC  = Mixer([0 0],3,'ENbondCC');
      for ipar = params
         m{ipar}.addENmodDiag(1,1,mixENdiagH);
         m{ipar}.addENmodDiag(6,1,mixENdiagC);
         m{ipar}.addENmodDiag(6,2,mixENdiagCp);
         m{ipar}.addENmodBonded(1,1,1,1,mixENbondHH);
         m{ipar}.addENmodBonded(1,6,1,1,mixENbondCH);
         m{ipar}.addENmodBonded(1,6,1,2,mixENbondCHp);
         m{ipar}.addENmodBonded(6,6,[1 2],[1 2],mixENbondCC);
      end
   end
   if (useDeltaCharges)
      for ipar = params
         for ienv = 1:m{ipar}.nenv
            m{ipar}.charges(:,ienv+1) = m{ipar}.charges(:,ienv+1) - m{ipar}.charges(:,1);
         end
          m{ipar}.charges(:,1) = m{ipar}.charges(:,1) - m{ipar}.charges(:,1);
      end
   end
end

if (debugModel)
   input junk;
   f1 = Fitme;
   for ipar = params
      f1.addFrag(m{ipar},HL{ipar,nhl});
   end
   f1.includeKE = includeKEmods;
   f1.includeEN = includeENmods * ones(1,6);
   f1.setEnvs(envs);
   npar = f1.npar;
   p = zeros(1,npar);
   grad = zeros(1,npar);
   y0 = f1.err(p);
   eps = 1.0e-3;
   for i=1:npar
      p1 = p;
      p1(i) = p1(i) + eps;
      y1 = f1.err(p1);
      grad(i) = norm(y1-y0)/eps;
      %m{1}.printMixers;
   end
end

if (handFit)
   ic = 0;
   for ipar = params
      ll = LL{ipar,1};
      hl = HL{ipar,nhl};
      t1 = ll.EKE;
      rr = (ic+1):(ic+size(t1,2));
      lke(rr) = t1;
      hke(rr) = hl.EKE;
      ic = ic + size(t1,2);
   end
   while 1
      
      p(1) = input('KE diag const H  ');
       p(2) = input('KE diag charge H ');
       p(3) = input('KE diag const C  ');
       p(4) = input('KE diag charge C ');
       p(5) = input('KE bond const HH ');
       p(6) = input('KE bond const CC ');
      %p(4) = input('EN diag const ');
      %p(5) = input('EN diag charge ');
      %p(6) = input('EN bond const ');
      ic = 0;
      for ipar = params
         disp(['starting calc on ipar ',num2str(ipar)]);
         m{ipar}.setPars(p);
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
      f1.addFrag(m{ipar},HL{ipar,nhl},plotNumber(ipar));
   end
   f1.includeKE = includeKEmods;
   f1.includeEN = includeENmods * ones(1,6);
   f1.setEnvs(envs);
   nfitpar = f1.npar;
   start = zeros(1,nfitpar);
   if (useStart)
      start = pstart;
   end
   % code from using matlabs optimization toolbox
   limits = [4 10 4 10 4 10 4 4 4 10 4 10 4 10 4 4 4 4];
   options = optimset('DiffMinChange',1.0e-5);
   pt = lsqnonlin(@f1.err, start,-limits,limits,options);
%    options = LMFnlsq;
%    options.Display =1;
%    options.FunTol = 1.0e-6;
%    options.XTol = 1.0e-5;
%    [pfit, Ssq, CNT, Res, XY] = LMFnlsq(@f1.err,start',options);
end
%%
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
      hl = HL{ipar,nhl};
      t1 = ll.EKE;
      rr = (ic+1):(ic+size(t1,2));
      lke(rr) = t1;
      hke(rr) = hl.EKE;
      le1(rr) = ll.Een(1);
      he1(rr) = hl.Een(1);
      disp(['starting calc on ipar ',num2str(ipar)]);
      m{ipar}.setPars(pt);
      m{ipar}.solveHF;
      mke(rr) = m{ipar}.EKE;
      me1(rr) = m{ipar}.Een(1);
      ic = ic + size(t1,2);
   end
   figure(100);
   hold off;
   plot(lke,hke,'r.');
   hold on;
   plot(lke,lke,'k.');
   plot(lke,mke,'b.');
   title('ke')
   figure(200);
   plot(mke,hke,'g.');
   title('ke')
   
   figure(101);
   hold off;
   plot(le1,he1,'r.');
   hold on;
   plot(le1,le1,'k.');
   plot(le1,me1,'b.');
   title('EN')
   figure(201);
   plot(me1,he1,'g.');
   title('EN')
   
   % all plots on one screen
   figure(500);
   subplot(2,2,1);
   hold off;
   plot(lke,hke,'r.');
   hold on;
   plot(lke,lke,'k.');
   plot(lke,mke,'b.');
   title('ke')
   subplot(2,2,2);
   hold off;
   plot(mke,hke,'g.');
   title('ke')
   
   subplot(2,2,3);
   hold off;
   plot(le1,he1,'r.');
   hold on;
   plot(le1,le1,'k.');
   plot(le1,me1,'b.');
   title('EN')
   subplot(2,2,4);
   hold off;
   plot(me1,he1,'g.');
   title('EN')

end


