%% Fitting multiple molecules, using makeFitme
clear classes;
myDir = 'T:\matdl\yaron\6-22-12\scaleconst\';
ftype = 2;
trainC{1}  = {'h2',[],'ch4',1:17,'envs',1:10};
testC{1} = {'h2',[],'ch4',1:17,'envs',20:30};
filePrefix{1} = 'ch4';

trainC{2}  = {'h2',[],'ethane',1:7,'envs',1:10};
testC{2} = {'h2',[],'ethane',1:7,'envs',20:30};
filePrefix{2} = 'c2h6';

trainC{3}  = {'h2',[],'ethylene',1:7,'envs',1:10};
testC{3} = {'h2',[],'ethylene',1:7,'envs',20:30};
filePrefix{3} = 'c2h4z2';

trainC{4}  = {'h2',[],'ch4',1:7,'ethane',1:7,'envs',1:10};
testC{4} = {'h2',[],'ch4',1:7,'ethane',1:7,'envs',20:30};
filePrefix{4} = 'ch4-c2h6';

trainC{5}  = {'h2',[],'ch4',1:7,'ethane',1:7,'ethylene',1:7,'envs',1:10};
testC{5} = {'h2',[],'ch4',1:7,'ethane',1:7,'ethylene',1:7,'envs',20:30};
filePrefix{5} = 'ch4-c2h6-c2h4';

trainC{6}  = {'h2',[],'ch4',1:19,'ethane',1:7,'envs',1:10};
testC{6} = {'h2',[],'ch4',1:19,'ethane',1:7,'envs',20:30};
filePrefix{6} = 'ch4f-c2h6';

trainC{7}  = {'h2',[],'ch4',1:19,'ethane',1:7,'ethylene',1:7,'envs',1:10};
testC{7} = {'h2',[],'ch4',1:19,'ethane',1:7,'ethylene',1:7,'envs',20:30};
filePrefix{7} = 'ch4f-c2h6-c2h4';

commonIn = {};
%
ofile = fopen('multi6-sum.txt','w');
molName = {'CH3' 'C2H6' 'C2H4'};
for iC = [1 2 3 4 6 7]
   for iPar = 1:3
      filePre = filePrefix{iC};
      dataDir = [myDir,filePre,'/fit-',num2str(iPar),'/'];
      load([dataDir,'all.mat']);
      % [pt,resnorm,residual,exitflag,output,lambda,jacobian] = ...
      %  lsqnonlin(@f1.err, start,lowLimits,highLimits,options);
      [err pnum etype] = f1.err(pt);
      err =err*627.509;
      eke = err(etype==1);
      eH  = err(etype==11);
      eC  = err(etype==16);
      e2  = err(etype==2);
      tot = norm(err)/sqrt(length(err));
      ke = norm(eke)/sqrt(length(eke));
      H = norm(eH)/sqrt(length(eH));
      C = norm(eC)/sqrt(length(eC));
      r2 = norm(e2)/sqrt(length(e2));
      fprintf(ofile,'%s fit %i  ',filePre,iPar);
      fprintf(ofile,'tot %4.2f ke %4.2f H %4.2f C %4.2f E2 %4.2f \n',...
         tot, ke, H, C, r2);
      if (length(unique(pnum)) > 1)
         for ip = unique(pnum)
            etot = err(pnum==ip);
            eke = err(etype==1 & pnum==ip);
            eH  = err(etype==11 & pnum==ip);
            eC  = err(etype==16 & pnum==ip);
            e2  = err(etype==2 & pnum==ip);
            tot = norm(etot)/sqrt(length(etot));
            ke = norm(eke)/sqrt(length(eke));
            H = norm(eH)/sqrt(length(eH));
            C = norm(eC)/sqrt(length(eC));
            r2 = norm(e2)/sqrt(length(e2));
            fprintf(ofile,'>> %s ',molName{ip-800});
            fprintf(ofile,'tot %4.2f ke %4.2f H %4.2f C %4.2f E2 %4.2f \n',...
               tot, ke, H, C, r2);
         end
      end
   end
   fprintf(ofile,'\n');
end
fclose(ofile);

%%


