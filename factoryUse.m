clear classes;
close all;

fct = cell(0,0);
ms = cell(0,0);
fitme = cell(0,0);
fname = cell(0,0);
% % dataExt = {'','-1c','-linrho','-diponly'};
% % for iext = 1:4
% %    load(['C:\matdl\yaron\11-29-12\factory\hybrid1\ch4r',dataExt{iext},'\all-3.mat'])
% %    fct{end+1} = fact;  fitme{end+1} = f1;  ms{end+1} = MSet.fromFitme(f1);
% % end
% % fname = {'orig','1c','lin','dip'};
% 
% dataExt = {'ethanerDat','propanerDat'};
% for iext = 1:length(dataExt)
%    load(['C:\matdl\yaron\dec12a\hybridslater\',dataExt{iext},'\all-3.mat'])
%    fct{end+1} = fact;  fitme{end+1} = f1;  ms{end+1} = MSet.fromFitme(f1);
% end
% fname = dataExt;
% 

dataExt = {''};%,'-1c','-linrho','-diponly'};
for iext = 1:length(dataExt)
   load(['C:\matdl\yaron\11-29-12\factory\hybrid1\ch4r',dataExt{iext},'\all-3.mat'])
   for imix = 1:length(f1.mixers)
      f1.mixers{imix}.bonded = 1;
   end
   f1.mixers{end}.bonded = 0;
   fct{end+1} = fact;  fitme{end+1} = f1;  ms{end+1} = MSet.fromFitme(f1);
end
fname = {'orig','1c','lin','dip'};

%%
nfact = length(fct);
ndata = length(fct);
esum = zeros(nfact,ndata);
for ifact = 1:nfact
   for idata = 1:ndata
      disp(['factory: ',fname{ifact},'   data: ',fname{idata}]);
      f2 = fct{ifact}.makeFitme(ms{idata},fitme{idata});
      f2.silent = 1;
%       res1 = f1.printEDetails;
%       disp(' ');
%       esum(ifact,idata) = res1{1}.etot;
      
         [e0 plotnum etype modelnum envnum] = f2.err(f2.getPars);
   f2.printEDetails;
   e0 =e0*627.509;
   eke = e0(etype==1);
   eH  = e0(etype==11);
   eC  = e0(etype==16);
   e2  = e0(etype==2);
   etot = e0(etype==3);
   disp(['mean ke ',num2str(mean(eke)),' H ',num2str(mean(eH)),' C ',...
      num2str(mean(eC)),' E2 ',num2str(mean(e2)),' Etot ', num2str(mean(etot))]);
   disp(['std  ke ',num2str(std(eke)),' H ',num2str(std(eH)),' C ',...
      num2str(std(eC)),' E2 ',num2str(std(e2)),' Etot ', num2str(std(etot))]);

      
   end
end
%%
ofile = 1;
for idata = 1:nfact
   fprintf(ofile,'\t%s ',fname{idata});
end
fprintf(ofile,'\n');

for ifact = 1:nfact
   fprintf(ofile,'%s ',fname{ifact});
   for idata = 1:ndata
      fprintf(ofile,'\t %5.2f ',esum(ifact,idata));
   end
   fprintf(ofile,'\n');
end

%% Error versus weight
clear classes
close all
% Make model set corresponding to methane
mch4 = MSet;
mch4.addData('datasets/ch4rDat.mat',11:20, 2:2:20 ,1,791);

dataroot = 'C:\matdl\yaron\dec12b\hybridslater\ethanerDat\all-3-weight';
lfiles = dir([dataroot,'/*.mat']);
for i = 1:length(lfiles)
   disp(lfiles(i).name);
   load([dataroot,'/',lfiles(i).name]);
   w(i) = f1.operWeights.Etot;
   [a b] = f1.printEDetails;
   etrain{i} = a{1};
   strain{i} = b{1};
   [a b] = ftest.printEDetails;
   etest{i} = a{1};
   stest{i} = b{1};
   fmeth = fact.makeFitme(mch4);
   [a b] = fmeth.printEDetails;
   emeth{i} = a{1};
   smeth{i} = b{1};
end
save([dataroot,'weightsum.mat'],'w','etrain','strain','etest', ...
   'stest','emeth','smeth');
%% 
clear all;
dataroot = 'C:\matdl\yaron\dec12b\hybridslater\ethanerDat\all-3-weight';
load([dataroot,'weightsum.mat']);
for iw = 1:length(w)
   errs{iw,1,1} = etrain{iw};
   errs{iw,1,2} = strain{iw};
   errs{iw,2,1} = etest{iw};
   errs{iw,2,2} = stest{iw};
   errs{iw,3,1} = emeth{iw};
   errs{iw,3,2} = smeth{iw};
end

close all
toplot = {'ke','H','C','e2','etot'};
psym = {'co','bo','b^','k^','ro'};
ltype = {'-','--'};
[ws,is] = sort(w);
etrain = cell(length(ws),1);
for id = 1:size(errs,2) % data set
   for il = 1:2 % err or standard deviation
      for it = 1:length(toplot) % data type
         for i = 1:length(w)
            x(i) = ws(i);
            er = errs{is(i),id,il};
            y(i) = getfield(er,toplot{it});
         end
         figure(70+id);
         subplot(2,1,il)
         hold on;
         plot(x,y,[psym{it},ltype{il}]);
      end
   end
   figure(70+id);
   ylabel('average error');
   xlabel('weight of Etotal');
   legend(toplot);
end


