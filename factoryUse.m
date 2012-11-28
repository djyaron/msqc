clear classes;
close all;

fct = cell(0,0);
ms = cell(0,0);
fitme = cell(0,0);
fname = cell(0,0);
dataExt = {'','-1c','-linrho','-diponly'};
for iext = 1:4
   load(['C:\matdl\yaron\11-29-12\factory\hybrid1\ch4r',dataExt{iext},'\all-3.mat'])
   fct{end+1} = fact;  fitme{end+1} = f1;  ms{end+1} = MSet.fromFitme(f1);
end
fname = {'orig','1c','lin','dip'};

nfact = length(fct);
ndata = length(fct);
esum = zeros(nfact,ndata);
for ifact = 1:nfact
   for idata = 1:ndata
      disp(['factory: ',fname{ifact},'   data: ',fname{idata}]);
      [f1 c1] = fct{ifact}.makeFitme(ms{idata},fitme{idata});
      res1 = f1.printEDetails;
      disp(' ');
      esum(ifact,idata) = res1{1}.etot;
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
