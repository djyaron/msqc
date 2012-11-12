commonDirName = 'T:\matdl\yaron\8-3-12\quadratic-sp\ch4-geom';
% bond lengths are:
% pars{2} = [0.97 0.97 0.97 0.97 109.47 109.47 109.47 120.0 -120.0];
% pars{20} = [1.02 1.02 1.02 1.02 109.47 109.47 109.47 120.0 -120.0];
% pars{21} = [1.07 1.07 1.07 1.07 109.47 109.47 109.47 120.0 -120.0];
% pars{1} = [1.12 1.12 1.12 1.12 109.47 109.47 109.47 120.0 -120.0];
% pars{22} = [1.17 1.17 1.17 1.17 109.47 109.47 109.47 120.0 -120.0];
% pars{23} = [1.22 1.22 1.22 1.22 109.47 109.47 109.47 120.0 -120.0];
% pars{3} = [1.27 1.27 1.27 1.27 109.47 109.47 109.47 120.0 -120.0];

extensions = {'2','20','21','1','22','23','3'};
postName = '\fit-9\';

for i = 1:length(extensions)
   dirName = [commonDirName,extensions{i},postName];
   disp(dirName);
   allFile = [dirName,'all.mat'];
   if (exist(allFile,'file'))
      load(allFile);
      for j = 1:length(f1.mixers)
         m{i,j} = f1.mixers{j};
      end
   else
      display(['could not find ',allFile]);
   end
end
%%
close all;
outfile = 1;
leg = {};
col = {'ko','bo','ro','co','go','yo','kx','bx','rx','cx','gx','yx'};
for j = 1:size(m,2)
   fprintf('%s ',m{1,j}.desc);
   for k = 1:length(m{1,j}.par)
      x = [];
      leg{end+1} = [m{1,j}.desc,' ',num2str(k)];
      for i = 1:size(m,1)
         fprintf('%6.3f ',m{i,j}.par(k));
         x(i) = m{i,j}.par(k);
      end
      figure(100 + k);
      hold on;
      plot(1:length(x),x,[col{j},'-']);
      fprintf('\n');
   end
end
figure(101);
legend(leg);