function [found,fileprefix] = findCalc(dataPath, targetConfig)
% looks to see if a file exists for this Ctarget
% if it does exist, it returns the file number
% if it does not exist, it returns the number that should be used to store
%  the new calculation.
% files begin with [config.template,#]
%  _cfg stores the config structure
%  _calc stores the data

found = 0;
fileprefix = '';
allFiles = dir(dataPath);
nfiles = length(allFiles);
ifile = 1;
while (~found && (ifile <= nfiles))
   fileName = allFiles(ifile).name;
   if (~isempty(regexp(fileName,'_cfg.mat','once')))
      %disp(['found ',fileName]);
      load([dataPath,filesep,fileName],'config');
      if targetConfig.compareConfigs( config )
        found = true;
        fullName = [dataPath,filesep,fileName];
        fileprefix = strrep(fullName,'_cfg.mat','');
      end
   end
   ifile = ifile + 1;   
end