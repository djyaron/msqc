function [found,fileprefix] = findCalc(dataPath,Ctarget)
% looks to see if a file exists for this Ctarget
% if it does exist, it returns the file number
% if it does not exist, it returns the number that should be used to store
%  the new calculation.
% files begin with [config.template,#]
%  _cfg stores the config structure
%  _calc stores the data


ifile = 0;
found = false;
failed = false;
while (~failed && ~found)
   ifile = ifile + 1;
   fileprefix = [dataPath,filesep,Ctarget.template,'_',int2str(ifile)];
   %disp(['looking for ',fileprefix]);
   cfgfilename = [fileprefix,'_cfg.mat'];
   try % try to open the file
      load(cfgfilename,'Cfile');
   catch
      failed = true;
      %disp('not found');
   end
   if (~failed)
      if (size(comp_struct(Cfile,Ctarget),1) == 0)
         found = true;
         %disp('found it');
      end
   end
end
