function loadAllEnv(obj)

% figure out how many enviroment files are available
ifile = 0;
failed = false;
while (~failed)
   ifile = ifile + 1;
   fileEnvPrefix = [obj.fileprefix,'\env_',int2str(ifile)];
   cfgfilename = [fileEnvPrefix,'_cfg.mat'];
   try % try to open the file
      load(cfgfilename,'envFile');
   catch
      failed = true;
   end
end

nfile = ifile-1;

if (nfile > 0)
   obj.setEnvSize(nfile);
   
   ifile = 0;
   failed = false;
   while ((~failed) && (ifile < nfile))
      ifile = ifile + 1;
      disp(['loading ',num2str(ifile)]);
      fileEnvPrefix = [obj.fileprefix,'_',int2str(ifile)];
      cfgfilename = [fileEnvPrefix,'_cfg.mat'];
      calcfilename = [fileEnvPrefix,'_calc.mat'];
      try % try to open the file
         load(cfgfilename,'envFile');
         load(calcfilename,'envResults');
      catch
         failed = true;
      end
      if (~failed)
         obj.nenv = obj.nenv + 1;
         disp(['loading env ',num2str(obj.nenv)]);
         obj.env(1,obj.nenv) = envFile;
         obj.H1Env(:,:,obj.nenv) = envResults.H1Env;
         obj.EhfEnv(1,obj.nenv)  = envResults.Ehf;
         obj.MP2Env(1,obj.nenv)  = envResults.MP2;
         obj.CorrEenv(1,obj.nenv)= envResults.CorrE;
         obj.EorbEnv(:,obj.nenv) = envResults.Eorb;
         obj.orbEnv(:,:,obj.nenv)= envResults.orb;
         obj.HnucEnv(:,obj.nenv) = envResults.Hnuc;
      end
   end
   
end