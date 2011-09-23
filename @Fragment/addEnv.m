function ifile = addEnv(obj, envTarget)

% The environment files will hold H1Env and EorbEnv for a particular
% environment. The file naming configuration is:
%   template_template#_env#
% obj.fileprefix contains template_template#

% first, look to see if a cfg file exists containing this environment
ifile = 0;
found = false;
failed = false;
while (~failed && ~found)
   ifile = ifile + 1;
   fileEnvPrefix = [obj.fileprefix,'\env_',int2str(ifile)];
   cfgfilename = [fileEnvPrefix,'_cfg.mat'];
   try % try to open the file
      load(cfgfilename,'envFile');
   catch
      failed = true;
      %disp('not found');
   end
   if (~failed)
      if (envFile.compare(envTarget))
         found = true;
         %disp('found it');
      end
   end
end
% at this point, calcfilename is either the file with the result
% or the next file that should be written
calcfilename = [fileEnvPrefix,'_calc.mat'];
% Either load envResults from calcfile, or generate it
if (found)
   disp(['found results, loading ',calcfilename]);
   load(calcfilename, 'envResults')
else
   disp(['not found, generating ',calcfilename]);

   % Do the calculation and read in data
   jobname = 'env';
   newline = char(10);
   
   % Need to put charge spec before environment, but gaussianFile has
   % header and environment.
   ctext = strrep(obj.gaussianFile, '!ENV', [envTarget.gaussianText()]);
   gjf_file = [jobname,'.gjf'];
   origdir = cd(obj.dataPath);
   fid1 = fopen(gjf_file,'w');
   fwrite(fid1, ctext, 'char');
   fclose(fid1);
   
   dataPath = obj.dataPath;
   system([obj.gaussianPath,'\',obj.gaussianExe,' ',gjf_file]);
   % convert checkpoint file to a formatted checkpoint file
   system([obj.gaussianPath,'\formchk.exe temp.chk temp.fch']);
   cd(origdir);
   % read in data from formatted checkpoint file
   try
      fid1 = fopen([dataPath,'\temp.fch'],'r');
      if (fid1 == -1)
         error('could not find fch file');
      end
   [Ehfe, Eorbe, orbe, ~,  ~, ~, ...
    dipolee, ~, ~, ~, ...
    ~, ~, ~] = ...
    Fragment.readfchk(fid1);
      fclose(fid1);
   catch
      disp('caught some stupid error');
      fclose(fid1);
      error('failed during env fchk read');
   end
   % read in data from the polyatom output file
   try
      fid1 = fopen([dataPath,'\fort.32'],'r');
      if (fid1 == -1)
         error(['could not find ',dataPath,'\fort.32']);
      end
      [~, H1e, ~, ~, Hnuce] = Fragment.readpolyatom(fid1);
      fclose(fid1);
   catch
      fclose(fid1);
      error('failed during polyatom read');
   end
   
   envResults.H1Env = H1e - obj.H1;
   envResults.Ehf   = Ehfe;
   envResults.Eorb  = Eorbe;
   envResults.orb   = orbe;
   envResults.Hnuc  = Hnuce;
   envResults.dipole = dipolee;
   % cleanup files
%    delete([dataPath,'\fort.32'], [dataPath,'\env.gjf'], ...
%       [dataPath,'\env.out'], [dataPath,'\temp.chk'], ...
%       [dataPath,'\temp.fch']);
   
   % save environment ot the cfg file
   envFile = envTarget;
   save(cfgfilename,'envFile');
   % save envResults to calcfile
   save(calcfilename,'envResults');
end % if (found)

% save data in the object
obj.nenv = obj.nenv + 1;
obj.env(1,obj.nenv) = envTarget;
obj.H1Env(:,:,obj.nenv) = envResults.H1Env;
obj.EhfEnv(1,obj.nenv)  = envResults.Ehf;
obj.EorbEnv(:,obj.nenv) = envResults.Eorb;
obj.HnucEnv(:,obj.nenv) = envResults.Hnuc;
obj.orbEnv(:,:,obj.nenv)  = envResults.orb;
obj.dipoleEnv(:,obj.nenv) = envResults.dipole;


