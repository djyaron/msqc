function ifile = addEnv(obj, envTarget)

% The environment files will hold H1Env and EorbEnv for a particular
% environment. The file naming configuration is:
%   template_template#_env#
% obj.fileprefix contains template_template#

% env files are stored in directory obj.fileprefix
% We check this directory for a *_cfg.mat that matches envTarget
% If found, we set found=1 and fileEnvPrefix to the full location of
% the *_cfg.mat file, with the _cfg.mat removed
% For instance, if the environment is in  env222_cfg.mat, then
% fileEnvPrefix is [obj.dataPath,filesep,obj.fileprefix,filesep,'env222'];
found = 0;
dirName = obj.fileprefix;
if (exist(dirName,'dir'))
   allFiles = dir(dirName);
   nfiles = length(allFiles);
   ifile = 1;
   while (~found && (ifile <= nfiles))
      fileName = allFiles(ifile).name;
      if (~isempty(regexp(fileName,'_cfg.mat','once')))
         %disp(['found ',fileName]);
         load([dirName,filesep,fileName],'envFile');
         if (envFile.compare(envTarget))
            found = true;
            fileEnvPrefix = [dirName,filesep,strrep(fileName,'_cfg.mat','')];
         end
      end
      ifile = ifile + 1;
   end
else
   mkdir(dirName);
end

% Backwards compatibility, load _calc.mat if present
if (found)
   calcfilename = [fileEnvPrefix,'_calc.mat'];
   if (exist(calcfilename,'file'))
      disp(['found env results, loading ',calcfilename]);
      loadMatFileFormat(obj,calcfilename,envTarget);
      return;
   end
end

if (found)
   envZipFileName = [fileEnvPrefix,'.zip'];
   disp(['found env results, loading ',envZipFileName]);
else
   temp1 = tempname('a'); % makes "a\uniquestring"
   uniqueStr = temp1(3:end);
   fileEnvPrefix = [dirName,filesep,'env_',uniqueStr];
   envZipFileName = [fileEnvPrefix,'.zip'];
   disp(['env results not found, generating ',envZipFileName]);
   envFile = envTarget;
   save([fileEnvPrefix,'_cfg.mat'], 'envFile');
   generateEnvZipData(obj,envZipFileName,envTarget);
end

loadEnvZipData(obj,envZipFileName,envTarget);

end

function generateEnvZipData(obj,envZipFileName,envTarget)
tempDir = tempname([obj.gaussianPath,filesep,'Scratch']);
mkdir(tempDir);
% Do the calculation and read in data
setenv('GAUSS_EXEDIR', obj.gaussianPath);
jobname = 'env';
newline = char(10);

ctext = obj.gaussianFile;

% Need to put charge spec before environment, but gaussianFile has
% header and environment.
if envTarget.ncharge > 0
   % add charge keyword, so calcs can be done in an environment
   ctext = strrep(ctext,'symm=noint','symm=noint charge');
   ctext = strrep(ctext, '!ENV', [envTarget.gaussianText()]);
end

% add the fields
% The way that Gaussian does this is really ugly and inconsistant
% There may be a better way to do this.
if envTarget.nfield > 0
   for ifield = 1:envTarget.nfield
      % turn into the format for Gaussian
      tmp = envTarget.fieldType(ifield, :);
      if (sum(tmp) < 4)
         dir = [repmat('X', 1, envTarget.fieldType(ifield, 1)), ...
            repmat('Y', 1, envTarget.fieldType(ifield, 2)), ...
            repmat('Z', 1, envTarget.fieldType(ifield, 3))];
      else
         dir = '';
         while max(tmp) > 0
            switch find(tmp == max(tmp), 1)
               case 1
                  dir = [dir, repmat('X', 1, tmp(1))];
                  tmp(1) = 0;
               case 2
                  dir = [dir, repmat('Y', 1, tmp(2))];
                  tmp(2) = 0;
               case 3
                  dir = [dir, repmat('Z', 1, tmp(3))];
                  tmp(3) = 0;
               otherwise
                  error('Invalid input');
            end
         end
      end
      mag = '';
      if envTarget.fieldMag(ifield) >= 0
         mag = '+';
      end
      mag = [mag, num2str(envTarget.fieldMag(ifield))];
      insert = ['Field=', dir, mag];
      % Add to the input file
      ctext = strrep(ctext, 'symm=noint', ['symm=noint ', insert]);
   end
end
gjf_file = [jobname,'.gjf'];
origdir = cd(tempDir);
fid1 = fopen(gjf_file,'w');
fwrite(fid1, ctext, 'char');
fclose(fid1);

resp = system([obj.gaussianPath,filesep,obj.gaussianExe,' ',gjf_file]);
if resp ~= 0
    disp('Retrying calculation with alternate convergence method.');
    ctext = strrep(ctext,'scf=conventional','scf=(conventional,qc)');
    fid1 = fopen(gjf_file,'w');
    fwrite(fid1, ctext, 'char');
    fclose(fid1);
    system([obj.gaussianPath,filesep,obj.gaussianExe,' ',gjf_file]);
end

toZip = {'temp.chk','fort.32','env.gjf','env.out'};
zip(envZipFileName,toZip);
cd(origdir);
status = rmdir(tempDir,'s');
while status ~= 1
    disp('  rmdir failed. Retrying...');
    pause(0.1);
    status = rmdir(tempDir,'s');
end
end

function loadEnvZipData(obj,envZipFileName,envTarget)

tempDir = tempname([obj.gaussianPath,filesep,'Scratch']);
mkdir(tempDir);
unzip(envZipFileName,tempDir);
origdir = cd(tempDir);
setenv('GAUSS_EXEDIR', obj.gaussianPath);
% convert checkpoint file to a formatted checkpoint file
system([obj.gaussianPath,'\formchk.exe temp.chk temp.fch']);
cd(origdir);
% read in data from formatted checkpoint file
try
   fid1 = fopen([tempDir,filesep,'temp.fch'],'r');
   if (fid1 == -1)
      error('could not find fch file');
   end
   
   [CorrEe, MP2e, Ehfe, Eorbe, orbe, junk,  junk2, junk3, ...
      dipolee, junk4, junk5, junk6, ...
      junk7, junk8, junk9] = ...
      Fragment.readfchk(fid1);
   fclose(fid1);
catch
   fclose(fid1);
   error('failed during env fchk read');
end
% read in data from the polyatom output file
try
   fid1 = fopen([tempDir,filesep,'fort.32'],'r');%,'b');
   if (fid1 == -1)
      error(['could not find env fort.32']);
   end
   [junk1, H1e, junk2, junk3, Hnuce] = Fragment.readpolyatom(fid1);
   fclose(fid1);
catch
   fclose(fid1);
   error('failed during env polyatom read');
end

try
    obj.nenv = obj.nenv + 1;
    obj.env(1,obj.nenv) = envTarget;
    obj.H1Env(:,:,obj.nenv) = H1e-obj.H1;
    obj.MP2Env(1,obj.nenv) = MP2e;
    obj.EhfEnv(1,obj.nenv)  = Ehfe;
    obj.EorbEnv(:,obj.nenv) = Eorbe;
    obj.HnucEnv(:,obj.nenv) = Hnuce;
    obj.orbEnv(:,:,obj.nenv)  = orbe;
    obj.dipoleEnv(:,obj.nenv) = dipolee;
catch
    error('failed to read in data.');
end

status = rmdir(tempDir,'s');
while status ~= 1
    disp('  rmdir failed. Retrying...');
    pause(0.1);
    status = rmdir(tempDir,'s');
end
end

function loadMatFileFormat(obj,calcfilename,envTarget)

load(calcfilename, 'envResults')
obj.nenv = obj.nenv + 1;
obj.env(1,obj.nenv) = envTarget;
obj.H1Env(:,:,obj.nenv) = envResults.H1Env;
if (isfield(envResults,'MP2'))
   obj.MP2Env(1,obj.nenv) = envResults.MP2;
end
obj.EhfEnv(1,obj.nenv)  = envResults.Ehf;
obj.EorbEnv(:,obj.nenv) = envResults.Eorb;
obj.HnucEnv(:,obj.nenv) = envResults.Hnuc;
obj.orbEnv(:,:,obj.nenv)  = envResults.orb;
obj.dipoleEnv(:,obj.nenv) = envResults.dipole;

end
