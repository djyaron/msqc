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
         tmp = envTarget.fieldType( ifield, : );
         if ( sum( tmp ) < 4 )
             dir = [ repmat( 'X', 1, envTarget.fieldType( ifield, 1 ) ), ...
                     repmat( 'Y', 1, envTarget.fieldType( ifield, 2 ) ), ...
                     repmat( 'Z', 1, envTarget.fieldType( ifield, 3 ) ) ]; 
         else
             dir = '';
             while max( tmp ) > 0
                switch find( tmp == max( tmp ), 1 )
                   case 1
                      dir = [ dir, repmat( 'X', 1, tmp( 1 ) ) ];
                      tmp( 1 ) = 0;
                   case 2
                      dir = [ dir, repmat( 'Y', 1, tmp( 2 ) ) ];
                      tmp( 2 ) = 0;
                   case 3
                      dir = [ dir, repmat( 'Z', 1, tmp( 3 ) ) ];
                      tmp( 3 ) = 0;
                   otherwise
                      error( 'Invalid input' );
                end
             end
         end
         mag = '';
         if envTarget.fieldMag( ifield ) >= 0
            mag = '+';
         end
         mag = [ mag, num2str( envTarget.fieldMag( ifield ) ) ];
         insert = [ 'Field=', dir, mag ];
         % Add to the input file
         ctext = strrep( ctext, 'symm=noint', [ 'symm=noint ', insert ] );
      end
   end
   
   gjf_file = [jobname,'.gjf'];
   origdir = cd(tempDir);
   fid1 = fopen(gjf_file,'w');
   fwrite(fid1, ctext, 'char');
   fclose(fid1);
   
   dataPath = tempDir;
   system([obj.gaussianPath,filesep,obj.gaussianExe,' ',gjf_file]);
   % convert checkpoint file to a formatted checkpoint file
   system([obj.gaussianPath,'\formchk.exe temp.chk temp.fch']);
   cd(origdir);
   % read in data from formatted checkpoint file
   try
      fid1 = fopen([dataPath,filesep,'temp.fch'],'r');
      if (fid1 == -1)
         error('could not find fch file');
      end

   [CorrEe, MP2e, Ehfe, Eorbe, orbe, junk,  junk2, junk3, ...
    dipolee, junk4, junk5, junk6, ...
    junk7, junk8, junk9] = ...
    Fragment.readfchk(fid1);
      fclose(fid1);
   catch
      disp('caught some stupid error');
      fclose(fid1);
      error('failed during env fchk read');
   end
   % read in data from the polyatom output file
   try
      fid1 = fopen([dataPath,filesep,'fort.32'],'r');%,'b');
      if (fid1 == -1)
         error(['could not find ',dataPath,filesep,'fort.32']);
      end
      [junk1, H1e, junk2, junk3, Hnuce] = Fragment.readpolyatom(fid1);
      fclose(fid1);
   catch
      fclose(fid1);
      error('failed during polyatom read');
   end
   
   envResults.H1Env = H1e - obj.H1;
   envResults.Ehf   = Ehfe;
   envResults.MP2   = MP2e;
   envResults.Eorb  = Eorbe;
   envResults.orb   = orbe;
   envResults.Hnuc  = Hnuce;
   envResults.dipole = dipolee;
   % cleanup files
   rmdir(tempDir,'s');
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
if (isfield(envResults,'MP2'))
   obj.MP2Env(1,obj.nenv) = envResults.MP2;
end
obj.EhfEnv(1,obj.nenv)  = envResults.Ehf;
obj.EorbEnv(:,obj.nenv) = envResults.Eorb;
obj.HnucEnv(:,obj.nenv) = envResults.Hnuc;
obj.orbEnv(:,:,obj.nenv)  = envResults.orb;
obj.dipoleEnv(:,obj.nenv) = envResults.dipole;
