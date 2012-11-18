function initializeZipData(obj,zipFileName)
% runs gaussian in scratch directories, and puts all relavant files
% into zipFileName
disp('initializing the data');
newline = char(10);

basisSet = obj.config.basisSet;
method   = obj.config.method;
charge   = obj.config.charge;
spin     = obj.config.spin;
par      = obj.config.par;
ctext    = obj.gaussianFile;
dataPath = obj.dataPath;
gaussianPath = obj.gaussianPath;
gaussianExe  = obj.gaussianExe;

% Do the calculation
jobname = 'full';
gjf_file = [jobname,'.gjf'];
tempDir = tempname([gaussianPath,'/','Scratch']);
mkdir(tempDir);
origdir = cd(tempDir); % save location so can move back
fid1 = fopen(gjf_file,'w');
fwrite(fid1, [ctext,newline,newline], 'char');
fclose(fid1);

setenv('GAUSS_EXEDIR', obj.gaussianPath);
resp1 = 1; 
while ( resp1 ~= 0)
    try
        resp1 = system([gaussianPath,'\',gaussianExe,' ',gjf_file]);
        %disp( resp1 );
        %disp( resp2 );
        if ( resp1 == 2057 )
            disp( '  removing temporary files' );
            delete( 'fort.6', 'gxx.d2e', 'gxx.inp', 'gxx.int', 'gxx.scr', ...
                    'temp.chk', 'temp.fch', 'temp.rwf' )
        end
    catch
        disp( 'Failed, retrying...' );
        resp1 = 1; 
    end
end
movefile('temp.chk','full.chk');
movefile('fort.32','full.f32');

toZip = {'full.gjf','full.chk','full.f32'};

% Do calculations with only one nucleus present at a time

% header for the Gaussian job file (input file)
header = ['%rwf=temp.rwf',newline,...
   '%nosave',newline,...
   '%chk=temp.chk',newline,...
   '# hf/',basisSet, newline...
   'nosymm int=noraff iop(99/6=1) ',...
   'scf=conventional',' symm=noint', newline, newline, ...
   'title', newline,newline];

% ______________________________________________

[n1,n2] = size(obj.H1);
natom = obj.natom;
obj.H1en = zeros(n1,n2,natom);

for iatom = 1:natom
   disp(['doing calc for atom ',num2str(iatom)]);
   ctext = header;
   % SCF can have trouble converging for these single-atom calculations
   % so we add a qc keyword, to help ensure convergence
   ctext = strrep(ctext,'scf=conventional','scf=(conventional,qc)');
   % charge and spin come next
   % To keep even number of electrons (and so spin 1), add an electron
   % by increasing the charge
   tt1 = obj.templateText;
   tt2 = strfind(tt1,['ATOM',num2str(iatom)]);
   atomSym = tt1(tt2-1);
   if (strcmp(atomSym,'H') || strcmp(atomSym,'F')) % should be made general, for any odd Z 
      tempCharge = charge -1;
   else
      tempCharge = charge;
   end
   ctext = [ctext, num2str(tempCharge), ' ', num2str(spin), newline];
   % For molecule specification, we first replace all ATOM# with spaces
   t1 = obj.templateText;
   % Iterate in reverse order, or replacements will not work properly 
   % with more than 10 atoms.
   for jatom = natom:-1:1
      if (jatom == iatom)
         t1 = strrep(t1, ['ATOM',num2str(jatom)],' ');
      else
         t1 = strrep(t1, ['ATOM',num2str(jatom)],'-Bq');
      end
   end
   % And replace all PAR# with the parameter values
   for ipar = obj.npar:-1:1
      t1 = strrep(t1, ['PAR',num2str(ipar)], num2str(par(ipar),'%23.12f'));
   end
   ctext = [ctext, t1];
   
   % Do the calculation and read in data
   jobname = ['atom',num2str(iatom)];
   gjf_file = [jobname,'.gjf'];
   fid1 = fopen(gjf_file,'w');
   fwrite(fid1, [ctext,newline,newline], 'char');
   fclose(fid1);
   resp1 = 1;
   while ( resp1 ~= 0 )
       try
           resp1 = system([gaussianPath,'\',gaussianExe,' ',gjf_file]);
           %disp( resp1 )
           if ( resp1 == 2057 )
               disp( '  removing temporary files' );
               delete( 'fort.6', 'gxx.d2e', 'gxx.inp', 'gxx.int', 'gxx.scr', ...
                   'temp.chk', 'temp.fch', 'temp.rwf' )
           end
       catch
           disp( 'Failed, retrying...' );
           resp1 = 1;
       end
   end
   movefile('fort.32',[jobname,'.f32']);
   % read in data from the polyatom output file

   toZip = {toZip{:},[jobname,'.f32']};
end
zip(zipFileName,toZip);

cd(origdir);
% cleanup files
%[status message messageid] = rmdir(tempDir,'s');
status = rmdir(tempDir,'s');
while status ~= 1
    disp('  rmdir failed. Retrying...');
    pause(0.1);
    status = rmdir(tempDir,'s');
end

