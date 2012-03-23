function initializeData(obj)

disp('initializing the data');
newline = char(10);

basisSet = obj.config.basisSet;
method   = obj.config.method;
charge   = obj.config.charge;
spin     = obj.config.spin;
par      = obj.config.par;
dataPath = obj.dataPath;
gaussianPath = obj.gaussianPath;
gaussianExe  = obj.gaussianExe;
% ___________________________________________________________
% header for the Gaussian job file (input file)
%  Note: for single atom calcs below, 'scf=conventional' is replaced
%        so if this keyword in header is changed, it needs to be changed
%        there as well
header = ['%rwf=temp.rwf',newline,...
   '%nosave',newline,...
   '%chk=temp.chk',newline,...
   '# ',method,'/',basisSet, newline...
   'nosymm int=noraff iop(99/6=1) ',...
   'scf=conventional',' symm=noint', newline, newline, ...
   'title', newline,newline];

% ____________________________________________________________
% Create Scratch directory within g09 to be used to store temp directories
% to store output files.

%if(~exist([gaussianPath, 'Scratch'], 'dir'));
%       mkdir(gaussianPath, 'Scratch');
%end
   
% Do calculation on entire fragment

% ctext will hold the Gaussian job file (input file)
% begin with the above header text
ctext = header;
% charge and spin come next
ctext = [ctext, num2str(charge), ' ', num2str(spin), newline];
% For molecule specification, we first replace all ATOM# with spaces
t1 = obj.templateText;
% Iterate in reverse order, or replacements will not work properly with 
% more than 10 atoms.
for iatom = obj.natom:-1:1
   t1 = strrep(t1, ['ATOM',num2str(iatom)], ' ');
end
% And replace all PAR# with the parameter values
for ipar = obj.npar:-1:1
   t1 = strrep(t1, ['PAR',num2str(ipar)], num2str(par(ipar),'%23.12f'));
end
ctext = [ctext, t1];

obj.gaussianFile = ctext;

% Do the calculation and read in data
jobname = 'full';
gjf_file = [jobname,'.gjf'];
tempDir = tempname([gaussianPath,'/','Scratch']);
mkdir(tempDir);
origdir = cd(tempDir); % should move into unique scratch directory
fid1 = fopen(gjf_file,'w');
fwrite(fid1, [ctext,newline,newline], 'char');
fclose(fid1);

setenv('GAUSS_EXEDIR', obj.gaussianPath);
system([gaussianPath,'\',gaussianExe,' ',gjf_file]);
% convert checkpoint file to a formatted checkpoint file
system([gaussianPath,'\formchk.exe temp.chk temp.fch']);
cd(origdir);
% read in data from formatted checkpoint file
try
   fid1 = fopen([tempDir,'/','temp.fch'],'r');
   if (fid1 == -1)
      error('could not find fch file');
   end
   [CorrE obj.MP2, obj.Ehf, obj.Eorb, obj.orb, obj.nelec,  obj.Z, obj.rcart, ...
    obj.dipole, obj.mulliken, obj.basisAtom, obj.basisType, ...
    obj.basisSubType, obj.basisNprims, obj.basisPrims ] = ...
    Fragment.readfchk(fid1);
   
   fclose(fid1);
catch
   fclose(fid1);
   error('failed during fchk read');
end
% read in data from the polyatom output file
try
%   fid1 = fopen([tempDir,'/','fort.32'],'r','b'); % only for mac version
   fid1 = fopen([tempDir,'/','fort.32'],'r');
   if (fid1 == -1)
      error(['could not find ',tempDir,'\','fort.32']);
   end
   [obj.S, obj.H1, obj.KE, obj.H2, obj.Hnuc] = Fragment.readpolyatom(fid1);
   fclose(fid1);
catch
   fclose(fid1);
   error('failed during polyatom read');
end
obj.nbasis = size(obj.H1,1);

% save files for debugging
system(['copy ', tempDir,filesep,'full.gjf ', tempDir,filesep,'debug.gjf']);
%system(['copy ', tempDir,filesep,'temp.fch ', tempDir,filesep,'debug.fch']);
%system(['copy ', tempDir,filesep,'full.out ', tempDir,filesep,'debug.out']);

% cleanup files
rmdir(tempDir,'s');

%delete([tempDir,filesep,'fort.32'], [tempDir,filesep,'full.gjf'], ...
%   [tempDir,filesep,'full.out'], [tempDir,filesep,'temp.chk'], ...
%   [tempDir,filesep,'temp.fch']);
% ____________________________________________________________
% Do calculation with only one nucleus present at a time

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
   if (rem(obj.Z(iatom),2) == 1) 
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
   tempDir = tempname([gaussianPath,'/','Scratch']);
   mkdir(tempDir);
   origdir = cd(tempDir);
   fid1 = fopen(gjf_file,'w');
   fwrite(fid1, [ctext,newline,newline], 'char');
   fclose(fid1);
   % run gaussian
%    for i = 1:3
%       if (exist(gfiles{i},'file'))
%          disp(['deleting ', gfiles{i}]);
%          delete(gfiles{i});
%       else
%          disp(['did not find ', gfiles{i}]);
%       end
%    end
%    input junk;
   
   system([gaussianPath,'\',gaussianExe,' ',gjf_file]);
   % check if the call to Gaussian failed.
   if ans ~= 0
       error( 'Gaussian failed. Check template, Gaussian input file, or Gaussian output file.' );
   end
   cd(origdir);
   % read in data from the polyatom output file
   try
      fid1 = fopen([tempDir,filesep,'fort.32'],'r');%,'b');
      if (fid1 == -1)
         error(['could not find ',tempDir,filesep,'fort.32']);
      end
      [junk, H1atom, KE, junk2, junk3] = Fragment.readpolyatom(fid1);
      fclose(fid1);
   catch
      %fclose(fid1);
      throw(['failed during polyatom read for atom ',num2str(iatom)]);
   end
   obj.H1en(:,:,iatom) = H1atom - KE;
   
% cleanup files
rmdir(tempDir,'s');
   
   %delete([tempDir,filesep,'fort.32'], [tempDir,filesep,jobname,'.gjf'], ...
   %   [tempDir,filesep,jobname,'.out'], [tempDir,filesep,'temp.chk']);
   
end


