function initializeData(obj)

disp('initializing the data');
newline = char(10);

basisSet = obj.config.basisSet;
charge   = obj.config.charge;
spin     = obj.config.spin;
par      = obj.config.par;
dataPath = obj.dataPath;
gaussianPath = obj.gaussianPath;
gaussianExe  = obj.gaussianExe;
% ___________________________________________________________
% header for the Gaussian job file (input file)
header = ['%rwf=temp.rwf',newline,...
   '%nosave',newline,...
   '%chk=temp.chk',newline,...
   '# hf/',basisSet, newline...
   'nosymm int=noraff iop(99/6=1) ',...
   'scf=conventional',' symm=noint', newline, newline, ...
   'title', newline,newline];

% ____________________________________________________________
% Do calculation on entire fragment

% ctext will hold the Gaussian job file (input file)
% begin with the above header text
ctext = header;
% charge and spin come next
ctext = [ctext, num2str(charge), ' ', num2str(spin), newline];
% For molecule specification, we first replace all ATOM# with spaces
t1 = obj.templateText;
for iatom = 1:obj.natom
   t1 = strrep(t1, ['ATOM',num2str(iatom)], ' ');
end
% And replace all PAR# with the parameter values
for ipar = 1:obj.npar
   t1 = strrep(t1, ['PAR',num2str(ipar)], num2str(par(ipar),'%23.12f'));
end
ctext = [ctext, t1];

% add charge keyword, so calcs can be done in an environment
obj.gaussianFile = strrep(ctext,'symm=noint','symm=noint charge');

% Do the calculation and read in data
jobname = 'full';
gjf_file = [jobname,'.gjf'];
origdir = cd(obj.dataPath);
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
   fid1 = fopen([dataPath,'\temp.fch'],'r');
   if (fid1 == -1)
      error('could not find fch file');
   end
   [obj.Ehf, obj.Eorb, obj.orb, obj.nelec,  obj.Z, obj.rcart, ...
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
   fid1 = fopen([dataPath,'\fort.32'],'r');
   if (fid1 == -1)
      error(['could not find ',dataPath,'\fort.32']);
   end
   [obj.S, obj.H1, obj.KE, obj.H2, obj.Hnuc] = Fragment.readpolyatom(fid1);
   fclose(fid1);
catch
   fclose(fid1);
   error('failed during polyatom read');
end
obj.nbasis = size(obj.H1,1);

% save files for debugging
% system(['copy ', dataPath,'\temp.fch ', dataPath,'\debug.fch']);
% system(['copy ', dataPath,'\full.out ', dataPath,'\debug.out']);
% cleanup files
delete([dataPath,'\fort.32'], [dataPath,'\full.gjf'], ...
   [dataPath,'\full.out'], [dataPath,'\temp.chk'], ...
   [dataPath,'\temp.fch']);
% ____________________________________________________________
% Do calculation with only one nucleus present at a time

[n1,n2] = size(obj.H1);
natom = obj.natom;
obj.H1en = zeros(n1,n2,natom);
for iatom = 1:natom
   disp(['doing calc for atom ',num2str(iatom)]);
   ctext = header;
   % charge and spin come next
   % To keep even number of electrons (and so spin 1), add an electron
   % by increasing the charge
   % TODO: this should only be done for certain atoms, should fix
   ctext = [ctext, num2str(charge-1), ' ', num2str(spin), newline];
   % For molecule specification, we first replace all ATOM# with spaces
   t1 = obj.templateText;
   for jatom = 1:natom
      if (jatom == iatom)
         t1 = strrep(t1, ['ATOM',num2str(jatom)],' ');
      else
         t1 = strrep(t1, ['ATOM',num2str(jatom)],'-Bq');
      end
   end
   % And replace all PAR# with the parameter values
   for ipar = 1:obj.npar
      t1 = strrep(t1, ['PAR',num2str(ipar)], num2str(par(ipar),'%23.12f'));
   end
   ctext = [ctext, t1];
   
   % Do the calculation and read in data
   jobname = ['atom',num2str(iatom)];
   gjf_file = [jobname,'.gjf'];
   origdir = cd(obj.dataPath);
   fid1 = fopen(gjf_file,'w');
   fwrite(fid1, ctext, 'char');
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
   cd(origdir);
   % read in data from the polyatom output file
   try
      fid1 = fopen([dataPath,'\fort.32'],'r');
      if (fid1 == -1)
         error(['could not find ',dataPath,'\fort.32']);
      end
      [~, H1atom, KE, ~, ~] = Fragment.readpolyatom(fid1);
      fclose(fid1);
   catch
      fclose(fid1);
      error(['failed during polyatom read for atom ',num2str(iatom)]);
   end
   obj.H1en(:,:,iatom) = H1atom - KE;
   
   % cleanup files
   delete([dataPath,'\fort.32'], [dataPath,'\',jobname,'.gjf'], ...
      [dataPath,'\',jobname,'.out'], [dataPath,'\temp.chk']);
   
end


