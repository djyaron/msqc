function loadZipData(obj,zipfile)

tempDir = tempname([obj.gaussianPath,filesep,'Scratch']);
mkdir(tempDir);
unzip(zipfile,tempDir);
origdir = cd(tempDir); % save location so can move back
% convert checkpoint file to a formatted checkpoint file
resp1 = system([obj.gaussianPath,'\formchk.exe full.chk full.fch']);
cd(origdir);
% read in data from formatted checkpoint file
try
   fid1 = fopen([tempDir,filesep,'full.fch'],'r');
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
   fid1 = fopen([tempDir,filesep,'full.f32'],'r');
   if (fid1 == -1)
      error(['could not find ',tempDir,'\','full.f32']);
   end
   [obj.S, obj.H1, obj.KE, obj.H2, obj.Hnuc] = Fragment.readpolyatom(fid1);
   fclose(fid1);
catch
   fclose(fid1);
   error('failed during polyatom read');
end
obj.nbasis = size(obj.H1,1);

[n1,n2] = size(obj.H1);
natom = obj.natom;
obj.H1en = zeros(n1,n2,natom);

for iatom = 1:natom
   jobname = ['atom',num2str(iatom)];
   try
      fid1 = fopen([tempDir,filesep,jobname,'.f32'],'r');%,'b');
      if (fid1 == -1)
         error(['could not find ',tempDir,filesep,jobname,'.f32']);
      end
      [junk, H1atom, KE, junk2, junk3] = Fragment.readpolyatom(fid1);
      fclose(fid1);
   catch
      %fclose(fid1);
      throw(['failed during polyatom read for atom ',num2str(iatom)]);
   end
   obj.H1en(:,:,iatom) = H1atom - KE;
end

cd(origdir);
% cleanup files
status = rmdir(tempDir,'s');
while status ~= 1
    disp('  rmdir failed. Retrying...');
    pause(0.1);
    status = rmdir(tempDir,'s');
end

