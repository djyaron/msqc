function loadZipData(frag,zipfile)

tempDir = tempname([frag.gaussianPath,filesep,'Scratch']);
mkdir(tempDir);
try
    unzip(zipfile,tempDir);
catch
    return
end
origdir = cd(tempDir); % save location so can move back
% convert checkpoint file to a formatted checkpoint file
resp1 = system([frag.gaussianPath,'\formchk.exe full.chk full.fch']);
cd(origdir);
% read in data from formatted checkpoint file
try
   fid1 = fopen([tempDir,filesep,'full.fch'],'r');
   if (fid1 == -1)
      error('could not find fch file');
   end
   [CorrE frag.MP2, frag.Ehf, frag.Eorb, frag.orb, frag.nelec,  frag.Z, frag.rcart, ...
    frag.dipole, frag.mulliken, frag.basisAtom, frag.basisType, ...
    frag.basisSubType, frag.basisNprims, frag.basisPrims ] = ...
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
   [frag.S, frag.H1, frag.KE, frag.H2, frag.Hnuc] = Fragment.readpolyatom(fid1);
   fclose(fid1);
catch
   fclose(fid1);
   error('failed during polyatom read');
end
frag.nbasis = size(frag.H1,1);

[n1,n2] = size(frag.H1);
natom = frag.natom;
frag.H1en = zeros(n1,n2,natom);

if frag.config.calcEn == 1 && frag.config.opt == 0
    for iatom = 1:natom
       jobname = ['atom',num2str(iatom)];
       try
          fid1 = fopen([tempDir,filesep,jobname,'.f32'],'r');%,'b');
          if (fid1 == -1)
             error(['could not find ',tempDir,filesep,jobname,'.f32']);
          end
          [~, H1atom, KE, ~, ~] = Fragment.readpolyatom(fid1);
          fclose(fid1);
       catch
          %fclose(fid1);
          throw(['failed during polyatom read for atom ',num2str(iatom)]);
       end
       frag.H1en(:,:,iatom) = H1atom - KE;
    end
end

cd(origdir);
% cleanup files
status = rmdir(tempDir,'s');
while status ~= 1
    disp('  rmdir failed. Retrying...');
    pause(0.1);
    status = rmdir(tempDir,'s');
end

