function runModelsParallel(dirName,nfiles)

parfor i = 1:nfiles
   infilename = [dirName,'todo',num2str(i),'.mat'];
   outfilename = [dirName,'done',num2str(i),'.mat'];
   load(infilename); % loads modFile and envsFile
   modFile.solveHF(envsFile);
   outFile = modFile;
   save(outfilename,'outFile');
   delete(infilename);
end

end