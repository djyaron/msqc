function runModelsParallel(dirName,nfiles)

mods = Model3.empty(nfiles,0);
envs = cell(nfiles,1);
for i = 1:nfiles
   infilename = [dirName,'todo',num2str(i),'.mat'];
   load(infilename); % 
   mods(i,1) = modFile1;
   envs{i,1} = envsFile1;
end

parfor i = 1:nfiles
   mods(i,1).solveHF(envs{i,1});
end

for i = 1:nfiles
   infilename = [dirName,'todo',num2str(i),'.mat'];
   outfilename = [dirName,'done',num2str(i),'.mat'];
   modFile2 = mods(i,1);
   save(outfilename,'modFile2');
   delete(infilename);
end

end