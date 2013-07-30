MSQC (Molecular Similarity in Quantum Chemistry)
================================================

TODO: Add big picture documentation.

Compilation of MEX C Code
-------------------------

There are a few MEX routines written in C that must be compiled before use. An attempt will be made to do the compilation automatically if the binaries do not exist. Alternatively, all of the code can be compiled by simply calling `mexCompile`. If you've never used `mex` before, you'll need to choose a compiler on the first run. If none are available, see http://www.mathworks.com/support/compilers/R2012b/win64.html. For 64 bit Windows and MATLAB R2012b, you probably want http://www.microsoft.com/en-us/download/details.aspx?id=8279.

Basic Flow of a Train and Test Run
----------------------------------

1. Make a set of models, by taking the frag's for the molecules you care about, and making Model3 classes out of them. At this point, those models are unmodifed (no parameters, so calling them will just do the bare LL method).  Make one set for the train, and then as many test sets as you like.
2. Create an MFactory, by establishing policies that say how mixers are created.
3. Create a Fitme object for the train MSet
4. Use some routine to optimize the parameters using that Fitme.
5. Now the MFactory has the mixers, with optimized parameters, and the policies about how to apply these mixers.
6. Create a Fitme object for each test data set, and when you call .err on this fitme, it will magically work.  You have to be careful that the ordering of parameters in this fitme object will not be the same as other fitme's created by the same MFactory.

Demo Code
---------

```matlab
%% Demonstration on 6/19/13, using MFactory, MixerC etc.
clear classes;
fileName = 'D:\dave\apoly\msqc\datasets\ch4rDat.mat';
ms = MSet;
ms.addData(fileName,1:2,1:2,1,101);
mtest = MSet;
mtest.addData(fileName,3:4,1:2,1,101);
mf = MFactory;
mf.setPolicies('hybridslater1');
mf.makeMixInfo([1 6]);
ftrain = mf.makeFitme(ms);
ftest = mf.makeFitme(mtest);

%%
maxIter = 100;
epsTest = 1.0e-5
updateContext = 0;
[err,pt, testErr, monitor] = contextFit4(ftrain,ftest,...
  maxIter,epsTest,updateContext);
```
