clear classes;
load('datasets/ch4rDat.mat');
f1 = HL{1,1};
%%
Eke1 = f1.EKE;
Ehf = Eke1;
for iatom = 1:f1.natom
   Een1{iatom} = f1.Een(iatom);
   Ehf = Ehf + Een1{iatom};
end
E21 = f1.E2;
for ienv=0:f1.nenv;
  E21test(ienv+1) = f1.partitionE2(ienv,f1.H2,{1:f1.nbasis});
end
Ehf = Ehf + E21test;
Eenv1 = f1.Eenv;
Ehf = Ehf + Eenv1;
Ehf = Ehf + f1.Hnuc;
Egaussian = [f1.Ehf, f1.EhfEnv];
disp(['diff in E2 from density trace ',num2str(max(abs(E21-E21test)))]);
disp(['diff in Ehf from gaussian ',num2str(max(abs(Egaussian-Ehf)))]);

%%
m1 = Model3(f1,f1,f1);
m1.solveHF;
Emodel3 = [m1.Ehf, m1.EhfEnv];
disp(['diff mod3 from stored ',num2str(max(abs(Egaussian-Emodel3)))]);
save('junk.mat','m1');
%%

% diagonalize a random symmetric matrix to get a orthogonal U
A = rand(f1.nbasis,f1.nbasis);
A = A' + A;
[U,d] = eig(A);

U = eye(9);
SC = f1.S(1:5,1:5);
U(1:5,1:5) = (SC)^(-1/2);

f1.transform(U);
%%
Eke2 = f1.EKE;
Ehf2 = Eke2;
for iatom = 1:f1.natom
   Een2{iatom} = f1.Een(iatom);
   Ehf2 = Ehf2 + Een2{iatom};
end
E22 = f1.E2;
for ienv=0:f1.nenv;
  E22test(ienv+1) = f1.partitionE2(ienv,f1.H2,{1:f1.nbasis});
end
Ehf2 = Ehf2 + E22test;
Eenv2 = f1.Eenv;
Ehf2 = Ehf2 + Eenv2;
Ehf2 = Ehf2 + f1.Hnuc;
Egaussian = [f1.Ehf, f1.EhfEnv];
disp('from transformed fragment');
disp(['diff in E2 from density trace ',num2str(max(abs(E22-E22test)))]);
disp(['diff in Ehf from gaussian ',num2str(max(abs(Egaussian-Ehf2)))]);

disp(['trans effect on Eke ',num2str(max(abs(Eke2-Eke1)))]);
disp(['trans effect on E2 ',num2str(max(abs(E22-E21)))]);
disp(['trans effect on E2test ',num2str(max(abs(E22test-E21test)))]);
disp(['trans effect on Eenv ',num2str(max(abs(Eenv2-Eenv1)))]);

