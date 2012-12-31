load('T:\matdl\yaron\dec12a\hybridslater\ch4rDat\start.mat');
f1.printMixers;
f1.printEDetails;
[e0 plotnum etype modelnum envnum] = ftest.err(f1.getPars);
e0 =e0*627.509;
%% Error from just shifting
[e0 plotnum etype modelnum envnum] = ftest.err(zeros(size(f1.getPars)));
e0 =e0*627.509;
%%
eke = e0(etype==1); eH  = e0(etype==11);
eC  = e0(etype==16); e2  = e0(etype==2);
etot = e0(etype==3);
disp(['ke ',num2str(std(eke)),' H ',num2str(std(eH)),' C ',...
   num2str(std(eC)),' E2 ',num2str(std(e2)),' Etot ', num2str(std(etot))]);