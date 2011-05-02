%parfor n_0=1:10
%  for p=1:5
for n_0=1:10
  for p=1:5
    n = n_0 * 10;
    name = ['n=',int2str(n),'p=',int2str(p)];
    experiment_gen(name, 'data2/fhydeFixed.tpl', '6-31G', 'data2/fhydeGen.tpl', n, p);
    experiment_comp(name);
  end
end
 