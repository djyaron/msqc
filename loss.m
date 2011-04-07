function ell = loss(frag, fraghat, experiments) 
  ell = 0;
  nenv = frag.nenv; % FIXME: memoization
  for i=1:size(experiments)
    [env0,env1,A,B] = experiments{i,:};  
    
    frag.addEnv(env0);
    frag.addEnv(env1);
    fraghat.addEnv(env0);
    fraghat.addEnv(env1);

    i0 = nenv + 1 + i; % FIXME: memoization
    i1 = nenv + 1 + i + 1;
    ell = ell + abs(Delta(frag, i0, i1, A, B) - Delta(fraghat, i0, i1, A, B));
  end
end
  
function d = Delta(frag, env0, env1, A, B)
  e0 = E(frag, env0, A, B);
  e1 = E(frag, env1, A, B);
  d = e0 - e1;     
end
  
function e = E(frag, env, A, B)
  rho = frag.density(env);
  gamma = frag.density2p(env);
  H1 = frag.H1 + frag.H1Env(:,:,env);
  H2 = frag.H2;
  e = sum(sum(rho(A,B).*H1(A,B))) + sum(sum(sum(sum(gamma(A,B,A,B).*H2(A,B,A,B)))));
end 