function ell = loss(frag_cfg, fraghat_cfg, tplpath, experiments)
  [frag, frag_dir] = tmpfrag(frag_cfg, tplpath);
  [fraghat, fraghat_dir] = tmpfrag(fraghat_cfg, tplpath);
  %onCleanup(@() rmdir(frag_dir, 's') & rmdir(fraghat_dir, 's'));
  
  ell = 0;
  for i=1:size(experiments)
    [env0,env1,A,B] = experiments{i,:};  
    
    frag.addEnv(env0);
    frag.addEnv(env1);
    fraghat.addEnv(env0);
    fraghat.addEnv(env1);
    
    i0 = i;
    i1 = i + 1;
    ell = ell + abs(Delta(frag, i0, i1, A, B) - Delta(fraghat, i0, i1, A, B));
  end
end
  
function [frag, dir] = tmpfrag(cfg, tplpath)
  dir = tempname;
  if (mkdir(dir) ~= 1), error('Could not create temporary directory.'), end;
  copyfile(tplpath, [dir filesep cfg.template '.tpl']);
  frag = Fragment(dir, cfg);
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