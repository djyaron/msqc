function d = Delta(frag, i0, i1, A, B)
  e0 = E(frag, i0, A, B);
  e1 = E(frag, i1, A, B);
  d = e0 - e1;     
end
  
function e = E(frag, ind, A, B)
  rho = frag.density(ind);
  gamma = frag.density2p(ind);
  H1 = frag.H1 + frag.H1Env(:,:,ind);
  H2 = frag.H2;
  e = sum(sum(rho(A,B).*H1(A,B))) + sum(sum(sum(sum(gamma(A,B,A,B).*H2(A,B,A,B)))));
end 