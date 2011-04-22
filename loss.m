function ell = loss(frag, fraghat, i0, i1, A, B) 
    ell = abs(Delta(frag, i0, i1, A, B) - Delta(fraghat, i0, i1, A, B));
end