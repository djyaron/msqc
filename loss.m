function ell = loss(frag, fraghat, i0, i1, A, B)
ell = abs(E1(frag, i0, i1, A, B) - E1(fraghat, i0, i1, A, B));
end

function e = E1(frag, ind, A, B)
% ind = environment number
% A = list of orbitals to include for "atom 1"
% B = list of orbitals to include for "atom 2"
% Example of getting above lists for atoms iatom and jatom
%  A = find(frag.basisAtom == iatom);
%  B = find(frag.basisAtom == jatom);
rho = frag.density(ind);
H1 = frag.H1 + frag.H1Env(:,:,ind);
e = sum(sum(rho(A,B).*H1(A,B)));
end

function e = E2(frag, ind, A, B, C, D)
% ind = environment number
% A,B,C,D = list of orbitals to include for each index
% if C,D are not passed, then equiv to: E2(frag,ind, A,A,B,B)
if (nargin == 4)
   C = B; 
   D = B;
   B = A;
end
p2hf = frag.density2p(ind);
O2   = frag.H2;
e = sum(sum(sum(sum( p2hf(A,B,C,D) .*O2(A,B,C,D) ))));
end