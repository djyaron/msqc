function transform(obj,U)
% Unitary transformation, via matrix U(oldBasis,newBasis)

% Transform each property that changes with basis set
% H1(nbasis,nbasis) full H1 operator of fragment
obj.H1 = U' * (obj.H1* U);
% H1en(nbasis,nbasis,natom) electron-nuclear interaction
for iatom = 1:obj.natom
   obj.H1en(:,:,iatom) = U' * (squeeze(obj.H1en(:,:,iatom))*U);
end
% KE(nbasis,nbasis) kinetic energy
obj.KE = U' * (obj.KE * U);
% S(nbasis,nbasis) overlap
obj.S = U' * (obj.S * U);
% H1Env(nbasis,nbasis,nenv) H1 due to environment
for ienv = 1:obj.nenv
   obj.H1Env(:,:,ienv) = U' * (squeeze(obj.H1Env(:,:,ienv)) * U);
end
% H2 (nbasis,nbasis,nbasis,nbasis) 2-elec interactions
N = obj.nbasis;
n2 = zeros(N,N,N,N);
for i=1:N
   for j=1:N
      for k=1:N
         for l=1:N
            for a=1:N
               for b=1:N
                  for c=1:N
                     for d=1:N
                        n2(i,j,k,l) = n2(i,j,k,l) + ...
                           obj.H2(a,b,c,d) * U(a,i) * U(b,j) ...
                           * U(c,k) * U(d,l);
                     end
                  end
               end
            end
         end
      end
   end
end
obj.H2 = n2;

% orb(nbasis,nbasis) molecular orbital coefficients
% orb is (atomic,molecular) so
obj.orb = U' * obj.orb;
% orbEnv(nbasis,nbasis,nenv) molecular orbitals in env
for ienv = 1:obj.nenv
   obj.orbEnv(:,:,ienv) = U' * squeeze(obj.orbEnv(:,:,ienv));
end


%
%       Ehf     % Hartree Fock energy
%       MP2     % MP2 Energy
%       %CorrE   % Correlation Energy (MP2-Ehf)
%       Eorb    % (nbasis,1)      molecular orbital energies
%       orb     % (nbasis,nbasis) molecular orbital coefficients
%       dipole  % (3,1)   dipole moment of molecule
%       mulliken % (1,natom)  mulliken charge on the atoms
% 
%       HnucEnv % (1,nenv)             Hnuc in environment
% 
%       EhfEnv   % (1,nenv)        Hartree-Fock energy in env
%       MP2Env   % (1,nenv)        MP2 energy in env
%       %CorrEenv % (1,nenv)        Correlation energy in env (MP2env-HFenv)
%       EorbEnv; % (nbasis,nenv)   molecular orbital energies in env
%       orbEnv;  % (nbasis,nbasis,nenv) molecular orbitals in env
%       dipoleEnv % (3,nenv) dipole moment in the environment
%       
%       basisAtom  % (nbasis,1) atom # on which the function is centered 
%       basisType  % (nbasis,1) l quantum number: 0=s 1=p 2=d 3=d etc
%       basisSubType % (nbasis,1) m quantum number: s=1 p=3 d=6 (cartesian)
%       basisNprims  % number of primitives in this function
%       basisPrims   % {nbasis,1} cell array of matrices of size (2,nprims)
%                    %    with (1,:) being contraction coefficients and
%                    %         (2,:) being primimitive exponents

