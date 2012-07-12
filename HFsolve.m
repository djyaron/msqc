function [P,orb,Eorb,Ehf,failed ] = ...
   HFsolve(H1, H2, S, Enuc, Nelec, guessDensity, ...
   eps,maxIter,minIter)
% Solve Hartree Fock equations
% Input:
%   H1: (nb2asis,nbasis) H1 matrix
%   H2: (nbasis,nbasis,nbasis,nbasis) H2 elements
%    S: (nbasis,nbasis) overlap matrix
%   Ennuc:  nuclear energy
%   Nelec:  number of electrons
%   guessDensity: Starting density matrix
%
% Output
%     P:   (nbasis,nbasis) density matrix
%   orb:   (nbasis,nbasis)
%            orbital coefficient matrix (atom basis #, mol orb #)
%   Eorb:  (nbasis,1)  molecular orbital energies
%   Ehf:    total Hartree-Fock energy
%   failed:  1 if failed to converge

if (nargin < 7)
   eps = 1.0e-10;
end
if (nargin < 8)
   maxIter = 5000;
end
if (nargin < 9)
   minIter = 5;
end

Nbasis = size(H1,1);
%  H2j {nbasis,nbasis} cell array of coulomb integrals
%  H2k {nbasis,nbasis} cell array of exchange integrals
H2j = cell(Nbasis,Nbasis);
H2k = cell(Nbasis,Nbasis);
for i=1:Nbasis
   for j=1:Nbasis
      H2j{i,j} = squeeze(H2(i,j,:,:));
      H2k{i,j} = squeeze(H2(i,:,:,j));
   end
end


%step 3 -- Calculate transformation matrix (eq. 3.167)
X = inv(sqrtm(S));
Pn = guessDensity;

Plast = Pn; % previous density matrix
iter = 0;
filled = 1:(Nelec/2);
finished = false;

%Begin iteration through
while (~finished) %step 11 -- Test convergence
   
   if (iter < maxIter/2)
      P = 0.5 * Pn + 0.5 * Plast;
   else
      P = Pn;
   end
   
   %step 5 -- Build 2-electron components of Fock matrix
   G = zeros(Nbasis);
   for i=1:Nbasis
      for j=1:Nbasis
         t1 = sum(sum( P'.* H2j{i,j} ));
         t2 = sum(sum( P'.* H2k{i,j} ));
         G(i,j) = G(i,j) + t1 - t2/2;
      end
   end
   %step 6 -- Obtain F (fock matrix)
   F = H1 + G;
   
   %step 7 -- Calculate the transformed F matrix
   Ft = X'*F*X; %#ok<MINV>
   
   %step 8 -- Find e and the transformed expansion coefficient matrices
   [Ct1,e1] = eig(Ft);
   e2 = diag(e1);
   [e, i1] = sort(e2);
   Ct = Ct1(:,i1);
   
   %step 9 -- Transform Ct back to C
   C = X*Ct; %#ok<MINV>
   
   %step 10 -- Calculate the new density matrix
   Plast = Pn;
   Pn = zeros(Nbasis);
   %Cj = conj(C);
   Pn = 2* C(:,filled)*( C(:,filled)');
   iter = iter + 1;
   %disp(['den change ',num2str( max(max(abs(P-Pn))))]);
   %disp(['diag Pn ',num2str(diag(Pn)')]);
   %Pn-Plast
   if (iter > maxIter)
      finished = true;
   elseif (iter > minIter)
      if (max(max(abs(P - Pn))) < eps)
         finished = true;
      end
   end
end
%End of iteration of steps 5-11

%if (iter < maxIter)
P = Pn; %for convenience

%Step 12: Output

%Total energy
%3.184: E0 = 1/2 Sum(i,j) {P(j,i)[H1(i,j) + F(i,j)]}
Ee = sum(sum(P.*(H1+F)));
Ehf = Ee/2 + Enuc;

%Orbital energies
Eorb = e;

%Molecular orbital components
orb = C;
failed = (iter+1 > maxIter);
if (failed)
   disp('convergence failure');
end
end

