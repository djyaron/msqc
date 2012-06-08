function [orb,Eorb,Ehf,status] = hartreeFock(obj,eps,maxIter,minIter)
% Solve Hartree Fock equations
% Input:
%   obj:  Holds hamiltonian information (H1,H2,S,nelec,Hnuc,H1env,HnucEnv)
%          [see Fragment class for definitions of these properties]
%   ienv:   environment number (uses H1env(:,:,env) and HnucEnv(env))
%          [defaults to 0, with 0 meaning isolated fragment]
%   eps:   convergence criteria
%          [defaults to 1e-8]
% Output
%   orb:   (nbasis,nbasis)
%            orbital coefficient matrix (atom basis #, mol orb #)
%   Eorb:  (nbasis,1)  molecular orbital energies
%   Ehf:    total Hartree-Fock energy
%   status: 0 if all is ok, 1 if HF iterations didn't converge

if (nargin < 2)
   eps = 1.0e-6;
end
if (nargin < 3)
   maxIter = 1000;
end
if (nargin < 4)
   minIter = 5;
end

status = 0; % initialize to ok
H1 = obj.H1;
S  = obj.S;
Enuc = obj.Hnuc(ienv);
Nelec = obj.frag.nelec;

Nbasis = size(H1,1); %Getting size of basis set

% Steps refer to the list below from the Szabo and Ostlund book
%step 3 -- Calculate transformation matrix (eq. 3.167)
X = inv(sqrtm(S));

%step 4 -- Guess at density matrix
Pn = obj.frag.density(ienv);

Plast = Pn; % previous density matrix
iter = 0;
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
         t1 = sum(sum( P'.* obj.H2j{i,j} ));
         t2 = sum(sum( P'.* obj.H2k{i,j} ));
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
   filled = 1:(Nelec/2);
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
obj.densitySave{ienv+1} = P;

%Step 12: Output

%Total energy
Ee = sum(sum(P.*(H1+F)));
Ehf = Ee/2 + Enuc;

%Orbital energies
Eorb = e;

%Molecular orbital components
orb = C;

if (iter+1 > maxIter)
   status = 1;
end

%{
Adapted from "Modern quantum chemistry", by Attila Szabó, Neil S. Ostlund
Numbered equations also adapted from here.
1. Specify a molecule
2. Calculate S(i,j), H^core (H1), and (i j|k l)(H2)
    -These first two steps are done by Gaussian
3. Diagonalize overlap matrix S and obtain X from 3.167
    3.167: X = S^(-1/2)
4. Guess the density matrix P (first guess is zeros here)
5. Calculate matrix G of 3.154 from P and H2
    G(i,j) = Sum(k, l){P(k,l)[(i j|l k)-1/2(i k|l j)]}
6. Add G to core-Hamiltonian  to get Fock matrix
    3.154: F(i,j) = H1(i,j) + G(i,j)
7. Calculate transformed Fock matrix F' = X'(t)FX
8. Diagonalize F' to obtain C' and epsilon
9. Calculate C = XC'
10. Form new density matrix P from C w/ 3.145
    3.145: P(i,j) = 2 Sum(1-Nelec/2){C(i,a) C*(j,a)}
11. Has P converged to within eps?
    No? -> Step 5 w/ new P from 10.
    Yes? -> Step 12
12. Use resultant solution, represented by C,P,F to calculate outputs
%}
