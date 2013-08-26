function [orb,Eorb,Ehf] = parallelHF(H1,H2,S,Enuc,Nelec,guessDensity,eps,maxIter,minIter)
% Solve Hartree Fock equations
% Input:
%   H1:    H1 matrix
%   H2:    H2 matrix
%   guessDensity:  starting density
%   eps:   convergence criteria
%          [defaults to 1e-8]
% Output:
%   orb:   (nbasis,nbasis)
%            orbital coefficient matrix (atom basis #, mol orb #)
%   Eorb:  (nbasis,1)  molecular orbital energies
%   Ehf:    total Hartree-Fock energy

if (nargin < 3)
    eps = 1.0e-10;
end
if (nargin < 4)
    maxIter = 5000;
end
if (nargin < 5)
    minIter = 5;
end

nbasis = size(H1,1);
% H2j {nbasis,nbasis} cell array of coulomb integrals
% H2k {nbasis,nbasis} cell array of exchange integrals
H2j = cell(nbasis,nbasis);
H2k = cell(nbasis,nbasis);
for i=1:nbasis
    for j=1:nbasis
        H2j{i,j} = reshape(H2(i,j,:,:), nbasis, nbasis);
        H2k{i,j} = reshape(H2(i,:,:,j), nbasis, nbasis);
    end
end


% Step 3 -- Calculate transformation matrix (eq. 3.167).
X = inv(sqrtm(S));

% Step 4 -- Guess at density matrix -- all zeros right now.
Pn = guessDensity;
Plast = Pn;  % previous density matrix
iter = 0;
finished = false;

%changeInDensityForCaution = sqrt(eps);
%changeInDensity = 10 * changeInDensityForCaution;
% Begin iteration through.
while (~finished)  % Step 11 -- Test convergence
    
%    if (changeInDensity > changeInDensityForCaution)
    if (iter < maxIter/2)
        P = 0.5 * Pn + 0.5 * Plast;
    else
        P = Pn;
    end
    
    % update charges
    %    if (iter > 0)
    %       Q = zeros(size(obj.Z));
    %       GAP = zeros(1,obj.natom);
    %       P1 = P.*obj.S;
    %       GOP = sum(P1,1);
    %       for i = 1:obj.natom
    %          GAP(i) = sum(GOP(1,arange{i}));
    %          Q(i) = obj.Z(i)-GAP(i);
    %       end
    %       obj.charges(:,ienv+1) = Q;
    %       H1 = obj.H1(ienv);
    %    end
    
    % Step 5 -- Build 2-electron components of Fock matrix.
    G = twoElecFock(P, H2j, H2k);
    
    
    % Step 6 -- Obtain F (fock matrix).
    F = H1 + G;
    
    % Step 7 -- Calculate the transformed F matrix.
    Ft = X'*F*X;
    
    % Step 8 -- Find e and the transformed expansion coefficient matrices.
    [Ct1,e1] = eig(Ft);
    e2 = diag(e1);
    [e, i1] = sort(e2);
    Ct = Ct1(:,i1);
    
    % Step 9 -- Transform Ct back to C.
    C = X*Ct; %#ok<MINV>
    
    % Step 10 -- Calculate the new density matrix.
    Plast = Pn;
    %Cj = conj(C);
    filled = 1:(Nelec/2);
    Pn = 2* C(:,filled)*( C(:,filled)');
    %     for i = 1:Nbasis
    %         for j = 1:Nbasis
    %             for a = 1:(Nelec/2)
    %                 Pn(i,j) = Pn(i,j) + (C(i,a)*Cj(j,a));
    %             end
    %             Pn(i,j) = Pn(i,j)*2;
    %         end
    %     end
    iter = iter + 1;
    %disp(['den change ',num2str( max(max(abs(P-Pn))))]);
    %disp(['diag Pn ',num2str(diag(Pn)')]);
    %Pn-Plast
    changeInDensity = max(max(abs(P - Pn)));
    %disp(['iter ',num2str(iter),' del ',num2str(changeInDensity)]);
    if (iter > maxIter)
        finished = true;
    elseif (iter > minIter)
        if (changeInDensity < eps)
            finished = true;
        end
    end
end
% End of iteration of steps 5-11.

P = Pn;  % For convenience.

% Step 12: Output.

%Total energy
%3.184: E0 = 1/2 Sum(i,j) {P(j,i)[H1(i,j) + F(i,j)]}
%    Ee = 0;
%    for i = 1:Nbasis
%       for j = 1:Nbasis
%          Ee = Ee + P(j,i)*(H1(i,j)+F(i,j));
%       end
%    end
Ee = sum(sum(P.*(H1+F)));
Ehf = Ee/2 + Enuc;

% Orbital energies.
Eorb = e;

% Molecular orbital components.
orb = C;

% if (iter+1 > maxIter)
%     disp('You are living on the edge.. hartree fock didn''t converge');
% end
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
